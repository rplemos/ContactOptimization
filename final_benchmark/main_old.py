import parser
import argparser
import contacts
import distances

from os import getpid
from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from psutil import Process
from numpy import zeros

def main():

    global_time_start = timer()
    file_list, fast, core, mode, cnum = argparser.cl_parse()
    
    if cnum:
        p = Process(getpid())
        p.cpu_affinity([cnum])
        print(f"Running on core {cnum}")
            
    # Get the list of pairs in the order they appear in the dictionary
    pairs = list(distances.distances.keys())

    distance_array = zeros(len(pairs))

    # Create a mapping from pairs to array indices
    pair_index = {pair: idx for idx, pair in enumerate(pairs)}
    
    # Fill the array with distances from the dictionary
    for (res1, res2), distance in distances.distances.items():
        pair = ''.join(sorted((res1, res2)))
        distance_array[pair_index[pair]] = distance
        
    if mode == "Single":
        
        for file in file_list:
            file_time_start = timer()
            if file.endswith(".pdb"):
                protein = parser.parse_pdb(file)
            else:
                protein = parser.parse_pdbx(file)
                
            if protein.true_count() > 25000:
                print(f"Skipping ID '{protein.id}'. Size: {protein.true_count()} residues")
                continue

            contacts_list, _ = contacts.contact_detection(protein, fast, pair_index, distance_array)
                
            file_time_end = timer()
            file_time = file_time_end - file_time_start
            print(protein.id, protein.true_count(), len(contacts_list), f"{file_time:.4f}")
                
    elif mode == "Multi":
            
        core = cpu_count() if core == 0 else core
        print(f"Starting processing with {core} cores\n")      
        
        with ProcessPoolExecutor(max_workers=core) as executor:
            future_to_file = {executor.submit(process_file, file_path, fast, pair_index, distance_array): file_path for file_path in file_list}
            
            for future in as_completed(future_to_file):
                try:
                    protein, contacts_list, process_time = future.result()
                    if protein is None: # skipped proteins
                        continue

                    print(protein.id, protein.true_count(), len(contacts_list), f"{process_time:.4f}")
                                                        
                except Exception as e:
                    print(f"Error: {e}")
                    
                finally:
                    del future_to_file[future] # cleans memory to avoid bloating
    
    global_time_end = timer()
    print(f"Total time elapsed: {global_time_end - global_time_start}\n")  

def process_file(file_path, fast, pair_index, distance_array):

    file_time_start = timer()
    
    try:
        if file_path.endswith(".pdb"):
            parsed_data = parser.parse_pdb(file_path)
        else:
            parsed_data = parser.parse_pdbx(file_path)
            
        if parsed_data.true_count() > 25000:
            print(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues")
            return None, None, None
        
        contacts_list, _ = contacts.contact_detection(parsed_data, fast, pair_index, distance_array)
        
        file_time_end = timer()
        file_time = file_time_end - file_time_start
        
        return parsed_data, contacts_list, file_time
    
    except KeyError as e:
        return (file_path, None, e)  # Return tuple with file_path, None result, and exception

if __name__ == "__main__":
    main()