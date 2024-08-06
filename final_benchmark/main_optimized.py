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
    file_list, core, mode, cnum = argparser.cl_parse()
    
    if cnum:
        p = Process(getpid())
        p.cpu_affinity([int(cnum)])
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
        single(file_list, pair_index, distance_array)
    elif mode == "Multi":
        multi(file_list, pair_index, distance_array, core)
    
    print(f"Total time elapsed: {timer() - global_time_start}\n")


def single(file_list, pair_index, distance_array):
    for file in file_list:
        try:
            result = process_file(file, pair_index, distance_array)
            if result:
                protein, contacts_list, process_time = result
                print(protein.id, protein.true_count(), len(contacts_list), f"{process_time:.4f}")
        except Exception as e:
            print(f"Error: {e}")
            
            
def multi(file_list, pair_index, distance_array, core):
    core = cpu_count() if core == 0 else core
    print(f"Starting processing with {core} cores") 
    
    with ProcessPoolExecutor(max_workers=core) as executor:
        futures = {executor.submit(process_file, file, pair_index, distance_array): file for file in file_list}
        
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    protein, contacts_list, process_time = result
                    print(protein.id, protein.true_count(), len(contacts_list), f"{process_time:.4f}")
            except Exception as e:
                print(f"Error: {e}")
            finally:
                del futures[future] # cleans memory to avoid bloating
           
                
def process_file(file_path, pair_index, distance_array):
    start_time = timer()
    
    try:
        parsed_data = parser.parse_pdb(file_path) if file_path.endswith(".pdb") else parser.parse_pdbx(file_path)
            
        if parsed_data.true_count() > 25000:
            print(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues")
            return None

        contacts_list = contacts.contact_detection(parsed_data, pair_index, distance_array)
        
        process_time = timer() - start_time
        return parsed_data, contacts_list, process_time
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None


if __name__ == "__main__":
    main()