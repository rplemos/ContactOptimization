import parser
import argparser
import contacts

from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count, Manager
import resource

def main():

    global_time_start = timer()
    pdb_files, fast, core, show_contacts, mode = argparser.cl_parse()
    file_list = pdb_files
    maximum_distances = {}
    
    if mode == "Single":
        for file in file_list:
            file_time_start = timer()
            if file.endswith(".pdb"):
                protein = parser.parse_pdb(file)
            else:
                protein = parser.parse_pdbx(file)
            contacts_list, _, maximum_distances = contacts.contact_detection(protein, fast, maximum_distances)
            
            if show_contacts:
                contacts.show_contacts(contacts_list)
                
            file_time_end = timer()
            file_time = file_time_end - file_time_start
            print(protein.id, protein.true_count(), len(contacts_list), f"{file_time:.4f}")
                
    elif mode == "Multi":
        if core == 0:
            core = cpu_count()
            
        print(f"Starting processing with {core} cores\n")      
        manager = Manager()
        maximum_distances = manager.dict()

        progress_list = manager.list()
        lock = manager.Lock()
        
        #total_files = len(file_list)  
        with ProcessPoolExecutor(max_workers=core) as executor:
            future_to_file = {executor.submit(process_file, file_path, fast, maximum_distances, lock, progress_list): file_path for file_path in file_list}
            
            for future in as_completed(future_to_file):
                
                try:
                    protein, contacts_list, process_time, maximum_distances = future.result()
                    
                    print(protein.id, protein.true_count(), len(contacts_list), f"{process_time:.4f}")

                    if show_contacts:
                        contacts.show_contacts(contacts_list)
                    
                    ################################################################
                    ##### debugging for # of files processed and memory usage ######
                    #
                    # with lock: 
                    #     processed_files = len(progress_list)
                    #     if processed_files % 100 == 0:
                    #         print()
                    #         print(processed_files,"/",total_files)
                    #         print_memory_usage()
                    #         print()
                    ################################################################
                                                        
                except Exception as e:
                    print(f"Error: {e}")
                    
                finally:
                    del future_to_file[future] # cleans memory to avoid bloating
    
    global_time_end = timer()
    print(f"Total time elapsed: {global_time_end - global_time_start}\n")

    maximum_distances = sorted(maximum_distances.items(), key=lambda x:x[1])
    print(maximum_distances)    

def process_file(file_path, fast, maximum_distances, lock, progress_list):

    file_time_start = timer()
    
    # Update progress list
    with lock:
        progress_list.append(1)
    try:
        if file_path.endswith(".pdb"):
            parsed_data = parser.parse_pdb(file_path)
        else:
            parsed_data = parser.parse_pdbx(file_path)
        contacts_list, _, maximum_distances = contacts.contact_detection(parsed_data, fast, maximum_distances)
        
        file_time_end = timer()
        file_time = file_time_end - file_time_start
        
        return parsed_data, contacts_list, file_time, maximum_distances
    
    except KeyError as e:
        return (file_path, None, e)  # Return tuple with file_path, None result, and exception

def print_memory_usage():
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.0  # Convert to MB
    print(f"Current memory usage: {usage} MB")

if __name__ == "__main__":
    main()