import parser
import argparser_core
import contacts_fast
import os

from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count, Manager
import resource
import psutil

def set_affinity(core_id):
    p = psutil.Process(os.getpid())
    p.cpu_affinity([core_id])


def main():

    global_time_start = timer()
    pdb_files, fast, core, show_contacts, mode, cnum = argparser_core.cl_parse()
    file_list = pdb_files
    
    core_id = cnum  
    set_affinity(int(core_id))
    print(f"Running on core {core_id}")
        
    if mode == "Single":
        for file in file_list:
            file_time_start = timer()
            if file.endswith(".pdb"):
                protein = parser.parse_pdb(file)
            else:
                protein = parser.parse_pdbx(file)

            contacts_list, _ = contacts_fast.contact_detection(protein, fast)
            
            if show_contacts:
                contacts_fast.show_contacts(contacts_list)
                
            file_time_end = timer()
            file_time = file_time_end - file_time_start
            print(protein.id, protein.true_count(), len(contacts_list), f"{file_time:.4f}")
                
    elif mode == "Multi":
        if core == 0:
            core = cpu_count()
            
        print(f"Starting processing with {core} cores\n")      
        manager = Manager()

        progress_list = manager.list()
        lock = manager.Lock()
        
        #total_files = len(file_list)  
        with ProcessPoolExecutor(max_workers=core) as executor:
            future_to_file = {executor.submit(process_file, file_path, fast, lock, progress_list): file_path for file_path in file_list}
            
            for future in as_completed(future_to_file):
                
                try:
                    protein, contacts_list, process_time = future.result()
                    
                    print(protein.id, protein.true_count(), len(contacts_list), f"{process_time:.4f}")

                    if show_contacts:
                        contacts_fast.show_contacts(contacts_list)
                                                        
                except Exception as e:
                    print(f"Error: {e}")
                    
                finally:
                    del future_to_file[future] # cleans memory to avoid bloating
    
    global_time_end = timer()
    print(f"Total time elapsed: {global_time_end - global_time_start}\n")  

def process_file(file_path, fast, lock, progress_list):

    file_time_start = timer()
    
    # Update progress list
    with lock:
        progress_list.append(1)
    try:
        if file_path.endswith(".pdb"):
            parsed_data = parser.parse_pdb(file_path)
        else:
            parsed_data = parser.parse_pdbx(file_path)
        contacts_list, _ = contacts_fast.contact_detection(parsed_data, fast)
        
        file_time_end = timer()
        file_time = file_time_end - file_time_start
        
        return parsed_data, contacts_list, file_time
    
    except KeyError as e:
        return (file_path, None, e)  # Return tuple with file_path, None result, and exception

def print_memory_usage():
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.0  # Convert to MB
    print(f"Current memory usage: {usage} MB")

if __name__ == "__main__":
    main()