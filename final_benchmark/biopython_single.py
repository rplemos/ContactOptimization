from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser, is_aa
import conditions_biop
from timeit import default_timer as timer
import os
import sys
import psutil

def set_affinity(core_id):
    p = psutil.Process(os.getpid())
    p.cpu_affinity([core_id])

core_id = sys.argv[2]   
set_affinity(int(core_id))
print(f"Running on core {core_id}")

global_start = timer()

folder = sys.argv[1]

files = []
for entry in os.scandir(folder):
    files.append(entry.name)

files = sorted(files)

for file in files:
    start = timer()

    try:
        if file.endswith(".pdb"):
            parser = PDBParser(QUIET=True) 
        elif file.endswith(".cif"):
            parser = MMCIFParser(QUIET=True) 
            
        file_path = f"{folder}{file}"
        full_structure = parser.get_structure('structure_id', file_path)
        structure = next(full_structure.get_models())
        
    except Exception as e:
        print(f"Error processing file {file}: {e}")
        continue
    
    protein_size = (len([_ for _ in structure.get_residues() if is_aa(_)]))
    if protein_size > 5000:
        print(f"Skipping ID '{file}'. Size: {protein_size} residues")
        continue

    # Filter out waters and heteroatoms
    atoms = [atom for atom in structure.get_atoms() if not atom.get_parent().get_resname() == 'HOH' and not atom.get_parent().id[0] == 'W']

    ns = NeighborSearch(atoms)
    radius = 6.0
    interactions = []

    for atom in atoms:
        center_coord = atom.coord
        neighbors = ns.search(center_coord, radius)
        residue_number = atom.get_parent().get_id()[1]
        chain1 = atom.get_parent().get_parent().get_id()
        for neighbor in neighbors:
            neighbor_number = neighbor.get_parent().get_id()[1]
            chain2 = neighbor.get_parent().get_parent().get_id()
            if atom != neighbor and residue_number < neighbor_number:
                atom1_name = f"{atom.get_parent().get_resname()}:{atom.get_name()}"
                atom2_name = f"{neighbor.get_parent().get_resname()}:{neighbor.get_name()}"
                distance = atom - neighbor
                
                if atom1_name in conditions_biop.contact_types and atom2_name in conditions_biop.contact_types:
                    for interaction, condition in conditions_biop.contact_conditions.items():
                        if interaction == 'hydrogen_bond' and (abs(neighbor_number - residue_number) <= 3):
                            continue
                        min_distance, max_distance = conditions_biop.categories[interaction]
                        if min_distance <= distance <= max_distance and condition(atom1_name, atom2_name):
                            interactions.append((atom1_name, residue_number, atom2_name, neighbor_number, interaction, distance, chain1, chain2))  
        
    # Print interactions
    interactions = sorted(interactions, key=lambda x:x[4])
    # for atom1_name, residue_number, atom2_name, neighbor_number, interaction, distance, chain1, chain2 in interactions:
    #     print(f"{chain1}-{residue_number}{atom1_name} and {chain2}-{neighbor_number}{atom2_name}: {distance:.2f} A. {interaction}")
    # for atom1_name, residue_number, atom2_name, neighbor_number, interaction, distance, chain1, chain2 in interactions:
    #     print(f"{chain1}-{residue_number+20}{atom1_name} and {chain2}-{neighbor_number+20}{atom2_name}: {distance:.2f} A. {interaction}")


    end = timer()
    time = end - start
    file = file.split(".")[0]
    print(file, protein_size, len(interactions), f"{time:.4f}")

global_end = timer()
print()
print(f"Total time: {global_end - global_start}")