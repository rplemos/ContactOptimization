from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser, is_aa
import conditions_biop
from timeit import default_timer as timer
import os
import sys
from psutil import Process

core_id = sys.argv[2]   
p = Process(os.getpid())
p.cpu_affinity([int(core_id)])
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
        parser = PDBParser(QUIET=True) if file.endswith(".pdb") else MMCIFParser(QUIET=True)
        structure = next(parser.get_structure('structure_id', os.path.join(folder, file)).get_models())
    except Exception as e:
        print(f"Error processing file {file}: {e}")
        continue
    
    protein_size = sum(1 for _ in structure.get_residues() if is_aa(_))
    if protein_size > 5000:
        print(f"Skipping ID '{file}'. Size: {protein_size} residues")
        continue

    # Filter out waters and heteroatoms
    atoms = [atom for atom in structure.get_atoms() if not atom.get_parent().get_resname() == 'HOH' and not atom.get_parent().id[0] == 'W']

    ns = NeighborSearch(atoms)
    radius = 6.0
    interactions = []

    for atom in atoms:
        residue_number = atom.get_parent().get_id()[1]
        chain1 = atom.get_parent().get_parent().get_id()
        atom1_name = f"{atom.get_parent().get_resname()}:{atom.get_name()}"
        
        if atom1_name not in conditions_biop.contact_types:
            continue

        for neighbor in ns.search(atom.coord, radius):
            neighbor_number = neighbor.get_parent().get_id()[1]
            chain2 = neighbor.get_parent().get_parent().get_id()
            
            if atom == neighbor or residue_number >= neighbor_number:
                continue

            atom2_name = f"{neighbor.get_parent().get_resname()}:{neighbor.get_name()}"
            if atom2_name not in conditions_biop.contact_types:
                continue

            distance = atom - neighbor
            
            if atom1_name in conditions_biop.contact_types and atom2_name in conditions_biop.contact_types:
                for interaction, condition in conditions_biop.contact_conditions.items():
                    if interaction == 'hydrogen_bond' and (abs(neighbor_number - residue_number) <= 3):
                        continue
                    min_distance, max_distance = conditions_biop.categories[interaction]
                    if min_distance <= distance <= max_distance and condition(atom1_name, atom2_name):
                        interactions.append((atom1_name, residue_number, atom2_name, neighbor_number, interaction, distance, chain1, chain2))  
        
    interactions = sorted(interactions, key=lambda x:x[4])

    file = file.split(".")[0]
    time = timer() - start
    print(file, protein_size, len(interactions), f"{time:.4f}")

print(f"Total time: {timer() - global_start}")