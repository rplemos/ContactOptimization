from classes import Protein, Chain, Residue, Atom

from numpy import mean, array
from numpy.linalg import svd
from math import dist

stacking = {
    'HIS':[10, 'CG','ND1','CE1','NE2','CD2'],
    'PHE':[11, 'CG','CD1','CE1','CZ','CE2','CD2'],
    'TRP':[14, 'CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
    'TYR':[12, 'CG','CD1','CE1','CZ','CE2','CD2'],
}


def parse_pdb(pdb_file):
    """
    Parses the PDB files and constructs Protein, Chain, Residue, and Atom objects.

    Args:
        pdb_files (list): List of PDB file names to parse.

    Returns:
        list: List of Protein objects representing the parsed proteins.
    """

    current_protein = Protein()
    current_chain = None
    current_residue = None
    
    valid_residues =  ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                        'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    with open(pdb_file) as f:
        for line in f:
            line = line.strip()
            
            if line == "ENDMDL":
                return current_protein

            if line.startswith("HEADER") and ("RNA" in line or "DNA" in line):
                current_protein.id = pdb_file.split("/")[-1][:4]
                current_protein.title = "DNA/RNA"
                break
            
            elif line.startswith("HEADER"):
                current_protein.id = line[62:]
                
            elif line.startswith("TITLE"):
                current_protein.set_title(line[10:])
                
            elif line.startswith("ATOM"):
                chain_id = line[21]
                resnum = int(line[22:26])
                if resnum <= 0:
                    continue
                resname = line[17:20]
                
                if resname == "HIE" or resname == "HID":  # alternative names for protonated histidines
                    resname = "HIS" 
                
                if resname not in valid_residues:
                    continue                       

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)

                if current_residue is None:  # new residue
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                    current_chain.residues.append(current_residue)
                
                if current_residue.resnum != resnum:
                    if len(current_residue.atoms) > 1:
                        current_chain.residues.append(current_residue) 
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                                                                
                atomname = line[12:16].replace(" ", "")
                if atomname == "OXT": # OXT is the C-terminal Oxygen atom. However, it exhibits the same properties of any Oxygen
                    atomname = atomname.replace("OXT","O")                
                
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occupancy = float(line[55:60])
                
                if occupancy >= 0.5: # ignores low quality atoms (arbitrary value!)
                    if current_residue.atoms:
                        distance = dist((current_residue.atoms[-1].x, current_residue.atoms[-1].y, current_residue.atoms[-1].z), (x, y, z))
                        if distance > 10:
                            # print(f"BIZARRE {distance:.3f}", atomname, occupancy, current_residue.resname, current_residue.resnum, current_chain.id, current_protein.id)
                            # print(current_residue.atoms[-1].x, current_residue.atoms[-1].y, current_residue.atoms[-1].z)
                            # print(x, y, z)
                            continue
                    atom = Atom(atomname, x, y, z, occupancy, current_residue) # creates atom
                    current_residue.atoms.append(atom)
                else:
                    continue
                    #print(atomname, occupancy, current_residue.resname, current_residue.resnum, current_chain.id, current_protein.id)

                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len(current_residue.atoms) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms[5:]]) # ignores [N, CA, C, O] and [RNG] atoms
                        centroid_atom = centroid(current_residue, ring_atoms)
                        current_residue.atoms.append(centroid_atom)
                        current_residue.ring = True

                        normal_vector = calc_normal_vector(ring_atoms)
                        current_residue.normal_vector = normal_vector

            elif line.startswith("END"):  
                # Handling cases where there is no ID
                if current_protein.id is None:
                    id = str(pdb_file).split("/")[-1]
                    id = id.split(".")[0]
                    current_protein.id = id  

    return current_protein


def parse_pdbx(pdbx_file):
            
    valid_residues =  ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                        'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
    valid_atoms = ['N', 'C', 'O', 'S']
    
    current_protein = Protein()
    current_chain = None
    current_residue = None
    atomsite_block = False # _atom_site. lines
    atominfo_block = False # ATOM        lines
    atom_lines = []

    with open(pdbx_file) as f:
        for line in f:
            line = line.strip()

            if line.startswith("_entry.id"):
                current_protein.id = line[-4:]
                
            elif line.startswith("TITLE"):
                current_protein.set_title(line[10:])
                            
            elif line.startswith("_atom_site.group_PDB"):
                atomsite_block = True
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and line.startswith("_atom_site"):
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and line.startswith("ATOM"):                                
                atomname_index = atom_lines.index("label_atom_id")
                resname_index = atom_lines.index("label_comp_id")
                chain_index = atom_lines.index("label_asym_id")
                resnum_index = atom_lines.index("label_seq_id")
                x_index = atom_lines.index("Cartn_x")
                y_index = atom_lines.index("Cartn_y")
                z_index = atom_lines.index("Cartn_z")
                occupancy_index = atom_lines.index("occupancy")
                model_index = atom_lines.index("pdbx_PDB_model_num")
                alt_occupancy_index = atom_lines.index("label_alt_id") # . if occupancy == 1, varies otherwise
                atom_element_index = atom_lines.index("type_symbol")
                                                               
                atomsite_block = False
                atominfo_block = True
                #print(atomname_index, resname_index, chain_index, resnum_index, x_index, y_index, z_index, occupancy_index, model_index)
                
            elif line.startswith("ATOM") and atominfo_block:
                line = line.split()
                
                element = line[atom_element_index]
                if element not in valid_atoms:
                    continue
                
                model = int(line[model_index])
                if model != 1: 
                    return current_protein
                
                chain_id = line[chain_index]
                
                resnum = int(line[resnum_index])
                if resnum <= 0:
                    continue
                resname = line[resname_index]

                if resname == "HIE" or resname == "HID":  # alternative names for protonated histidines
                    resname = "HIS"  

                if resname not in valid_residues:
                    continue                            

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)
                    current_residue = None

                if current_residue is None:  # new residue
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                    current_chain.residues.append(current_residue)
                
                if current_residue.resnum != resnum:
                    if len(current_residue.atoms) > 1:
                        current_chain.residues.append(current_residue) 
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                    #current_chain.residues.append(current_residue)
                                                                
                atomname = line[atomname_index]
                if atomname == "OXT": # OXT is the C-terminal Oxygen atom. However, it exhibits the same properties of any Oxygen
                    continue
                    atomname = atomname.replace("OXT","O")    
                    
                x, y, z = float(line[x_index]), float(line[y_index]), float(line[z_index])
                occupancy = float(line[occupancy_index])
                    
                if occupancy >= 0.5 and (line[alt_occupancy_index] == "." or line[alt_occupancy_index] == "A"): # ignores low quality atoms (arbitrary value!)
                    if current_residue.atoms:
                        distance = dist((current_residue.atoms[-1].x, current_residue.atoms[-1].y, current_residue.atoms[-1].z), (x, y, z))
                        if distance > 5:
                            # print(f"\tBIZARRE {distance:.3f}", atomname, occupancy, current_residue.resname, current_residue.resnum, current_chain.id, current_protein.id)
                            # print("\tAtom1:", current_residue.atoms[-1].x, current_residue.atoms[-1].y, current_residue.atoms[-1].z, current_residue.atoms[-1].atomname)
                            # print("\tAtom2:", x, y, z, atomname)
                            continue
                    atom = Atom(atomname, x, y, z, occupancy, current_residue) # creates atom
                    current_residue.atoms.append(atom)
                else:
                    continue
                                
                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len(current_residue.atoms) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms[5:]]) # ignores [N, CA, C, O] and [RNG] atoms
                        centroid_atom = centroid(current_residue, ring_atoms)
                        current_residue.atoms.append(centroid_atom)
                        current_residue.ring = True

                        normal_vector = calc_normal_vector(ring_atoms)
                        current_residue.normal_vector = normal_vector

            elif atominfo_block and line == "#":
                if resname in valid_residues:
                    current_chain.residues.append(current_residue) # appends the last residue
                atominfo_block = False 
    
    return current_protein


def centroid(residue, ring_atoms):
    centroid = mean(ring_atoms, axis = 0)
    centroid_atom = Atom("RNG", centroid[0], centroid[1], centroid[2], 1, residue)
    
    return centroid_atom


def calc_normal_vector(ring_atoms):
    centroid = mean(ring_atoms, axis = 0) # axis=0 -> mean through the columns
    centered_ring_atoms = ring_atoms - centroid # normalizes to origin
    
    # Use singular value decomposition (SVD) to calculate the plane
    _, _, vh = svd(centered_ring_atoms) # vh = V^T
    normal_vector = vh[2]  # The normal vector is the last row of the V^T matrix
    
    return normal_vector