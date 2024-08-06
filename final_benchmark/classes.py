class Protein:
    def __init__(self):
        self.title = None
        self.id = None
        self.chains = []

    def set_title(self, title):
        if self.title is None:
            self.title = title
        else:
            self.title += " " + title.strip()

    def get_chains(self):
        for chain in self.chains:
            yield chain

    def get_residues(self):
        for chain in self.get_chains():
            for residue in chain.residues:
                yield residue
                
    def true_count(self):
        return sum(1 for _ in self.get_residues())
    
    def full_count(self):
        chain_residues = {}
        current_size = 0
        total_size = 0
        
        for chain in self.chains:
            residue_count = chain.count_residues()
            chain_residues[chain.id] = current_size
            current_size += residue_count
            total_size += residue_count

        return chain_residues, total_size
            

class Chain:
    def __init__(self, id, residues):
        self.id = id
        self.residues = residues
    
    def count_residues(self):
        return self.residues[-1].resnum if len(self.residues) > 1 else 0


class Residue:
    def __init__(self, resnum, resname, atoms, chain, ring, normal_vector):
        self.resnum = resnum
        self.resname = resname
        self.atoms = atoms
        self.chain = chain
        self.ring = ring
        self.normal_vector = normal_vector
        
class Atom:
    def __init__(self, atomname, x, y, z, occupancy, residue):
        self.atomname = atomname
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.residue = residue


class Contact: 
    def __init__(self, id1, chain1, residue_num1, residue_name1, atom1, 
                 id2, chain2, residue_num2, residue_name2, atom2, 
                 distance, type, atom_object1, atom_object2):
        self.id1 = id1
        self.chain1 = chain1
        self.residue_num1 = residue_num1
        self.residue_name1 = residue_name1
        self.atom1 = atom1
        self.id2 = id2
        self.chain2 = chain2
        self.residue_num2 = residue_num2
        self.residue_name2 = residue_name2
        self.atom2 = atom2
        self.distance = distance
        self.type = type
        self.atom_object1 = atom_object1
        self.atom_object2 = atom_object2
    
    def print_values(self):
        all_values = list(self.__dict__.values())
        return [f"{all_values[0]}:{all_values[1]}", f"{all_values[2]}{all_values[3]}:{all_values[4]}",
                f"{all_values[5]}:{all_values[6]}", f"{all_values[7]}{all_values[8]}:{all_values[9]}",
                all_values[10], all_values[11]]
    
    def print_text(self):
        all_values = list(self.__dict__.values())
        return f"{all_values[1]}-{all_values[2]}{all_values[3]}:{all_values[4]} and {all_values[6]}-{all_values[7]}{all_values[8]}:{all_values[9]}: {all_values[10]} A. {all_values[11].capitalize()}"

