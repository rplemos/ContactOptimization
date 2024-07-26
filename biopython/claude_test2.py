from Bio.PDB import *
from Bio.PDB.NeighborSearch import NeighborSearch
import itertools

# Define contact types and their distance ranges
CONTACT_TYPES = {
    'salt_bridge': (0, 3.9),
    'hydrophobic': (2, 4.5),
    'hydrogen_bond': (0, 3.9),
    'repulsive': (2, 6),
    'attractive': (3.9, 6),
    'disulfide_bond': (0, 2.8),
}

def categorize_contact(atom1, atom2, distance):
    """Categorize the contact between two atoms based on their properties and distance."""
    element1 = atom1.element
    element2 = atom2.element
    residue1 = atom1.get_parent()
    residue2 = atom2.get_parent()

    for contact_type, (min_dist, max_dist) in CONTACT_TYPES.items():
        if min_dist <= distance <= max_dist:
            if contact_type == 'salt_bridge' and ((element1 in ['N', 'O'] and element2 in ['N', 'O']) or
                                                  (element1 in ['C', 'S'] and element2 in ['C', 'S'])):
                return contact_type
            elif contact_type == 'hydrophobic' and element1 == 'C' and element2 == 'C':
                return contact_type
            elif contact_type == 'hydrogen_bond' and ((element1 in ['O', 'N'] and element2 in ['O', 'N', 'H']) or
                                                      (element1 == 'H' and element2 in ['O', 'N'])):
                return contact_type
            elif contact_type == 'repulsive' and ((element1 in ['N', 'O'] and element2 in ['N', 'O']) or
                                                  (element1 in ['C', 'S'] and element2 in ['C', 'S'])):
                return contact_type
            elif contact_type == 'attractive' and ((element1 in ['N', 'O'] and element2 in ['C', 'S']) or
                                                   (element1 in ['C', 'S'] and element2 in ['N', 'O'])):
                return contact_type
            elif contact_type == 'disulfide_bond' and (residue1.resname == 'CYS' and residue2.resname == 'CYS' and
                                                       element1 == 'S' and element2 == 'S'):
                return contact_type
    return None

def find_contacts(structure, distance_cutoff=6.0):
    """Find and categorize contacts in the structure."""
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    contacts = {}

    for atom1, atom2 in itertools.combinations(atoms, 2):
        # Skip if atoms are from the same residue
        if atom1.get_parent() == atom2.get_parent():
            continue

        distance = atom1 - atom2
        if distance <= distance_cutoff:
            contact_type = categorize_contact(atom1, atom2, distance)
            if contact_type:
                key = tuple(sorted([atom1.get_full_id(), atom2.get_full_id()]))
                contacts[key] = (contact_type, distance)

    return contacts

# Load the structure
parser = MMCIFParser()
structure = parser.get_structure("protein", "../8uw4.cif")

# Find contacts
contacts = find_contacts(structure)
