from Bio.PDB import MMCIFParser, NeighborSearch
import conditions
from timeit import default_timer as timer

start = timer()

def process_contacts(neighbor_search, categories, contact_conditions, contact_types):
    contacts = {}
    processed_pairs = set()

    for atom1, atom2 in neighbor_search.search_all(6.0):  # Use the maximum distance from categories
        # Skip if atoms are from the same residue
        if atom1.parent is atom2.parent:
            continue

        # Create a unique identifier for the atom pair
        pair_id = tuple(sorted((atom1.full_id, atom2.full_id)))
        
        # Skip if this pair has already been processed
        if pair_id in processed_pairs:
            continue
        
        processed_pairs.add(pair_id)

        # Get atom names in the format "RESNAME:ATOMNAME"
        name1 = f"{atom1.parent.resname}:{atom1.name}"
        name2 = f"{atom2.parent.resname}:{atom2.name}"

        # Skip if either atom is not in contact_types
        if name1 not in contact_types or name2 not in contact_types:
            continue

        distance = atom1 - atom2

        for category, (min_dist, max_dist) in categories.items():
            if min_dist <= distance <= max_dist:
                condition = contact_conditions[category]
                if condition(name1, name2):
                    if category not in contacts:
                        contacts[category] = []
                    if category == 'hydrogen_bond' and (abs(atom2.get_parent().get_id()[1] - atom1.get_parent().get_id()[1]) <= 3):
                        continue
                    contacts[category].append((atom1, atom2, distance))

    return contacts

# Usage example:
parser = MMCIFParser()
structure = parser.get_structure("protein", "../8uw4.cif")
neighbor_search = NeighborSearch(list(structure.get_atoms()))

contacts = process_contacts(neighbor_search, conditions.categories, conditions.contact_conditions, conditions.contact_types)

print(len(contacts))
for key in contacts.keys():
    print(key, len(contacts[key]))

# for category, contact_list in contacts.items():
#     print(f"{category}:")
#     for atom1, atom2, distance in contact_list:
#         print(f"  {atom1.full_id} - {atom2.full_id}: {distance:.2f} Å")


end = timer()
print(f"Total time elapsed: {end - start}\n")