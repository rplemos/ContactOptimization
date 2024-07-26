from math import dist
from timeit import default_timer as timer
from numpy import dot, arccos, degrees
from numpy.linalg import norm

from classes import Contact
import conditions
import distances

def contact_detection(protein, fast, maximum_distances):
    start = timer()
    
    residues = list(protein.get_residues())
    contacts = []
    
    for i, residue1 in enumerate(residues[1:]):
        for j, residue2 in enumerate(residues[i+1:], start=i+1):
            
            if residue1.resnum == residue2.resnum and residue1.chain.id == residue2.chain.id: # ignores same residue
                continue
            
            if len(residue1.atoms) > 1 and len(residue2.atoms) > 1:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons

                distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
            
                # pair = tuple(sorted((residue1.resname, residue2.resname)))
                # if distance_ca > (distances.distances[pair] + 0.01):
                #    continue

                if distance_ca > 21:
                    continue
                
            else:
                continue              
            
            # CHECKING FOR AROMATIC STACKINGS
            if residue1.ring and residue2.ring:
                ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                if distance >= 2 and distance <= 5: # within aromatic stacking limits
                    if (160 <= angle < 180) or (0 <= angle < 20):
                        stack_type = "-parallel"
                    elif (80 <= angle < 100):
                        stack_type = "-perpendicular"
                    else:
                        stack_type = "-other"

                    contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, ring1.atomname, 
                                    protein.id, residue2.chain.id, residue2.resnum, residue2.resname, ring2.atomname, 
                                    float(f"{distance:.2f}"), "stacking"+stack_type, ring1, ring2)
                    
                    contacts.append(contact)
                    
            for atom1 in residue1.atoms:
                for atom2 in residue2.atoms:
                    name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from contacts dictionary
                    name2 = f"{atom2.residue.resname}:{atom2.atomname}"
                    
                    if name1 in conditions.contact_types and name2 in conditions.contact_types: # excludes the RNG atom and any different other
                        
                        distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                        
                        if distance <= 6: # max distance for contacts
                            for contact_type, distance_range in conditions.categories.items():
                                
                                if not fast:
                                    if contact_type == 'hydrogen_bond' or contact_type == 'hydrophobic':
                                        continue                                
                                
                                if contact_type == 'hydrogen_bond' and (abs(residue2.resnum - residue1.resnum) <= 3): # skips alpha-helix for h-bonds
                                    continue
                                
                                if distance_range[0] <= distance <= distance_range[1]: # fits the range
                                    if conditions.contact_conditions[contact_type](name1, name2): # fits the type of contact
                                                                                                
                                        contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, atom1.atomname, 
                                                        protein.id, residue2.chain.id, residue2.resnum, residue2.resname, atom2.atomname, 
                                                        float(f"{distance:.2f}"), contact_type, atom1, atom2)

                                        contacts.append(contact)
                                    
                                        # ####################
                                        # # BLOCK FOR CONSTRUCTING MAXIMUM DISTANCES LIST                   
                                        if (residue1.resname, residue2.resname) not in maximum_distances:
                                            maximum_distances[residue1.resname, residue2.resname] = [float(f"{distance_ca:.2f}"), float(f"{distance:.2f}"), residue1.resnum, residue2.resnum, protein.id, atom1.atomname, atom2.atomname, residue1.chain.id, residue2.chain.id, contact_type]
                                            #print(f"Creating {residue1.resname, residue2.resname} : {distance_ca: .2f} | {residue1.resnum}, {residue2.resnum}, {protein.id}, {atom1.atomname}, {atom2.atomname}, {residue1.chain.id}, {residue2.chain.id}, {contact_types[0]}")
                                        elif distance_ca > maximum_distances[residue1.resname, residue2.resname][0]:
                                            #print(f"Changing {residue1.resname, residue2.resname} from {maximum_distances[residue1.resname, residue2.resname][0]} to {distance_ca: .2f}  | {residue1.resnum}, {residue2.resnum}, {protein.id}, {atom1.atomname}, {atom2.atomname}, {residue1.chain.id}, {residue2.chain.id}, {contact_types[0]}")
                                            maximum_distances[residue1.resname, residue2.resname] = [float(f"{distance_ca:.2f}"), float(f"{distance:.2f}"), residue1.resnum, residue2.resnum, protein.id, atom1.atomname, atom2.atomname, residue1.chain.id, residue2.chain.id, contact_type]
                                        # # ###################
                                                                                                                
    end = timer()
    current_time = end - start

    return contacts, current_time, maximum_distances


def show_contacts(contacts):
    category_counts = {}

    for contact in contacts:
        category = contact.type
        category_counts[category] = category_counts.get(category, 0) + 1
    
    sorted_categories = sorted(category_counts.items(), key=lambda x: x[1])

    for category, count in sorted_categories:
        print(f"\nNumber of {category} occurrences:", count)
        if count <= 9999:
            print(f"All entries for {category}:")
            for entry in contacts:
                if entry.type == category:
                    print("\t",entry.print_text())

def calc_angle(vector1, vector2):
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product
    angle = arccos(dot_product / magnitude_product) # angle in radians   
    
    return degrees(angle)
