from sys import exit
from argparse import ArgumentParser, ArgumentError, ArgumentTypeError

def cl_parse():
    try:
        parser = ArgumentParser(description='PDB/mmcif parser and fast contact detection')
        parser.add_argument('-pdb', nargs='+', required=True, type=validate_file, help='List of PDB files (at least one required)')
        parser.add_argument('-fast', action='store_false', required=False, help='Set if Hydrogen Bond and Hydrophobic contacts are calculated')
        parser.add_argument('-show', '--show_contacts', required=False, action='store_true', help='Shows a list of all contacts')
        parser.add_argument('-core', type=int, required=False, default=0, help='Number of cores to use (only needed on Multi mode)')
        parser.add_argument('-mode', required=False, default='Single', help='Select "SingleCore" or "MultiCore" mode')
        parser.add_argument('-cnum', required=False, help='Select specific core')

        args = parser.parse_args()

        pdb_files = args.pdb
        fast = args.fast
        core = args.core
        show_contacts = args.show_contacts
        mode = args.mode
        modes = ["Single", "Multi"]
        if mode not in modes:
            raise ValueError("Invalid Mode!")
        cnum = args.cnum
        
    except ArgumentError as e:
        print(f"Argument Error: {str(e)}")
        exit(1)

    except ValueError as e:
        print(f"Error: {str(e)}")
        exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        exit(1)
    
    return pdb_files, fast, core, show_contacts, mode, cnum
        
def validate_file(value):
    if value.endswith('.pdb') or value.endswith('.pdbx') or value.endswith('.cif'):
        return value
    else:
        raise ArgumentTypeError(f"{value} is not a valid file. File must end with '.pdb', '.pdbx', or '.cif'")
