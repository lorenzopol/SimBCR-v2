import os
import gzip
from Bio.PDB import PDBParser
from io import StringIO


def handle_gzip(look_up_directory):
    files = [f for f in os.listdir(look_up_directory) if
             os.path.isfile(os.path.join(look_up_directory, f)) and f.endswith(".gz")]
    files_sorted_by_time = sorted(files, key=lambda x: os.path.getmtime(os.path.join(look_up_directory, x)),
                                  reverse=True)
    file = files_sorted_by_time[0]

    # dispatch file to cwd
    input_file_path = os.path.join(look_up_directory, file)
    os.rename(input_file_path, input_file_path.replace(" ", ""))
    input_file_path = input_file_path.replace(" ", "")
    output_filename = ".".join(file.split(r"\ ")[-1].split(".")[:-1])

    # dump content
    with gzip.open(input_file_path, 'rb') as f_in:
        with open(output_filename, 'wb') as f_out:
            # Read and decompress the content from the input .gz file and write it to the output file
            f_out.write(f_in.read())
    return output_filename


def is_valid_pdb(input_string):
    try:
        parser = PDBParser()
        pdb_file = StringIO(input_string)
        _ = parser.get_structure('temp_structure', pdb_file)
        return True  # Valid PDB format
    except Exception as e:
        print(f"Not a valid PDB file: {e}")
        return False  # Not a valid PDB format


def show_3D_from_pdb(filename, mode=None, cdr3_range=None):
    """expects a .pdb files and renders it in html-ipynb"""
    assert ".pdb" in filename, f"ERROR: pdb_render expected a .pdb file. {filename} was given"
    os.system(f"python show_3D.py {filename} --mode {mode} --cdr3 {cdr3_range}")

def parse_pdb_content(pdb_file_path):
    # Create a PDB parser
    parser = PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure('pdb_structure', pdb_file_path)
    print(f"{structure.center_of_mass() = }")

    for model in structure:
        print(f"{model.child_dict = }")

        for chain in model:
            print(f"{chain.get_unpacked_list() = }")

            for residue in chain:
                print(f"{residue}")
                print(f"{residue.child_dict = }")

                for atom in residue:
                    print(f"    {atom.name = }")
                    if hasattr(atom, "charge"):
                        print(f"    {atom.charge = }")
                    print(f"    {atom.coord = }")
                    print(f"    {atom.radius = }")

                    print(f"    {atom.altloc = }")

                    if hasattr(atom, "anisou"):
                        print(f"    {atom.anisou = }")

                    if hasattr(atom, "sigatm"):
                        print(f"    {atom.sigatm = }")

                    if hasattr(atom, "siguij"):
                        print(f"    {atom.siguij = }")

                    print(f"    {atom.bfactor = }")

                    print(f"    {atom.level = }")
                    print(f"    {atom.occupancy = }")
                    print(f"    {atom.parent = }")
                    if hasattr(atom, "vector"):
                        print(f"    {atom.vector = }")

                    print(f"    {atom.is_disordered() = }")
                    print("")
                print("")
