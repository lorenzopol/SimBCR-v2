import os
import gzip
from Bio.PDB import PDBParser
from BiologicalConstants import *


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


def pdb_to_obj(pdb_file, output_file, return_pos=False,
               add_relative_radius=False) -> list[tuple[int, int, int]] | None:
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    atom_coord = []

    with open(output_file, 'w') as obj_file:
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        if return_pos:
                            atom_coord.append((x, y, z))
                        obj_file.write(f"v {x} {y} {z}\n")

        # Write faces data to OBJ file
        obj_file.write("g Protein\n")
        obj_file.write("s off\n")
        # backbone edge
        # for i in range(len(back_bone_atom_id_list) - 1):
        #     obj_file.write(f"f {back_bone_atom_id_list[i]} {back_bone_atom_id_list[i+1]}\n")

        for model in structure:
            for chain in model:
                for residue_idx, residue in enumerate(chain):
                    bond_table = residue_bond_table[residue.resname]

                    # compute inter-aa bonds
                    for atom1, binds_to_list in bond_table.items():
                        for atom2 in binds_to_list:
                            print(residue.child_dict[atom2].radius)
                            obj_file.write(
                                f"f {residue.child_dict[atom1].serial_number} {residue.child_dict[atom2].serial_number}\n")

                    if residue_idx == len(chain) - 1:
                        break
        # compute peptide bonds
        residues = tuple(structure.get_residues())
        for idx in range(1, len(residues)):
            prev_aa = residues[idx - 1].child_dict
            current_aa = residues[idx].child_dict
            c_atom = prev_aa.get("C")
            n_atom = current_aa.get("N")
            assert c_atom is not None and n_atom is not None, f"{residues[idx - 1]} at index {idx - 1} does not have a C " \
                                                              f"atom or {residues[idx]} at index {idx} does not have a " \
                                                              f"n atom. Check pdb file"

            obj_file.write(f"f {c_atom.serial_number} {n_atom.serial_number}\n")

    if return_pos:
        return atom_coord
