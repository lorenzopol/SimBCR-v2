import os
import gzip
from Bio.PDB import PDBParser
from BiologicalConstants import *
import math


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


def pdb_to_obj(pdb_file, output_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    back_bone_atoms = ("N", "CA", "C")
    atom_counter = 1

    with open(output_file, 'w') as obj_file:
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        if atom.name in back_bone_atoms:
                            atom.set_serial_number(atom_counter)
                        obj_file.write(f"v {x} {y} {z}\n")
                        atom_counter += 1

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

                    # compute peptide bonds
                    if residue_idx == len(chain) - 1:
                        break
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


def mine_generate_uv_sphere_obj(position, radius, n_slices, n_stack):
    position = [position[0], position[2], -position[1]]

    # add top vertex
    vertex = [(position[0], position[1]+radius, position[2]), ]
    for i in range(n_stack - 1):
        phi = math.pi * (i + 1) / n_stack

        for j in range(n_slices):
            theta = 2 * math.pi * j / n_slices
            x = position[0] + radius * math.sin(phi) * math.cos(theta)
            z = position[2] + radius * math.sin(phi) * math.sin(theta)
            y = position[1] + radius * math.cos(phi)
            vertex.append((x, y, z))

    vertex.append([position[0], position[1]-radius, position[2]])
    quads = []
    for j in range(n_stack-2):
        j0 = j * n_slices + 1
        j1 = (j + 1) * n_slices + 1
        for i in range(n_slices):
            i0 = j0 + i
            i1 = j0 + (i + 1) % n_slices
            i2 = j1 + (i + 1) % n_slices
            i3 = j1 + i
            quads.append((i0, i1, i2, i3))

    with open("obj_files/mine_uvSphere.obj", "w") as file:
        for point in vertex:
            file.write(f"v {point[0]} {point[1]} {point[2]}\n")

        for a, b, c, d in quads:
            file.write(f"f {a+1} {b+1} {c+1}\n")
            file.write(f"f {a+1} {c+1} {d+1}\n")

        # add tris to top vertex
        for i in range(2, n_slices+1):
            file.write(f"f {1} {i} {i+1}\n")
        file.write(f"f {1} {2} {n_slices+1}\n")

        # add tris to last vertex
        last_vertx = n_slices * (n_stack - 1) + 2
        for i in range(1, n_slices):
            file.write(f"f {last_vertx} {last_vertx-i} {last_vertx-i-1}\n")
        file.write(f"f {last_vertx} {last_vertx-1} {last_vertx-n_slices}\n")


if __name__ == "__main__":
    mine_generate_uv_sphere_obj([-1, 2, -3], 1, 10, 10)
    # pdb_to_obj(r'C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb',
    #            r'C:\Users\loren\PycharmProjects\SimBCR-v2\obj_files\MineOutputFaces.obj')
