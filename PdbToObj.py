import os
from Bio.PDB import PDBParser
import BiologicalConstants as bc
from dataclasses import dataclass
import math


@dataclass
class PdbToObjConverter:
    pdb_file_path: str

    def __post_init__(self):
        # consider making parser static if not used elsewhere
        self.parser = PDBParser()
        self.structure = self.parser.get_structure('protein', self.pdb_file_path)
        self.all_atom_coords = self.get_all_atom_coords()
        self.back_bone_atom_coords = self.get_only_back_bone_atom_coords()

    def get_all_atom_coords(self):
        all_atom_coord = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x, y, z = atom.coord
                        all_atom_coord.append((x, y, z))
        return all_atom_coord
    
    def get_only_back_bone_atom_coords(self):
        back_bone_atoms = ("N", "CA", "C")
        back_bone_atom_coord = []

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.name in back_bone_atoms:
                            x, y, z = atom.coord
                            back_bone_atom_coord.append((x, y, z))
        return back_bone_atom_coord

    def get_only_residue_atom_coords(self):
        back_bone_atoms = ("N", "CA", "C")
        residue_atom_coord = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.name not in back_bone_atoms:
                            x, y, z = atom.coord
                            residue_atom_coord.append((x, y, z))
        return residue_atom_coord

    def compute_inter_aa_bonds_relationship(self):
        """uses atom.serial_number(ie second col in .pdb files) to reference bonds relationship between atoms inside
        an AA"""
        inter_aa_bonds_relationship_container = []
        for model in self.structure:
            for chain in model:
                for residue_idx, residue in enumerate(chain):
                    bond_table = bc.residue_bond_table[residue.resname]

                    # compute inter-aa bonds
                    for atom1, binds_to_list in bond_table.items():
                        for atom2 in binds_to_list:
                            inter_aa_bonds_relationship_container.append(
                                (residue.child_dict[atom1].serial_number, residue.child_dict[atom2].serial_number))

                    # do we really need to optimise this lol
                    if residue_idx == len(chain) - 1:
                        break
        return inter_aa_bonds_relationship_container

    def compute_peptide_bonds_relationship(self):
        peptide_bonds_relationship_container = []
        residues = tuple(self.structure.get_residues())
        for idx in range(1, len(residues)):
            prev_aa = residues[idx - 1].child_dict
            current_aa = residues[idx].child_dict
            c_atom = prev_aa.get("C")
            n_atom = current_aa.get("N")
            assert c_atom is not None and n_atom is not None, f"{residues[idx - 1]} at index {idx - 1} does not have a C " \
                                                              f"atom or {residues[idx]} at index {idx} does not have a " \
                                                              f"n atom. Check pdb file"

            peptide_bonds_relationship_container.append((c_atom.serial_number, n_atom.serial_number))
        return peptide_bonds_relationship_container

    @staticmethod
    def calculate_sphere_on_atom_coord(atom_coords: list | tuple, radius: int | float,
                                       n_slices: int, n_stack: int,
                                       iteration: int = 0,
                                       v_container: list | None = None, f_container: list | None = None):
        """calculate vertex and faces based on atom_pos.
        > single use
        calculate_sphere_on_atom_coord(atom_coords=[x0, y0, z0], radius=0.2,
                                       n_slices=7, n_stack=7)

        > iterative use
        calculate_sphere_on_atom_coord(atom_coords=[x_i, y_i, z_i], radius=0.2,
                                       n_slices=7, n_stack=7,
                                       iteration=i,
                                       v_container=[], f_container=[])
        returns v_container and f_container with the correctly populated data for vertex and faces value for obj dump"""
        if v_container is None or f_container is None:
            v_container = []
            f_container = []
        last_vertx = n_slices * (n_stack - 1) + 2
        position = [atom_coords[0], atom_coords[2], -atom_coords[1]]
        starting_idx = last_vertx * iteration

        # add top vertex
        vertex = [(position[0], position[1] + radius, position[2]), ]
        for i in range(n_stack - 1):
            phi = math.pi * (i + 1) / n_stack

            for j in range(n_slices):
                theta = 2 * math.pi * j / n_slices
                x = position[0] + radius * math.sin(phi) * math.cos(theta)
                z = position[2] + radius * math.sin(phi) * math.sin(theta)
                y = position[1] + radius * math.cos(phi)
                vertex.append((x, y, z))

        vertex.append([position[0], position[1] - radius, position[2]])
        quads = []
        for j in range(n_stack - 2):
            j0 = j * n_slices + 1
            j1 = (j + 1) * n_slices + 1
            for i in range(n_slices):
                i0 = starting_idx + (j0 + i)
                i1 = starting_idx + (j0 + (i + 1) % n_slices)
                i2 = starting_idx + (j1 + (i + 1) % n_slices)
                i3 = starting_idx + (j1 + i)
                quads.append((i0, i1, i2, i3))
        for point in vertex:
            v_container.append((point[0], point[1], point[2]))

        for a, b, c, d in quads:
            f_container.append((a + 1, b + 1, c + 1))
            f_container.append((a + 1, c + 1, d + 1))

        for i in range(starting_idx + 2, starting_idx + n_slices + 1):
            f_container.append((starting_idx + 1, i, i + 1))
        f_container.append((starting_idx + 1, starting_idx + 2, starting_idx + n_slices + 1))

        for i in range(1, n_slices):
            f_container.append((starting_idx + last_vertx,
                                starting_idx + last_vertx - i,
                                starting_idx + last_vertx - i - 1))
        f_container.append((starting_idx + last_vertx, starting_idx + last_vertx, starting_idx + last_vertx - n_slices))

        return v_container, f_container

    @staticmethod
    def dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container):
        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        for face in f_container:
            obj_file_handle.write(f"f {face[0]} {face[1]} {face[2]}\n")

    def convert(self, obj_output_filepath):
        v_container = []
        f_container = []
        for index, atom_pos in enumerate(self.all_atom_coords):
            v_container, f_container = self.calculate_sphere_on_atom_coord(
                atom_pos, 0.2, 9, 9, index, v_container, f_container)
        with open(obj_output_filepath, "w") as obj_file_handle:
            self.dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container)


if __name__ == "__main__":
    conv = PdbToObjConverter(r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb")
    conv.convert("obj_files/first_try.obj")
