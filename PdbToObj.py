from Bio.PDB import PDBParser
import BiologicalConstants as bc
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class PdbToObjConverter:
    pdb_file_path: str

    def __post_init__(self):
        # consider making parser static if not used elsewhere
        self.parser = PDBParser()
        self.structure = self.parser.get_structure('protein', self.pdb_file_path)
        self.all_atom_coords = self.get_all_atom_coords()
        self.all_atoms = self.get_all_atoms()
        self.back_bone_atom_coords = self.get_only_back_bone_atom_coords()

    def get_all_atoms(self):
        all_atom = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        all_atom.append(atom)
        return all_atom

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
                                       iteration=0, v_container=None, f_container=None):
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

    def calculate_sphere_on_bulk_atom_coords_list(self, atom_pos_container: tuple | list, radius,
                                                  n_slices, n_stack):
        """calculate vertex and faces based on atom_pos.
                > single use
                calculate_sphere_on_atom_coord(atom_coords=[[x0, y0, z0]], radius=0.2,
                                               n_slices=7, n_stack=7)

                > iterative use
                calculate_sphere_on_atom_coord(atom_coords=[[x_i, y_i, z_i], ...], radius=0.2,
                                               n_slices=7, n_stack=7,
                                               iteration=i,
                                               v_container=[], f_container=[])
                returns v_container and f_container with the correctly populated data for vertex and faces value for obj dump"""

        # there is only 1 atom to draw
        if len(atom_pos_container) == 1:
            v_container, f_container = self.calculate_sphere_on_atom_coord(
                atom_pos_container, radius, n_slices, n_stack, 0)
        else:
            v_container = []
            f_container = []
            for index, atom_pos in enumerate(atom_pos_container):
                v_container, f_container = self.calculate_sphere_on_atom_coord(
                    atom_pos, radius, n_slices, n_stack, index, v_container, f_container)
        return v_container, f_container

    def calculate_sphere_on_bulk_atom_list(self):
        """uses raw input from .pdb files and not manually parsed coords"""
        v_container = []
        f_container = []
        for index, atom in enumerate(self.all_atoms):
            v_container, f_container = self.calculate_sphere_on_atom_coord(
                atom.coord, bc.relative_atom_radius.get(atom.element) / 4, 7, 7, index, v_container, f_container)
        return v_container, f_container

    @staticmethod
    def dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container):
        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        for face in f_container:
            obj_file_handle.write(f"f {face[0]} {face[1]} {face[2]}\n")

    def dump_bond_containers_to_obj_file(self, obj_file_handle, inter_aa_bonds_rel, peptide_bonds_rel):
        for point in self.all_atom_coords:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")

        for serial_number_atom_a, serial_number_atom_b in inter_aa_bonds_rel:
            obj_file_handle.write(f"f {serial_number_atom_a} {serial_number_atom_b}\n")
        for serial_number_atom_a, serial_number_atom_b in peptide_bonds_rel:
            obj_file_handle.write(f"f {serial_number_atom_a} {serial_number_atom_b}\n")

    def populate_v_f_bond_containers_to_cylinder_file(self, inter_aa_bonds_rel, peptide_bonds_rel):
        v_container = []
        f_container = []
        for index, (atom1_idx, atom2_idx) in enumerate(inter_aa_bonds_rel):
            pos1 = self.all_atom_coords[atom1_idx - 1]
            pos2 = self.all_atom_coords[atom2_idx - 1]
            v_container, f_container = my_create_cylinder_obj(pos1, pos2, 0.1, 4, index, v_container, f_container)
            last_idx = index
        for index, (atom1_idx, atom2_idx) in enumerate(peptide_bonds_rel):
            pos1 = self.all_atom_coords[atom1_idx - 1]
            pos2 = self.all_atom_coords[atom2_idx - 1]
            v_container, f_container = my_create_cylinder_obj(pos1, pos2, 0.1, 4, index + last_idx, v_container,
                                                              f_container)
        return v_container, f_container

    @staticmethod
    def dump_bond_containers_to_cylinder_file(obj_file_handle, v_container, f_container):
        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        for face in f_container:
            obj_file_handle.write(f"f {face[0]} {face[1]} {face[2]}\n")

    # todo: should separate point clouds dumpers from sphere/cyl dumper
    def convert_atom_pos_from_coords(self, atom_coords_obj_output_filepath):
        v_container, f_container = self.calculate_sphere_on_bulk_atom_coords_list(
            self.all_atom_coords, 0.2, 11, 5)
        with open(atom_coords_obj_output_filepath, "w") as obj_file_handle:
            self.dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container)

    def convert_atom_pos(self, atom_coords_obj_output_filepath):
        v_container, f_container = self.calculate_sphere_on_bulk_atom_list()
        with open(atom_coords_obj_output_filepath, "w") as obj_file_handle:
            self.dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container)

    def convert_bond_pos(self, bond_coords_obj_output_filepath):
        inter_aa_bonds_rel = self.compute_inter_aa_bonds_relationship()
        peptide_bonds_rel = self.compute_peptide_bonds_relationship()
        with open(bond_coords_obj_output_filepath, "w") as obj_file_handle:
            self.dump_bond_containers_to_obj_file(obj_file_handle, inter_aa_bonds_rel, peptide_bonds_rel)

    def convert_bond_pos_to_cylinder(self, bond_coords_obj_output_filepath):
        inter_aa_bonds_rel = self.compute_inter_aa_bonds_relationship()
        peptide_bonds_rel = self.compute_peptide_bonds_relationship()
        v_container, f_container = self.populate_v_f_bond_containers_to_cylinder_file(inter_aa_bonds_rel,
                                                                                      peptide_bonds_rel)
        with open(bond_coords_obj_output_filepath, "w") as obj_file_handle:
            self.dump_bond_containers_to_cylinder_file(obj_file_handle, v_container, f_container)


def my_create_cylinder_obj(point1, point2, radius, num_segments, iteration, v_container, f_container):
    point1 = np.array(point1)
    point2 = np.array(point2)
    starting_idx = num_segments * 2 * iteration
    theta = np.linspace(0, 2 * np.pi, num_segments)

    m = len(radius) if isinstance(radius, list) else 1  # Check number of radius values

    if m == 1:
        radius = [radius, radius]
        m = 2

    X = np.zeros((m, num_segments))
    Y = np.zeros((m, num_segments))
    Z = np.zeros((m, num_segments))

    v = (point2 - point1) / np.linalg.norm(point2 - point1)
    R2 = np.random.rand(3)
    x2 = v - R2 / np.dot(R2, v)
    x2 /= np.linalg.norm(x2)
    x3 = np.cross(v, x2)
    x3 /= np.linalg.norm(x3)

    r1x, r1y, r1z = point1
    r2x, r2y, r2z = point2
    x2x, x2y, x2z = x2
    x3x, x3y, x3z = x3

    time = np.linspace(0, 1, m)

    for j in range(m):
        t = time[j]
        X[j, :] = r1x + (r2x - r1x) * t + radius[j] * np.cos(theta) * x2x + radius[j] * np.sin(theta) * x3x
        Y[j, :] = r1y + (r2y - r1y) * t + radius[j] * np.cos(theta) * x2y + radius[j] * np.sin(theta) * x3y
        Z[j, :] = r1z + (r2z - r1z) * t + radius[j] * np.cos(theta) * x2z + radius[j] * np.sin(theta) * x3z
    vertex = []
    for ext in range(len(X)):
        for idx in range(len(X[ext])):
            v_container.append((X[ext][idx], Y[ext][idx], Z[ext][idx]))
    for i in range(1, num_segments):
        f_container.append((starting_idx + i, starting_idx + i + 1, starting_idx + num_segments + i))
        f_container.append((starting_idx + i + 1, starting_idx + num_segments + i, starting_idx + num_segments + i + 1))
    f_container.append((starting_idx + 1, starting_idx + num_segments, starting_idx + 2 * num_segments))
    f_container.append((starting_idx + 1, starting_idx + num_segments + 1, starting_idx + 2 * num_segments))
    return v_container, f_container


if __name__ == "__main__":
    # g_point1 = np.array([-0.5, 3, 0])  # Replace with your first center point
    # g_point2 = np.array([0, 1, 2])  # Replace with your second center point
#
    # g_v_container, g_f_container = my_create_cylinder_obj(g_point1, g_point2, 1, 32, 0, [], [])
    # with open(r"C:\Users\loren\PycharmProjects\SimBCR-v2\obj_files\cylinder.obj", "w") as file:
    #     for g_point in g_v_container:
    #         file.write(f"v {g_point[0]} {g_point[1]} {g_point[2]}\n")
    #     for g_face in g_f_container:
    #         file.write(f"f {g_face[0]} {g_face[1]} {g_face[2]}\n")

    conv = PdbToObjConverter(r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb")
    conv.convert_atom_pos("obj_files/atom_coords.obj")
    conv.convert_bond_pos_to_cylinder("obj_files/bond_coords.obj")
