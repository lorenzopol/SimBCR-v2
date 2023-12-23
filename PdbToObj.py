from Bio.PDB import PDBParser
import BiologicalConstants as bc
from dataclasses import dataclass
import math
import numpy as np


class GeometryBuilder:
    @staticmethod
    def calculate_cylinder_from_caps_pos(cap_pos1: list | tuple, cap_pos2: list | tuple,
                                         radius: float | int | list,
                                         num_segments: int, iteration: int,
                                         v_container, f_container):
        """
        taken from https://it.mathworks.com/matlabcentral/fileexchange/5468-cylinder-between-2-points
        cap_pos1 and cap_pos2 indicates the center of the circles that enclose the cylinder
        if radius is a list, each value will be the radius of a slice"""
        cap_pos1 = np.array(cap_pos1)
        cap_pos2 = np.array(cap_pos2)
        starting_idx = num_segments * 2 * iteration
        theta = np.linspace(0, 2 * np.pi, num_segments)

        m = len(radius) if isinstance(radius, list) else 1  # Check number of radius values

        if m == 1:
            radius = [radius, radius]
            m = 2

        X = np.zeros((m, num_segments))
        Y = np.zeros((m, num_segments))
        Z = np.zeros((m, num_segments))

        direction = (cap_pos2 - cap_pos1) / np.linalg.norm(cap_pos2 - cap_pos1)
        random_3v = np.random.rand(3)

        x2 = direction - random_3v / np.dot(random_3v, direction)
        x2 /= np.linalg.norm(x2)
        x3 = np.cross(direction, x2)
        x3 /= np.linalg.norm(x3)

        p1x, p1y, p1z = cap_pos1
        p2x, p2y, p2z = cap_pos2
        x2x, x2y, x2z = x2
        x3x, x3y, x3z = x3

        time = np.linspace(0, 1, m)

        for j in range(m):
            t = time[j]
            X[j, :] = p1x + (p2x - p1x) * t + radius[j] * np.cos(theta) * x2x + radius[j] * np.sin(theta) * x3x
            Y[j, :] = p1y + (p2y - p1y) * t + radius[j] * np.cos(theta) * x2y + radius[j] * np.sin(theta) * x3y
            Z[j, :] = p1z + (p2z - p1z) * t + radius[j] * np.cos(theta) * x2z + radius[j] * np.sin(theta) * x3z

        for ext in range(len(X)):
            for idx in range(len(X[ext])):
                v_container.append((X[ext][idx], Y[ext][idx], Z[ext][idx]))
        for i in range(1, num_segments):
            f_container.append((starting_idx + i, starting_idx + i + 1, starting_idx + num_segments + i))
            f_container.append(
                (starting_idx + i + 1, starting_idx + num_segments + i, starting_idx + num_segments + i + 1))
        f_container.append((starting_idx + 1, starting_idx + num_segments, starting_idx + 2 * num_segments))
        f_container.append((starting_idx + 1, starting_idx + num_segments + 1, starting_idx + 2 * num_segments))
        return v_container, f_container

    @staticmethod
    def calculate_sphere_on_coord(center_coord: list | tuple, radius: int | float,
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
        position = [center_coord[0], center_coord[2], -center_coord[1]]
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


class ObjWriters:
    @staticmethod
    def dump_v_f_container_to_obj_file(obj_file_handle, v_container: list | tuple, f_container):
        """general dumper for obj file.
        Expected input:
            :arg obj_file_handle, file handle from with statement
            :arg v_container, iterable that stores the XYZ position of the vertices.
                ex -> v_container = [[0.0, 1.0, -2.0], [1.2, -3.4, 1.0], ...]
            :arg f_container, iterable that stores the indexes of the vertices in v_container that form a face
                ex -> f_container = [[1, 2, 3], [2, 3, 4]]"""
        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        for face in f_container:
            obj_file_handle.write(f"f {face[0]} {face[1]} {face[2]}\n")

    @staticmethod
    def dump_atom_pos_point_cloud(obj_file_handle, parser3d):
        """dump point cloud representation of the atom's position of a 3dparser"""
        for point in parser3d.all_atom_coords:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")

    @staticmethod
    def dump_bond_containers_to_obj_file(obj_file_handle, parser3d):
        """dump the stick representation of the bonds in a 3dparser"""
        for serial_number_atom_a, serial_number_atom_b in parser3d.inter_aa_bonds_rel:
            obj_file_handle.write(f"f {serial_number_atom_a} {serial_number_atom_b}\n")
        for serial_number_atom_a, serial_number_atom_b in parser3d.peptide_bonds_rel:
            obj_file_handle.write(f"f {serial_number_atom_a} {serial_number_atom_b}\n")


@dataclass
class PdbParser3D:
    pdb_file_path: str

    def __post_init__(self):
        # consider making parser static if not used elsewhere
        self.parser = PDBParser()
        self.structure = self.parser.get_structure('protein', self.pdb_file_path)
        self.all_atom_coords = self.get_all_atom_coords()
        self.all_atoms = self.get_all_atoms()
        self.back_bone_atom_coords = self.get_only_back_bone_atom_coords()
        self.inter_aa_bonds_rel = self.compute_inter_aa_bonds_relationship()
        self.peptide_bonds_rel = self.compute_peptide_bonds_relationship()

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


class PdbToObjConverter:
    def __init__(self, parser3d: PdbParser3D):
        self.parser3d = parser3d

    @staticmethod
    def instance_spheres_on_atom_coords_list(atom_pos_container: tuple | list, radius,
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
            v_container, f_container = GeometryBuilder.calculate_sphere_on_coord(
                atom_pos_container, radius, n_slices, n_stack, 0)
        else:
            v_container = []
            f_container = []
            for index, atom_pos in enumerate(atom_pos_container):
                v_container, f_container = GeometryBuilder.calculate_sphere_on_coord(
                    atom_pos, radius, n_slices, n_stack, index, v_container, f_container)
        return v_container, f_container

    def convert_atom_pos_from_coords(self, atom_coords_obj_output_filepath):
        v_container, f_container = self.instance_spheres_on_atom_coords_list(
            self.parser3d.all_atom_coords, 0.2, 11, 5)
        with open(atom_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container)

    # =========================================== BONDS ========================================
    def instance_cylinder_on_bond_coords_list(self, inter_aa_bonds_rel, peptide_bonds_rel):
        v_container = []
        f_container = []
        for index, (atom1_idx, atom2_idx) in enumerate(inter_aa_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container = GeometryBuilder.calculate_cylinder_from_caps_pos(pos1, pos2, 0.1, 4, index,
                                                                                        v_container, f_container)
            last_idx = index
        for index, (atom1_idx, atom2_idx) in enumerate(peptide_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container = GeometryBuilder.calculate_cylinder_from_caps_pos(pos1, pos2, 0.1, 4,
                                                                                        index + last_idx, v_container,
                                                                                        f_container)
        return v_container, f_container

    def convert_bond_pos_to_cylinder(self, bond_coords_obj_output_filepath):
        v_container, f_container = self.instance_cylinder_on_bond_coords_list(self.parser3d.inter_aa_bonds_rel,
                                                                              self.parser3d.peptide_bonds_rel)
        with open(bond_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_v_f_container_to_obj_file(obj_file_handle, v_container, f_container)


if __name__ == "__main__":
    parser = PdbParser3D(r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb")
    conv = PdbToObjConverter(parser)
    conv.convert_atom_pos_from_coords("obj_files/atom_coords.obj")
    conv.convert_bond_pos_to_cylinder("obj_files/bond_coords.obj")
