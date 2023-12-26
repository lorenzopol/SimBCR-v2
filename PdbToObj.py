from Bio.PDB import PDBParser
import BiologicalConstants as bc
from dataclasses import dataclass
import numpy as np


class GeometryBuilder:
    @staticmethod
    def calculate_face_normal(v0, v1, v2):
        dir0 = np.array(v1) - np.array(v0)
        dir1 = np.array(v2) - np.array(v0)
        face_normal = np.cross(dir0, dir1)
        magnitude = np.linalg.norm(face_normal)
        return face_normal / magnitude if magnitude else np.zeros_like(face_normal)

    @staticmethod
    def is_vector_facing_point(vector, point):
        vector_to_point = point - vector
        return np.dot(vector, vector_to_point) > 0

    @staticmethod
    def calculate_cylinder_from_caps_pos(cap_pos1: list | tuple, cap_pos2: list | tuple,
                                         radius: float | int | list,
                                         num_segments: int, iteration: int,
                                         v_container: list | tuple, f_container: list | tuple,
                                         n_container: list | tuple, t_container: list | tuple,
                                         fake_normals: bool, fake_texture: bool):
        """
        taken from https://it.mathworks.com/matlabcentral/fileexchange/5468-cylinder-between-2-points
        cap_pos1 and cap_pos2 indicates the center of the circles that enclose the cylinder
        if radius is a list, each value will be the radius of a slice"""
        cap_pos1 = np.array(cap_pos1)
        cap_pos2 = np.array(cap_pos2)

        starting_idx = num_segments * 2 * iteration
        theta = np.linspace(0, 2 * np.pi, num_segments)
        neo_f_container = []

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
            neo_f_container.append((starting_idx + i + 1, starting_idx + num_segments + i + 1,
                                    starting_idx + num_segments + i, starting_idx + i))
        neo_f_container.append((starting_idx + num_segments, starting_idx + 2 * num_segments,
                                starting_idx + num_segments + 1, starting_idx + 1))
        if fake_normals:
            n_container.extend([[0.00, 0.00, 1.00] for _ in neo_f_container])
        else:
            # we can take three arbitrary vertices because each cylinder face is coplanar by construction
            for face in neo_f_container:
                v0, v1, v2 = v_container[face[0] - 1], v_container[face[1] - 1], v_container[face[2] - 1]
                face_normal = GeometryBuilder.calculate_face_normal(v0, v1, v2)
                n_container.append(face_normal)

        if fake_texture:
            t_container.append([0.00, 0.00])
        f_container.extend(neo_f_container)
        return v_container, f_container, n_container, t_container

    @staticmethod
    def vertex(x, y, z, radius):
        length = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        return [(i / length) * radius for idx, i in enumerate([x, y, z])]

    @staticmethod
    def make_std_icosphere(radius):
        PHI = (1 + np.sqrt(5)) / 2
        verts = [
            GeometryBuilder.vertex(-1, PHI, 0, radius),
            GeometryBuilder.vertex(1, PHI, 0, radius),
            GeometryBuilder.vertex(-1, -PHI, 0, radius),
            GeometryBuilder.vertex(1, -PHI, 0, radius),
            GeometryBuilder.vertex(0, -1, PHI, radius),
            GeometryBuilder.vertex(0, 1, PHI, radius),
            GeometryBuilder.vertex(0, -1, -PHI, radius),
            GeometryBuilder.vertex(0, 1, -PHI, radius),
            GeometryBuilder.vertex(PHI, 0, -1, radius),
            GeometryBuilder.vertex(PHI, 0, 1, radius),
            GeometryBuilder.vertex(-PHI, 0, -1, radius),
            GeometryBuilder.vertex(-PHI, 0, 1, radius),
        ]

        faces = [
            # 5 faces around point 0
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            # Adjacent faces
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            # 5 faces around 3
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            # Adjacent faces
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ]
        return verts, faces

    @staticmethod
    def middle_point(point_1, point_2, radius, verts, middle_point_cache):
        smaller_index = min(point_1, point_2)
        greater_index = max(point_1, point_2)
        key = f"{smaller_index}-{greater_index}"

        if key in middle_point_cache:
            return middle_point_cache[key], verts
            # If it's not in cache, then we can cut it
        vert_1 = verts[point_1]
        vert_2 = verts[point_2]
        middle = [sum(i) / 2 for i in zip(vert_1, vert_2)]
        verts.append(GeometryBuilder.vertex(*middle, radius))
        index = len(verts) - 1
        middle_point_cache[key] = index
        return index, verts

    @staticmethod
    def calculate_sphere_on_coord(position: list[int | float, ...] | tuple[int | float, ...],
                                  radius: int | float,
                                  subdiv: int, iteration: int,
                                  fake_normals: bool, fake_texture: bool):
        # from https://sinestesia.co/blog/tutorials/python-icospheres/
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
        start_idx = (10 * (4 ** subdiv) + 2) * iteration
        middle_point_cache = {}
        position = [position[0], position[2], -position[1]]

        v_container, f_container = GeometryBuilder.make_std_icosphere(radius)
        n_container = []
        t_container = []
        for i in range(subdiv):

            # populate f_container. !!! pre_v_container gets modified in middle_point
            temp_f_container = []
            for tri in f_container:
                v1, v_container = GeometryBuilder.middle_point(tri[0], tri[1], radius,
                                                               v_container, middle_point_cache)
                v2, v_container = GeometryBuilder.middle_point(tri[1], tri[2], radius,
                                                               v_container, middle_point_cache)
                v3, v_container = GeometryBuilder.middle_point(tri[2], tri[0], radius,
                                                               v_container, middle_point_cache)

                temp_f_container.append([tri[0], v1, v3])
                temp_f_container.append([tri[1], v2, v1])
                temp_f_container.append([tri[2], v3, v2])
                temp_f_container.append([v1, v2, v3])
            f_container = temp_f_container
        # populate n_container
        if fake_normals:
            n_container.extend([[0.00, 0.00, 1.00] for _ in f_container])
        else:
            for face_idx, face in enumerate(f_container):
                face_idx_0 = face[0]
                face_idx_1 = face[1]
                face_idx_2 = face[2]
                v0, v1, v2 = v_container[face_idx_0], v_container[face_idx_1], v_container[face_idx_2]
                face_normal = GeometryBuilder.calculate_face_normal(v0, v1, v2)
                n_container.append(face_normal.tolist())

        # populate t_container
        if fake_texture:
            t_container = [[0.00, 0.00, 1.00] for _ in f_container]
        for idx_face in range(len(f_container)):
            f_container[idx_face][0] = f_container[idx_face][0] + 1 + start_idx
            f_container[idx_face][1] = f_container[idx_face][1] + 1 + start_idx
            f_container[idx_face][2] = f_container[idx_face][2] + 1 + start_idx

        for idx_vert in range(len(v_container)):
            v_container[idx_vert][0] = v_container[idx_vert][0] + position[0]
            v_container[idx_vert][1] = v_container[idx_vert][1] + position[1]
            v_container[idx_vert][2] = v_container[idx_vert][2] + position[2]
        return v_container, f_container, n_container, t_container


class ObjWriters:
    @staticmethod
    def dump_TRIS_containers(obj_file_handle,
                             v_container: list | tuple, f_container: list | tuple,
                             n_container: list | tuple, t_container: list | tuple):
        """general TRIS dumper for obj file.
            :param obj_file_handle: file handle from with statement
            :param v_container: iterable that stores the XYZ position of the vertices.
                ex -> v_container = [[0.0, 1.0, -2.0], [1.2, -3.4, 1.0], ...]
            :param f_container: iterable that stores the indexes of the vertices in v_container that form a face
                ex -> f_container = [[1, 2, 3], [2, 3, 4]]
            :param n_container: iterable that stores the normal of the face
                ex -> n_container = [[0.00, 0.00, 1.00], [0.50, -0.50, 0.00]]
            :param t_container: iterable that stores the vertex texture in the uv space
                ex -> t_container = [[0.25, 0.00], [0.00, 0.00]]
                """

        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        if True:
            fake_texture_idx = 1
            obj_file_handle.write(f"vt 0.00 0.00\n")

        for normal in n_container:
            obj_file_handle.write(f"vn {normal[0]} {normal[1]} {normal[2]}\n")

        for face_idx, face in enumerate(f_container):
            obj_file_handle.write(f"f {face[0]}/{fake_texture_idx}/{face_idx + 1}"
                                  f" {face[1]}/{fake_texture_idx}/{face_idx + 1}"
                                  f" {face[2]}/{fake_texture_idx}/{face_idx + 1}\n")

    @staticmethod
    def dump_QUADS_containers(obj_file_handle, v_container: list | tuple, f_container: list | tuple,
                              n_container, t_container):
        """general QUADS dumper for obj file.
            :param obj_file_handle: file handle from with statement
            :param v_container: iterable that stores the XYZ position of the vertices.
                ex -> v_container = [[0.0, 1.0, -2.0], [1.2, -3.4, 1.0], ...]
            :param f_container: iterable that stores the indexes of the vertices in v_container that form a face
                ex -> f_container = [[1, 2, 3], [2, 3, 4]]
            :param n_container: iterable that stores the normal of the face
                ex -> n_container = [[0.00, 0.00, 1.00], [0.50, -0.50, 0.00]]
            :param t_container: iterable that stores the vertex texture in the uv space
                ex -> t_container = [[0.25, 0.00], [0.00, 0.00]]
                """

        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        if True:
            fake_texture_idx = 1
            obj_file_handle.write(f"vt 0.00 0.00\n")

        for normal in n_container:
            obj_file_handle.write(f"vn {normal[0]} {normal[1]} {normal[2]}\n")

        for face_idx, face in enumerate(f_container):
            obj_file_handle.write(f"f {face[0]}/{fake_texture_idx}/{face_idx + 1}"
                                  f" {face[1]}/{fake_texture_idx}/{face_idx + 1}"
                                  f" {face[2]}/{fake_texture_idx}/{face_idx + 1}"
                                  f" {face[3]}/{fake_texture_idx}/{face_idx + 1}\n")

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
    def instance_spheres_on_atom_coords_list(atom_pos_container: tuple | list, radius, subdiv,
                                             fake_normals, fake_texture):
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
            v_container, f_container, n_container, t_container = GeometryBuilder.calculate_sphere_on_coord(
                atom_pos_container, radius, subdiv, 0, fake_normals, fake_texture)

        else:
            v_container = []
            f_container = []
            n_container = []
            t_container = []
            for index, atom_pos in enumerate(atom_pos_container):
                nth_v_container, nth_f_container, nth_n_container, nth_t_container = \
                    GeometryBuilder.calculate_sphere_on_coord(
                        atom_pos, radius, subdiv, index, fake_normals, fake_texture)

                v_container.extend(nth_v_container)
                f_container.extend(nth_f_container)
                n_container.extend(nth_n_container)
                t_container.extend(nth_t_container)
        return v_container, f_container, n_container, t_container

    def convert_atom_pos_from_coords(self, atom_coords_obj_output_filepath, radius, subdiv,
                                     fake_normals, fake_texture):
        v_container, f_container, n_container, t_container = self.instance_spheres_on_atom_coords_list(
            self.parser3d.all_atom_coords, radius, subdiv, fake_normals, fake_texture)
        with open(atom_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_TRIS_containers(obj_file_handle, v_container, f_container,
                                            n_container, t_container)

    # =========================================== BONDS ========================================
    def instance_cylinder_on_bond_coords_list(self, radius, num_segments,
                                              inter_aa_bonds_rel, peptide_bonds_rel,
                                              fake_normals, fake_texture):
        v_container = []
        f_container = []
        n_container = []
        t_container = []
        for index, (atom1_idx, atom2_idx) in enumerate(inter_aa_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container, n_container, t_container = GeometryBuilder.calculate_cylinder_from_caps_pos(
                pos1, pos2, radius, num_segments, index,
                v_container, f_container, n_container, t_container, fake_normals, fake_texture)
            last_idx = index
        for index, (atom1_idx, atom2_idx) in enumerate(peptide_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container, n_container, t_container = GeometryBuilder.calculate_cylinder_from_caps_pos(
                pos1, pos2, radius, num_segments, index + last_idx,
                v_container, f_container, n_container, t_container, fake_normals, fake_texture)
        return v_container, f_container, n_container, t_container

    def convert_bond_pos_to_cylinder(self, bond_coords_obj_output_filepath, radius, num_segments,
                                     fake_normals, fake_texture):
        v_container, f_container, n_container, t_container = self.instance_cylinder_on_bond_coords_list(
            radius, num_segments + 1,
            self.parser3d.inter_aa_bonds_rel, self.parser3d.peptide_bonds_rel, fake_normals, fake_texture)
        with open(bond_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_QUADS_containers(obj_file_handle, v_container,
                                             f_container, n_container, t_container)


if __name__ == "__main__":
    parser = PdbParser3D(r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb")
    conv = PdbToObjConverter(parser)
    conv.convert_atom_pos_from_coords("obj_files/atom_coords.obj", radius=.5, subdiv=1,
                                      fake_normals=False, fake_texture=True)
    conv.convert_bond_pos_to_cylinder("obj_files/bond_coords.obj", radius=0.25, num_segments=6,
                                      fake_normals=False, fake_texture=True)

    # from random import randrange

    # v_cont = []
    # f_cont = []
    # n_cont = []
    # t_cont = []
    # for idx in range(10):
    #     nth_v_cont, nth_f_cont, nth_n_cont, nth_t_cont = GeometryBuilder.calculate_sphere_on_coord(
    #         [randrange(-10, 10), randrange(-10, 10), randrange(-10, 10)],
    #         radius=1, subdiv=1, iteration=idx,
    #         fake_normals=False, fake_texture=True)
    #     v_cont.extend(nth_v_cont)
    #     f_cont.extend(nth_f_cont)
    #     n_cont.extend(nth_n_cont)
    #     t_cont.extend(nth_t_cont)
    # with open("obj_files/icosphere.obj", "w") as file:
    #     ObjWriters.dump_TRIS_containers(file, v_cont, f_cont, n_cont, t_cont)
