from Bio.PDB import PDBParser
import BiologicalConstants as bc
from dataclasses import dataclass
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
            f_container.append((starting_idx + i, starting_idx + num_segments + i,
                                starting_idx + num_segments + i + 1, starting_idx + i + 1))
        f_container.append((starting_idx + 1, starting_idx + 2 * num_segments,
                            starting_idx + num_segments + 1, starting_idx + num_segments))
        return v_container, f_container

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
    def middle_point(verts, point_1, point_2, radius, middle_point_cache):
        # We check if we have already cut this edge first # to avoid duplicated verts
        smaller_index = min(point_1, point_2)
        greater_index = max(point_1, point_2)
        key = f"{smaller_index}-{greater_index}"
        if key in middle_point_cache:
            return middle_point_cache[key], middle_point_cache
            # If it's not in cache, then we can cut it
        vert_1 = verts[point_1]
        vert_2 = verts[point_2]
        middle = [sum(i) / 2 for i in zip(vert_1, vert_2)]
        verts.append(GeometryBuilder.vertex(*middle, radius))
        index = len(verts) - 1
        middle_point_cache[key] = index
        return index, middle_point_cache

    @staticmethod
    def calculate_sphere_on_coord(position, radius, subdiv, iteration, v_container, f_container):
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
        start_idx = (40 * subdiv + 2) * iteration + 1
        position = [position[0], position[2], -position[1]]
        middle_point_cache = {}

        pre_v_container, pre_f_container = GeometryBuilder.make_std_icosphere(radius)
        for i in range(subdiv):
            for tri in pre_f_container:
                v1, middle_point_cache = GeometryBuilder.middle_point(pre_v_container, tri[0], tri[1], radius,
                                                                      middle_point_cache)
                v2, middle_point_cache = GeometryBuilder.middle_point(pre_v_container, tri[1], tri[2], radius,
                                                                      middle_point_cache)
                v3, middle_point_cache = GeometryBuilder.middle_point(pre_v_container, tri[2], tri[0], radius,
                                                                      middle_point_cache)
                f_container.append([start_idx + tri[0], start_idx + v1, start_idx + v3])
                f_container.append([start_idx + tri[1], start_idx + v2, start_idx + v1])
                f_container.append([start_idx + tri[2], start_idx + v3, start_idx + v2])
                f_container.append([start_idx + v1, start_idx + v2, start_idx + v3])
        for v in pre_v_container:
            v_container.append([(v[0] + position[0]),
                                (v[1] + position[1]),
                                (v[2] + position[2])])

        return v_container, f_container


class ObjWriters:
    @staticmethod
    def dump_TRIS_containers(obj_file_handle,
                             v_container: list | tuple, f_container: list | tuple,
                             fake_normals: bool, fake_texture: bool):
        """general dumper for obj file.
        Expected input:
            :param obj_file_handle: file handle from with statement
            :param v_container: iterable that stores the XYZ position of the vertices.
                ex -> v_container = [[0.0, 1.0, -2.0], [1.2, -3.4, 1.0], ...]
            :param f_container: iterable that stores the indexes of the vertices in v_container that form a face
                ex -> f_container = [[1, 2, 3], [2, 3, 4]]
            :param n_container: iterable that stores the normal of the face
                ex -> n_container = [[0.00, 0.00, 1.00], [0.50, -0.50, 0.00]]
            :param t_container: iterable that stores the vertex texture in the uv space
                ex -> t_container = [[0.25, 0.00], [0.00, 0.00]]
            :param fake_normals: insert fake normal vector for each face
            :param fake_texture: insert fake UV coords for each vertex
                """

        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        if fake_texture:
            fake_texture_idx = 1
            obj_file_handle.write(f"vt 0.00 0.00\n")

        if fake_normals:
            fake_normal_idx = 1
            obj_file_handle.write(f"vn 0.00 0.00 0.00\n")

        for face in f_container:
            if fake_normals and fake_texture:
                obj_file_handle.write(f"f {face[0]}/{fake_texture_idx}/{fake_normal_idx}"
                                      f" {face[1]}/{fake_texture_idx}/{fake_normal_idx}"
                                      f" {face[2]}/{fake_texture_idx}/{fake_normal_idx}\n")

    @staticmethod
    def dump_v_f_QUADS_containers(obj_file_handle, v_container: list | tuple, f_container: list | tuple,
                                  fake_normals: bool, fake_texture: bool):
        """general dumper for obj file.
        Expected input:
            :param obj_file_handle: file handle from with statement
            :param v_container: iterable that stores the XYZ position of the vertices.
                ex -> v_container = [[0.0, 1.0, -2.0], [1.2, -3.4, 1.0], ...]
            :param f_container: iterable that stores the indexes of the vertices in v_container that form a face
                ex -> f_container = [[1, 2, 3], [2, 3, 4]]
            :param fake_normals: insert fake normal vector for each face
            :param fake_texture: insert fake UV coords for each vertex
                """

        for point in v_container:
            obj_file_handle.write(f"v {point[0]} {point[1]} {point[2]}\n")
        if fake_texture:
            fake_texture_idx = 1
            obj_file_handle.write(f"vt 0.00 0.00\n")

        if fake_normals:
            fake_normal_idx = 1
            obj_file_handle.write(f"vn 0.00 0.00 0.00\n")

        for face in f_container:
            if fake_normals and fake_texture:
                obj_file_handle.write(f"f {face[0]}/{fake_texture_idx}/{fake_normal_idx}"
                                      f" {face[1]}/{fake_texture_idx}/{fake_normal_idx}"
                                      f" {face[2]}/{fake_texture_idx}/{fake_normal_idx}"
                                      f" {face[3]}/{fake_texture_idx}/{fake_normal_idx}\n")

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
    def instance_spheres_on_atom_coords_list(atom_pos_container: tuple | list, radius, subdiv, fake_normals, fake_texture):
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
                atom_pos_container[0], radius, subdiv, 0, [], [])

        else:
            v_container = []
            f_container = []
            for index, atom_pos in enumerate(atom_pos_container):
                v_container, f_container = GeometryBuilder.calculate_sphere_on_coord(
                    atom_pos, radius, subdiv, index, v_container, f_container)
        return v_container, f_container

    def convert_atom_pos_from_coords(self, atom_coords_obj_output_filepath, radius, subdiv,
                                     fake_normals, fake_texture):
        v_container, f_container = self.instance_spheres_on_atom_coords_list(
            self.parser3d.all_atom_coords, radius, subdiv, fake_normals, fake_texture)
        with open(atom_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_TRIS_containers(obj_file_handle, v_container, f_container,
                                            fake_normals, fake_texture)

    # =========================================== BONDS ========================================
    def instance_cylinder_on_bond_coords_list(self, radius, num_segments,
                                              inter_aa_bonds_rel, peptide_bonds_rel):
        v_container = []
        f_container = []
        for index, (atom1_idx, atom2_idx) in enumerate(inter_aa_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container = GeometryBuilder.calculate_cylinder_from_caps_pos(pos1, pos2, radius,
                                                                                        num_segments, index,
                                                                                        v_container, f_container)
            last_idx = index
        for index, (atom1_idx, atom2_idx) in enumerate(peptide_bonds_rel):
            pos1 = self.parser3d.all_atom_coords[atom1_idx - 1]
            pos2 = self.parser3d.all_atom_coords[atom2_idx - 1]
            v_container, f_container = GeometryBuilder.calculate_cylinder_from_caps_pos(pos1, pos2, radius,
                                                                                        num_segments,
                                                                                        index + last_idx, v_container,
                                                                                        f_container)
        return v_container, f_container

    def convert_bond_pos_to_cylinder(self, bond_coords_obj_output_filepath, radius, num_segments,
                                     fake_normals, fake_texture):
        v_container, f_container = self.instance_cylinder_on_bond_coords_list(
            radius, num_segments + 1,
            self.parser3d.inter_aa_bonds_rel, self.parser3d.peptide_bonds_rel)
        with open(bond_coords_obj_output_filepath, "w") as obj_file_handle:
            ObjWriters.dump_v_f_QUADS_containers(obj_file_handle, v_container, f_container,
                                                 fake_normals, fake_texture)


if __name__ == "__main__":
    parser = PdbParser3D(r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\tiny.pdb")
    conv = PdbToObjConverter(parser)
    conv.convert_atom_pos_from_coords("obj_files/atom_coords.obj", radius=.5, subdiv=1,
                                      fake_normals=True, fake_texture=True)
    conv.convert_bond_pos_to_cylinder("obj_files/bond_coords.obj", radius=0.25, num_segments=4,
                                      fake_normals=True, fake_texture=True)

    v_cont, f_cont = GeometryBuilder.calculate_sphere_on_coord([0, 0, 0], 1, 1, 0, [], [])
    with open("obj_files/icosphere.obj", "w") as file:
        ObjWriters.dump_TRIS_containers(file, v_cont, f_cont, fake_normals=True, fake_texture=True)
