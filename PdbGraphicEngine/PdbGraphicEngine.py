from pyray import *
from PIL import Image

from BiologicalConstants import relative_atom_radius
from PdbToObj import PdbParser3D


class Globals:
    ATOM_RADIUS_FACTOR = 1.5
    BOND_RADIUS_FACTOR = 1.5
    Y_GLOBAL_ROT_MATRIX = matrix_rotate(Vector3(0.0, 1.0, 0.0), 90)
    X_GLOBAL_ROT_MATRIX = matrix_rotate(Vector3(0.0, 0.0, 1.0), 90)

    @staticmethod
    def get_default_camera() -> Camera3D:
        return Camera3D(
            Vector3(5.6, 6.3, 63.7),
            Vector3(0.0, 0.0, 0.0),
            Vector3(0.0, 1.0, 0.0),
            45.0, CameraProjection.CAMERA_PERSPECTIVE)


def in_place_array_mat_mul(array, matrix):
    for i in range(len(array)):
        array[i] = matrix_multiply(array[i], matrix)
    return array


class GraphicEngine:
    def __init__(self, pdb_parser: PdbParser3D):
        self.pdb_parser = pdb_parser
        self.SCREEN_WIDTH = 1920
        self.SCREEN_HEIGHT = 1080
        self.velocity = 0.2
        self.is_first_draw_call = True
        self.camera = Globals.get_default_camera()
        self.init_engine()

        # draw calls attributes
        self.gen_sphere = gen_mesh_sphere(0.3, 12, 12)
        self.gen_cylinder = gen_mesh_cylinder(0.125, 1, 12)

    def init_engine(self):
        init_window(self.SCREEN_WIDTH, self.SCREEN_HEIGHT, "raylib [core] example - 3d camera mode")
        disable_cursor()
        set_target_fps(60)

    def compute_bonds_transforms(self, bonds_container):
        transforms_container = []
        up = Vector3(0, 1, 0)
        for bond in bonds_container:
            prev_atom_pos, next_atom_pos = Vector3(*self.pdb_parser.all_atom_coords[bond[0] - 1]), \
                Vector3(*self.pdb_parser.all_atom_coords[bond[1] - 1])
            bond_axis = vector3_subtract(next_atom_pos, prev_atom_pos)
            bond_length = vector3_length(bond_axis)
            perp = vector3_cross_product(up, bond_axis)
            angle = vector3_angle(up, bond_axis)
            transMat = matrix_translate(prev_atom_pos.x, prev_atom_pos.y, prev_atom_pos.z)
            rotMat = matrix_rotate(perp, angle)
            scaleMat = matrix_scale(Globals.BOND_RADIUS_FACTOR, bond_length, Globals.BOND_RADIUS_FACTOR)
            std_transform = matrix_multiply(
                matrix_multiply(scaleMat, rotMat),
                transMat)
            transforms_container.append(
                matrix_multiply(
                    matrix_multiply(std_transform, Globals.X_GLOBAL_ROT_MATRIX),
                    Globals.Y_GLOBAL_ROT_MATRIX)
            )
        return transforms_container


    def run(self):
        shader = load_shader("shaders/default.vert", "shaders/default.frag")

        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MVP] = get_shader_location(shader, "mvp")
        shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW] = get_shader_location(shader, "viewPos")
        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MODEL] = get_shader_location_attrib(shader,
                                                                                              "instanceTransform")

        # todo: add lights proper implementation
        set_shader_value(shader, get_shader_location(shader, "ambient"), Vector4(2, 2, 2, 1.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC4)
        set_shader_value(shader, get_shader_location(shader, "ambient_intensity"), Vector3(.1, .1, .1),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "diffuse_intensity"), Vector3(.8, .8, .8),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "specular_intensity"), Vector3(1.0, 1.0, 1.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "light_pos"), Vector3(-10.0, 0.0, 24.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)

        # materials
        base_atoms_material: Material = load_material_default()
        base_atoms_material.shader = shader
        base_atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = GRAY

        cdr1_atoms_material: Material = load_material_default()
        cdr1_atoms_material.shader = shader
        cdr1_atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = RED

        cdr2_atoms_material: Material = load_material_default()
        cdr2_atoms_material.shader = shader
        cdr2_atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = GREEN

        cdr3_atoms_material: Material = load_material_default()
        cdr3_atoms_material.shader = shader
        cdr3_atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = BLUE

        bonds_material: Material = load_material_default()
        bonds_material.shader = shader
        bonds_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = BLACK

        # compute atoms transforms
        base_atoms_transforms = []
        cdr1_atoms_transforms = []
        cdr2_atoms_transforms = []
        cdr3_atoms_transforms = []
        max_x, max_y = 0, 0
        for atom in self.pdb_parser.all_atoms:
            x, y, z = self.pdb_parser.all_atom_coords[atom.serial_number - 1]
            max_x = max(max_x, abs(x))
            max_y = max(max_y, abs(y))
            atom_relative_radius = relative_atom_radius[atom.element] * Globals.ATOM_RADIUS_FACTOR
            std_transform = matrix_multiply(matrix_multiply(
                matrix_scale(atom_relative_radius, atom_relative_radius, atom_relative_radius),
                matrix_rotate(Vector3(1.0, 1.0, 1.0), 0)
            ), matrix_translate(x, y, z))
            transform = matrix_multiply(
                matrix_multiply(std_transform, Globals.X_GLOBAL_ROT_MATRIX),
                Globals.Y_GLOBAL_ROT_MATRIX)

            if atom.serial_number in self.pdb_parser.cdr1_atoms:
                cdr1_atoms_transforms.append(transform)
            elif atom.serial_number in self.pdb_parser.cdr2_atoms:
                cdr2_atoms_transforms.append(transform)
            elif atom.serial_number in self.pdb_parser.cdr3_atoms:
                cdr3_atoms_transforms.append(transform)
            else:
                base_atoms_transforms.append(transform)

        # compute bonds transforms
        base_bonds_transform = self.compute_bonds_transforms(self.pdb_parser.base_bonds)
        cdr1_bonds_transform = self.compute_bonds_transforms(self.pdb_parser.cdr1_bonds)
        cdr2_bonds_transform = self.compute_bonds_transforms(self.pdb_parser.cdr2_bonds)
        cdr3_bonds_transform = self.compute_bonds_transforms(self.pdb_parser.cdr3_bonds)

        # mainloop
        is_camera_orbit_control = False
        show_grid = True
        is_rendering = False
        spin_counter = 0
        spin = 0
        cdr_show = 5
        while not window_should_close():
            # camera mode controller
            if is_key_pressed(KeyboardKey.KEY_O):
                is_camera_orbit_control = not is_camera_orbit_control
            if is_key_pressed(KeyboardKey.KEY_H):
                show_grid = not show_grid
            if is_camera_orbit_control:
                self.camera.target = Vector3(0.0, 0.0, 0.0)
                update_camera(self.camera, CameraMode.CAMERA_ORBITAL)
            else:
                update_camera(self.camera, CameraMode.CAMERA_FREE)
            if is_key_pressed(KeyboardKey.KEY_Z):
                is_rendering = not is_rendering
            if is_key_pressed(KeyboardKey.KEY_X):
                spin_counter += 1
                spin = (spin_counter % 3) - 1
            # reset camera
            if is_key_pressed(KeyboardKey.KEY_R):
                self.camera.target = Vector3(0.0, 0.0, 0.0)
                self.camera.up = Vector3(0.0, 1.0, 0.0)

            if is_key_down(KeyboardKey.KEY_RIGHT_BRACKET):  # equals minus keyboard key
                self.camera.fovy -= 1
            if is_key_down(KeyboardKey.KEY_SLASH):  # equals plus keyboard key
                self.camera.fovy += 1

            if is_key_pressed(KeyboardKey.KEY_ONE):
                cdr_show = 1
            if is_key_pressed(KeyboardKey.KEY_TWO):
                cdr_show = 2
            if is_key_pressed(KeyboardKey.KEY_THREE):
                cdr_show = 3
            if is_key_pressed(KeyboardKey.KEY_FOUR):
                cdr_show = 4
            if is_key_pressed(KeyboardKey.KEY_FIVE):
                cdr_show = 5

            if spin != 0:
                spin_matrix = matrix_rotate(Vector3(0.0, 1.0, 0.0), get_frame_time() * spin)

                base_atoms_transforms = in_place_array_mat_mul(base_atoms_transforms, spin_matrix)
                cdr1_atoms_transforms = in_place_array_mat_mul(cdr1_atoms_transforms, spin_matrix)
                cdr2_atoms_transforms = in_place_array_mat_mul(cdr2_atoms_transforms, spin_matrix)
                cdr3_atoms_transforms = in_place_array_mat_mul(cdr3_atoms_transforms, spin_matrix)

                base_bonds_transform = in_place_array_mat_mul(base_bonds_transform, spin_matrix)
                cdr1_bonds_transform = in_place_array_mat_mul(cdr1_bonds_transform, spin_matrix)
                cdr2_bonds_transform = in_place_array_mat_mul(cdr2_bonds_transform, spin_matrix)
                cdr3_bonds_transform = in_place_array_mat_mul(cdr3_bonds_transform, spin_matrix)

            set_shader_value(shader, shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW], self.camera.position,
                             ShaderUniformDataType.SHADER_UNIFORM_VEC3)
            # Draw
            begin_drawing()
            clear_background(RAYWHITE)
            begin_mode_3d(self.camera)

            if cdr_show in (1, 4, 5):
                draw_mesh_instanced(self.gen_sphere, cdr1_atoms_material, cdr1_atoms_transforms,
                                    len(cdr1_atoms_transforms))
                draw_mesh_instanced(self.gen_cylinder, bonds_material, cdr1_bonds_transform, len(cdr1_bonds_transform))
            if cdr_show in (2, 4, 5):
                draw_mesh_instanced(self.gen_sphere, cdr2_atoms_material, cdr2_atoms_transforms,
                                    len(cdr2_atoms_transforms))
                draw_mesh_instanced(self.gen_cylinder, bonds_material, cdr2_bonds_transform, len(cdr2_bonds_transform))
            if cdr_show in (3, 4, 5):
                draw_mesh_instanced(self.gen_sphere, cdr3_atoms_material, cdr3_atoms_transforms,
                                    len(cdr3_atoms_transforms))
                draw_mesh_instanced(self.gen_cylinder, bonds_material, cdr3_bonds_transform, len(cdr3_bonds_transform))
            if cdr_show == 5:
                draw_mesh_instanced(self.gen_sphere, base_atoms_material, base_atoms_transforms,
                                    len(base_atoms_transforms))
                draw_mesh_instanced(self.gen_cylinder, bonds_material, base_bonds_transform, len(base_bonds_transform))

            if show_grid:
                draw_grid(int(max(max_x, max_y)) * 2, 1.0)

            if is_rendering:
                export_image(load_image_from_screen(), "export.png")
                is_rendering = not is_rendering
            end_mode_3d()

            # Draw GUI
            draw_rectangle(10, 40, 240, 200, fade(SKYBLUE, 0.5))
            draw_rectangle_lines(10, 40, 240, 200, BLUE)

            draw_text("Free camera default controls:", 20, 45, 10, BLACK)
            draw_text("- WASD to move", 40, 60, 10, DARKGRAY)
            draw_text("- SPACE/L-CTRL to move up/down ", 40, 80, 10, DARKGRAY)
            draw_text("- Z to render", 40, 100, 10, DARKGRAY)
            draw_text("- X to make the receptor spin", 40, 120, 10, DARKGRAY)
            draw_text("- R to reset camera", 40, 140, 10, DARKGRAY)
            draw_text("- H to toggle the grid", 40, 160, 10, DARKGRAY)
            draw_text("- O to toggle the orbit view", 40, 180, 10, DARKGRAY)
            draw_text("- Mouse Wheel or +/- to Zoom in-out", 40, 200, 10, DARKGRAY)
            draw_text("- Mouse Wheel Pressed to Pan", 40, 220, 10, DARKGRAY)

            draw_rectangle(10, 240, 100, 100, fade(LIGHTGRAY, 0.5))
            draw_rectangle_lines(10, 240, 100, 100, LIGHTGRAY)
            draw_text("CDR INFO", 20, 260, 10, BLACK)
            draw_text(f"CDR1: {self.pdb_parser.cdr1_range[0]}>{self.pdb_parser.cdr1_range[1]}", 20, 280, 10, DARKGRAY)
            draw_text(f"CDR2: {self.pdb_parser.cdr2_range[0]}>{self.pdb_parser.cdr2_range[1]}", 20, 300, 10, DARKGRAY)
            draw_text(f"CDR3: {self.pdb_parser.cdr3_range[0]}>{self.pdb_parser.cdr3_range[1]}", 20, 320, 10, DARKGRAY)

            draw_fps(10, 10)

            end_drawing()

        # De-Initialization
        close_window()


if __name__ == "__main__":
    pdb_path = r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb"
    _3d_parser = PdbParser3D(pdb_path, "94-112")

    engine = GraphicEngine(_3d_parser)
    engine.run()
