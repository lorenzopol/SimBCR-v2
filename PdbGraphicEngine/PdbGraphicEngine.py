from pyray import *

from BiologicalConstants import relative_atom_radius
from rlights import Light, LightType
from PdbToObj import PdbParser3D


class Globals:
    ATOM_RADIUS_FACTOR = 1.5
    BOND_RADIUS_FACTOR = 1.5

    @staticmethod
    def get_default_camera() -> Camera3D:
        return Camera3D(
            Vector3(5.6, 6.3, 63.7),
            Vector3(0.0, 0.0, 0.0),
            Vector3(0.0, 1.0, 0.0),
            45.0, CameraProjection.CAMERA_PERSPECTIVE)


class GraphicEngine:
    def __init__(self, pdb_parser: PdbParser3D):
        self.pdb_parser = pdb_parser
        self.SCREEN_WIDTH = 1920
        self.SCREEN_HEIGHT = 1080
        self.velocity = 0.2
        self.camera = Globals.get_default_camera()
        self.init_engine()

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
            transforms_container.append(
                matrix_multiply(
                    matrix_multiply(scaleMat, rotMat),
                    transMat)
            )
        return transforms_container

    def run(self):
        # check https://www.reddit.com/r/raylib/comments/v70krp/how_does_rlglh_batching_work/
        # batch = rl_load_render_batch(1, 8192)

        shader = load_shader("shaders/default.vert", "shaders/default.frag")
        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MVP] = get_shader_location(shader, "mvp")
        shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW] = get_shader_location(shader, "viewPos")
        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MODEL] = get_shader_location_attrib(shader,
                                                                                              "instanceTransform")

        set_shader_value(shader, get_shader_location(shader, "ambient"), Vector4(2, 2, 2, 1.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC4)
        set_shader_value(shader, get_shader_location(shader, "ambient_intensity"), Vector3(.1, .1, .1),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "diffuse_intensity"), Vector3(.8, .8, .8),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "specular_intensity"), Vector3(1.0, 1.0, 1.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)
        set_shader_value(shader, get_shader_location(shader, "light_pos"), Vector3(0.0, 0.0, 24.0),
                         ShaderUniformDataType.SHADER_UNIFORM_VEC3)

        # todo: add lights proper implementation
        atoms_material: Material = load_material_default()
        atoms_material.shader = shader
        atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = RED
        bonds_material: Material = load_material_default()
        bonds_material.shader = shader
        bonds_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = BLUE

        atoms_transforms = []
        gen_sphere = gen_mesh_sphere(0.3, 12, 12)
        max_x, max_y = 0, 0
        for atom in self.pdb_parser.all_atoms:
            x, y, z = self.pdb_parser.all_atom_coords[atom.serial_number-1]
            max_x = max(max_x, abs(x))
            max_y = max(max_y, abs(y))
            atom_relative_radius = relative_atom_radius[atom.element] * Globals.ATOM_RADIUS_FACTOR
            atoms_transforms.append(matrix_multiply(
                matrix_multiply(
                    matrix_scale(atom_relative_radius, atom_relative_radius, atom_relative_radius),
                    matrix_rotate(Vector3(1.0, 1.0, 1.0), 0)
                ), matrix_translate(x, y, z))
            )
        inter_aa_bonds_transforms = self.compute_bonds_transforms(self.pdb_parser.inter_aa_bonds_rel)
        peptide_bonds_transforms = self.compute_bonds_transforms(self.pdb_parser.peptide_bonds_rel)
        gen_cylinder = gen_mesh_cylinder(0.125, 1, 12)

        is_camera_orbit_control = False
        show_grid = True
        # Main game loop
        while not window_should_close():
            time = Vector2(get_time(), 0)
            set_shader_value(shader, get_shader_location(shader, "u_time"), time, ShaderUniformDataType.SHADER_UNIFORM_VEC2)

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

            # reset camera
            if is_key_pressed(KeyboardKey.KEY_R):
                del self.camera
                self.camera = Globals.get_default_camera()

            set_shader_value(shader, shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW], self.camera.position,
                             ShaderUniformDataType.SHADER_UNIFORM_VEC3)
            # Draw
            begin_drawing()
            clear_background(RAYWHITE)
            begin_mode_3d(self.camera)
            draw_mesh_instanced(gen_sphere, atoms_material, atoms_transforms, len(atoms_transforms))
            draw_mesh_instanced(gen_cylinder, bonds_material, inter_aa_bonds_transforms, len(inter_aa_bonds_transforms))
            draw_mesh_instanced(gen_cylinder, bonds_material, peptide_bonds_transforms, len(peptide_bonds_transforms))

            if show_grid:
                draw_grid(int(max(max_x, max_y))*2, 1.0)
            end_mode_3d()

            # Draw GUI
            draw_rectangle(10, 40, 240, 200, fade(SKYBLUE, 0.5))
            draw_rectangle_lines(10, 40, 240, 200, BLUE)

            draw_text("Free camera default controls:", 20, 45, 10, BLACK)
            draw_text("- WASD to move", 40, 60, 10, DARKGRAY)
            draw_text("- SPACE/L-CTRL to move up/down ", 40, 80, 10, DARKGRAY)
            draw_text("- Z to look at center", 40, 100, 10, DARKGRAY)
            draw_text("- move the mouse to free look", 40, 120, 10, DARKGRAY)
            draw_text("- R to reset camera", 40, 140, 10, DARKGRAY)
            draw_text("- H to toggle the grid", 40, 160, 10, DARKGRAY)
            draw_text("- O to toggle the orbit view", 40, 180, 10, DARKGRAY)
            draw_text("- Mouse Wheel to Zoom in-out", 40, 200, 10, DARKGRAY)
            draw_text("- Mouse Wheel Pressed to Pan", 40, 220, 10, DARKGRAY)
            draw_text(
                f"Camera position: {self.camera.position.x} | {self.camera.position.y} | {self.camera.position.z}", 40,
                240, 10, RED)
            draw_fps(10, 10)

            end_drawing()

        # De-Initialization
        close_window()


if __name__ == "__main__":
    pdb_path = r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb"
    _3d_parser = PdbParser3D(pdb_path, None)

    engine = GraphicEngine(_3d_parser)
    engine.run()
