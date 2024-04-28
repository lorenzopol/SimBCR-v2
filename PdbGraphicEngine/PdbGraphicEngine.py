from pyray import *

from rlights import Light, LightType
from PdbToObj import PdbParser3D


class GraphicEngine:
    def __init__(self, pdb_parser: PdbParser3D):
        self.pdb_parser = pdb_parser
        self.SCREEN_WIDTH = 1600
        self.SCREEN_HEIGHT = 900
        self.velocity = 0.2
        self.camera = Camera3D(
            Vector3(10.0, 10.0, 10.0),
            Vector3(0.0, 0.0, 0.0),
            Vector3(0.0, 1.0, 0.0),
            45.0,
            CameraProjection.CAMERA_PERSPECTIVE
        )
        self.init_engine()
        print(f"atoms_coords: {_3d_parser.all_atom_coords}")
        print(f"inter_aa_bonds_rel: {_3d_parser.inter_aa_bonds_rel}")
        print(f"peptide_bonds_rel: {_3d_parser.peptide_bonds_rel}")

    def init_engine(self):
        init_window(self.SCREEN_WIDTH, self.SCREEN_HEIGHT, "raylib [core] example - 3d camera mode")
        disable_cursor()
        set_target_fps(60)

    def run(self):
        shader = load_shader("shaders/lighting_instancing.vs", "shaders/lighting.fs")
        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MVP] = get_shader_location(shader, "mvp")
        shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW] = get_shader_location(shader, "viewPos")
        shader.locs[ShaderLocationIndex.SHADER_LOC_MATRIX_MODEL] = get_shader_location_attrib(shader, "instanceTransform")

        set_shader_value(shader, get_shader_location(shader, "ambient"), Vector4(2, 2, 2, 1.0), ShaderUniformDataType.SHADER_UNIFORM_VEC4)
        Light(LightType.POINT, (5.0, 10.0, 10.0), (0.0, 0.0, 0.0), (255.0, 255.0, 255.0, 1.0), shader)

        atoms_material: Material = load_material_default()
        atoms_material.shader = shader
        atoms_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = RED

        bonds_material: Material = load_material_default()
        bonds_material.shader = shader
        bonds_material.maps[MaterialMapIndex.MATERIAL_MAP_ALBEDO].color = BLUE

        sphere_transforms = []
        gen_sphere = gen_mesh_sphere(0.5, 12, 12)
        for atom_coord in self.pdb_parser.all_atom_coords:
            sphere_transforms.append(
                matrix_multiply(
                    matrix_rotate(Vector3(1.0, 1.0, 1.0), 0),
                    matrix_translate(*atom_coord)
                )
            )
        cylinder_transforms = []
        gen_cylinder = gen_mesh_cylinder(0.25, 1, 12)
        for bond in self.pdb_parser.inter_aa_bonds_rel:
            prev_atom_pos, next_atom_pos = self.pdb_parser.all_atom_coords[bond[0]-1], self.pdb_parser.all_atom_coords[bond[1]-1]
            bond_axis = vector3_normalize(vector3_subtract(next_atom_pos, prev_atom_pos))
            cylinder_transforms.append(
                matrix_multiply(
                    matrix_rotate(bond_axis, 0),
                    matrix_translate(*prev_atom_pos)
                )
            )

        is_camera_orbit_control = False
        show_grid = True
        # Main game loop
        while not window_should_close():
            if is_key_pressed(KeyboardKey.KEY_Z):
                self.camera.target = Vector3(0.0, 0.0, 0.0)
            if is_key_down(KeyboardKey.KEY_X):
                self.camera.position.y -= (1.0 * self.velocity)
            if is_key_down(KeyboardKey.KEY_C):
                self.camera.position.y += (1.0 * self.velocity)

            # camera mode controller
            if is_key_pressed(KeyboardKey.KEY_O):
                is_camera_orbit_control = not is_camera_orbit_control
            if is_key_pressed(KeyboardKey.KEY_H):
                show_grid = not show_grid
            if is_camera_orbit_control:
                update_camera(self.camera, CameraMode.CAMERA_ORBITAL)
            else:
                update_camera(self.camera, CameraMode.CAMERA_FIRST_PERSON)

            # reset camera
            if is_key_pressed(KeyboardKey.KEY_R):
                del self.camera
                self.camera = Camera3D(
                    Vector3(10.0, 10.0, 10.0),
                    Vector3(0.0, 0.0, 0.0),
                    Vector3(0.0, 1.0, 0.0),
                    45.0,
                    CameraProjection.CAMERA_PERSPECTIVE
                )

            set_shader_value(shader, shader.locs[ShaderLocationIndex.SHADER_LOC_VECTOR_VIEW], self.camera.position, ShaderUniformDataType.SHADER_UNIFORM_VEC3)
            # Draw
            begin_drawing()
            clear_background(RAYWHITE)
            begin_mode_3d(self.camera)
            draw_mesh_instanced(gen_sphere, atoms_material, sphere_transforms, len(sphere_transforms))
            draw_mesh_instanced(gen_cylinder, bonds_material, cylinder_transforms, len(cylinder_transforms))

            if show_grid:
                draw_grid(100, 1.0)
            end_mode_3d()
            draw_rectangle(10, 40, 240, 200, fade(SKYBLUE, 0.5))
            draw_rectangle_lines(10, 40, 240, 200, BLUE)

            draw_text("Free camera default controls:", 20, 45, 10, BLACK)
            draw_text("- WASD to move", 40, 60, 10, DARKGRAY)
            draw_text("- C/X to move up/down ", 40, 80, 10, DARKGRAY)
            draw_text("- Z to look at center", 40, 100, 10, DARKGRAY)
            draw_text("- move the mouse to free look", 40, 120, 10, DARKGRAY)
            draw_text("- R to reset camera", 40, 140, 10, DARKGRAY)
            draw_text("- H to toggle the grid", 40, 160, 10, DARKGRAY)
            draw_text("- O to toggle the orbit view", 40, 180, 10, DARKGRAY)
            draw_text("- Mouse Wheel to Zoom in-out", 40, 200, 10, DARKGRAY)
            draw_text("- Mouse Wheel Pressed to Pan", 40, 220, 10, DARKGRAY)

            draw_fps(10, 10)

            end_drawing()

        # De-Initialization
        close_window()


if __name__ == "__main__":
    pdb_path = r"C:\Users\loren\PycharmProjects\SimBCR-v2\pdb_files\first_try.pdb"
    _3d_parser = PdbParser3D(pdb_path, None)

    engine = GraphicEngine(_3d_parser)
    engine.run()
