from pyray import *

MAX_LIGHTS = 4
LIGHTS_COUNT = 0


class LightType:
    DIRECTIONAL = 0
    POINT = 1


class Light:
    def __init__(self, light_type: int, position: list | tuple, target: list | tuple, color: list | tuple,
                 shader):
        self.light_type = light_type
        self.position = position
        self.target = target
        self.color = color
        self.enabled = True

        self.enabledLoc = None
        self.typeLoc = None
        self.positionLoc = None
        self.targetLoc = None
        self.colorLoc = None

        self.create_light(shader)

    def create_light(self, shader):
        global LIGHTS_COUNT
        if LIGHTS_COUNT < MAX_LIGHTS:
            self.enabledLoc = get_shader_location(shader, f"lights{LIGHTS_COUNT}enabled")
            self.typeLoc = get_shader_location(shader, f"lights{LIGHTS_COUNT}type")
            self.positionLoc = get_shader_location(shader, f"lights{LIGHTS_COUNT}position")
            self.targetLoc = get_shader_location(shader, f"lights{LIGHTS_COUNT}target")
            self.colorLoc = get_shader_location(shader, f"lights{LIGHTS_COUNT}color")

            set_shader_value(shader, self.positionLoc, Vector3(*self.position), ShaderUniformDataType.SHADER_UNIFORM_VEC3)
            set_shader_value(shader, self.targetLoc, Vector3(*self.target), ShaderUniformDataType.SHADER_UNIFORM_VEC3)
            set_shader_value(shader, self.colorLoc, Vector4(
                self.color[0] / 255,
                self.color[1] / 255,
                self.color[2] / 255,
                self.color[3]), ShaderUniformDataType.SHADER_UNIFORM_VEC4)
            set_shader_value(shader, self.typeLoc, hex(id(self.light_type)), ShaderUniformDataType.SHADER_UNIFORM_INT)
            set_shader_value(shader, self.enabledLoc, hex(id(self.enabled)), ShaderUniformDataType.SHADER_UNIFORM_INT)
            LIGHTS_COUNT += 1
        else:
            print(f"INFO: SHADER: [rlights] Max lights reached: {MAX_LIGHTS}")
