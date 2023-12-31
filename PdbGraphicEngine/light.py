import glm


class Light:
    def __init__(self, position=(4, 12, 13), color=(1, 1, 1)):
        self.position = glm.vec3(position)
        self.color = glm.vec3(color)
        # intensities
        self.Ia = 0.1 * self.color  # ambient color
        self.Id = 0.8 * self.color  # diffuse color
        self.Is = 1.0 * self.color  # specular color
