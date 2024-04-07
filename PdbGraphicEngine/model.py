import glm


class BaseModel:
    def __init__(self, app, vao_name, tex_id, pos=(0, 0, 0), rot=(0, 0, 0), scale=(1, 1, 1)):
        self.app = app
        self.pos = pos
        self.rot = glm.vec3([glm.radians(a) for a in rot])
        self.scale = scale
        self.m_model = self.get_model_matrix()
        self.tex_id = tex_id
        self.vao = app.mesh.vao.vaos[vao_name]
        self.program = self.vao.program
        self.camera = self.app.camera

    def update(self, rotate):
        ...

    def get_model_matrix(self):
        m_model = glm.mat4()

        # translation
        m_model = glm.translate(m_model, self.pos)

        # rotations
        m_model = glm.rotate(m_model, self.rot.x, glm.vec3(1, 0, 0))
        m_model = glm.rotate(m_model, self.rot.y, glm.vec3(0, 1, 0))
        m_model = glm.rotate(m_model, self.rot.z, glm.vec3(0, 0, 1))

        # scaling
        m_model = glm.scale(m_model, self.scale)
        return m_model

    def render(self, rotate):
        self.update(rotate)
        self.vao.render()


class Cube(BaseModel):
    def __init__(self, app, vao_name="cube", tex_id=0, pos=(0, 0, 0), rot=(0, 0, 0), scale=(1, 1, 1)):
        super().__init__(app, vao_name, tex_id, pos, rot, scale)
        self.texture = self.app.mesh.texture.textures[self.tex_id]
        self.on_init()

    def update(self, rotate):
        self.texture.use()
        self.program["camPos"].write(self.app.camera.position)
        self.program["m_view"].write(self.app.camera.m_view)
        self.program["m_model"].write(self.m_model)

    def on_init(self):
        # texture
        self.program["u_texture_0"] = 0
        self.texture.use()

        # mvp mats
        self.program["m_model"].write(self.m_model)
        self.program["m_view"].write(self.app.camera.m_view)
        self.program["m_proj"].write(self.app.camera.m_proj)

        # light
        self.program["light.position"].write(self.app.light.position)
        self.program["light.Ia"].write(self.app.light.Ia)
        self.program["light.Id"].write(self.app.light.Id)
        self.program["light.Is"].write(self.app.light.Is)


class Obj(BaseModel):
    def __init__(self, app, vao_name, tex_id="",
                 pos=(0, 0, 0), rot=(0, 0, 0), scale=(1, 1, 1)):
        super().__init__(app, vao_name, tex_id, pos, rot, scale)
        self.texture = self.app.mesh.texture.textures[self.tex_id]
        self.on_init()

    def update(self, rotate):
        self.texture.use()
        if rotate:
            up = (0, 1, 0) if round(self.rot[0], 2) == -3.14 else (0, 0, -1)
            self.m_model = glm.rotate(self.m_model, self.app.time/1200, glm.vec3(up))
        self.program["camPos"].write(self.app.camera.position)
        self.program["m_view"].write(self.app.camera.m_view)
        self.program["m_model"].write(self.m_model)

        self.program["m_model"].write(self.m_model)

    def on_init(self):
        # texture
        self.program["u_texture_0"] = 0
        self.texture.use()

        # mvp mats
        self.program["m_model"].write(self.m_model)
        self.program["m_view"].write(self.app.camera.m_view)
        self.program["m_proj"].write(self.app.camera.m_proj)

        # light
        self.program["light.position"].write(self.app.light.position)
        self.program["light.Ia"].write(self.app.light.Ia)
        self.program["light.Id"].write(self.app.light.Id)
        self.program["light.Is"].write(self.app.light.Is)
