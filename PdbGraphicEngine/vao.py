from .vbo import VBO
from .shader_program import ShaderProgram


class VAO:
    def __init__(self, ctx):
        self.ctx = ctx
        self.vbo = VBO(ctx)
        self.program = ShaderProgram(ctx)
        self.vaos = {
            "atoms": self.get_vao(
                program=self.program.programs["default"],
                vbo=self.vbo.vbos["atoms"]),
            "bonds": self.get_vao(
                program=self.program.programs["default"],
                vbo=self.vbo.vbos["bonds"]),
        }

    def get_vao(self, program, vbo):
        vao = self.ctx.vertex_array(program, [(vbo.vbo, vbo.format, *vbo.attrib)])
        return vao

    def destroy(self):
        self.vbo.destroy()
        self.program.destroy()
