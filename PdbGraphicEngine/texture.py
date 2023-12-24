import pygame as pg
import moderngl as mgl
import os


class Texture:
    def __init__(self, ctx):
        self.ctx = ctx
        self.textures = {
            "": self.get_texture(path=f"{os.path.dirname(__file__)}/textures/default.png")
        }

    def get_texture(self, path):
        texture = pg.image.load(path).convert()
        texture = pg.transform.flip(texture, flip_x=False, flip_y=True)
        texture = self.ctx.texture(size=texture.get_size(), components=3,
                                   data=pg.image.tostring(texture, "RGB"))
        texture.filter = (mgl.LINEAR_MIPMAP_LINEAR, mgl.LINEAR)
        texture.build_mipmaps()
        texture.anisotropy = 32.0
        return texture

    def destroy(self):
        [tex.release() for tex in self.textures.values()]
