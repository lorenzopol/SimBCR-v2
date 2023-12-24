from .model import *


class Scene:
    def __init__(self, app):
        self.app = app
        self.objects = []
        self.load()

    def add_object(self, obj):
        self.objects.append(obj)

    def load(self):
        app = self.app
        add = self.add_object
        # default addition
        # add(Cube(app))

        # "custom" add
        # add(Cube(app, tex_id=1, pos=(-2.5, 0, 0), rot=(45, 0, 0), scale=(1, 2, 1)))

        add(Obj(app, vao_name="atoms", pos=(-14.712911152061071, -2, -20)))
        add(Obj(app, vao_name="bonds", pos=(-14.712911152061071, -2, -20), rot=(-90, 0, 0)))

    def render(self):
        for obj in self.objects:
            obj.render()
