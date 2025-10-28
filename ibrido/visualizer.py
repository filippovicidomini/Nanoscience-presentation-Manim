# show_orbital_obj.py
from manim import *
import numpy as np
import trimesh
from manim import config
#config.renderer = "opengl"   # evita problemi Cairo/gradienti

MESH_FILE = "/Users/filippovicidomini/Library/CloudStorage/GoogleDrive-filippo.vicidomini@alumni.unive.it/My Drive/University /Nanotechnology & Nanomaterials/Nanoscience-presentation-Manim/tetrahedral_touch_calibrated.obj"  # cambia se hai generato altri

class ShowOrbitalOBJ(ThreeDScene):
    def construct(self):
        self.camera.background_color = "#0b0e10"
        self.set_camera_orientation(phi=60*DEGREES, theta=30*DEGREES)

        tm = trimesh.load(MESH_FILE)
        V = np.asarray(tm.vertices, dtype=float)
        F = np.asarray(tm.faces, dtype=int)

        # centra e scala la mesh
        V -= V.mean(axis=0)                   # centra sull’origine
        max_r = np.max(np.linalg.norm(V, axis=1))
        scale = 1                           # <-- fattore di riduzione (0.5 = metà)
        if max_r > 1e-9:
            V *= (scale * 3.5 / max_r)        # porta il raggio ~3.5 e poi riduci con 'scale'

        # Ulteriore riduzione (opzionale)
        max_faces = 3500
        if len(F) > max_faces:
            keep = np.random.choice(len(F), max_faces, replace=False)
            F = F[keep]

        tris = VGroup()
        for a, b, c in F:
            pa, pb, pc = V[a], V[b], V[c]
            poly = Polygon(
                [pa[0], pa[1], pa[2]],
                [pb[0], pb[1], pb[2]],
                [pc[0], pc[1], pc[2]],
                stroke_width=0
            ).set_fill(BLUE_C, opacity=0.6)
            tris.add(poly)
            
        def color_from_density(p):
            r = np.linalg.norm(p)
            # più vicino al centro → colore più brillante
            return interpolate_color(BLUE_E, WHITE, np.exp(-r))

        tris = VGroup()
        for a, b, c in F:
            pa, pb, pc = V[a], V[b], V[c]
            color = color_from_density((pa+pb+pc)/3)
            poly = Polygon(
                [pa[0], pa[1], pa[2]],
                [pb[0], pb[1], pb[2]],
                [pc[0], pc[1], pc[2]],
                stroke_width=0
            ).set_fill(color, opacity=0.6)
            tris.add(poly)

        self.play(FadeIn(tris))
        # rotation in x e y
        self.play(Rotate(tris, angle=2*PI, axis=[1,1,0], run_time=5))