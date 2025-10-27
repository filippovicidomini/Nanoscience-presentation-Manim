from manim import *
import numpy as np
config.renderer = "opengl" 

POS_COLOR = BLUE_C     # fase +
NEG_COLOR = RED_C      # fase -
NUCLEUS_COLOR = GRAY_B

def sphere_surface(center=ORIGIN, R=1.0):
    S = Surface(
        lambda u, v: center + R * np.array([
            np.sin(u) * np.cos(v),
            np.sin(u) * np.sin(v),
            np.cos(u)
        ]),
        u_range=[0, PI],
        v_range=[0, TAU],
        resolution=(24, 48),
    )
    S.set_fill(GRAY_C, opacity=0.18).set_stroke(width=0)
    return S

def p_orbital(axis="z", R=1.3, opacity=0.35):
    # funzione base r(u,v)
    base = lambda u, v: np.array([
        np.sin(u)*np.cos(v),
        np.sin(u)*np.sin(v),
        np.cos(u)
    ])
    if axis == "x":
        amp = lambda u, v: abs(np.sin(u)*np.cos(v))
        sign = lambda u, v: np.sign(np.sin(u)*np.cos(v))
    elif axis == "y":
        amp = lambda u, v: abs(np.sin(u)*np.sin(v))
        sign = lambda u, v: np.sign(np.sin(u)*np.sin(v))
    else:  # z
        amp = lambda u, v: abs(np.cos(u))
        sign = lambda u, v: np.sign(np.cos(u))

    # lobo “positivo” (mostra tutto ma col colore POS)
    pos = Surface(
        lambda u, v: (R * amp(u, v)) * base(u, v),
        u_range=[0, PI],
        v_range=[0, TAU],
        resolution=(24, 48),
    )
    pos.set_fill(POS_COLOR, opacity=opacity).set_stroke(width=0)

    # lobo “negativo”: mascheriamo la metà con segno < 0
    neg = Surface(
        lambda u, v: (R * amp(u, v)) * base(u, v) if sign(u, v) < 0 else np.array([np.nan, np.nan, np.nan]),
        u_range=[0, PI],
        v_range=[0, TAU],
        resolution=(24, 48),
    )
    neg.set_fill(NEG_COLOR, opacity=opacity).set_stroke(width=0)

    return VGroup(pos, neg)

def sp3_lobe(direction=np.array([1,1,1]), R=1.6, opacity=0.38):
    d = direction / (np.linalg.norm(direction) + 1e-9)
    def base(u, v):
        return np.array([np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u)])
    def amp(u, v):
        vec = base(u, v)
        return max(0.0, np.dot(vec, d))
    L = Surface(
        lambda u, v: (R * amp(u, v)) * base(u, v),
        u_range=[0, PI],
        v_range=[0, TAU],
        resolution=(24, 48),
    )
    L.set_fill(POS_COLOR, opacity=opacity).set_stroke(width=0)
    return L

class SiliconOrbitalsAndSP3_Safe(ThreeDScene):
    def construct(self):
        self.camera.background_color = "#0b0e10"
        self.set_camera_orientation(phi=65*DEGREES, theta=30*DEGREES)

        # Nucleo
        nucleus = Sphere(radius=0.18)
        nucleus.set_fill(NUCLEUS_COLOR, opacity=0.35).set_stroke(width=0)

        # 3s + 3p
        s3 = sphere_surface(R=1.0)
        px = p_orbital("x", R=1.35, opacity=0.35)
        py = p_orbital("y", R=1.35, opacity=0.35)
        pz = p_orbital("z", R=1.35, opacity=0.35)

        self.play(FadeIn(nucleus))
        self.play(FadeIn(s3))
        self.play(FadeIn(px), FadeIn(py), FadeIn(pz))
        self.wait(0.6)

        # sp3 (niente gradienti)
        dirs = [
            np.array([+1, +1, +1]),
            np.array([+1, -1, -1]),
            np.array([-1, +1, -1]),
            np.array([-1, -1, +1]),
        ]
        sp3 = VGroup(*[sp3_lobe(d, R=1.65, opacity=0.40) for d in dirs])

        self.play(FadeOut(px), FadeOut(py), FadeOut(pz), FadeOut(s3), run_time=0.8)
        self.play(FadeIn(sp3), run_time=0.8)

        # “Respiro” via sola opacità (no gradienti)
        t = ValueTracker(0.0)
        def breathe_opacity(m, dt):
            t.increment_value(dt)
            alpha = 0.40 + 0.08*np.sin(TAU*0.7*t.get_value())
            for l in sp3:
                l.set_fill(POS_COLOR, opacity=alpha)
        sp3.add_updater(breathe_opacity)

        # “Vicini” come sferette semplici (niente Dot3D)
        neighbors = VGroup()
        for d in dirs:
            p = d / (np.linalg.norm(d)+1e-9) * 1.9
            s = Sphere(radius=0.06).move_to(p)
            s.set_fill(WHITE, opacity=0.7).set_stroke(width=0)
            neighbors.add(s)

        self.play(FadeIn(neighbors, scale=0.9))
        self.wait(1.0)

        self.play(FadeOut(neighbors), FadeOut(sp3), FadeIn(s3), FadeIn(px), FadeIn(py), FadeIn(pz), run_time=1.0)
        self.wait(0.5)