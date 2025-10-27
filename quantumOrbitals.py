from manim import *
import numpy as np

# Palette
POS = BLUE_C
NEG = RED_C
NUC = GRAY_B

def to3d(v):
    v = np.array(v, dtype=float)
    if v.shape == (2,):
        v = np.array([v[0], v[1], 0.0])
    elif v.shape != (3,):
        # fallback prudente
        v = np.array([1.0, 0.0, 0.0])
    return v

def lobo_ellittico(dir_vec=RIGHT, length=2.2, width=1.1, color=POS, alpha=0.42):
    """
    Lobo 'tipo orbitale' super leggero: Ellipse allungata orientata lungo dir_vec.
    dir_vec viene sempre trattato come vettore 3D (z=0 se mancante).
    """
    d3 = to3d(dir_vec)
    n = np.linalg.norm(d3)
    if n < 1e-9:
        d3 = np.array([1.0, 0.0, 0.0])
        n = 1.0
    d3 /= n

    center = 0.9 * d3  # 3D
    ell = Ellipse(width=length, height=width, color=color, stroke_width=0)
    ell.set_fill(color, opacity=alpha)
    # orienta nel piano xy
    ang = np.arctan2(d3[1], d3[0])
    ell.rotate(ang).move_to(center)
    return ell

def doppio_lobo(dir_vec=RIGHT, length=2.4, width=1.1, color_pos=POS, color_neg=NEG, alpha=0.42):
    d3 = to3d(dir_vec)
    return VGroup(
        lobo_ellittico(d3,  length, width, color_pos, alpha),
        lobo_ellittico(-d3, length, width, color_neg, alpha),
    )

def s_orbital(levels=1, base_radius=0.85):
    g = VGroup()
    for k in range(levels):
        R = base_radius * (1 + 0.45*k)
        col = POS if k % 2 == 0 else NEG
        disk = Circle(radius=R, color=col, stroke_width=0).set_fill(col, opacity=0.18 if k == 0 else 0.10)
        g.add(disk)
    nuc = Circle(radius=0.10, color=NUC, stroke_width=0).set_fill(NUC, opacity=0.6)
    g.add(nuc)
    return g

def sp3_quattro_lobi(scale=1.0, alpha=0.42):
    # Direzioni tetraedriche proiettate (in piano: z=0)
    dirs = [
        np.array([+1, +1, 0.0]),
        np.array([+1, -1, 0.0]),
        np.array([-1, +1, 0.0]),
        np.array([-1, -1, 0.0]),
    ]
    g = VGroup()
    for i, d in enumerate(dirs):
        col = POS if i % 2 == 0 else NEG
        g.add(lobo_ellittico(d, length=2.3*scale, width=1.1*scale, color=col, alpha=alpha))
    nuc = Circle(radius=0.10*scale, color=NUC, stroke_width=0).set_fill(NUC, opacity=0.6)
    g.add(nuc)
    return g

class QMOrbitalsUltraLite2D(Scene):
    def construct(self):
        self.camera.background_color = "#0b0e10"

        col_x = [-4.5, 0.0, +4.5]
        row_y = [+1.8, -2.0]

        s1 = s_orbital(levels=1, base_radius=0.95).move_to(np.array([col_x[0], row_y[0], 0]))
        s2 = s_orbital(levels=2, base_radius=0.80).move_to(np.array([col_x[0], row_y[1], 0]))

        px = doppio_lobo(RIGHT,      length=2.6, width=1.1).move_to(np.array([col_x[1], row_y[0], 0]))
        py = doppio_lobo(UP,         length=2.6, width=1.1).move_to(np.array([col_x[1], row_y[1], 0]))
        pz = doppio_lobo(RIGHT+UP,   length=2.6, width=1.1).move_to(np.array([col_x[2], row_y[0], 0]))

        sp3 = sp3_quattro_lobi(scale=1.0, alpha=0.42).move_to(np.array([col_x[2], row_y[1], 0]))

        self.play(FadeIn(s1, shift=0.2*UP), FadeIn(s2, shift=0.2*DOWN), run_time=0.6)
        self.play(FadeIn(px, shift=0.2*LEFT), FadeIn(py, shift=0.2*DOWN), run_time=0.6)
        self.play(FadeIn(pz, shift=0.2*UP), FadeIn(sp3, shift=0.2*RIGHT), run_time=0.6)

        # opzionale: leggerissima pulsazione dei lobi
        # opzionale: leggerissima pulsazione dei lobi
        amp = 0.03
        t = ValueTracker(0.0)

        def breathe_opacity(_mobj, dt):   # updater di Mobject: (mobj, dt)
            t.increment_value(dt)
            k = 0.40 + amp*np.sin(TAU*0.5*t.get_value())
            for grp in [px, py, pz, sp3]:
                for l in grp:
                    if isinstance(l, VMobject):
                        l.set_fill(l.get_color(), opacity=k)

        # crea un “ticker” invisibile a cui attaccare l’updater
        ticker = VMobject()
        self.add(ticker)
        ticker.add_updater(breathe_opacity)

        self.wait(1.8)

        # (facoltativo) pulisci alla fine
        ticker.remove_updater(breathe_opacity)