from manim import *
from manim import config
import numpy as np
import random

# Forza OpenGL: evita problemi di gradienti Cairo
config.renderer = "opengl"

POS = BLUE_C      # fase +
NEG = RED_C       # fase -
NUCLEUS = GRAY_B

# --- Sferiche armoniche reali: l=0,1 (s e p) ---
def Y_real(l, m, theta, phi):
    if l == 0:
        return 0.5/np.sqrt(np.pi) * np.ones_like(theta)
    if l == 1:
        c = np.sqrt(3/(4*np.pi))
        if m == 0:   # p_z ~ cos(theta)
            return c * np.cos(theta)
        if m == +1:  # p_x ~ sin(theta) cos(phi)
            return c * np.sin(theta) * np.cos(phi)
        if m == -1:  # p_y ~ sin(theta) sin(phi)
            return c * np.sin(theta) * np.sin(phi)
    return np.zeros_like(theta)

# --- Radiali idrogeno-like (qualitativo): 3s, 3p ---
def R_nl(n, l, r, a0=1.0):
    x = r/a0
    if n == 3 and l == 0:   # 3s: due nodi radiali
        return (1 - (2*x)/3 + (2*x**2)/27) * np.exp(-x/3)
    if n == 3 and l == 1:   # 3p: un nodo radiale
        return (x) * (1 - x/6) * np.exp(-x/3)
    return np.exp(-r/(n*a0))

# --- Superficie angolare r(θ,φ) ∝ |Y_lm|: mostra lobi e piani nodali (niente gradient) ---
def angular_lobe_surface(l=1, m=0, scale=1.6, opacity=0.40):
    def base(u, v):
        th, ph = u, v
        Y = Y_real(l, m, th, ph)
        r = scale * np.abs(Y)
        x = r * np.sin(th) * np.cos(ph)
        y = r * np.sin(th) * np.sin(ph)
        z = r * np.cos(th)
        return np.array([x, y, z])

    res = (28, 56)  # moderato per stabilità

    pos = Surface(
        lambda u, v: base(u, v) if Y_real(l, m, u, v) >= 0 else np.array([np.nan, np.nan, np.nan]),
        u_range=[0, PI], v_range=[0, TAU], resolution=res
    ).set_fill(POS, opacity=opacity).set_stroke(width=0)

    neg = Surface(
        lambda u, v: base(u, v) if Y_real(l, m, u, v) < 0 else np.array([np.nan, np.nan, np.nan]),
        u_range=[0, PI], v_range=[0, TAU], resolution=res
    ).set_fill(NEG, opacity=opacity).set_stroke(width=0)

    return VGroup(pos, neg)

# --- Nuvola punti da |ψ|^2 via rejection sampling ---
def sample_orbital_points(n=3, l=1, m=0, N=1800, a0=1.0, r_max=8.0, max_trials=180000):
    pts, signs = [], []
    def pdf(r, th, ph):
        Y = Y_real(l, m, th, ph)
        return (R_nl(n, l, r, a0)**2) * (Y**2) * (r**2) * np.sin(th)

    M = 0.1  # bound iniziale
    trials = 0
    while len(pts) < N and trials < max_trials:
        trials += 1
        r  = np.random.rand() * r_max
        th = np.arccos(1 - 2*np.random.rand())   # uniforme sulla sfera
        ph = np.random.rand() * TAU
        val = pdf(r, th, ph)
        if val > M:
            M = val
        if np.random.rand() * M < val:
            x = r * np.sin(th) * np.cos(ph)
            y = r * np.sin(th) * np.sin(ph)
            z = r * np.cos(th)
            pts.append([x, y, z])
            signs.append(1 if Y_real(l, m, th, ph) >= 0 else -1)

    return np.array(pts, dtype=float), np.array(signs, dtype=int)

def point_cloud_from_orbital(n=3, l=1, m=0, N=1800, scale=0.45, dot_r=0.03):
    P, S = sample_orbital_points(n, l, m, N=N)
    P = P * scale
    grp = VGroup()
    for i in range(P.shape[0]):
        d = Sphere(radius=dot_r).move_to(P[i])
        d.set_fill(POS if S[i] > 0 else NEG, opacity=0.95).set_stroke(width=0)
        grp.add(d)
    return grp

class QMSiliconLikeOrbitalsSafe(ThreeDScene):
    def construct(self):
        self.camera.background_color = "#0b0e10"
        self.set_camera_orientation(phi=65*DEGREES, theta=35*DEGREES)

        # ROOT per "zoom" robusto (scala il gruppo invece della camera)
        ROOT = VGroup()
        self.add(ROOT)

        # Nucleo
        nucleus = Sphere(radius=0.12).set_fill(NUCLEUS, opacity=0.5).set_stroke(width=0)
        ROOT.add(nucleus)

        # 3s: nuvola doppio guscio (nodo radiale), superficie angolare isotropa tenue
        cloud_3s_inner = point_cloud_from_orbital(n=3, l=0, m=0, N=900,  scale=0.22, dot_r=0.028)
        cloud_3s_outer = point_cloud_from_orbital(n=3, l=0, m=0, N=1100, scale=0.42, dot_r=0.028)
        s_shell = angular_lobe_surface(l=0, m=0, scale=1.0, opacity=0.12)
        group_3s = VGroup(s_shell, cloud_3s_inner, cloud_3s_outer)
        ROOT.add(group_3s)

        # 3p_z: superficie angolare + nuvola con piano nodale
        pz_surface = angular_lobe_surface(l=1, m=0, scale=1.8, opacity=0.40)
        cloud_3p   = point_cloud_from_orbital(n=3, l=1, m=0, N=2200, scale=0.60, dot_r=0.028)

        # Mostra 3s
        self.play(FadeIn(nucleus))
        self.play(FadeIn(group_3s))
        self.wait(0.8)

        # “Zoom” dolce via ROOT (compatibile ovunque)
        self.play(ROOT.animate.scale(1.1), run_time=0.6)

        # Passa a 3p_z
        self.play(FadeOut(group_3s, run_time=0.7))
        ROOT.add(pz_surface, cloud_3p)
        self.play(FadeIn(pz_surface, run_time=0.7))
        self.play(FadeIn(cloud_3p, run_time=0.5))
        self.wait(1.0)

        # Mostra rapidamente p_x e p_y (rotazioni tramite ReplacementTransform)
        px_surface = angular_lobe_surface(l=1, m=+1, scale=1.8, opacity=0.40)
        py_surface = angular_lobe_surface(l=1, m=-1, scale=1.8, opacity=0.40)
        self.play(ReplacementTransform(pz_surface.copy(), px_surface), run_time=0.6)
        self.play(ReplacementTransform(px_surface.copy(), py_surface), run_time=0.6)
        self.play(ReplacementTransform(py_surface.copy(), pz_surface), run_time=0.6)
        self.wait(0.5)

        # Uscita pulita
        self.play(FadeOut(cloud_3p), FadeOut(pz_surface), FadeOut(nucleus))
        self.wait(0.3)