from manim import *
import numpy as np
import math

# Palette
NUCLEUS_COLOR = GRAY_B
SHELL_COLOR   = GRAY_D
E_CORE_COLOR  = BLUE_C     # K, L
E_VAL_COLOR   = YELLOW     # M (valenza)
POS_COLOR     = BLUE_C     # fase +
NEG_COLOR     = RED_C      # fase -

# ---------- Helpers “orbitali” ultra-light (2D) ----------
def to3d(v):
    v = np.array(v, dtype=float)
    if v.shape == (2,):
        v = np.array([v[0], v[1], 0.0])
    elif v.shape != (3,):
        v = np.array([1.0, 0.0, 0.0])
    return v

def lobo_ellittico(dir_vec=RIGHT, length=2.2, width=1.1, color=POS_COLOR, alpha=0.42):
    d3 = to3d(dir_vec); n = np.linalg.norm(d3) or 1.0; d3 = d3 / n
    center = 0.9 * d3
    ell = Ellipse(width=length, height=width, color=color, stroke_width=0)
    ell.set_fill(color, opacity=alpha)
    ang = math.atan2(d3[1], d3[0])
    return ell.rotate(ang).move_to(center)

def doppio_lobo(dir_vec=RIGHT, length=2.4, width=1.15, cpos=POS_COLOR, cneg=NEG_COLOR, alpha=0.42):
    d3 = to3d(dir_vec)
    return VGroup(
        lobo_ellittico(d3,  length, width, cpos, alpha),
        lobo_ellittico(-d3, length, width, cneg, alpha),
    )

def orbitale_s_stilizzato(levels=2, base_radius=0.9):
    g = VGroup()
    for k in range(levels):
        R = base_radius * (1 + 0.45*k)
        col = POS_COLOR if k % 2 == 0 else NEG_COLOR
        disk = Circle(radius=R, color=col, stroke_width=0).set_fill(col, opacity=0.14 if k == 0 else 0.10)
        g.add(disk)
    nuc = Circle(radius=0.12, color=NUCLEUS_COLOR, stroke_width=0).set_fill(NUCLEUS_COLOR, opacity=0.55)
    g.add(nuc)
    return g

# ---------- Elettroni sui gusci (2D) ----------
def electron_on_shell(center, R, speed=0.7, color=E_CORE_COLOR, radius=0.06, phase=0.0, start_active=False):
    dot = Dot(radius=radius, color=color)
    theta = ValueTracker(phase)
    running = [start_active]
    def upd(m, dt):
        if not running[0]:
            return
        dt = min(max(dt, 0.0), 1/60)       # evita “salti” sul primo frame
        theta.increment_value(speed * dt * TAU)
        ang = theta.get_value()
        pos = center + R * np.array([np.cos(ang), np.sin(ang), 0.0])
        m.move_to(pos)
    dot.add_updater(upd)
    dot.move_to(center + R * np.array([np.cos(phase), np.sin(phase), 0.0]))
    dot.running = running
    return dot

class SiliconAtomFull2D(Scene):
    def construct(self):
        self.camera.background_color = "#0b0e10"

        # ---- Fase A: “orbitali” qualitativi (s + p_x, p_y, p_z) ----
        # s centrato
        s_vis = orbitale_s_stilizzato(levels=2, base_radius=0.92)

        # tre p ortogonali (doppio lobo con fasi ±)
        px = doppio_lobo(RIGHT,       length=2.6, width=1.15)
        py = doppio_lobo(UP,          length=2.6, width=1.15)
        pz = doppio_lobo(RIGHT+UP,    length=2.6, width=1.15)

        # Nucleo sobrio (rimane attraverso le fasi)
        nucleus = Circle(radius=0.12, color=NUCLEUS_COLOR, stroke_width=0).set_fill(NUCLEUS_COLOR, opacity=0.55)

        # Intro orbitali
        self.play(FadeIn(nucleus), run_time=0.4)
        self.play(FadeIn(s_vis), run_time=0.4)
        #self.play(FadeIn(pz), run_time=0.4)
        #self.play(ReplacementTransform(pz.copy(), px), run_time=0.4)
        #self.play(ReplacementTransform(px.copy(), py), run_time=0.4)
        #self.play(ReplacementTransform(py.copy(), pz), run_time=0.4)
        self.play(FadeIn(px), runtime =0.4, lag_ratio=0.2)
        self.play(FadeIn(py), runtime =0.4, lag_ratio=0.2)
        self.play(FadeIn(pz), runtime =0.4, lag_ratio=0.2)
        # Pulizia verso la fase “gusci”
        self.play(FadeOut(px), FadeOut(py), FadeOut(pz), FadeOut(s_vis), run_time=0.5)

        # ---- Fase B: gusci K/L/M con 2,8,4 elettroni ----
        C = np.array([0.0, 0.0, 0.0])
        R_K, R_L, R_M = 0.55, 1.05, 1.55

        shellK = Circle(radius=R_K, color=SHELL_COLOR, stroke_width=2).set_stroke(opacity=0.6).move_to(C)
        shellL = Circle(radius=R_L, color=SHELL_COLOR, stroke_width=2).set_stroke(opacity=0.6).move_to(C)
        shellM = Circle(radius=R_M, color=SHELL_COLOR, stroke_width=2).set_stroke(opacity=0.75).move_to(C)

        eK_ph = [0, np.pi]
        eL_ph = [k*np.pi/4 for k in range(8)]
        eM_ph = [0, np.pi/2, np.pi, 3*np.pi/2]

        eK = VGroup(*[electron_on_shell(C, R_K, speed=0.85, color=E_CORE_COLOR, phase=ph) for ph in eK_ph])
        eL = VGroup(*[electron_on_shell(C, R_L, speed=0.70, color=E_CORE_COLOR, phase=ph) for ph in eL_ph])
        eM = VGroup(*[electron_on_shell(C, R_M, speed=0.58, color=E_VAL_COLOR,  phase=ph) for ph in eM_ph])

        # Disegna gusci + elettroni
        self.play(FadeIn(shellK), FadeIn(shellL), FadeIn(shellM), run_time=0.5)
        self.play(FadeIn(eK), FadeIn(eL), FadeIn(eM), run_time=0.5)

        # Avvia il moto dopo l’ingresso
        for d in list(eK) + list(eL) + list(eM):
            d.running[0] = True

        self.wait(1.4)

        # ---- Fase C: focus valenza (opzionale) ----
        for d in list(eK) + list(eL):
            d.running[0] = False
        self.play(
            FadeOut(shellK, rate_func=linear),
            FadeOut(shellL, rate_func=linear),
            FadeOut(eK,     rate_func=linear),
            FadeOut(eL,     rate_func=linear),
            run_time=0.7
        )
        # Lascia nucleus, shell M ed eM in scena per 1 s
        self.wait(1.0)