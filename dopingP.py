from manim import *
import numpy as np

# Palette
SI_COLOR = GRAY_B
P_COLOR  = BLUE_B
E_COLOR  = YELLOW
RAIL_COLOR = WHITE

def make_atom(center, label="Si", color=SI_COLOR, r=0.35):
    ring = Circle(radius=r, color=color, stroke_width=4).move_to(center)
    fill = Circle(radius=r, color=color, stroke_width=0).move_to(center)
    fill.set_opacity(0.15)
    txt  = Text(label, font_size=22).move_to(center)
    return VGroup(fill, ring, txt)

def bond(p, q):
    seg = Line(p, q, color=WHITE, stroke_width=5)
    seg.set_opacity(0.25)
    return seg

def electron_on_bond(p, q, speed=0.8, radius=0.06, color=E_COLOR):
    """Elettrone che si muove avanti/indietro lungo il legame p<->q."""
    d = Dot(radius=radius, color=color)
    t = ValueTracker(np.random.rand())
    direction = {"d": 1}
    def upd(m, dt):
        new_t = t.get_value() + direction["d"] * speed * dt
        if new_t > 1.0:
            new_t = 2.0 - new_t
            direction["d"] *= -1
        if new_t < 0.0:
            new_t = -new_t
            direction["d"] *= -1
        t.set_value(new_t)
        pos = (1-new_t) * p + new_t * q
        m.move_to(pos)
    d.add_updater(upd)
    d.move_to((1-t.get_value())*p + t.get_value()*q)
    return d

class NTypeDopingPhosphorus2D(Scene):
    # Gioca animazioni solo se ce ne sono
    def play_if_any(self, anims, **kwargs):
        anims = [a for a in anims if isinstance(a, Animation)]
        if len(anims) == 0:
            return self.wait(0.001)
        return self.play(*anims, **kwargs)

    def construct(self):
        self.camera.background_color = "#111111"

        # --- Layout e “zone” --------------------------------------------------
        # Frame tipico: width ~14, height ~8
        TITLE_Y = 3.5
        RAIL_Y  = 2.2
        GRID_Y0 = -0.2  # abbassiamo la griglia
        LEGEND_Y = -3.4

        # Zona didascalie a destra (pannello fisso)
        panel_w, panel_h = 3.4, 3.2
        panel_center = np.array([5.0, 0.8, 0])  # x>0 per stare a destra
        caption_panel = RoundedRectangle(corner_radius=0.15, width=panel_w, height=panel_h, color=WHITE)
        caption_panel.set_fill(WHITE, opacity=0.06).set_stroke(opacity=0.6)
        caption_panel.move_to(panel_center)

        caption = Text("", font_size=24, line_spacing=0.9).move_to(panel_center)
        def set_caption(t):
            caption.become(Text(t, font_size=24, line_spacing=0.9).move_to(panel_center))

        # Titolo compatto
        title = Text("n-type (P): elettrone donato e conduzione", font_size=32).move_to(np.array([0, TITLE_Y, 0]))

        # --- Griglia atomi 2x3 ------------------------------------------------
        rows, cols = 2, 3
        dx, dy = 2.6, 2.0
        x0, y0 = -dx, GRID_Y0 + 0.6
        centers = [[np.array([x0 + c*dx, y0 - r*dy, 0]) for c in range(cols)] for r in range(rows)]

        atoms = VGroup()
        bonds = VGroup()
        edges = []
        for r in range(rows):
            for c in range(cols):
                atoms.add(make_atom(centers[r][c]))
        for r in range(rows):
            for c in range(cols):
                p = centers[r][c]
                if c+1 < cols:
                    q = centers[r][c+1]
                    bonds.add(bond(p, q)); edges.append((p, q))
                if r+1 < rows:
                    q = centers[r+1][c]
                    bonds.add(bond(p, q)); edges.append((p, q))

        electrons = VGroup()
        for (p, q) in edges:
            electrons.add(
                electron_on_bond(p, q, speed=0.9, radius=0.06),
                electron_on_bond(q, p, speed=0.9, radius=0.06)
            )

        # Legenda compatta in basso
        legend = VGroup(
            VGroup(Dot(radius=0.09, color=E_COLOR), Text("e⁻ (elettrone)", font_size=20, color=E_COLOR)).arrange(RIGHT, buff=0.2),
            VGroup(Line(ORIGIN, 0.8*RIGHT, stroke_width=6, color=RAIL_COLOR), Text("corsia di conduzione", font_size=20, color=RAIL_COLOR)).arrange(RIGHT, buff=0.2),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to(np.array([-6.0, LEGEND_Y, 0]))
        legend_panel = SurroundingRectangle(legend, color=WHITE, buff=0.25, corner_radius=0.15).set_fill(WHITE, opacity=0.06)

        # --- Corsia di conduzione (limitata a sinistra per non invadere il pannello) ---
        xL, xR = -6.0, 3.4  # ATTENZIONE: xR < pannello (che inizia ~3.3)
        rail = Line(np.array([xL, RAIL_Y, 0]), np.array([xR, RAIL_Y, 0]), stroke_width=6, color=RAIL_COLOR).set_opacity(0.35)
        efield = Arrow(np.array([xL, RAIL_Y+0.5, 0]), np.array([xR, RAIL_Y+0.5, 0]), buff=0).set_stroke(width=6)
        efield_lbl = Text("Campo elettrico", font_size=22).next_to(efield, DOWN, buff=0.12)

        # --- Intro ------------------------------------------------------------
        self.play(FadeIn(bonds), *[FadeIn(a) for a in atoms], FadeIn(title))
        self.play(FadeIn(electrons))
        self.play(FadeIn(rail), FadeIn(efield), FadeIn(efield_lbl))
        self.play(FadeIn(caption_panel))
        set_caption("Silicio intrinseco:\nElettroni di legame si muovono lungo i legami.")
        self.play(FadeIn(caption))
        self.play(FadeIn(legend_panel), FadeIn(legend))
        self.wait(0.6)

        # --- Inseriamo Fosforo (donor) in alto-centro (0,1) -------------------
        rP, cP = 0, 1
        phosphorus = make_atom(centers[rP][cP], label="P", color=P_COLOR)
        self.play(Transform(atoms[rP*cols + cP], phosphorus))
        set_caption("Dopiamo con Fosforo (donor):\n5 e⁻ di valenza → 1 elettrone in eccesso.")
        self.play(caption.animate.set_opacity(1.0))
        self.wait(0.2)

        # Elettrone extra vicino a P
        extra_e = Dot(radius=0.08, color=E_COLOR).move_to(centers[rP][cP] + 0.5*UP)
        self.play(FadeIn(extra_e, scale=0.7))

        # --- Ionizzazione del donor → elettrone sale sulla corsia -------------
        set_caption("Il donor si ionizza facilmente:\nL’elettrone extra diventa (quasi) libero\n→ banda di conduzione.")
        self.play(caption.animate.set_opacity(1.0))

        up_point = centers[rP][cP] + 1.1*UP
        self.play(extra_e.animate.move_to(up_point), run_time=0.6)
        # proiezione sulla rail (limitata a xR per non entrare nel pannello)
        proj_x = np.clip(centers[rP][cP][0], xL+0.5, xR-0.5)
        on_rail = np.array([proj_x, RAIL_Y, 0])
        self.play(extra_e.animate.move_to(on_rail), run_time=0.6)

        # Drift lungo rail (ciclo tra xL e xR senza invadere la zona testi)
        t = ValueTracker(0.15)
        def update_free(m, dt):
            v = 1.6
            new_t = (t.get_value() + v*dt) % 1.0
            t.set_value(new_t)
            x = (1-new_t)*xL + new_t*xR
            m.move_to(np.array([x, RAIL_Y, 0]))
        extra_e.add_updater(update_free)
        self.wait(1.0)

        # --- Nota finale nel pannello e piccolo drift dei legami --------------
        set_caption("Conduzione n-type:\nGli elettroni di legame restano;\nquello donato scorre lungo la corsia.")
        self.play(caption.animate.set_opacity(1.0))

        drift = 0.15*RIGHT
        movers = [e for e in electrons if e in self.mobjects]
        self.play_if_any([e.animate.shift(drift) for e in movers], run_time=1.0)

        # Chiusura
        final = Text("P (donor) → elettrone extra quasi libero → corrente n-type", font_size=26).to_edge(DOWN)
        self.play(FadeIn(final))
        self.wait(1.0)