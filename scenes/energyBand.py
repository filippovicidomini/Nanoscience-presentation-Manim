from manim import *
import numpy as np

# =========================================================
# SCENA: ENERGY BANDS DELLA GIUNZIONE PN
# =========================================================

class PNJunctionEnergyBands(Scene):
    def construct(self):
        # ---- Stile coerente con le tue scene
        self.camera.background_color = "#0b0f14"

        # Palette (coerente con reticolo: elettroni blu, highlight giallo)
        COL_COND = "#ffe680"   # conduction band (giallo caldo)
        COL_VAL  = "#ffcc66"   # valence band (giallo più aranciato)
        COL_FERMI = "#6EA8FE"  # blu elettroni
        COL_AXIS = "#E6E6E6"
        COL_ACCENT = YELLOW_B  # per barrier highlight
        COL_BLOCK = RED        # reverse bias "blocco" (solo accento)

        def boxed_label(text, color, font_size=26, padding=0.16, opacity=0.80):
            t = Text(text, font_size=font_size, color=color)
            bg = RoundedRectangle(
                corner_radius=0.12,
                width=t.width + 2 * padding,
                height=t.height + 2 * padding,
                stroke_width=0,
                fill_color="#0b0f14",
                fill_opacity=opacity,
            )
            grp = VGroup(bg, t)
            grp.text = t
            grp.bg = bg
            return grp

        def pointer(from_point, to_point, color=WHITE, width=3):
            return Arrow(from_point, to_point, buff=0.05, stroke_width=width, color=color)

        # -----------------------------
        # Assi energia vs posizione
        # -----------------------------
        axes = Axes(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            x_length=11,
            y_length=6.6,
            tips=False,
            axis_config={"stroke_color": COL_AXIS, "stroke_width": 2, "include_numbers": False},
        ).to_edge(UP, buff=0.75)

        # Axis labels (boxed and moved out of busy area)
        x_label = boxed_label("position", color=COL_AXIS, font_size=22, opacity=0.65)
        x_label.next_to(axes, DOWN, buff=0.28)

        # Keep 'energy' horizontal (cleaner) and place it outside to the left
        y_label = boxed_label("energy", color=COL_AXIS, font_size=22, opacity=0.65)
        y_label.next_to(axes, LEFT, buff=0.32).shift(UP * 2.2)

        # vertical junction marker (x=0)
        junction_line = DashedLine(
            axes.c2p(0, -5), axes.c2p(0, 5),
            dash_length=0.12, color=COL_AXIS, stroke_width=2
        ).set_opacity(0.55)

        # Side labels (N / P) - y changed for more headroom
        label_N = Text("N", font_size=34, color=COL_AXIS).move_to(axes.c2p(-4.2, 4.25))
        label_P = Text("P", font_size=34, color=COL_AXIS).move_to(axes.c2p( 4.2, 4.25))

        def make_curve(xs, ys, color, width=4):
            pts = [axes.c2p(x, y) for x, y in zip(xs, ys)]
            return VMobject().set_points_smoothly(pts).set_stroke(color=color, width=width)

        stage_title = boxed_label("Before contact", color=COL_AXIS, font_size=28, opacity=0.72)
        # Place above the plot area, centered, so it doesn't collide with N/P labels
        stage_title.next_to(axes, UP, buff=0.18)

        self.play(
            Create(axes),
            FadeIn(x_label),
            FadeIn(y_label),
            FadeIn(junction_line),
            FadeIn(label_N),
            FadeIn(label_P),
            FadeIn(stage_title),
            run_time=1.6
        )
        self.wait(0.3)

        # =========================================================
        # 1) PRIMA DEL CONTATTO: due regioni separate con Fermi diversi
        # =========================================================

        # Diagrammi "separati": lasciamo un piccolo gap attorno a x=0
        xs_left  = np.linspace(-5, -0.5, 80)
        xs_right = np.linspace(0.5,  5, 80)

        # Bande piatte (prima del contatto)
        Ec_left_y  = np.full_like(xs_left,  1.6)  # conduction N
        Ev_left_y  = np.full_like(xs_left, -1.4)  # valence N
        Ec_right_y = np.full_like(xs_right, 1.9)  # conduction P leggermente diverso
        Ev_right_y = np.full_like(xs_right,-1.1)  # valence P

        # Fermi: N alto, P basso
        Ef_left_y  = 0.6
        Ef_right_y = -0.2

        Ec_left  = make_curve(xs_left,  Ec_left_y,  COL_COND, width=5)
        Ev_left  = make_curve(xs_left,  Ev_left_y,  COL_VAL,  width=5)
        Ef_left  = DashedLine(axes.c2p(-5, Ef_left_y),  axes.c2p(-0.5, Ef_left_y),
                              dash_length=0.10, color=COL_FERMI, stroke_width=3)

        Ec_right = make_curve(xs_right, Ec_right_y, COL_COND, width=5)
        Ev_right = make_curve(xs_right, Ev_right_y, COL_VAL,  width=5)
        Ef_right = DashedLine(axes.c2p(0.5, Ef_right_y), axes.c2p(5, Ef_right_y),
                              dash_length=0.10, color=COL_FERMI, stroke_width=3)

        Ec_before = VGroup(Ec_left, Ec_right)
        Ev_before = VGroup(Ev_left, Ev_right)
        Ef_before = VGroup(Ef_left, Ef_right)


        # Intro fade in
        self.play(
            Create(Ec_left), Create(Ev_left), Create(Ef_left),
            Create(Ec_right), Create(Ev_right), Create(Ef_right),
            run_time=2.0
        )
        # --- Helpers for dynamic curve-following labels ---
        def point_on_curve_near_x(curve: VMobject, x_world: float):
            pts = curve.get_points()
            if pts is None or len(pts) == 0:
                return ORIGIN
            idx = int(np.argmin(np.abs(pts[:, 0] - x_world)))
            return pts[idx]

        def curve_label(text: str, color, curve: VMobject, x_world: float, offset=RIGHT * 0.35):
            # Small boxed label that tracks a curve point near a given x
            lbl = boxed_label(text, color=color, font_size=22, opacity=0.70)
            p = point_on_curve_near_x(curve, x_world)
            lbl.move_to(p + offset)
            return lbl
        self.wait(1.0)  # <-- Voice: "Prima del contatto, ... N Fermi alto, P basso"

        # =========================================================
        # 2) DOPO IL CONTATTO (EQUILIBRIO): Fermi unico + band bending
        # =========================================================

        # Fermi unico: linea tratteggiata orizzontale su tutta la giunzione
        Ef_eq_y = 0.25
        Ef_eq = DashedLine(
            axes.c2p(-5, Ef_eq_y), axes.c2p(5, Ef_eq_y),
            dash_length=0.10, color=COL_FERMI, stroke_width=3
        )

        # Costruiamo curve con band bending (smooth step)
        xs = np.linspace(-5, 5, 220)

        # funzione di bending (sigmoide)
        def bend(x, x0=0.0, k=1.3):
            return np.tanh(k*(x-x0))

        # Ampiezza piegatura: crea "barriera" energetica vicino alla giunzione
        bend_amp = 0.75
        bend_profile = bend(xs, 0.0, 1.0)

        # conduction: sale verso destra (tipico band bending)
        Ec_eq_y = 1.75 + bend_amp * bend_profile
        Ev_eq_y = -1.25 + bend_amp * bend_profile

        Ec_eq = make_curve(xs, Ec_eq_y, COL_COND, width=5)
        Ev_eq = make_curve(xs, Ev_eq_y, COL_VAL,  width=5)

        # Evidenziazione barriera (una brace verticale in zona giunzione)
        barrier_brace = BraceBetweenPoints(
            axes.c2p(0.15, Ev_eq_y[np.argmin(np.abs(xs-0.15))]),
            axes.c2p(0.15, Ec_eq_y[np.argmin(np.abs(xs-0.15))]),
            direction=RIGHT
        ).set_color(COL_ACCENT).set_opacity(0.9)

        barrier_tag = boxed_label("Potential barrier", color=COL_ACCENT, font_size=24, opacity=0.70)
        barrier_tag.next_to(barrier_brace, RIGHT, buff=0.18)
        # Keep barrier label inside the plot area (avoid going out of frame)
        right_limit = axes.c2p(4.6, 0)[0]
        if barrier_tag.get_right()[0] > right_limit:
            barrier_tag.shift(LEFT * (barrier_tag.get_right()[0] - right_limit))

        bending_callout = boxed_label("Band bending", color=COL_ACCENT, font_size=24, opacity=0.70)
        bending_callout.move_to(axes.c2p(-3.55, 3.2))
        bending_arrow = pointer(
            bending_callout.get_right() + RIGHT * 0.05,
            axes.c2p(-0.2, 2.1),
            color=COL_ACCENT,
            width=3
        )

        stage_eq = boxed_label("After contact (equilibrium)", color=COL_AXIS, font_size=28, opacity=0.72)
        stage_eq.next_to(axes, UP, buff=0.18)
        # --- Dynamic labels that FOLLOW the current curves ---
        # Keep them inside the axes: choose an x close to the right edge of the plot area
        x_label_world = axes.c2p(4.15, 0)[0]

        Ec_lbl = always_redraw(lambda: curve_label("Ec", COL_COND, Ec_eq, x_label_world, offset=RIGHT * 0.35 + UP * 0.15))
        Ev_lbl = always_redraw(lambda: curve_label("Ev", COL_VAL,  Ev_eq, x_label_world, offset=RIGHT * 0.35 + DOWN * 0.15))
        Ef_lbl = always_redraw(lambda: curve_label("Ef", COL_FERMI, Ef_eq, x_label_world, offset=RIGHT * 0.35 + UP * 0.05))

        # Push Ef slightly up so it doesn't overlap the band labels when bending is small
        Ef_lbl.add_updater(lambda m: m.shift(UP * 0.25))

        self.add(Ec_lbl, Ev_lbl, Ef_lbl)

        # Transizione: fondi i due Fermi in uno, e piega le bande
        self.play(
            Transform(stage_title, stage_eq),
            ReplacementTransform(Ef_before, Ef_eq),
            ReplacementTransform(Ec_before, Ec_eq),
            ReplacementTransform(Ev_before, Ev_eq),
            run_time=2.2,
            rate_func=rate_functions.ease_in_out_cubic
        )

        self.play(
            FadeIn(barrier_brace),
            FadeIn(barrier_tag),
            FadeIn(bending_callout),
            GrowArrow(bending_arrow),
            run_time=0.9
        )
        self.wait(1.3)  # <-- Voice: "Quando mettiamo a contatto... Fermi unico... band bending... campo interno..."

        # Frase chiave (senza testo in scena, ma lasciamo tempo)
        self.wait(1.2)  # <-- Voice: "La barriera di potenziale nello spazio reale è la stessa che vediamo nello spazio energetico."

        # =========================================================
        # 3) FORWARD BIAS: barriera energetica si abbassa + crossing carriers
        # =========================================================
        # --- carrier helpers (DEVONO stare prima dell'uso) ---
        def electron_dot(x, y):
            return Dot(axes.c2p(x, y), radius=0.06, color=COL_FERMI)

        def hole_dot(x, y):
            return Circle(
                radius=0.075,
                stroke_width=3,
                color=YELLOW_C
            ).move_to(axes.c2p(x, y))
            
        # In forward: bending ridotto (barriera più bassa)
        bend_amp_fwd = 0.12
        Ec_fwd_y = 1.75 + bend_amp_fwd * bend_profile
        Ev_fwd_y = -1.25 + bend_amp_fwd * bend_profile
        Ec_fwd = make_curve(xs, Ec_fwd_y, COL_COND, width=5)
        Ev_fwd = make_curve(xs, Ev_fwd_y, COL_VAL,  width=5)

        # brace aggiornata più corta
        barrier_brace_fwd = BraceBetweenPoints(
            axes.c2p(0.15, Ev_fwd_y[np.argmin(np.abs(xs-0.15))]),
            axes.c2p(0.15, Ec_fwd_y[np.argmin(np.abs(xs-0.15))]),
            direction=RIGHT
        ).set_color(COL_ACCENT).set_opacity(0.9)

        stage_fwd = boxed_label("Forward bias", color=COL_AXIS, font_size=28, opacity=0.72)
        stage_fwd.next_to(axes, UP, buff=0.18)

        fwd_note = boxed_label("Barrier ↓", color=COL_ACCENT, font_size=26)
        fwd_note.move_to(axes.c2p(-3.7, -3.2))

        # Animazione di riduzione barriera
        self.play(
            Transform(stage_title, stage_fwd),
            Transform(Ec_eq, Ec_fwd),
            Transform(Ev_eq, Ev_fwd),
            Transform(barrier_brace, barrier_brace_fwd),
            FadeIn(fwd_note),
            FadeOut(bending_callout),
            FadeOut(bending_arrow),
            run_time=2.0,
            rate_func=rate_functions.ease_in_out_cubic
        )
        self.wait(0.6)  # <-- Voice: "In polarizzazione diretta... abbassa la barriera energetica..."

        # Streaming current (particle-like): single carriers emitted one-after-another
        n_particles = 100
        lanes_e = [1.75, 1.55, 1.35, 1.15, 0.95, 0.75]
        lanes_h = [-0.65, -0.85, -1.05, -1.25, -1.45, -1.65]

        # Pre-create particles (opacity 0) so we can run a clean LaggedStart
        electrons = VGroup(*[
            electron_dot(-4.35, float(np.random.choice(lanes_e))) for _ in range(n_particles)
        ])
        holes = VGroup(*[
            hole_dot(4.35, float(np.random.choice(lanes_h))) for _ in range(n_particles)
        ])

        electrons.set_opacity(0)
        holes.set_opacity(0)
        self.add(electrons, holes)

        def particle_trip(mobj, shift_vec):
            # appear -> travel -> disappear (symmetric: same fade & scale)
            t_fade = 0.22
            t_travel = 0.90
            s = 0.92
            return Succession(
                FadeIn(mobj, scale=s, run_time=t_fade, rate_func=rush_into),
                mobj.animate.shift(shift_vec).set_opacity(1).set_run_time(t_travel).set_rate_func(rush_into),
                FadeOut(mobj, scale=s, run_time=t_fade, rate_func=rush_into),
            )

        e_anims = [particle_trip(e, RIGHT * 6.9) for e in electrons]
        h_anims = [particle_trip(h, LEFT * 6.9) for h in holes]

        # Staggered start times -> continuous flow (not a block)
        self.play(
            AnimationGroup(
                LaggedStart(*e_anims, lag_ratio=0.1),
                LaggedStart(*h_anims, lag_ratio=0.1),
                lag_ratio=0.0,
            ),
            rate_func=rush_into,
            run_time=7.4,
        )

        # Remove particles (clean scene)
        self.remove(electrons, holes)

        self.wait(0.8)

        # =========================================================
        # 4) REVERSE BIAS: barriera aumenta, bande più inclinate, flusso bloccato
        # =========================================================

        bend_amp_rev = 1.55
        Ec_rev_y = 1.75 + bend_amp_rev * bend_profile
        Ev_rev_y = -1.25 + bend_amp_rev * bend_profile
        Ec_rev = make_curve(xs, Ec_rev_y, COL_COND, width=5)
        Ev_rev = make_curve(xs, Ev_rev_y, COL_VAL,  width=5)

        barrier_brace_rev = BraceBetweenPoints(
            axes.c2p(0.15, Ev_rev_y[np.argmin(np.abs(xs-0.15))]),
            axes.c2p(0.15, Ec_rev_y[np.argmin(np.abs(xs-0.15))]),
            direction=RIGHT
        ).set_color(COL_ACCENT).set_opacity(0.95)

        stage_rev = boxed_label("Reverse bias", color=COL_AXIS, font_size=28, opacity=0.72)
        stage_rev.next_to(axes, UP, buff=0.18)

        rev_note = boxed_label("Barrier ↑", color=COL_ACCENT, font_size=26)
        rev_note.move_to(axes.c2p(-3.7, -3.2))

        self.play(
            Transform(stage_title, stage_rev),
            Transform(fwd_note, rev_note),
            Transform(Ec_eq, Ec_rev),
            Transform(Ev_eq, Ev_rev),
            Transform(barrier_brace, barrier_brace_rev),
            run_time=2.0,
            rate_func=rate_functions.ease_in_out_cubic
        )

        # effetto "blocco": una X tenue sulla giunzione
        block = boxed_label("Current blocked", color=COL_BLOCK, font_size=26, opacity=0.65)
        block.move_to(axes.c2p(0.0, -3.25))

        x_mark = VGroup(
            Line(axes.c2p(-0.35, -2.9), axes.c2p(0.35, -2.2), color=COL_BLOCK, stroke_width=5),
            Line(axes.c2p(-0.35, -2.2), axes.c2p(0.35, -2.9), color=COL_BLOCK, stroke_width=5),
        ).set_opacity(0.55)

        self.play(FadeIn(block), FadeIn(x_mark), run_time=0.6)
        self.wait(1.6)  # <-- Voice: "In polarizzazione inversa... barriera aumenta... flusso bloccato."

        # Outro clean
        self.play(
            FadeOut(block),
            FadeOut(x_mark),
            FadeOut(barrier_tag),
            FadeOut(barrier_brace),
            FadeOut(fwd_note),
            run_time=0.7
        )
        self.wait(0.6)

        # Cleanup updaters (se mai)
        for m in self.mobjects:
            m.clear_updaters()