from manim import *
import numpy as np

# =========================================================
# SCENE: PN JUNCTION – RECOMBINATION & PHOTON (LED)
# =========================================================

class PNJunctionLED(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # -------------------------------------------------
        # Color palette (coherent with previous scenes)
        # -------------------------------------------------
        COL_ELECTRON = "#6EA8FE"
        COL_HOLE = YELLOW_C
        COL_PHOTON_RED = RED_C
        COL_PHOTON_GREEN = GREEN_C
        COL_PHOTON_BLUE = BLUE_C
        COL_AXIS = "#E6E6E6"
        COL_BOX = "#0b0f14"

        # -------------------------------------------------
        # Helpers
        # -------------------------------------------------
        def boxed_label(text, color, font_size=26, padding=0.18, opacity=0.85):
            t = Text(text, font_size=font_size, color=color)
            bg = RoundedRectangle(
                corner_radius=0.12,
                width=t.width + 2 * padding,
                height=t.height + 2 * padding,
                stroke_width=0,
                fill_color=COL_BOX,
                fill_opacity=opacity,
            )
            return VGroup(bg, t)

        def electron(pos):
            return Dot(pos, radius=0.07, color=COL_ELECTRON)

        def hole(pos):
            return Circle(radius=0.09, stroke_width=3, color=COL_HOLE).move_to(pos)

        def photon_wave(color, radius=0.22):
            # Larger initial radius so emission clearly reads as light
            return Circle(radius=radius, stroke_width=5, color=color)

        # -------------------------------------------------
        # Axes: position vs energy (simplified)
        # -------------------------------------------------
        axes = Axes(
            x_range=[-5, 5, 1],
            y_range=[-3, 3, 1],
            x_length=10,
            y_length=5,
            tips=False,
            axis_config={"stroke_color": COL_AXIS, "stroke_width": 2},
        ).to_edge(UP, buff=0.8)

        x_label = boxed_label("position", COL_AXIS, font_size=22)
        y_label = boxed_label("energy", COL_AXIS, font_size=22)
        x_label.next_to(axes, DOWN, buff=0.3)
        y_label.next_to(axes, LEFT, buff=0.3).shift(UP * 1.5)

        self.play(Create(axes), FadeIn(x_label), FadeIn(y_label), run_time=1.4)

        # -------------------------------------------------
        # Bands (slightly curved like the previous scene)
        # -------------------------------------------------
        xs = np.linspace(-5, 5, 240)

        # Soft bending profile (small in forward bias, just for continuity)
        k_bend = 0.9
        bend_amp = 0.22
        bend_profile = np.tanh(k_bend * xs)

        def band_y(x, level):
            # band value at a single x (for placing carriers/pointers)
            return level + bend_amp * np.tanh(k_bend * x)

        def make_band(level, color):
            ys = level + bend_amp * bend_profile
            pts = [axes.c2p(x, y) for x, y in zip(xs, ys)]
            return VMobject().set_points_smoothly(pts).set_stroke(color=color, width=5)

        Ec_level = 1.4
        Ev_level = -1.4
        Ec = make_band(Ec_level, "#ffe680")
        Ev = make_band(Ev_level, "#ffcc66")

        Ec_lbl = boxed_label("Ec", "#ffe680", font_size=22)
        Ev_lbl = boxed_label("Ev", "#ffcc66", font_size=22)

        # Place labels near the right edge, aligned with the levels
        Ec_lbl.move_to(axes.c2p(4.25, 1.85))
        Ev_lbl.move_to(axes.c2p(4.25, -1.85))

        # Small pointers so it is obvious which line is which
        Ec_ptr = Arrow(Ec_lbl.get_left(), axes.c2p(4.6, band_y(4.6, Ec_level)), buff=0.08, stroke_width=3, color="#ffe680")
        Ev_ptr = Arrow(Ev_lbl.get_left(), axes.c2p(4.6, band_y(4.6, Ev_level)), buff=0.08, stroke_width=3, color="#ffcc66")

        self.play(Create(Ec), Create(Ev), FadeIn(Ec_lbl), FadeIn(Ev_lbl), GrowArrow(Ec_ptr), GrowArrow(Ev_ptr), run_time=1.6)

        # -------------------------------------------------
        # Title
        # -------------------------------------------------
        title = boxed_label("Recombination and photon emission (LED)", COL_AXIS, font_size=28)
        title.next_to(axes, UP, buff=0.25)
        self.play(FadeIn(title))
        self.wait(0.5)

        # -------------------------------------------------
        # Forward bias: repeated carrier injection -> recombination -> red photon
        # -------------------------------------------------
        e_start = axes.c2p(-4.2, band_y(-4.2, Ec_level))
        h_start = axes.c2p(4.2,  band_y(4.2,  Ev_level))
        # Use the axes' true origin to avoid any subtle coordinate mismatch
        recomb_point = axes.get_origin()

        def photon_pulse(color, center, max_scale=6.8):
            # Bigger expanding wave -> clear light emission
            p = photon_wave(color).move_to(center)
            # Scale about the emission point to avoid any drift
            return Succession(
                Create(p),
                p.animate.scale(max_scale, about_point=center).set_opacity(0).set_run_time(0.9),
                FadeOut(p, run_time=0.01),
            )

        # Electron "fall" arrow (reused)
        fall_arrow = Arrow(
            recomb_point + UP * 0.6,
            recomb_point + DOWN * 0.6,
            buff=0,
            stroke_width=4,
            color=COL_ELECTRON,
        ).set_opacity(0.85)

        self.play(FadeIn(fall_arrow), run_time=0.35)

        # --- Flow ramp: from sparse recombinations to near-continuous emission ---
        lanes = np.linspace(-0.35, 0.35, 7)  # vertical lanes around the active region

        def recomb_event(dy, photon_color):
            e = electron(e_start + UP * dy)
            h = hole(h_start + UP * dy)
            return Succession(
                AnimationGroup(FadeIn(e, run_time=0.10), FadeIn(h, run_time=0.10)),
                AnimationGroup(
                    e.animate.move_to(recomb_point + UP * (0.0 + dy)).set_rate_func(rate_functions.ease_in_out_cubic),
                    h.animate.move_to(recomb_point + DOWN * (0.0 - dy)).set_rate_func(rate_functions.ease_in_out_cubic),
                    run_time=0.55,
                ),
                AnimationGroup(
                    FadeOut(e, run_time=0.10),
                    FadeOut(h, run_time=0.10),
                    photon_pulse(photon_color, recomb_point, max_scale=15.8),
                    run_time=0.70,
                ),
            )

        def run_stream(photon_color, n_events, lag_ratio, run_time):
            events = [recomb_event(float(np.random.choice(lanes)), photon_color) for _ in range(n_events)]
            self.play(LaggedStart(*events, lag_ratio=lag_ratio), run_time=run_time)

        # Phase 1: a few separated events (viewer understands the mechanism)
        for dy in [-0.25, 0.0, 0.25]:
            self.play(recomb_event(dy, COL_PHOTON_RED), run_time=1.10)

        # Phase 2: ramp up (more frequent)
        run_stream(COL_PHOTON_RED, n_events=10, lag_ratio=0.22, run_time=3.2)

        # Phase 3: near-continuous stream (overlapping events -> continuous red light)
        run_stream(COL_PHOTON_RED, n_events=50, lag_ratio=0.08, run_time=8.0)

        self.wait(0.6)

        # -------------------------------------------------
        # Bandgap → photon color: animate the gap changing (red → green → blue)
        # -------------------------------------------------
    

        def set_gap(Ec_level, Ev_level, color, caption):
            # Build new CURVED bands (same bending profile)
            Ec_new = make_band(Ec_level, "#ffe680")
            Ev_new = make_band(Ev_level, "#ffcc66")

            # Update pointers to match the new curved levels at x=4.6
            Ec_ptr_new = Arrow(Ec_lbl.get_left(), axes.c2p(4.6, band_y(4.6, Ec_level)), buff=0.08, stroke_width=3, color="#ffe680")
            Ev_ptr_new = Arrow(Ev_lbl.get_left(), axes.c2p(4.6, band_y(4.6, Ev_level)), buff=0.08, stroke_width=3, color="#ffcc66")

            cap = boxed_label(caption, color, font_size=22)
            cap.next_to(axes, DOWN, buff=0.4)

            return Ec_new, Ev_new, Ec_ptr_new, Ev_ptr_new, cap

        # Start from small gap (red)
        cap_now = None

        gap_steps = [
            (1.4, -1.4, COL_PHOTON_RED,   "Small gap → red"),
            (1.8, -1.8, COL_PHOTON_GREEN, "Medium gap → green"),
            (2.2, -2.2, COL_PHOTON_BLUE,  "Large gap → blue"),
        ]

        for Ec_level, Ev_level, c, txt in gap_steps:
            Ec_new, Ev_new, Ec_ptr_new, Ev_ptr_new, cap = set_gap(Ec_level, Ev_level, c, txt)

            anims = [
                Transform(Ec, Ec_new),
                Transform(Ev, Ev_new),
                Transform(Ec_ptr, Ec_ptr_new),
                Transform(Ev_ptr, Ev_ptr_new),
            ]
            if cap_now is None:
                anims.append(FadeIn(cap))
            else:
                anims.append(Transform(cap_now, cap))

            self.play(*anims, run_time=1.2, rate_func=rate_functions.ease_in_out_cubic)
            cap_now = cap

            # Continuous emission at the new bandgap (short stream)
            run_stream(c, n_events=14, lag_ratio=0.08, run_time=2.4)
            self.wait(0.25)

        self.wait(1.0)

        self.wait(3)