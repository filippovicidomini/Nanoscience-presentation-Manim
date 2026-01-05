from manim import *
import numpy as np

# =========================================================
# SCENE: DIODE I–V CURVE (macroscopic signature)
# Duration target: ~30–40s (tweak run_time values if needed)
# =========================================================

class DiodeIVCurve(ZoomedScene):
    def __init__(self, **kwargs):
        super().__init__(
            zoom_factor=0.38,
            zoomed_display_height=2.9,
            zoomed_display_width=4.6,
            image_frame_stroke_width=3,
            zoomed_camera_config={"default_frame_stroke_width": 2},
            **kwargs,
        )

    def construct(self):
        self.camera.background_color = "#0b0f14"

        # --- Palette (coerente con scene precedenti) ---
        COL_AXIS = "#E6E6E6"
        COL_BOX  = "#0b0f14"
        COL_FWD  = YELLOW_B
        COL_REV  = BLUE_C
        COL_CURVE = "#ffe680"  # warm yellow like bands
        COL_KNEE  = RED_C

        # -----------------------------
        # Helpers
        # -----------------------------
        def boxed_label(text, color=COL_AXIS, font_size=26, padding=0.18, opacity=0.85):
            t = Text(text, font_size=font_size, color=color)
            bg = RoundedRectangle(
                corner_radius=0.14,
                width=t.width + 2 * padding,
                height=t.height + 2 * padding,
                stroke_width=0,
                fill_color=COL_BOX,
                fill_opacity=opacity,
            )
            g = VGroup(bg, t)
            g.text = t
            g.bg = bg
            return g

        # -----------------------------
        # Axes
        # -----------------------------
        # Voltage in V, Current in arbitrary units (or mA).
        axes = Axes(
            x_range=[-3, 3, 1],
            y_range=[-0.6, 3.0, 0.5],
            x_length=10.8,
            y_length=6.2,
            tips=False,
            axis_config={"stroke_color": COL_AXIS, "stroke_width": 2},
        ).to_edge(UP, buff=0.75)

        x_label = boxed_label("V (voltage)", COL_AXIS, font_size=22, opacity=0.65)
        y_label = boxed_label("I (current)", COL_AXIS, font_size=22, opacity=0.65)
        x_label.next_to(axes, DOWN, buff=0.28)
        y_label.next_to(axes, LEFT, buff=0.32).shift(UP * 2.0)

        title = boxed_label("Diode I–V characteristic", COL_AXIS, font_size=28, opacity=0.75)
        title.next_to(axes, UP, buff=0.18)

        self.play(Create(axes), FadeIn(x_label), FadeIn(y_label), FadeIn(title), run_time=1.6)

        # -----------------------------
        # Diode-like I–V model (simple, visually clear)
        # -----------------------------
        # Knee voltage around ~0.7 V (silicon-like). Change if you want LED-ish.
        V_knee = 0.7

        def I_of_V(V):
            """Smooth diode-like I–V (no jump at V=0).
            - Forward: Shockley-like exponential
            - Reverse: tiny leakage that saturates (almost 0)
            - Optional: soft reverse breakdown beyond V_br (conceptual)
            """
            Is = 0.03          # scale factor (visual)
            nVt = 0.26         # sets the exponential steepness (visual)

            # Base Shockley (continuous at 0): I = Is*(exp(V/nVt) - 1)
            I = Is * (np.exp(V / nVt) - 1.0)

            # In reverse, clamp to a small leakage (keeps it near 0 and readable)
            I_leak = -0.045
            if I < I_leak:
                I = I_leak

            # Soft reverse breakdown ramp (conceptual): more negative current after V_br
            V_br = -2.4
            if V < V_br:
                # smooth extra current that grows after breakdown
                I -= 0.22 * (1.0 - np.exp((V - V_br) / 0.25))

            return I

        xs = np.linspace(-3, 3, 800)
        ys = np.array([I_of_V(v) for v in xs])

        # clamp to axes y_range max to avoid going off-frame
        ys = np.clip(ys, axes.y_range[0], axes.y_range[1])

        curve = axes.plot_line_graph(
            x_values=xs,
            y_values=ys,
            line_color=COL_CURVE,
            add_vertex_dots=False,
            stroke_width=5,
        )

        # Animate curve drawing
        self.play(Create(curve), run_time=2.2, rate_func=rate_functions.ease_in_out_cubic)

        # -----------------------------
        # Highlight regions: Reverse & Forward
        # -----------------------------
        # Reverse region shading (left half)
        rev_rect = Rectangle(
            width=axes.x_axis.get_length() / 2,
            height=axes.y_axis.get_length(),
            stroke_width=0,
            fill_color=COL_REV,
            fill_opacity=0.08,
        )
        rev_rect.align_to(axes.c2p(-3, 0), LEFT).align_to(axes, DOWN)

        # Forward region shading (right half)
        fwd_rect = Rectangle(
            width=axes.x_axis.get_length() / 2,
            height=axes.y_axis.get_length(),
            stroke_width=0,
            fill_color=COL_FWD,
            fill_opacity=0.06,
        )
        fwd_rect.align_to(axes.c2p(0, 0), LEFT).align_to(axes, DOWN)

        rev_tag = boxed_label("Reverse bias\n(I ≈ 0)", COL_REV, font_size=22, opacity=0.75)
        fwd_tag = boxed_label("Forward bias\n(I rises fast)", COL_FWD, font_size=22, opacity=0.75)

        rev_tag.move_to(axes.c2p(-2.05, 2.25))
        fwd_tag.move_to(axes.c2p( 2.05, 2.25))

        self.play(FadeIn(rev_rect), FadeIn(fwd_rect), run_time=0.8)
        self.play(FadeIn(rev_tag), FadeIn(fwd_tag), run_time=0.8)
        self.wait(1.0)

        # -----------------------------
        # Clean zoom (magnifier) around 0 V using ZoomedScene
        # -----------------------------
        # Define the zoom window in data coordinates
        x0, x1 = -0.55, 0.95
        y0, y1 = -0.06, 0.55

        # Zoom frame size in scene units (match axes scaling)
        frame_w = axes.c2p(x1, 0)[0] - axes.c2p(x0, 0)[0]
        frame_h = axes.c2p(0, y1)[1] - axes.c2p(0, y0)[1]

        zoom_frame = self.zoomed_camera.frame
        zoom_frame.set_width(frame_w)
        zoom_frame.set_height(frame_h)
        zoom_frame.set_stroke(COL_AXIS)
        zoom_frame.set_stroke_width(3)
        zoom_frame.set_opacity(0.55)
        zoom_frame.move_to(axes.get_origin())

        # Place the zoom display in the bottom-right
        self.zoomed_display.to_corner(DR, buff=0.50)
        self.zoomed_display.set_stroke(COL_AXIS)
        self.zoomed_display.set_stroke_width(2)
        self.zoomed_display.set_opacity(0.35)

        zoom_title = boxed_label("Zoom around 0 V", COL_AXIS, font_size=20, opacity=0.70)
        zoom_title.next_to(self.zoomed_display, UP, buff=0.12)

        # A tiny origin marker to make it obvious what we are zooming
        origin_marker = Dot(axes.get_origin(), radius=0.045, color=COL_AXIS)

        self.add(origin_marker)
        self.play(FadeIn(zoom_frame), run_time=0.5)

        # Activate zooming and pop out the display
        self.activate_zooming()
        self.play(
            FadeIn(self.zoomed_display),
            FadeIn(zoom_title),
            run_time=0.8,
        )

        # Optional: subtle connector line (very light)
        connector = DashedLine(
            zoom_frame.get_corner(DR),
            self.zoomed_display.get_corner(UL),
            dash_length=0.10,
            stroke_width=2,
            color=COL_AXIS,
        ).set_opacity(0.18)
        self.play(Create(connector), run_time=0.5)

        self.wait(1.1)

        # -----------------------------
        # Knee voltage annotation
        # -----------------------------
        knee_point = axes.c2p(V_knee, I_of_V(V_knee))
        knee_dot = Dot(knee_point, radius=0.07, color=COL_KNEE)

        knee_line = DashedLine(
            axes.c2p(V_knee, 0),
            knee_point,
            dash_length=0.10,
            stroke_width=2,
            color=COL_KNEE,
        ).set_opacity(0.9)

        knee_lbl = boxed_label("knee voltage", COL_KNEE, font_size=22, opacity=0.8)
        knee_lbl.next_to(knee_point, UP + RIGHT, buff=0.15)

        knee_arrow = Arrow(
            knee_lbl.get_left(),
            knee_point,
            buff=0.08,
            stroke_width=3,
            color=COL_KNEE,
        )

        self.play(FadeIn(knee_dot), Create(knee_line), run_time=0.7)
        self.play(FadeIn(knee_lbl), GrowArrow(knee_arrow), run_time=0.8)
        self.wait(0.6)

        # -----------------------------
        # Reverse breakdown (conceptual marker)
        # -----------------------------
        V_br = -2.4
        br_line = DashedLine(
            axes.c2p(V_br, axes.y_range[0]),
            axes.c2p(V_br, axes.y_range[1]),
            dash_length=0.12,
            stroke_width=2,
            color=COL_REV,
        ).set_opacity(0.55)

        br_lbl = boxed_label("reverse breakdown", COL_REV, font_size=22, opacity=0.75)
        br_lbl.next_to(axes.c2p(V_br, 2.7), UP, buff=0.12)

        self.play(Create(br_line), FadeIn(br_lbl), run_time=1.0)
        self.wait(1.0)

        # -----------------------------
        # Animated sweep ("measuring" the I–V curve) — placed AFTER annotations
        # -----------------------------
        v_tracker = ValueTracker(-3.0)

        def clamp_I(v):
            return float(np.clip(I_of_V(v), axes.y_range[0], axes.y_range[1]))

        sweep_dot = always_redraw(
            lambda: Dot(
                axes.c2p(v_tracker.get_value(), clamp_I(v_tracker.get_value())),
                radius=0.055,
                color=COL_AXIS,
            )
        )

        trail = TracedPath(
            sweep_dot.get_center,
            dissipating_time=1.6,
            stroke_width=4,
            stroke_color=COL_CURVE,
            stroke_opacity=0.55,
        )

        sweep_tag = boxed_label("Voltage sweep", COL_AXIS, font_size=22, opacity=0.70)
        sweep_tag.move_to(axes.c2p(-2.35, -0.45))

        self.add(trail, sweep_dot)
        self.play(FadeIn(sweep_tag), run_time=0.5)

        # Sweep reverse -> forward (slower and more readable)
        self.play(v_tracker.animate.set_value(3.0), run_time=6.0, rate_func=linear)
        self.play(FadeOut(sweep_tag), run_time=0.3)

        # Keep the dot/trail briefly, then remove for a clean final key message
        self.wait(0.6)
        self.remove(trail, sweep_dot)

        # -----------------------------
        # Key message (short, powerful)
        # -----------------------------
        key = boxed_label("I–V curve = macroscopic signature\nof the potential barrier", COL_AXIS, font_size=24, opacity=0.75)
        key.to_edge(DOWN, buff=0.6)

        self.play(FadeIn(key), run_time=0.9)
        self.wait(3.6)

        # Clean up updaters if any (none here, but safe)
        for m in self.mobjects:
            m.clear_updaters()