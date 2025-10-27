from manim import *
import numpy as np

# ==============================================================
#  Quantum Confinement Showcase — NO LaTeX VERSION
#  (uses Text/DecimalNumber only; safe on machines without TeX)
# ==============================================================

class QuantumConfinementNoTex(Scene):
    def construct(self):
        # Parameters
        L0, Lmin = 8.0, 3.0
        n_levels = 3
        m_scale = 1.0
        L = ValueTracker(L0)

        # Left panel axes
        left_panel = Axes(
            x_range=[0, L0, L0/8],
            y_range=[-1.2, 1.2, 0.4],
            x_length=7,
            y_length=3.2,
            tips=False,
        ).to_edge(LEFT, buff=0.7).shift(UP*0.5)

        title = VGroup(
            Text("Quantum Confinement", weight=BOLD),
            Text("Particle in an Infinite Well", slant=ITALIC, font_size=28)
        ).arrange(DOWN, buff=0.1).to_edge(UP)

        # Walls
        left_wall = always_redraw(lambda: left_panel.get_vertical_line(left_panel.c2p(0,0)).set_color(BLUE))
        right_wall = always_redraw(lambda: left_panel.get_vertical_line(left_panel.c2p(L.get_value(),0)).set_color(BLUE))

        # Brace and label for L
        L_brace = always_redraw(lambda: BraceBetweenPoints(
            left_panel.c2p(0, 1.15), left_panel.c2p(L.get_value(), 1.15), buff=0.06
        ))
        L_value = always_redraw(lambda: DecimalNumber(L.get_value(), num_decimal_places=1))
        L_text = Text("L = ")
        L_label = always_redraw(lambda: VGroup(L_text, L_value).arrange(RIGHT, buff=0.08).next_to(L_brace, UP, buff=0.08))

        # Wavefunctions (no LaTeX labels; use unicode ψ)
        def psi_n(n):
            def f(x):
                Lval = L.get_value()
                if 0 <= x <= Lval:
                    return np.sqrt(2/Lval) * np.sin(n*np.pi*x/Lval)
                return 0.0
            return f

        psi_colors = [YELLOW, ORANGE, TEAL]
        psi_graphs = [
            always_redraw(lambda n=n: left_panel.plot(
                lambda x: m_scale*psi_n(n)(x),
                x_range=[0, L.get_value()],
                use_smoothing=True,
                color=psi_colors[n-1]
            )) for n in range(1, 3)
        ]

        psi1_label = always_redraw(lambda: Text("ψ1").scale(0.6).next_to(psi_graphs[0], UP, buff=0.1))
        psi2_label = always_redraw(lambda: Text("ψ2").scale(0.6).next_to(psi_graphs[1], DOWN, buff=0.1))

        # Probability density area for n=1
        def prob_density(x):
            Lval = L.get_value()
            if 0 <= x <= Lval:
                return (2/Lval)*np.sin(np.pi*x/Lval)**2
            return 0.0
        pd_graph = always_redraw(lambda: left_panel.plot(
            lambda x: 0.6*prob_density(x), x_range=[0, L.get_value()], color=YELLOW
        ))
        pd_area = always_redraw(lambda: left_panel.get_area(pd_graph, x_range=[0, L.get_value()], opacity=0.2, color=YELLOW))

        # Energy ladder (right panel)
        E_axes = Axes(
            x_range=[0, 1, 1], y_range=[0, 6, 1], x_length=1.2, y_length=4.2, tips=False
        ).to_edge(RIGHT, buff=1.3).shift(UP*0.5)
        E_title = Text("Energy Levels").next_to(E_axes, UP)

        def En(n):
            return 0.6 + 2.2*(n*n)/(L.get_value()**2)  # visual mapping

        E_lines = [
            always_redraw(lambda n=n: E_axes.get_horizontal_line(E_axes.c2p(0, En(n))).set_color(WHITE))
            for n in range(1, n_levels+1)
        ]
        E_labels = [
            always_redraw(lambda n=n: Text(f"E{n}").scale(0.6).next_to(E_axes.c2p(0.9, En(n)), RIGHT, buff=0.06))
            for n in range(1, n_levels+1)
        ]

        # Info box (plain text)
        info_lines = VGroup(
            Text("En ∝ n² / L²", font_size=28),
            Text("k = nπ/L   →   λ = 2L/n", font_size=28)
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        info_panel = RoundedRectangle(corner_radius=0.2, height=1.6, width=5.4).set_opacity(0.15)
        info_group = VGroup(info_panel, info_lines)
        info_group.to_corner(DOWN+RIGHT, buff=0.6)
        info_lines.move_to(info_group)

        # Build scene
        self.play(FadeIn(title, shift=DOWN))
        self.play(Create(left_panel), FadeIn(left_wall), FadeIn(right_wall))
        self.play(FadeIn(L_brace), FadeIn(L_label))
        self.play(Create(pd_graph), FadeIn(pd_area))
        self.play(Create(psi_graphs[0]), FadeIn(psi1_label))
        self.play(Create(psi_graphs[1]), FadeIn(psi2_label))
        self.play(Create(E_axes), FadeIn(E_title))
        self.play(*[FadeIn(line) for line in E_lines], *[FadeIn(lbl) for lbl in E_labels])
        self.play(FadeIn(info_group))
        self.wait(0.4)

        # Animate shrinking L
        self.play(L.animate.set_value(Lmin), run_time=4, rate_func=smooth)
        self.wait(0.3)

        # Photon emission (color encodes ΔE)
        def photon_color():
            Lval = L.get_value()
            dE = (4 - 1)/(Lval**2)
            dE_min = (4 - 1)/(L0**2)
            dE_max = (4 - 1)/(Lmin**2)
            alpha = (dE - dE_min)/(dE_max - dE_min + 1e-6)
            return interpolate_color(RED, BLUE, alpha)

        photon = always_redraw(lambda: Dot(color=photon_color()).scale(0.8).move_to(E_axes.c2p(0.6, En(2)+0.1)))
        self.play(FadeIn(photon, scale=0.5), run_time=0.5)
        self.play(photon.animate.move_to(E_axes.c2p(0.6, En(1)+0.08)), run_time=0.7)
        self.play(photon.animate.shift(LEFT*5 + DOWN*1.5).set_opacity(0.6), run_time=1.4)
        self.play(FadeOut(photon))

        # Breathing
        for new_L in [3.6, 3.2, 3.8, 3.4]:
            self.play(L.animate.set_value(new_L), run_time=0.7)
        self.wait(0.5)

        # Outro
        highlight = SurroundingRectangle(psi_graphs[0], buff=0.15, color=YELLOW).set_opacity(0)
        self.play(highlight.animate.set_opacity(0.5), run_time=0.5)
        self.play(FadeOut(highlight), FadeOut(info_group))

        credit = Text("ΔE grows as L shrinks → bluer photons")
        credit.next_to(E_axes, DOWN, buff=0.3)
        self.play(Write(credit))
        self.wait(0.8)

        self.play(*map(FadeOut, [
            credit, title, left_panel, left_wall, right_wall, L_brace, L_label,
            pd_graph, pd_area, psi_graphs[0], psi_graphs[1], psi1_label, psi2_label,
            E_axes, E_title, *E_lines, *E_labels
        ]))
        self.wait(0.2)
