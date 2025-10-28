from manim import *
import numpy as np

class QuantumConfinement(Scene):
    def construct(self):
        L = ValueTracker(6.0)  # larghezza dinamica
        axes = Axes(x_range=[0,6,1], y_range=[-1.2,1.2,0.5], x_length=7, y_length=3)
        self.add(axes)

        # pareti della buca
        left_wall = always_redraw(lambda: axes.get_vertical_line(axes.c2p(0,0)))
        right_wall = always_redraw(lambda: axes.get_vertical_line(axes.c2p(L.get_value(),0)))
        box_brace = always_redraw(lambda: BraceBetweenPoints(axes.c2p(0,1.1), axes.c2p(L.get_value(),1.1), buff=0.1))
        L_label = always_redraw(lambda: box_brace.get_text(f"L = {L.get_value():.1f}"))

        self.add(left_wall, right_wall, box_brace, L_label)

        def psi(n):
            def f(x):
                Lval = L.get_value()
                if x<0 or x>Lval: return 0
                return np.sqrt(2/Lval)*np.sin(n*np.pi*x/Lval)
            return f

        psi1_graph = always_redraw(lambda: axes.plot(psi(1), x_range=[0, L.get_value()], use_smoothing=True))
        psi2_graph = always_redraw(lambda: axes.plot(psi(2), x_range=[0, L.get_value()], use_smoothing=True))
        psi1_label = MathTex(r"\psi_1(x)").next_to(psi1_graph, UP)
        psi2_label = MathTex(r"\psi_2(x)").next_to(psi2_graph, DOWN)

        self.play(Create(psi1_graph), FadeIn(psi1_label))
        self.play(Create(psi2_graph), FadeIn(psi2_label))

        # scala dei livelli a destra
        E_axes = Axes(x_range=[0,1,1], y_range=[0,5,1], x_length=1, y_length=4).to_edge(RIGHT)
        self.play(FadeIn(E_axes))
        E_line_1 = always_redraw(lambda: E_axes.get_horizontal_line(E_axes.c2p(0, 1/(L.get_value()**2)+0.5)))
        E_line_2 = always_redraw(lambda: E_axes.get_horizontal_line(E_axes.c2p(0, 4/(L.get_value()**2)+0.5)))
        self.play(FadeIn(E_line_1), FadeIn(E_line_2))

        # “stringi” la scatola
        self.play(L.animate.set_value(3.0), run_time=3)

        # “emissione” E2 -> E1 (solo effetto visivo)
        photon = Dot(color=YELLOW).move_to(E_axes.c2p(0.5, (4-1)/(L.get_value()**2)+0.5))
        self.play(FadeIn(photon), photon.animate.shift(LEFT*4+DOWN*2), run_time=2)
        self.wait()