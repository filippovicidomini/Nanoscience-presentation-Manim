from manim import *

class MeissnerEffect(Scene):
    def construct(self):
        # Campo applicato: frecce verticali
        arrows = VGroup(*[
            Arrow(start=UP*3, end=DOWN*3, buff=0, stroke_width=2)
            .shift(LEFT*5 + RIGHT*i*1.0) for i in range(11)
        ])
        self.play(LaggedStart(*[Create(a) for a in arrows], lag_ratio=0.05))
        
        # Campione: cilindro visto di fronte -> disco
        sc = Circle(radius=1.2, color=BLUE_E, fill_opacity=0.15, stroke_width=4)
        self.play(FadeIn(sc))
        self.wait(0.5)

        # Fase "normale": le frecce lo attraversano (nessun cambiamento)
        self.wait(0.5)

        # Cool down sotto Tc -> Meissner: devia le linee ai lati
        def bend_arrow(a: Arrow):
            # se l'asse X della freccia cade "dentro" il disco, piegala ai lati
            x = a.get_center()[0]
            if abs(x) < 1.2:
                # ricrea l'arrow deviandola a lato del disco
                sign = -1 if x < 0 else 1
                new = CubicBezier(
                    a.get_start(),
                    a.get_start() + RIGHT*sign*1.5 + DOWN*1.0,
                    a.get_end()   + RIGHT*sign*1.5 + UP*1.0,
                    a.get_end()
                ).set_stroke(width=2)
                return new
            else:
                return a.copy()

        bent = VGroup(*[bend_arrow(a) for a in arrows])
        self.play(
            *[Transform(arrows[i], bent[i]) for i in range(len(arrows))],
            sc.animate.set_fill(BLUE, opacity=0.25),
            run_time=2
        )
        # Etichette (senza LaTeX)
        title = Text("Meissner effect: perfect diamagnet", weight=BOLD, font_size=36)
        title.to_edge(UP)
        self.play(Write(title))
        self.wait(1.5)
        
from manim import *
import numpy as np

class LondonPenetration(Scene):
    def construct(self):
        # Lastra SC: rettangolo, superficie a x=0, materiale per x>0
        plate = Rectangle(height=4, width=6).set_stroke(WHITE, 2)
        plate.shift(RIGHT*1.0)
        self.play(Create(plate))

        # Colormap "manuale": rettangoli verticali con opacità che decresce ~ exp(-x/lambda)
        lam = 1.0
        bars = VGroup()
        for i, x in enumerate(np.linspace(0.1, 2.8, 16)):
            h = 3.6
            bar = Rectangle(height=h, width=0.15)
            bar.set_fill(color=BLUE, opacity=float(np.exp(-x/lam)))
            bar.set_stroke(opacity=0)
            bar.move_to(plate.get_left() + RIGHT*(x+0.2))
            bars.add(bar)
        self.play(FadeIn(bars, lag_ratio=0.05))

        # Frecce del campo che "entrano" di poco
        arrows = VGroup()
        for y in np.linspace(-1.5, 1.5, 7):
            a = Arrow(start=LEFT*4 + UP*y, end=plate.get_left()+LEFT*0.2+UP*y, buff=0)
            a.set_stroke(width=2)
            arrows.add(a)
        self.play(LaggedStart(*[Create(a) for a in arrows], lag_ratio=0.05))

        # Testo esplicativo (niente LaTeX)
        t1 = Text("Magnetic field decays inside over λ (London depth)", font_size=32)
        t1.to_edge(UP)
        t2 = Text("B(x) ≈ B0 · exp(-x/λ)", font_size=30, color=YELLOW).next_to(t1, DOWN)
        self.play(Write(t1), FadeIn(t2))
        self.wait(1.0)

        # Slider λ: diminuisci λ e mostra una penetrazione più corta
        lam_tracker = ValueTracker(1.0)

        def update_bars(group):
            lam = lam_tracker.get_value()
            for i, bar in enumerate(group):
                x = bar.get_center()[0] - plate.get_left()[0]
                bar.set_fill(opacity=float(np.exp(-max(x,0)/max(lam, 1e-3))))
            return group

        bars.add_updater(update_bars)
        self.play(lam_tracker.animate.set_value(0.5), run_time=2)
        self.play(lam_tracker.animate.set_value(1.5), run_time=2)
        bars.remove_updater(update_bars)
        self.wait(0.5)