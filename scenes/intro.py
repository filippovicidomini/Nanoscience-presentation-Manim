from manim import *

class IntroScene(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # -----------------------------
        # TESTI
        # -----------------------------
        title = Text(
            "From Silicon to Light",
            font_size=72,
            color=WHITE
        )

        subtitle = Text(
            "PN Junctions and Light Emitting Diodes",
            font_size=36,
            color=BLUE_B
        )

        university = Text(
            "University of Venice Ca' Foscari",
            font_size=28,
            color=GREY_B
        )

        course = Text(
            "Nanotechnology & Nanomaterials",
            font_size=26,
            color=GREY_B
        )
        author = Text("Filippo Vicidomini", font_size=24, color=GREY_C)

        # -----------------------------
        # POSIZIONAMENTO
        # -----------------------------
        title.to_edge(UP, buff=1.6)
        subtitle.next_to(title, DOWN, buff=0.4)

        line = Line(
            LEFT * 4.5,
            RIGHT * 4.5,
            stroke_width=2,
            color=GREY_D
        )
        line.next_to(subtitle, DOWN, buff=0.6)

        university.next_to(line, DOWN, buff=0.6)
        course.next_to(university, DOWN, buff=0.25)
        author.next_to(course, DOWN, buff=0.25)

        # -----------------------------
        # ANIMAZIONI
        # -----------------------------
        self.play(
            FadeIn(title),
            run_time=3,
            lag_ratio=0.9
        )

        self.play(
            FadeIn(subtitle, shift=UP * 0.3),
            run_time=2
        )

        self.play(
            Create(line),
            run_time=1.5
        )

        self.play(
            FadeIn(university, shift=UP * 0.2),
            FadeIn(course, shift=UP * 0.2),
            run_time=2
        )
        self.play(
            FadeIn(author),
            run_time=2,
            lag_ratio=0.4
        )

        self.wait(1.8)

        # -----------------------------
        # USCITA SCENA
        # -----------------------------
        self.play(
            FadeOut(VGroup(
                title, subtitle, line, university, course, author
            )),
            run_time=0.8
        )