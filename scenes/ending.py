from manim import *

class OutroThanks(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # ---------- Title ----------
        title = Text("Thank you for watching!", font_size=64, color=WHITE)
        title.to_edge(UP, buff=1.2)

        line = Line(LEFT * 5.0, RIGHT * 5.0, stroke_width=2, color=GREY_D)
        line.next_to(title, DOWN, buff=0.45)

        # ---------- GitHub ----------
        github_label = Text("Code on GitHub:", font_size=28, color=GREY_B)
        github_link = Text(
            "github.com/filippovicidomini/Nanoscience-presentation-Manim",
            font_size=26,
            color=BLUE_B
        )
        github_group = VGroup(github_label, github_link).arrange(DOWN, buff=0.18)
        github_group.next_to(line, DOWN, buff=0.7)

        # ---------- References ----------
        ref_title = Text("References", font_size=34, color=WHITE)

        refs = VGroup(
            Text("• IUPAC Gold Book — Light-emitting diode (LT07414), doi:10.1351/goldbook.LT07414", font_size=20, color=GREY_B),
            Text("• Wikipedia — Doping (semiconductor) (rev. 2025-11-22)", font_size=20, color=GREY_B),
            Text("• Wikipedia — Intrinsic semiconductor (rev. 2025-05-27)", font_size=20, color=GREY_B),
            Text("• Jha & Kumar (2024) — The Role of Doping in Modifying Semiconductor Properties, doi:10.48175/IJARSCT-19475", font_size=20, color=GREY_B),
            Text("• Electronics-Tutorials — PN Junction Diode and diode characteristics (2013)", font_size=20, color=GREY_B),
            Text("• Wayback Machine PDF — 100 years of optoelectronics (archived 2012)", font_size=20, color=GREY_B),
            Text("• WaferPro — What is doping in semiconductors (2024)", font_size=20, color=GREY_B),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22)

        ref_block = VGroup(ref_title, refs).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
        ref_block.next_to(github_group, DOWN, buff=0.85).to_edge(LEFT, buff=1.0)

        # Slightly scale down if your frame is crowded
        ref_block.scale(0.92)

        # ---------- Animations ----------
        self.play(FadeIn(title, shift=UP * 0.3), run_time=3.0)
        self.play(Create(line), run_time=0.5)
        self.play(FadeIn(github_group, shift=UP * 0.2), run_time=2)
        self.play(FadeIn(ref_title, shift=UP * 0.15), run_time=2)
        self.play(LaggedStart(*[FadeIn(r, shift=UP * 0.05) for r in refs], lag_ratio=0.08), run_time=2)

        self.wait(10.0)

        # Outro fade
        self.play(FadeOut(VGroup(title, line, github_group, ref_block)), run_time=2)