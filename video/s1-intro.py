from manim import *

# Scene 1 — Intro
# Manim Community Edition ≥ 0.18
# Render:  manim -pqh s1_intro.py S1Intro

class S1Intro(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f17"

        title = Text(
            "From silicon to LED",
            font_size=64,
            weight=BOLD,
            slant=ITALIC,
            color="#E6EDF3",
        ).to_edge(UP, buff=0.6)

        # Piccolo motivo grafico: un “atomo” stilizzato con 4 elettroni di valenza
        atom_group = self.make_atom_motif(level_counts=[2,8,4], base_radius=0.7, radius_step=0.55)
        atom_group.scale(1.0).to_edge(RIGHT, buff=1.0).shift(0.0*UP)
        atom_group.set_opacity(0.9)

        # Punti elenco agenda
        agenda_items = [
            "From orbitals to shell model",
            "Covalent bonding in silicon bulk",
            "Doping: phosphorus (n-type) and boron (p-type)",
            "p–n junction: depletion and electric fields",
            "Quick LED simulation",
        ]
        bullets = VGroup(*[
            self.bullet_line(txt, i) for i, txt in enumerate(agenda_items)
        ])
        bullets.arrange(DOWN, aligned_edge=LEFT, buff=0.25)
        bullets.to_edge(LEFT, buff=1.0)

        
        # Intro animations
        self.play(
            LaggedStart(
                Write(title),
                #FadeIn(meta, shift=0.15*DOWN),
                FadeIn(atom_group, shift=0.2*UP, scale=0.95),
                run_time=2.0,
                lag_ratio=0.2,
            )
        )
        
        self.play(LaggedStart(*[Write(b) for b in bullets], lag_ratio=0.15, run_time=3))

        # Rotazione morbida: ruota i gruppi di elettroni in senso alternato
        spin = [Rotate(g, angle=(2*PI if i%2==0 else -2*PI)) for i, g in enumerate(atom_group.electron_groups)]
        self.play(*spin, run_time=6, rate_func=linear)

        # Outro della scena (senza tagline)
        self.wait(0.8)
        self.play(*[FadeOut(m) for m in [title, bullets]], FadeOut(atom_group, scale=0.9))
        self.wait(0.2)

    # ——————————————————————————————————————————————————————————————
    def bullet_line(self, text: str, idx: int) -> VGroup:
        dot = Dot(radius=0.05, color="#6EA8FE")
        line_text = Text(text, font_size=30, color="#D7E2EE")
        group = VGroup(dot, line_text)
        group.arrange(RIGHT, buff=0.3)
        # leggera sfalsatura per un look "timeline"
        group.shift(0.04*DOWN*idx)
        return group

    def make_atom_motif(self, level_counts=[2,8,4], base_radius=0.7, radius_step=0.55) -> VGroup:
        """Motivo atomico del silicio con tutti gli elettroni: 2, 8, 4 su tre gusci.
        level_counts controlla la popolazione per guscio; base_radius e radius_step il raggio dei gusci.
        """
        nucleus = Circle(radius=0.18, color="#66E0A3", fill_opacity=1.0).set_stroke(width=0)
        shell_circles = VGroup()
        electron_groups = []  # lista di VGroup, uno per guscio

        for i, cnt in enumerate(level_counts):
            r = base_radius + i * radius_step
            shell = Circle(radius=r, color="#2A3B55", stroke_width=2)
            shell_circles.add(shell)

            eg = VGroup()
            for k in range(cnt):
                angle = k * TAU / cnt
                pos = r * np.array([np.cos(angle), np.sin(angle), 0])
                # Puntini leggermente più piccoli quando il guscio è affollato
                dot_r = 0.06 if cnt >= 6 else 0.07
                e = Dot(point=pos, radius=dot_r, color="#6EA8FE")
                eg.add(e)
            electron_groups.append(eg)

        all_electrons = VGroup(*electron_groups)
        atom = VGroup(shell_circles, nucleus, all_electrons)
        # attribuiti utili per animazioni successive
        atom.shell_circles = shell_circles
        atom.electron_groups = electron_groups
        return atom

    def start_spin(self, atom: VGroup, omega: float = 0.9):
        """Attiva uno spin continuo sugli elettroni con velocità alternate per guscio."""
        for i, g in enumerate(atom.electron_groups):
            sign = 1 if (i % 2 == 0) else -1
            w = omega * (1 - 0.12 * i)
            g.add_updater(lambda m, dt, s=sign, ww=w, a=atom: m.rotate(s*ww*dt, about_point=a.get_center()))

    def stop_spin(self, atom: VGroup):
        """Ferma lo spin rimuovendo gli updaters."""
        for g in atom.electron_groups:
            g.clear_updaters()
