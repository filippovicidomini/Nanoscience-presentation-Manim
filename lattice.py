from manim import *
import numpy as np

# --- SCENA 4: RETICOLO DI SILICIO CON DOPANTE BORO (schematico 2D) ---
class SiliconLatticeBoron(Scene):
    def construct(self):
        # Parametri reticolo (proiezione 2D stilizzata, non struttura cristallografica reale)
        rows, cols = 6, 9
        a = 0.9  # passo (distanza tra nodi)
        origin = ORIGIN

        title = VGroup(
            Text("Silicon lattice (2D schematic)", font_size=34, weight=BOLD),
            Text("Insert one Boron atom → acceptor / hole", font_size=26, slant=ITALIC)
        ).arrange(DOWN, buff=0.15).to_edge(UP)

        # Genera posizioni su griglia romboidale (effetto diamond-like)
        def pos(i, j):
            # offset alternato per creare pattern a rombo
            offset = (j % 2) * (a * 0.5)
            return origin + np.array([(j - cols/2) * a, (i - rows/2) * a + offset, 0])

        # Atomi di silicio
        atoms = {}
        atom_group = VGroup()
        for i in range(rows):
            for j in range(cols):
                d = Dot(pos(i, j), radius=0.08, color=GREY_B, fill_opacity=1.0)
                atoms[(i, j)] = d
                atom_group.add(d)

        # Legami (collega a destra e in alto-destra per una maglia pulita)
        bonds = VGroup()
        for i in range(rows):
            for j in range(cols):
                p = pos(i, j)
                # destra
                if j + 1 < cols:
                    q = pos(i, j + 1)
                    bonds.add(Line(p, q, stroke_width=2, color=GREY_D, z_index=-1))
                # alto-destra
                if i - 1 >= 0 and j + 1 < cols:
                    q2 = pos(i - 1, j + 1)
                    bonds.add(Line(p, q2, stroke_width=2, color=GREY_D, z_index=-1))

        # Costruisci scena base
        self.play(FadeIn(title, shift=DOWN))
        self.play(LaggedStartMap(Create, bonds, lag_ratio=0.01, run_time=1.2))
        self.play(LaggedStart(*[FadeIn(a, scale=0.7) for a in atom_group], lag_ratio=0.01, run_time=1.0))
        self.wait(0.2)

        # Scegli un sito centrale per sostituzione con Boro
        ci, cj = rows // 2, cols // 2
        si_dot = atoms[(ci, cj)]
        b_color = PURPLE_A

        # Evidenzia l'intorno del sito (bagliore + cerchio)
        glow = Circle(radius=0.22, color=YELLOW).move_to(si_dot.get_center()).set_stroke(width=4)
        self.play(Create(glow), Flash(si_dot, flash_radius=0.3, color=YELLOW))

        # Sostituzione: Si (grigio) -> B (viola/rosa): atomo accettore (p-type)
        b_dot = Dot(si_dot.get_center(), radius=0.1, color=b_color)
        b_label = Text("B", font_size=28, weight=BOLD).move_to(b_dot.get_center() + 0.3*UP)
        self.play(Transform(si_dot, b_dot), FadeIn(b_label, shift=UP*0.2))
        self.wait(0.2)

        # Visualizza una "lacuna" (hole) come cerchietto semi-trasparente con segno +
        hole = Circle(radius=0.10, color=YELLOW, stroke_width=3).move_to(b_dot.get_center() + 0.25*RIGHT)
        plus = Text("+", font_size=22).move_to(hole.get_center())
        hole_group = VGroup(hole, plus).set_opacity(0.9)
        self.play(FadeIn(hole_group, scale=0.6))

        # Anima la delocalizzazione della lacuna su legami vicini (drift breve tra vicini)
        nbr_offsets = [0.25*RIGHT, 0.22*UR, 0.25*UP, 0.22*UL, 0.25*LEFT, 0.22*DL, 0.25*DOWN, 0.22*DR]
        for off in nbr_offsets:
            self.play(hole_group.animate.move_to(b_dot.get_center() + off), run_time=0.25)
        self.play(hole_group.animate.move_to(b_dot.get_center() + 0.25*RIGHT), run_time=0.25)

        # Piccola legenda in basso a destra
        legend = VGroup(
            Dot(radius=0.08, color=GREY_B), Text("Si", font_size=22),
            Dot(radius=0.1, color=b_color), Text("B (acceptor)", font_size=22),
            Circle(radius=0.08, color=YELLOW, stroke_width=3), Text("hole", font_size=22)
        ).arrange_in_grid(rows=3, cols=2, buff=0.15, cell_alignment=LEFT)
        legend_panel = RoundedRectangle(width=4.4, height=2.2, corner_radius=0.15).set_opacity(0.12)
        legend_group = VGroup(legend_panel, legend).to_corner(DOWN+RIGHT, buff=0.4)
        legend.move_to(legend_group)
        self.play(FadeIn(legend_group))

        # Sfumatura finale
        self.wait(0.8)
        self.play(FadeOut(glow))
        self.wait(0.2)
        
        
class SiliconLatticeBoronPro(Scene):
    def construct(self):
        # --------- Parameters ---------
        rows, cols = 8, 12
        a = 0.8             # lattice pitch
        origin = ORIGIN + LEFT*2
        si_color = GREY_B
        bond_color = GREY_D
        b_color = PURPLE_A
        hole_color = YELLOW

        # --------- Title ---------
        title = VGroup(
            Text("Silicon lattice — Boron acceptor", font_size=36, weight=BOLD),
            Text("Substitutional B creates a hole (p-type)", font_size=26, slant=ITALIC)
        ).arrange(DOWN, buff=0.1).to_edge(UP)
        self.play(FadeIn(title, shift=DOWN))

        # --------- Lattice positions (rhombic pattern) ---------
        def pos(i, j):
            offset = (j % 2) * (a * 0.5)
            return origin + np.array([(j - cols/2) * a, (i - rows/2) * a + offset, 0])

        # Silicon atoms
        atoms = {}
        atom_group = VGroup()
        for i in range(rows):
            for j in range(cols):
                d = Dot(pos(i, j), radius=0.07, color=si_color)
                atoms[(i, j)] = d
                atom_group.add(d)

        # Bonds: right and up-right neighbors
        bonds = VGroup()
        for i in range(rows):
            for j in range(cols):
                p = pos(i, j)
                if j + 1 < cols:
                    bonds.add(Line(p, pos(i, j+1), stroke_width=2, color=bond_color, z_index=-1))
                if i - 1 >= 0 and j + 1 < cols:
                    bonds.add(Line(p, pos(i-1, j+1), stroke_width=2, color=bond_color, z_index=-1))

        self.play(LaggedStartMap(Create, bonds, lag_ratio=0.006, run_time=1.2))
        self.play(LaggedStart(*[FadeIn(a, scale=0.7) for a in atom_group], lag_ratio=0.004, run_time=1.0))
        self.wait(0.2)

        # --------- Choose substitution site (central) ---------
        ci, cj = rows//2, cols//2
        si_dot = atoms[(ci, cj)]

        # Glow focus
        glow = Circle(radius=0.22, color=hole_color).set_stroke(width=5).move_to(si_dot)
        pulse = glow.copy().set_opacity(0.5)
        self.play(Create(glow), FadeIn(pulse, scale=1.2))
        self.play(pulse.animate.scale(1.8).set_opacity(0), run_time=0.8)

        # Substitute with Boron
        b_dot = Dot(si_dot.get_center(), radius=0.1, color=b_color)
        b_label = Text("B", font_size=26, weight=BOLD).move_to(b_dot.get_center()+0.28*UP)
        self.play(Transform(si_dot, b_dot), FadeIn(b_label, shift=UP*0.2))

        # --------- Create a hole near B and animate delocalization ---------
        hole = Circle(radius=0.11, color=hole_color, stroke_width=3).move_to(b_dot.get_center()+0.24*RIGHT)
        plus = Text("+", font_size=20).move_to(hole.get_center())
        hole_g = VGroup(hole, plus).set_opacity(0.95)
        self.play(FadeIn(hole_g, scale=0.6))

        # Neighbor positions around B
        nbr_offsets = [0.24*RIGHT, 0.2*UR, 0.24*UP, 0.2*UL, 0.24*LEFT, 0.2*DL, 0.24*DOWN, 0.2*DR]
        for _ in range(2):
            for off in nbr_offsets:
                self.play(hole_g.animate.move_to(b_dot.get_center()+off), run_time=0.18)
        self.play(hole_g.animate.move_to(b_dot.get_center()+0.24*RIGHT), run_time=0.2)

        # --------- Right inset: simple band diagram (no LaTeX) ---------
        inset_origin = RIGHT*4.4 + DOWN*0.3
        width = 3.6
        # Conduction band Ec and valence band Ev
        Ec_y = 1.2
        Ev_y = -1.2
        Ec = Line(inset_origin + LEFT*width/2 + UP*Ec_y,
                  inset_origin + RIGHT*width/2 + UP*Ec_y,
                  color=BLUE_E, stroke_width=4)
        Ev = Line(inset_origin + LEFT*width/2 + UP*Ev_y,
                  inset_origin + RIGHT*width/2 + UP*Ev_y,
                  color=RED_E, stroke_width=4)
        # Acceptor level EA slightly above Ev
        EA_y = Ev_y + 0.35
        EA = DashedLine(inset_origin + LEFT*width/2 + UP*EA_y,
                        inset_origin + RIGHT*width/2 + UP*EA_y,
                        color=b_color, stroke_width=3, dash_length=0.12, dashed_ratio=0.6)
        # Labels
        lab_Ec = Text("Ec", font_size=26, color=BLUE_E).next_to(Ec, LEFT, buff=0.1)
        lab_Ev = Text("Ev", font_size=26, color=RED_E).next_to(Ev, LEFT, buff=0.1)
        lab_EA = Text("EA (acceptor)", font_size=24, color=b_color).next_to(EA, LEFT, buff=0.1)
        panel = RoundedRectangle(width=width+0.6, height=3.2, corner_radius=0.2).set_opacity(0.10)
        panel.move_to(inset_origin)

        self.play(FadeIn(panel))
        self.play(Create(Ec), Create(Ev))
        self.play(Create(EA))
        self.play(FadeIn(lab_Ec), FadeIn(lab_Ev), FadeIn(lab_EA))

        # Drop a small dot from EA to Ev to suggest ionization (hole creation)
        carrier = Dot(color=b_color).move_to(inset_origin + LEFT*0.2 + UP*(EA_y))
        self.play(FadeIn(carrier, scale=0.6))
        self.play(carrier.animate.shift(DOWN*(EA_y - Ev_y - 0.08)), run_time=0.8)
        self.play(Flash(hole_g, flash_radius=0.25, color=hole_color))
        self.play(FadeOut(carrier), run_time=0.3)

        # Subtle field arrows (qualitative) pointing from B^- toward hole path
        arrows = VGroup()
        for k, off in enumerate(nbr_offsets[::2]):
            arr = Arrow(b_dot.get_center()+off*0.4, b_dot.get_center()+off*0.8, stroke_width=2, buff=0)
            arr.set_color(YELLOW_D).set_opacity(0.6)
            arrows.add(arr)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.2, run_time=1.0))

        # Outro pause
        self.wait(0.8)