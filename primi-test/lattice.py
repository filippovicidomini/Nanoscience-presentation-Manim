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
        
        

    
# --- SCENA 4C (pro): RETICOLO Si + B con elettroni, lacuna e elettrone libero ---
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
        e_color = BLUE_E

        # --------- Title ---------
        title = VGroup(
            Text("Silicon lattice — Boron acceptor", font_size=36, weight=BOLD),
            Text("Electrons on bonds, hole creation, free electron", font_size=24, slant=ITALIC)
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

        # Bonds: right and up-right neighbors (store for electron placement)
        bonds = VGroup()
        bond_segments = []
        for i in range(rows):
            for j in range(cols):
                p = pos(i, j)
                if j + 1 < cols:
                    q = pos(i, j+1)
                    seg = Line(p, q, stroke_width=2, color=bond_color, z_index=-1)
                    bonds.add(seg)
                    bond_segments.append((p, q, seg))
                if i - 1 >= 0 and j + 1 < cols:
                    q2 = pos(i-1, j+1)
                    seg2 = Line(p, q2, stroke_width=2, color=bond_color, z_index=-1)
                    bonds.add(seg2)
                    bond_segments.append((p, q2, seg2))

        self.play(LaggedStartMap(Create, bonds, lag_ratio=0.006, run_time=1.2))
        self.play(LaggedStart(*[FadeIn(a, scale=0.7) for a in atom_group], lag_ratio=0.004, run_time=1.0))
        self.wait(0.2)

        # --------- Choose substitution site (central) ---------
        ci, cj = rows//2, cols//2
        si_dot = atoms[(ci, cj)]
        center = si_dot.get_center()

        # Glow focus
        glow = Circle(radius=0.22, color=hole_color).set_stroke(width=5).move_to(center)
        pulse = glow.copy().set_opacity(0.5)
        self.play(Create(glow), FadeIn(pulse, scale=1.2))
        self.play(pulse.animate.scale(1.8).set_opacity(0), run_time=0.8)

        # --------- Place electrons on bonds (pairs near midpoints around center) ---------
        electron_group = VGroup()
        for (p, q, seg) in bond_segments:
            mid = (p + q)/2
            if np.linalg.norm(mid - center) < a*2.2:
                dir_vec = q - p
                perp = np.array([-dir_vec[1], dir_vec[0], 0])
                if np.linalg.norm(perp) == 0:
                    perp = np.array([0,1,0])
                perp = perp/np.linalg.norm(perp)*0.06
                e1 = Dot(mid + perp, radius=0.045, color=e_color)
                e2 = Dot(mid - perp, radius=0.045, color=e_color)
                electron_group.add(e1, e2)
        self.play(LaggedStart(*[FadeIn(e, scale=0.6) for e in electron_group], lag_ratio=0.01, run_time=0.6))

        # --------- Substitute with Boron ---------
        b_dot = Dot(center, radius=0.1, color=b_color)
        b_label = Text("B", font_size=26, weight=BOLD).move_to(center+0.28*UP)
        self.play(Transform(si_dot, b_dot), FadeIn(b_label, shift=UP*0.2))

        # --------- Capture one electron by B (acceptor) and create a hole on that bond ---------
        if len(electron_group) > 0:
            closest_e = min(electron_group, key=lambda d: np.linalg.norm(d.get_center()-center))
            prev_pos = closest_e.get_center()
            self.play(closest_e.animate.move_to(center), run_time=0.6)
            self.play(FadeOut(closest_e, scale=0.5), run_time=0.3)
        else:
            prev_pos = center + RIGHT*0.24

        hole = Circle(radius=0.11, color=hole_color, stroke_width=3).move_to(prev_pos)
        plus = Text("+", font_size=20).move_to(prev_pos)
        hole_g = VGroup(hole, plus).set_opacity(0.95)
        self.play(FadeIn(hole_g, scale=0.6))

        # --------- Delocalize the hole around B ---------
        nbr_offsets = [0.24*RIGHT, 0.2*UR, 0.24*UP, 0.2*UL, 0.24*LEFT, 0.2*DL, 0.24*DOWN, 0.2*DR]
        for _ in range(2):
            for off in nbr_offsets:
                self.play(hole_g.animate.move_to(center+off), run_time=0.16)
        self.play(hole_g.animate.move_to(center+0.24*RIGHT), run_time=0.2)

        # --------- Right inset: band diagram + electron capture + free electron ---------
        inset_origin = RIGHT*4.6 + DOWN*0.3
        width = 3.8
        Ec_y = 1.2
        Ev_y = -1.2
        EA_y = Ev_y + 0.35

        panel = RoundedRectangle(width=width+0.6, height=3.2, corner_radius=0.2).set_opacity(0.10).move_to(inset_origin)
        Ec = Line(inset_origin + LEFT*width/2 + UP*Ec_y, inset_origin + RIGHT*width/2 + UP*Ec_y, color=BLUE_E, stroke_width=4)
        Ev = Line(inset_origin + LEFT*width/2 + UP*Ev_y, inset_origin + RIGHT*width/2 + UP*Ev_y, color=RED_E, stroke_width=4)
        EA = DashedLine(inset_origin + LEFT*width/2 + UP*EA_y, inset_origin + RIGHT*width/2 + UP*EA_y, color=b_color, stroke_width=3, dash_length=0.12, dashed_ratio=0.6)
        lab_Ec = Text("Ec", font_size=26, color=BLUE_E).next_to(Ec, LEFT, buff=0.1)
        lab_Ev = Text("Ev", font_size=26, color=RED_E).next_to(Ev, LEFT, buff=0.1)
        lab_EA = Text("EA (acceptor)", font_size=24, color=b_color).next_to(EA, LEFT, buff=0.1)

        self.play(FadeIn(panel), Create(Ec), Create(Ev), Create(EA))
        self.play(FadeIn(lab_Ec), FadeIn(lab_Ev), FadeIn(lab_EA))

        # Electron captured from Ev to EA (leaves a hole in Ev)
        e_band = Dot(color=e_color).move_to(inset_origin + LEFT*0.9 + UP*(Ev_y))
        self.play(FadeIn(e_band, scale=0.6))
        self.play(e_band.animate.move_to(inset_origin + LEFT*0.4 + UP*(EA_y)), run_time=0.7)
        self.play(Flash(hole_g, flash_radius=0.25, color=hole_color))
        self.play(FadeOut(e_band), run_time=0.3)

        # Minority free electron in conduction band (animated motion)
        free_e = Dot(color=e_color).move_to(inset_origin + LEFT*width/2 + RIGHT*0.25 + UP*(Ec_y))
        self.play(FadeIn(free_e, scale=0.6))
        self.play(free_e.animate.move_to(inset_origin + RIGHT*width/2 + LEFT*0.25 + UP*(Ec_y)), run_time=1.2)
        self.play(free_e.animate.move_to(inset_origin + LEFT*width/2 + RIGHT*0.9 + UP*(Ec_y)), run_time=1.0)

        # Subtle field arrows around B
        arrows = VGroup()
        for off in nbr_offsets[::2]:
            arr = Arrow(center+off*0.4, center+off*0.8, stroke_width=2, buff=0).set_color(YELLOW_D).set_opacity(0.6)
            arrows.add(arr)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.2, run_time=1.0))

        # Outro pause
        self.wait(0.8)
