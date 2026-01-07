from manim import *
import numpy as np
import random

config.frame_rate = 60

class PNJunctionReverseBiasWithDopants(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # -----------------------------
        # PARAMETRI
        # -----------------------------
        rows = 5
        cols = 4
        spacing = 1.8
        
        origin_n = LEFT * 1.8 
        origin_p = RIGHT * 1.8

        # -----------------------------
        # 1) SETUP AMBIENTE
        # -----------------------------
        atom_grid_n, bonds_n, bond_group_n, lattice_n = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=origin_n
        )
        atom_grid_p, bonds_p, bond_group_p, lattice_p = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=origin_p
        )

        interface_pairs = VGroup()
        for i in range(rows):
            left_atom = atom_grid_n[i][cols - 1]
            right_atom = atom_grid_p[i][0]
            b = self.make_bond(left_atom, right_atom, (i, cols-1), (i, 0), sep=0.18)
            interface_pairs.add(b["pair"])

        # --- DROGANTI ---
        dopants_n_coords = [(1, 1), (3, 2), (0, 3)]
        dopants_p_coords = [(1, 2), (3, 1), (4, 0)]

        # Donatori su N (arancione)
        for (i, j) in dopants_n_coords:
            atom = atom_grid_n[i][j]
            p_nucleus = Circle(radius=0.28, color=ORANGE, fill_opacity=1.0).set_stroke(width=0)
            p_nucleus.move_to(atom.nucleus.get_center())
            atom.nucleus.become(p_nucleus)

        # Accettori su P (viola)
        for (i, j) in dopants_p_coords:
            atom = atom_grid_p[i][j]
            b_nucleus = Circle(radius=0.22, color=PURPLE_C, fill_opacity=1.0).set_stroke(width=0)
            b_nucleus.move_to(atom.nucleus.get_center())
            atom.nucleus.become(b_nucleus)

        semiconductor = VGroup(lattice_n, lattice_p, interface_pairs)
        self.add(semiconductor)

        # -----------------------------
        # 2) DEPLETION REGION INIZIALE
        # -----------------------------
        depletion_width_initial = 2.5
        depletion = RoundedRectangle(
            corner_radius=0.25,
            width=depletion_width_initial,
            height=(rows - 1) * spacing + 1.5,
            stroke_width=2,
            stroke_color=YELLOW_B,
            fill_color=YELLOW,
            fill_opacity=0.12,
        ).move_to(ORIGIN)

        self.play(FadeIn(depletion), run_time=1.5)
        self.wait(2)

        # -----------------------------
        # 3) APPLICAZIONE TENSIONE (REVERSE BIAS)
        # -----------------------------
        wire_y = -4.0
        contact_left = lattice_n.get_left() + LEFT * 0.5   # lato N
        contact_right = lattice_p.get_right() + RIGHT * 0.5  # lato P
        
        circuit = VMobject(color=WHITE, stroke_width=3)
        circuit.set_points_as_corners([
            contact_left + UP*1.5,
            contact_left,
            [contact_left[0], wire_y, 0],
            [contact_right[0], wire_y, 0],
            contact_right,
            contact_right + UP*1.5
        ])

        battery = self.get_battery_symbol().move_to([0, wire_y, 0])

        # Reverse bias: N al +, P al -
        plus_label  = MathTex("+").set_color(WHITE).scale(0.9)
        minus_label = MathTex("-").set_color(WHITE).scale(0.9)

        plus_label.next_to(contact_left + UP*1.5, UP, buff=0.15)    # N side = +
        minus_label.next_to(contact_right + UP*1.5, UP, buff=0.15) # P side = -

        self.play(
            Create(circuit),
            FadeIn(battery),
            FadeIn(plus_label),
            FadeIn(minus_label),
            run_time=4.0,
            lag_ratio=0.8
        )

        # -----------------------------
        # 4) EFFETTO: AUMENTO BARRIERA (DEPLETION SI ALLARGA)
        # -----------------------------
        widened_width = 3.8  # prova 3.4 / 3.8 / 4.2 a gusto

        self.play(
            depletion.animate.stretch_to_fit_width(widened_width).set_fill(opacity=0.18),
            run_time=2.6,
            rate_func=rate_functions.ease_in_out_cubic
        )
        self.wait(2)

        # -----------------------------
        # 5) CORRENTE QUASI NULLA (NO INIEZIONE)
        # -----------------------------
        # Opzione A (super minimal): non far muovere nulla.
        # self.wait(2.0)

        # Opzione B: mostra un *piccolissimo leakage* (1 elettrone e 1 lacuna che si muovono poco e svaniscono)
        for _ in range(5):
            y = random.choice([-1.2, 0.0, 1.2])

            # Minority-like drift molto piccolo (giusto per suggerire "leakage")
            e = self.make_free_electron(lattice_n.get_left() + RIGHT*0.6 + UP*y)
            h = self.make_hole(lattice_p.get_right() + LEFT*0.6 + UP*y)

            self.add(e, h)
            self.play(
                e.animate.shift(RIGHT*0.9),
                h.animate.shift(LEFT*0.9),
                run_time=1.2,
                rate_func=linear
            )
            self.play(
                FadeOut(e, scale=0.9),
                FadeOut(h, scale=0.9),
                run_time=2
            )

        # Freccia corrente: piccola e tenue (reverse -> ~0)
        tiny_current = Arrow(LEFT*1.0, RIGHT*1.0, color=YELLOW, buff=0, stroke_width=3)\
            .next_to(battery, DOWN, buff=0.5)\
            .set_opacity(0.35)
        self.play(GrowArrow(tiny_current), run_time=1.2)
        self.wait(3)

        for m in self.mobjects:
            m.clear_updaters()

    # =========================================================
    # FUNZIONI DI SUPPORTO (IDENTICHE ALLE TUE)
    # =========================================================

    def get_battery_symbol(self):
        plate_neg = Line(UP*0.15, DOWN*0.15, color=WHITE, stroke_width=4).move_to(LEFT*0.15)
        plate_pos = Line(UP*0.3, DOWN*0.3, color=WHITE, stroke_width=4).move_to(RIGHT*0.15)
        wire_l = Line(LEFT*0.5, plate_neg.get_center(), color=WHITE)
        wire_r = Line(plate_pos.get_center(), RIGHT*0.5, color=WHITE)
        return VGroup(wire_l, wire_r, plate_neg, plate_pos)

    def build_lattice_with_bonds(self, rows, cols, spacing, origin=ORIGIN):
        atom_grid = [[None for _ in range(cols)] for _ in range(rows)]
        bond_group = VGroup()
        atoms = VGroup()
        bonds = []

        for i in range(rows):
            for j in range(cols):
                x = (j - (cols - 1) / 2) * spacing
                y = (i - (rows - 1) / 2) * spacing
                pos = np.array([x, y, 0.0]) + origin

                atom = self.make_valence_atom()
                atom.move_to(pos)
                atom.electrons.set_opacity(0)

                atom_grid[i][j] = atom
                atoms.add(atom)

        for i in range(rows):
            for j in range(cols - 1):
                a_atom = atom_grid[i][j]
                b_atom = atom_grid[i][j + 1]
                b = self.make_bond(a_atom, b_atom, (i, j), (i, j + 1))
                bonds.append(b)
                bond_group.add(b["pair"])

        for i in range(rows - 1):
            for j in range(cols):
                a_atom = atom_grid[i][j]
                b_atom = atom_grid[i + 1][j]
                b = self.make_bond(a_atom, b_atom, (i, j), (i + 1, j))
                bonds.append(b)
                bond_group.add(b["pair"])

        lattice = VGroup(atoms, bond_group)
        return atom_grid, bonds, bond_group, lattice

    def make_bond(self, atom_a, atom_b, a_idx, b_idx, sep=0.18):
        p1 = atom_a.get_center()
        p2 = atom_b.get_center()
        mid0 = 0.5 * (p1 + p2)

        e1 = Dot(mid0, radius=0.07, color="#6EA8FE")
        e2 = Dot(mid0, radius=0.07, color="#6EA8FE")

        self.attach_bond_osc_dynamic(e1, atom_a, atom_b, side=+1, sep=sep, phase_shift=0.0)
        self.attach_bond_osc_dynamic(e2, atom_a, atom_b, side=-1, sep=sep, phase_shift=PI)

        pair = VGroup(e1, e2)
        return {"a": a_idx, "b": b_idx, "atom_a": atom_a, "atom_b": atom_b,
                "e1": e1, "e2": e2, "pair": pair}

    def make_valence_atom(self, orbit_radius=0.8) -> VGroup:
        nucleus = Circle(radius=0.22, color="#66E0A3", fill_opacity=1.0).set_stroke(width=0)
        orbit = Circle(radius=orbit_radius, stroke_width=2.0, stroke_opacity=0.85)
        orbit.set_stroke(color="#ffe680")

        directions = {"up": UP, "right": RIGHT, "down": DOWN, "left": LEFT}
        electrons = VGroup()
        for _, d in directions.items():
            e = Dot(point=orbit_radius * d, radius=0.07, color="#6EA8FE")
            electrons.add(e)

        atom = VGroup(orbit, nucleus, electrons)
        atom.orbit = orbit
        atom.nucleus = nucleus
        atom.electrons = electrons
        atom.core = VGroup(orbit, nucleus)
        return atom

    def attach_bond_osc_dynamic(self, dot, atom_a, atom_b, side=+1, sep=0.18,
                                amplitude=0.03, omega=3.0, phase_shift=0.0, noise_amp=0.004):
        dot.clear_updaters()
        t = [0.0]
        def updater(mobj, dt):
            t[0] += dt
            p1 = atom_a.get_center()
            p2 = atom_b.get_center()
            mid = 0.5 * (p1 + p2)
            d = p2 - p1
            d[2] = 0.0
            nrm = np.linalg.norm(d[:2]) or 1.0
            d = d / nrm
            perp = np.array([-d[1], d[0], 0.0])
            base_offset = side * sep * perp
            osc = amplitude * np.sin(omega * t[0] + phase_shift)
            jitter = noise_amp * np.random.normal(size=3)
            jitter[2] = 0.0
            mobj.move_to(mid + base_offset + osc * perp + jitter)
        dot.add_updater(updater)

    def make_free_electron(self, pos):
        return Dot(radius=0.085, color=RED).move_to(pos)

    def make_hole(self, pos):
        return Circle(radius=0.11, stroke_width=4, color=YELLOW_C).move_to(pos)