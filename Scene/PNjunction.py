from manim import *
import numpy as np
import random


class PNJunctionFormation(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # -----------------------------
        # PARAMETRI
        # -----------------------------
        rows = 5
        cols = 4
        spacing = 1.8

        # dopanti (pochi, sparsi) – indici coerenti con griglia rows x cols
        dopants_n = [(1, 0), (3, 2), (2, 3), (0, 1)]  # fosforo (donatori)
        dopants_p = [(1, 1), (3, 2), (2, 0), (4, 3)]  # boro (accettori)

        # -----------------------------
        # 1) COSTRUISCI DUE RETICOLI (stesso stile delle altre scene)
        # -----------------------------
        atom_grid_n, bonds_n, bond_group_n, lattice_n = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=LEFT * 6.0
        )
        atom_grid_p, bonds_p, bond_group_p, lattice_p = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=RIGHT * 6.0
        )

        # Mostra: core + bond electrons
        # (nelle altre scene aggiungi atoms e poi bond_group; qui facciamo uguale)
        self.add(lattice_n)
        self.add(lattice_p)
        self.wait(1 / 60)

        label_n = Text("n-type", font_size=34, color=BLUE_B).next_to(lattice_n, UP, buff=0.35)
        label_p = Text("p-type", font_size=34, color=RED_B).next_to(lattice_p, UP, buff=0.35)
        self.play(FadeIn(label_n), FadeIn(label_p), run_time=0.6)

        # -----------------------------
        # 2) DOPING VISIVO (P e B) – stesso stile trasform nuclei + label
        # -----------------------------

        for (i, j) in dopants_n:
            atom = atom_grid_n[i][j]
            p_nucleus = Circle(radius=0.28, color=ORANGE, fill_opacity=1.0).set_stroke(width=0)
            p_nucleus.move_to(atom.nucleus.get_center())
            self.play(Transform(atom.nucleus, p_nucleus), run_time=0.25)

        for (i, j) in dopants_p:
            atom = atom_grid_p[i][j]
            b_nucleus = Circle(radius=0.22, color=PURPLE_C, fill_opacity=1.0).set_stroke(width=0)
            b_nucleus.move_to(atom.nucleus.get_center())
            self.play(Transform(atom.nucleus, b_nucleus), run_time=0.25)

        self.wait(0.2)

        # -----------------------------
        # 3) AVVICINAMENTO (formazione giunzione)
        # -----------------------------
        self.play(
            lattice_n.animate.shift(RIGHT * 2.4),
            label_n.animate.shift(RIGHT * 2.4),
            lattice_p.animate.shift(LEFT * 2.4),
            label_p.animate.shift(LEFT * 2.4),
            run_time=1.6,
            rate_func=smooth,
        )

        # -----------------------------
        # 3b) CREA LEGAMI ALL'INTERFACCIA (continuazione reticolo)
        #     Due elettroni condivisi tra atomi di bordo (n|p)
        # -----------------------------
        interface_pairs = VGroup()
        interface_bonds = []

        for i in range(rows):
            left_atom = atom_grid_n[i][cols - 1]   # colonna di bordo lato n
            right_atom = atom_grid_p[i][0]         # colonna di bordo lato p
            b = self.make_bond(left_atom, right_atom, (i, cols - 1), (i, 0), sep=0.18)
            interface_bonds.append(b)
            interface_pairs.add(b["pair"])

        # aggiungi i nuovi elettroni di legame sopra tutto, e falli comparire
        self.add(interface_pairs)
        self.play(FadeIn(interface_pairs, shift=UP * 0.05), run_time=0.5)

        # -----------------------------
        # 4) DEPLETION REGION (highlight) – coerente con palette
        # -----------------------------
        depletion = RoundedRectangle(
            corner_radius=0.25,
            width=2.0,
            height=(rows - 1) * spacing + 2.0,
            stroke_width=2,
            stroke_color=YELLOW_B,
            fill_color=YELLOW,
            fill_opacity=0.08,
        ).move_to(0.5 * (atom_grid_n[rows // 2][cols - 1].get_center() + atom_grid_p[rows // 2][0].get_center()))
        self.play(FadeIn(depletion), run_time=0.5)

        # -----------------------------
        # 5) DIFFUSIONE PORTATORI + RICOMBINAZIONE
        #    (solo vicino all’interfaccia)
        # -----------------------------
        carriers = VGroup()
        electrons = []
        holes = []

        for k in range(4):
            y = (k - 1.5) * 0.65
            e = self.make_free_electron(LEFT * 1.2 + UP * y)
            h = self.make_hole(RIGHT * 1.2 + UP * y)
            carriers.add(e, h)
            electrons.append(e)
            holes.append(h)

        self.play(FadeIn(carriers), run_time=0.35)

        # diffusione: e- n->p, h+ p->n
        self.play(
            *[e.animate.shift(RIGHT * 1.8) for e in electrons],
            *[h.animate.shift(LEFT * 1.8) for h in holes],
            run_time=1.0,
            rate_func=rate_functions.ease_in_out_sine,
        )

        # ricombinazione: effetto anello + fade
        for e, h in zip(electrons, holes):
            self.recombine(e, h)

        self.wait(0.2)

        # -----------------------------
        # 6) IONI FISSI + CAMPO E
        # -----------------------------
        fixed = VGroup()
        for k in range(4):
            y = (k - 1.5) * 0.65
            dplus = MathTex("D^+", font_size=30, color=BLUE_B).move_to(LEFT * 0.65 + UP * y)
            aminus = MathTex("A^-", font_size=30, color=RED_B).move_to(RIGHT * 0.65 + UP * y)
            fixed.add(dplus, aminus)
        self.play(FadeIn(fixed), run_time=0.6)

        arrows = VGroup()
        for k in range(3):
            y = (k - 1) * 0.8
            a = Arrow(RIGHT * 0.85 + UP * y, LEFT * 0.85 + UP * y, buff=0.0, stroke_width=4)
            arrows.add(a)
        self.play(*[Create(a) for a in arrows], run_time=0.6)

        eq = MathTex("J_{diff}+J_{drift}=0", font_size=44, color=WHITE).to_edge(DOWN)
        self.play(FadeIn(eq), run_time=0.6)
        self.wait(1.2)

        # pulizia updater (stesso approccio robusto)
        for mob in lattice_n.family_members_with_points():
            mob.clear_updaters()
        for mob in lattice_p.family_members_with_points():
            mob.clear_updaters()

    # =========================================================
    # COSTRUZIONE RETICOLO + LEGAMI (coppie vibranti) – COPIATO DAL TUO STILE
    # =========================================================
    def build_lattice_with_bonds(self, rows, cols, spacing, origin=ORIGIN):
        atom_grid = [[None for _ in range(cols)] for _ in range(rows)]
        bond_group = VGroup()
        atoms = VGroup()
        bonds = []  # lista di dict come in dopingHole

        # atomi (core + elettroni di valenza “spenti”)
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

        # bond orizzontali
        for i in range(rows):
            for j in range(cols - 1):
                a_atom = atom_grid[i][j]
                b_atom = atom_grid[i][j + 1]
                b = self.make_bond(a_atom, b_atom, (i, j), (i, j + 1))
                bonds.append(b)
                bond_group.add(b["pair"])

        # bond verticali
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
        """Create a covalent bond as a vibrating electron pair.

        IMPORTANT: electrons are updated dynamically from atom_a/atom_b positions,
        so when the lattice shifts, the bond electrons shift coherently too.
        """
        # initialize near the correct midpoint
        p1 = atom_a.get_center()
        p2 = atom_b.get_center()
        mid0 = 0.5 * (p1 + p2)

        # placeholder dots; their true position is controlled by updaters
        e1 = Dot(mid0, radius=0.07, color="#6EA8FE")
        e2 = Dot(mid0, radius=0.07, color="#6EA8FE")

        self.attach_bond_osc_dynamic(e1, atom_a, atom_b, side=+1, sep=sep, phase_shift=0.0)
        self.attach_bond_osc_dynamic(e2, atom_a, atom_b, side=-1, sep=sep, phase_shift=PI)

        pair = VGroup(e1, e2)
        return {
            "a": a_idx,
            "b": b_idx,
            "atom_a": atom_a,
            "atom_b": atom_b,
            "e1": e1,
            "e2": e2,
            "pair": pair,
        }

    # =========================================================
    # ATOMO BASE (uguale alle altre scene)
    # =========================================================
    def make_valence_atom(self, orbit_radius=0.8) -> VGroup:
        nucleus = Circle(radius=0.22, color="#66E0A3", fill_opacity=1.0).set_stroke(width=0)
        orbit = Circle(radius=orbit_radius, stroke_width=2.0, stroke_opacity=0.85)
        orbit.set_stroke(color="#ffe680")

        directions = {"up": UP, "right": RIGHT, "down": DOWN, "left": LEFT}
        electrons = VGroup()
        e_by_dir = {}
        for name, d in directions.items():
            e = Dot(point=orbit_radius * d, radius=0.07, color="#6EA8FE")
            electrons.add(e)
            e_by_dir[name] = e

        atom = VGroup(orbit, nucleus, electrons)
        atom.orbit = orbit
        atom.nucleus = nucleus
        atom.electrons = electrons
        atom.e_by_dir = e_by_dir
        atom.core = VGroup(orbit, nucleus)
        return atom

    # =========================================================
    # OSCILLAZIONE ELETTRONI DI LEGAME DINAMICA (aggiornatori legati a due atomi)
    # =========================================================
    def attach_bond_osc_dynamic(
        self,
        dot,
        atom_a,
        atom_b,
        side=+1,
        sep=0.18,
        amplitude=0.03,
        omega=3.0,
        phase_shift=0.0,
        noise_amp=0.004,
    ):
        """Attach an updater that keeps a bond-electron tied to (atom_a, atom_b)."""
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

            # oscillate along perp direction (same visual style as your other scenes)
            mobj.move_to(mid + base_offset + osc * perp + jitter)

        dot.add_updater(updater)

    # =========================================================
    # PORTATORI MOBILI + RICOMBINAZIONE (robusta)
    # =========================================================
    def make_free_electron(self, pos):
        # Evidenzia elettrone libero: dot rosso + doppio glow pulsante
        e = Dot(radius=0.085, color=RED).move_to(pos)

        glow1 = Circle(radius=0.18, stroke_width=2, color=BLUE_C, stroke_opacity=0.8)
        glow2 = Circle(radius=0.28, stroke_width=2, color=BLUE_C, stroke_opacity=0.5)
        glow1.move_to(pos)
        glow2.move_to(pos)

        def pulse(m, dt):
            if not hasattr(m, "t"):
                m.t = 0
            m.t += dt
            s = 1.0 + 0.15 * np.sin(3 * m.t)
            m.set(width=m.width * s)

        glow2.add_updater(pulse)
        glow1.add_updater(lambda m, dt: m.move_to(e.get_center()))
        glow2.add_updater(lambda m, dt: m.move_to(e.get_center()))

        label = MathTex("e^-", font_size=26, color=BLUE_C)
        label.add_updater(lambda m: m.next_to(e, UP, buff=0.12))

        return VGroup(e, glow1, glow2, label)

    def make_hole(self, pos):
        # Evidenzia lacuna: anello giallo + label
        h = Circle(radius=0.11, stroke_width=4, color=YELLOW_C).move_to(pos)
        label = MathTex("h^+", font_size=26, color=YELLOW_C)
        label.add_updater(lambda m: m.next_to(h, UP, buff=0.12))
        return VGroup(h, label)

    def recombine(self, e_group, h_mob):
        # e_group è VGroup(e_dot, glow, label)
        e_dot = e_group[0]
        center = e_dot.get_center()
        ring = Circle(radius=0.10, stroke_width=4, color=YELLOW).move_to(center)
        self.play(
            AnimationGroup(
                ring.animate.scale(3.0).set_stroke(opacity=0.0),
                FadeOut(e_group),
                FadeOut(h_mob),
                lag_ratio=0.0,
            ),
            run_time=0.45,
        )