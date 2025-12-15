from manim import *
import numpy as np
import random

class SiliconLatticeWithFreeElectron(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f17"

        # -----------------------------
        # PARAMETRI RETICOLO
        # -----------------------------
        rows = 5
        cols = 9
        spacing = 2.0

        # -----------------------------
        # 1) COSTRUZIONE RETICOLO (CORE ONLY)
        # -----------------------------
        atom_grid = [[None for _ in range(cols)] for _ in range(rows)]
        atoms = VGroup()

        for i in range(rows):
            for j in range(cols):
                x = (j - (cols - 1) / 2) * spacing
                y = (i - (rows - 1) / 2) * spacing
                pos = np.array([x, y, 0.0])

                atom = self.make_valence_atom()
                atom.move_to(pos)

                # NON mostriamo gli elettroni di valenza
                atom.electrons.set_opacity(0)

                atom_grid[i][j] = atom
                atoms.add(atom)
                
        self.add(atoms)
        self.wait(0.3)

        # -----------------------------
        # 2) LEGAMI COVALENTI (COPPIE e-)
        # -----------------------------
        bond_electrons = VGroup()

        for i in range(rows):
            for j in range(cols - 1):
                pair = self.make_bond_pair(
                    atom_grid[i][j].get_center(),
                    atom_grid[i][j + 1].get_center()
                )
                bond_electrons.add(pair)

        for i in range(rows - 1):
            for j in range(cols):
                pair = self.make_bond_pair(
                    atom_grid[i][j].get_center(),
                    atom_grid[i + 1][j].get_center()
                )
                bond_electrons.add(pair)
                
        self.add(bond_electrons)
        self.wait(0.5)

        # -----------------------------
        # 3) DOPING: FOSFORO
        # -----------------------------
        dop_i, dop_j = 2, 4
        dop_atom = atom_grid[dop_i][dop_j]

        self.play(Indicate(dop_atom.core, scale_factor=1.05), run_time=0.7)

        p_nucleus = Circle(
            radius=0.28,
            color=ORANGE,
            fill_opacity=1.0
        ).set_stroke(width=0).move_to(dop_atom.nucleus.get_center())

        p_label = Text("P", font_size=30, color=ORANGE).move_to(dop_atom.get_center())

        self.play(
            Transform(dop_atom.nucleus, p_nucleus),
            FadeIn(p_label, shift=UP * 0.1),
            run_time=0.9
        )
        self.wait(0.3)

        # -----------------------------
        # 4) ELETTRONE LIBERO (QUASI BALISTICO)
        # -----------------------------
        e_free = Dot(radius=0.075, color=RED)
        e_free.move_to(dop_atom.get_center() + 0.6 * RIGHT + 0.25 * UP)

        glow = Circle(radius=0.16).set_stroke(BLUE_C, width=2, opacity=0.55)
        glow.move_to(e_free.get_center())

        e_label = Text("e⁻", font_size=28, color=BLUE_C).next_to(e_free, UP, buff=0.12)
        e_label.add_updater(lambda m: m.next_to(e_free, UP, buff=0.12))
        glow.add_updater(lambda m, dt: m.move_to(e_free.get_center()))

        self.play(
            FadeIn(e_free, scale=0.6),
            FadeIn(glow, scale=0.6),
            FadeIn(e_label, shift=UP * 0.1),
            run_time=0.7
        )

        self.add_free_electron_motion_potential(e_free, atom_grid)
        self.wait(20)

    # =========================================================
    # ATOMO BASE
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
    # COPPIA DI LEGAME (CON VIBRAZIONE)
    # =========================================================
    def make_bond_pair(self, p1, p2, sep=0.18):
        mid = 0.5 * (p1 + p2)

        d = p2 - p1
        d[2] = 0
        d /= np.linalg.norm(d[:2]) or 1

        perp = np.array([-d[1], d[0], 0.0])

        e1 = Dot(mid + sep * perp, radius=0.07, color="#6EA8FE")
        e2 = Dot(mid - sep * perp, radius=0.07, color="#6EA8FE")

        self.add_oscillation_to_bond_electron(e1, mid, sep * perp, phase_shift=0)
        self.add_oscillation_to_bond_electron(e2, mid, -sep * perp, phase_shift=PI)

        return VGroup(e1, e2)

    # =========================================================
    # VIBRAZIONE ELETTRONI DI LEGAME
    # =========================================================
    def add_oscillation_to_bond_electron(
        self,
        dot,
        base_center,
        base_offset,
        amplitude=0.03,
        omega=3.0,
        phase_shift=0.0,
        noise_amp=0.004,
    ):
        t = [0.0]
        perp = base_offset / (np.linalg.norm(base_offset[:2]) or 1)

        def updater(m, dt):
            t[0] += dt
            osc = amplitude * np.sin(omega * t[0] + phase_shift)
            jitter = noise_amp * np.random.normal(size=3)
            jitter[2] = 0
            m.move_to(base_center + base_offset + osc * perp + jitter)

        dot.add_updater(updater)

    # =========================================================
    # ELETTRONE LIBERO — PRESET A
    # =========================================================
    def add_free_electron_motion_potential(
        self,
        e_dot,
        atom_grid,
        D=0.98,
        gamma=0.15,
        g_vec=np.array([0.0, 0.0, 0.0]),
        k_rep=3.2,
        a=0.30,
        p=6.0,
        max_speed=5.5,
        bound_margin=0.95,
    ):
        nuclei = np.array([atom_grid[i][j].get_center()
                           for i in range(len(atom_grid))
                           for j in range(len(atom_grid[0]))])

        min_x = nuclei[:, 0].min() + bound_margin
        max_x = nuclei[:, 0].max() - bound_margin
        min_y = nuclei[:, 1].min() + bound_margin
        max_y = nuclei[:, 1].max() - bound_margin

        v = np.zeros(3)

        def repulsion(x):
            r = x - nuclei
            r[:, 2] = 0
            r2 = np.sum(r[:, :2]**2, axis=1) + 1e-8
            denom = (r2 + a * a) ** ((p + 2) / 2)
            return (k_rep * (a**p) * r / denom[:, None]).sum(axis=0)

        def updater(m, dt):
            nonlocal v
            dt = max(dt, 1e-3)
            x = m.get_center()

            F = g_vec + repulsion(x)
            noise = np.sqrt(2 * D * dt) * np.array([np.random.normal(), np.random.normal(), 0])
            v[:] = v + (F - gamma * v) * dt + noise

            speed = np.linalg.norm(v[:2])
            if speed > max_speed:
                v[:] *= max_speed / speed

            x_new = x + v * dt
            if x_new[0] < min_x or x_new[0] > max_x:
                v[0] *= -0.7
            if x_new[1] < min_y or x_new[1] > max_y:
                v[1] *= -0.7

            m.move_to(x + v * dt)

        e_dot.add_updater(updater)