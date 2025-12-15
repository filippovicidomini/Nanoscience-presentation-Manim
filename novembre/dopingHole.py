from manim import *
import numpy as np
import random

class BoronLatticeWithHole(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f17"

        # -----------------------------
        # PARAMETRI
        # -----------------------------
        rows = 5
        cols = 9
        spacing = 2.0

        # -----------------------------
        # 1) RETICOLO: nuclei+orbite (già presenti)
        #    + elettroni di legame (già presenti e vibranti)
        # -----------------------------
        atom_grid, bonds, bond_group = self.build_lattice_with_bonds(rows, cols, spacing)

        self.add(VGroup(*[atom_grid[i][j].core for i in range(rows) for j in range(cols)]))
        self.add(bond_group)
        self.wait(1/60)  # stabilizza updater al primo frame

        # -----------------------------
        # 2) DOPING: sostituisci un atomo con BORO
        # -----------------------------
        dop_i, dop_j = 2, 4
        dop_atom = atom_grid[dop_i][dop_j]

        self.play(Indicate(dop_atom.core, scale_factor=1.05), run_time=0.7)

        b_nucleus = Circle(radius=0.22, color=PURPLE_C, fill_opacity=1.0).set_stroke(width=0)
        b_nucleus.move_to(dop_atom.nucleus.get_center())

        b_label = Text("B", font_size=30, color=PURPLE_C).move_to(dop_atom.get_center())

        self.play(
            Transform(dop_atom.nucleus, b_nucleus),
            FadeIn(b_label, shift=UP * 0.1),
            run_time=0.9
        )
        self.wait(0.2)

        # -----------------------------
        # 3) CREA UNA LACUNA: togli 1 elettrone da un legame vicino al dopante
        # -----------------------------
        # scegli un bond che tocca il dopante
        near = self.bonds_touching_atom(bonds, (dop_i, dop_j))
        hole_bond = random.choice(near)

        # togli uno dei due elettroni (per esempio e2)
        self.make_hole_in_bond(hole_bond, missing="e2")

        # marker visivo della lacuna
        hole_marker = Circle(radius=0.10, color=YELLOW_C).set_stroke(width=4).move_to(hole_bond["mid"])
        hole_label = Text("h⁺", font_size=28, color=YELLOW_C).next_to(hole_marker, UP, buff=0.12)

        hole_label.add_updater(lambda m: m.next_to(hole_marker, UP, buff=0.12))
        self.add(hole_marker, hole_label)

        self.wait(0.4)

        # -----------------------------
        # 4) RANDOM WALK della lacuna sui legami
        # -----------------------------
        steps = 16
        current = hole_bond
        prev = None

        for _ in range(steps):
            # scegli un bond adiacente (condivide un atomo) per far “migrare” la lacuna
            neighbors = self.neighbor_bonds(bonds, current)
            if prev is not None and len(neighbors) > 1:
                neighbors = [b for b in neighbors if b is not prev] or neighbors
            nxt = random.choice(neighbors)

            # “salto” elettronico: un elettrone dal bond nxt riempie il buco in current
            self.hop_hole(current, nxt, hole_marker)

            prev, current = current, nxt

        self.wait(2.0)

    # =========================================================
    # COSTRUZIONE RETICOLO + LEGAMI (coppie vibranti)
    # =========================================================
    def build_lattice_with_bonds(self, rows, cols, spacing):
        atom_grid = [[None for _ in range(cols)] for _ in range(rows)]
        bonds = []
        bond_group = VGroup()

        # atomi (core + elettroni di valenza "spenti")
        for i in range(rows):
            for j in range(cols):
                x = (j - (cols - 1)/2) * spacing
                y = (i - (rows - 1)/2) * spacing
                pos = np.array([x, y, 0.0])

                atom = self.make_valence_atom()
                atom.move_to(pos)
                atom.electrons.set_opacity(0)  # non li usiamo in questa scena

                atom_grid[i][j] = atom

        # crea bond orizzontali
        for i in range(rows):
            for j in range(cols - 1):
                p1 = atom_grid[i][j].get_center()
                p2 = atom_grid[i][j+1].get_center()
                b = self.make_bond(p1, p2, (i, j), (i, j+1))
                bonds.append(b)
                bond_group.add(b["pair"])

        # crea bond verticali
        for i in range(rows - 1):
            for j in range(cols):
                p1 = atom_grid[i][j].get_center()
                p2 = atom_grid[i+1][j].get_center()
                b = self.make_bond(p1, p2, (i, j), (i+1, j))
                bonds.append(b)
                bond_group.add(b["pair"])

        return atom_grid, bonds, bond_group

    def make_bond(self, p1, p2, a_idx, b_idx, sep=0.18):
        mid = 0.5 * (p1 + p2)

        d = p2 - p1
        d[2] = 0
        d /= (np.linalg.norm(d[:2]) or 1.0)
        perp = np.array([-d[1], d[0], 0.0])

        off1 = sep * perp
        off2 = -sep * perp

        e1 = Dot(mid + off1, radius=0.07, color="#6EA8FE")
        e2 = Dot(mid + off2, radius=0.07, color="#6EA8FE")

        # vibrazioni (come il tuo stile)
        self.attach_bond_osc(e1, mid, off1, phase_shift=0.0)
        self.attach_bond_osc(e2, mid, off2, phase_shift=PI)

        pair = VGroup(e1, e2)

        return {
            "a": a_idx, "b": b_idx,
            "mid": mid,
            "off1": off1, "off2": off2,
            "e1": e1, "e2": e2,
            "pair": pair,
            "missing": None  # "e1" o "e2"
        }

    # =========================================================
    # LACUNA: rimuovi 1 elettrone da un bond
    # =========================================================
    def make_hole_in_bond(self, bond, missing="e2"):
        bond["missing"] = missing
        dot = bond[missing]
        dot.clear_updaters()
        dot.set_fill(opacity=0)
        dot.set_stroke(opacity=0)

    # =========================================================
    # ADIACENZE TRA LEGAMI
    # =========================================================
    def bonds_touching_atom(self, bonds, atom_idx):
        out = []
        for b in bonds:
            if b["a"] == atom_idx or b["b"] == atom_idx:
                out.append(b)
        return out

    def neighbor_bonds(self, bonds, bond):
        # due bond sono vicini se condividono almeno un atomo
        a1, b1 = bond["a"], bond["b"]
        out = []
        for bb in bonds:
            if bb is bond:
                continue
            if bb["a"] in (a1, b1) or bb["b"] in (a1, b1):
                out.append(bb)
        return out

    # =========================================================
    # HOP: un elettrone da nxt riempie la lacuna in current
    #      e la lacuna si sposta su nxt
    # =========================================================
    def hop_hole(self, current, nxt, hole_marker):
        # current DEVE avere una lacuna
        if current["missing"] is None:
            return

        # scegli quale elettrone "donare" da nxt: preferisci uno visibile
        donor_name = "e1" if nxt["e1"].get_fill_opacity() > 0.5 else "e2"
        donor = nxt[donor_name]

        # posizione da riempire nel bond current
        missing_name = current["missing"]
        target_pos = current["mid"] + (current["off1"] if missing_name == "e1" else current["off2"])

        # crea una copia temporanea che si muove (più controllabile)
        temp = donor.copy()
        temp.clear_updaters()
        temp.set_fill(opacity=1.0)
        temp.set_stroke(opacity=1.0)
        self.add(temp)

        # "svuota" il donor nel bond nxt -> nuova lacuna lì
        donor.clear_updaters()
        donor.set_fill(opacity=0)
        donor.set_stroke(opacity=0)
        nxt["missing"] = donor_name

        # animazione: temp va a riempire la lacuna
        self.play(
            temp.animate.move_to(target_pos),
            run_time=0.25,
            rate_func=smooth
        )

        # ripristina nel bond current: riaccendi il dot mancante (quello originale)
        fill_dot = current[missing_name]
        fill_dot.move_to(target_pos)
        fill_dot.set_fill(opacity=1.0)
        fill_dot.set_stroke(opacity=1.0)
        self.attach_bond_osc(
            fill_dot,
            current["mid"],
            current["off1"] if missing_name == "e1" else current["off2"],
            phase_shift=0.0 if missing_name == "e1" else PI
        )
        current["missing"] = None

        # rimuovi la particella temp (era solo “trasferimento visivo”)
        self.remove(temp)

        # sposta marker della lacuna sul bond nxt
        self.play(hole_marker.animate.move_to(nxt["mid"]), run_time=0.2, rate_func=smooth)

    # =========================================================
    # ATOMO BASE (uguale al tuo)
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
    # OSCILLAZIONE ELETTRONE DI LEGAME (tua logica, ripulita)
    # =========================================================
    def attach_bond_osc(self, dot, base_center, base_offset,
                        amplitude=0.03, omega=3.0, phase_shift=0.0, noise_amp=0.004):
        # evita doppio updater
        dot.clear_updaters()

        base_center = np.array(base_center)
        base_offset = np.array(base_offset)
        perp = base_offset.copy()
        if np.linalg.norm(perp[:2]) == 0:
            perp = np.array([1.0, 0.0, 0.0])
        perp = perp / np.linalg.norm(perp[:2])

        t = [0.0]

        def updater(mobj, dt):
            t[0] += dt
            osc = amplitude * np.sin(omega * t[0] + phase_shift)
            jitter = noise_amp * np.random.normal(size=3)
            jitter[2] = 0.0
            mobj.move_to(base_center + base_offset + osc * perp + jitter)

        dot.add_updater(updater)
        