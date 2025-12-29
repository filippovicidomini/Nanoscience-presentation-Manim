from manim import *
import numpy as np
import random

class PNJunctionForwardBiasWithDopants(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f14"

        # -----------------------------
        # PARAMETRI
        # -----------------------------
        rows = 5
        cols = 4
        spacing = 1.8
        
        # Centri dei due reticoli (già vicini)
        origin_n = LEFT * 1.8 
        origin_p = RIGHT * 1.8

        # -----------------------------
        # 1) SETUP AMBIENTE (Reticoli Base + Legami Interfaccia)
        # -----------------------------
        atom_grid_n, bonds_n, bond_group_n, lattice_n = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=origin_n
        )
        atom_grid_p, bonds_p, bond_group_p, lattice_p = self.build_lattice_with_bonds(
            rows, cols, spacing, origin=origin_p
        )

        # Creazione legami di interfaccia (ponte tra N e P)
        interface_pairs = VGroup()
        for i in range(rows):
            left_atom = atom_grid_n[i][cols - 1]
            right_atom = atom_grid_p[i][0]
            # Usiamo la tua funzione make_bond dinamica
            b = self.make_bond(left_atom, right_atom, (i, cols-1), (i, 0), sep=0.18)
            interface_pairs.add(b["pair"])

        # --- MODIFICA: AGGIUNTA ATOMI DROGANTI VISIBILI ---
        # Coordinate scelte per spargere i droganti
        dopants_n_coords = [(1, 1), (2, 2), (0, 2)]
        dopants_p_coords = [(3, 2), (3, 1), (4, 2)]

        # Applica stile donatori (Fosforo - Arancione) su N
        for (i, j) in dopants_n_coords:
            atom = atom_grid_n[i][j]
            p_nucleus = Circle(radius=0.28, color=ORANGE, fill_opacity=1.0).set_stroke(width=0)
            p_nucleus.move_to(atom.nucleus.get_center())
            # Sostituiamo il nucleo verde con quello arancione istantaneamente
            atom.nucleus.become(p_nucleus)
            # Aggiungo etichetta "P" per chiarezza
            p_label = Text("P", font_size=22, color=WHITE).move_to(atom.nucleus)
            atom.add(p_label)

        # Applica stile accettori (Boro - Viola) su P
        for (i, j) in dopants_p_coords:
            atom = atom_grid_p[i][j]
            b_nucleus = Circle(radius=0.22, color=PURPLE_C, fill_opacity=1.0).set_stroke(width=0)
            b_nucleus.move_to(atom.nucleus.get_center())
            # Sostituiamo nucleo
            atom.nucleus.become(b_nucleus)
            # Aggiungo etichetta "B"
            b_label = Text("B", font_size=22, color=WHITE).move_to(atom.nucleus)
            atom.add(b_label)
        # --------------------------------------------------


        # Aggiungiamo tutto alla scena staticamente (stato iniziale)
        semiconductor = VGroup(lattice_n, lattice_p, interface_pairs)
        self.add(semiconductor)

        label_n = Text("n-type (Maggioritari: e⁻)", font_size=24, color=BLUE_B).next_to(lattice_n, UP, buff=0.5)
        label_p = Text("p-type (Maggioritari: h⁺)", font_size=24, color=RED_B).next_to(lattice_p, UP, buff=0.5)
        self.add(label_n, label_p)

        # -----------------------------
        # 2) VISUALIZZAZIONE DEPLETION REGION INIZIALE
        # -----------------------------
        # All'equilibrio è larga
        depletion_width_initial = 2.5
        depletion = RoundedRectangle(
            corner_radius=0.25,
            width=depletion_width_initial,
            height=(rows - 1) * spacing + 1.5,
            stroke_width=2,
            stroke_color=YELLOW_B,
            fill_color=YELLOW,
            fill_opacity=0.15,
        ).move_to(ORIGIN)
        
        dep_label = Text("Barriera di Potenziale (Equilibrio)", font_size=24, color=YELLOW).next_to(depletion, DOWN, buff=0.25)
        
        self.play(FadeIn(depletion), FadeIn(dep_label), run_time=1.0)
        self.wait(0.5)

        # -----------------------------
        # 3) APPLICAZIONE TENSIONE (BATTERIA)
        # -----------------------------
        wire_y = -4.0
        contact_left = lattice_n.get_left() + LEFT * 0.5
        contact_right = lattice_p.get_right() + RIGHT * 0.5
        
        circuit = VMobject(color=WHITE, stroke_width=3)
        circuit.set_points_as_corners([
            contact_left + UP*1.5,
            contact_left,
            [contact_left[0], wire_y, 0],
            [contact_right[0], wire_y, 0],
            contact_right,
            contact_right + UP*1.5
        ])

        # Simbolo Batteria: FORWARD BIAS (+ a destra su P, - a sinistra su N)
        battery = self.get_battery_symbol().move_to([0, wire_y, 0])
        voltage_text = MathTex("V_{bias} > 0", color=GREEN).next_to(battery, UP, buff=0.2)
        
        self.play(
            Create(circuit),
            FadeIn(battery),
            Write(voltage_text),
            run_time=1.5
        )

        # -----------------------------
        # 4) EFFETTO: RIDUZIONE BARRIERA
        # -----------------------------
        new_width = 0.6
        
        self.play(
            depletion.animate.stretch_to_fit_width(new_width).set_fill(opacity=0.05),
            dep_label.animate.become(Text("Barriera Ridotta!", font_size=24, color=GREEN).move_to(dep_label.get_center())),
            run_time=2.0,
            rate_func=rate_functions.ease_in_out_cubic
        )
        self.wait(0.5)

        # -----------------------------
        # 5) INIEZIONE PORTATORI E CORRENTE (Spiegazione)
        # -----------------------------
        explanation_title = Text("Iniezione di Portatori Maggioritari", font_size=32, color=WHITE).to_edge(UP, buff=1)
        explanation_sub = Text("La barriera bassa permette il flusso massiccio di N→P e P→N", font_size=24, color=GRAY).next_to(explanation_title, DOWN)
        self.play(Write(explanation_title), Write(explanation_sub))

        # Etichette per il flusso
        e_flow_label = Text("Elettroni (N→P)", font_size=20, color=BLUE_C).move_to(UP*2 + LEFT*3)
        h_flow_label = Text("Lacune (P→N)", font_size=20, color=YELLOW_C).move_to(UP*2 + RIGHT*3)
        self.play(FadeIn(e_flow_label), FadeIn(h_flow_label))
        
        # Ciclo di animazione flusso
        for wave in range(4):
            electrons = VGroup()
            holes = VGroup()
            
            # Generazione ai lati opposti
            for k in range(3):
                y = (k - 1) * 1.2
                # Elettrone parte da N (sinistra)
                start_e = lattice_n.get_left() + RIGHT*0.5 + UP*y
                e_mob = self.make_free_electron(start_e)
                electrons.add(e_mob)
                
                # Lacuna parte da P (destra)
                start_h = lattice_p.get_right() + LEFT*0.5 + UP*y
                h_mob = self.make_hole(start_h)
                holes.add(h_mob)

            self.add(electrons, holes)
            
            # Movimento incrociato attraverso la giunzione
            self.play(
                electrons.animate.shift(RIGHT * 6.0), # Elettroni vanno a destra
                holes.animate.shift(LEFT * 6.0),      # Lacune vanno a sinistra
                run_time=2.0,
                rate_func=linear
            )
            
            # Ricombinazione (sparizione nella regione opposta)
            self.play(
                FadeOut(electrons),
                FadeOut(holes),
                run_time=0.2
            )

        # Conclusione
        current_arrow = Arrow(LEFT*2, RIGHT*2, color=YELLOW, buff=0).next_to(battery, DOWN, buff=0.5)
        current_text = Text("Grande Corrente Diretta (I)", font_size=28, color=YELLOW).next_to(current_arrow, DOWN)
        self.play(GrowArrow(current_arrow), Write(current_text))
        self.wait(3)

        # Pulizia finale
        for m in self.mobjects:
            m.clear_updaters()

    # =========================================================
    # FUNZIONI DI SUPPORTO (TUE, ORIGINALI)
    # =========================================================

    def get_battery_symbol(self):
        # Polo lungo (+) a destra, corto (-) a sinistra
        plate_neg = Line(UP*0.15, DOWN*0.15, color=WHITE, stroke_width=4).move_to(LEFT*0.15)
        plate_pos = Line(UP*0.3, DOWN*0.3, color=WHITE, stroke_width=4).move_to(RIGHT*0.15)
        wire_l = Line(LEFT*0.5, plate_neg.get_center(), color=WHITE)
        wire_r = Line(plate_pos.get_center(), RIGHT*0.5, color=WHITE)
        
        plus = Tex("+", color=RED, font_size=30).next_to(plate_pos, UP, buff=0.1)
        minus = Tex("-", color=BLUE, font_size=30).next_to(plate_neg, UP, buff=0.1)
        return VGroup(wire_l, wire_r, plate_neg, plate_pos, plus, minus)

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
        p1 = atom_a.get_center()
        p2 = atom_b.get_center()
        mid0 = 0.5 * (p1 + p2)

        e1 = Dot(mid0, radius=0.07, color="#6EA8FE")
        e2 = Dot(mid0, radius=0.07, color="#6EA8FE")

        self.attach_bond_osc_dynamic(e1, atom_a, atom_b, side=+1, sep=sep, phase_shift=0.0)
        self.attach_bond_osc_dynamic(e2, atom_a, atom_b, side=-1, sep=sep, phase_shift=PI)

        pair = VGroup(e1, e2)
        return {
            "a": a_idx, "b": b_idx, "atom_a": atom_a, "atom_b": atom_b,
            "e1": e1, "e2": e2, "pair": pair,
        }

    def make_valence_atom(self, orbit_radius=0.8) -> VGroup:
        nucleus = Circle(radius=0.22, color="#66E0A3", fill_opacity=1.0).set_stroke(width=0)
        orbit = Circle(radius=orbit_radius, stroke_width=2.0, stroke_opacity=0.85)
        orbit.set_stroke(color="#ffe680")

        directions = {"up": UP, "right": RIGHT, "down": DOWN, "left": LEFT}
        electrons = VGroup()
        for name, d in directions.items():
            e = Dot(point=orbit_radius * d, radius=0.07, color="#6EA8FE")
            electrons.add(e)

        atom = VGroup(orbit, nucleus, electrons)
        atom.orbit = orbit
        atom.nucleus = nucleus
        atom.electrons = electrons
        atom.core = VGroup(orbit, nucleus)
        return atom

    def attach_bond_osc_dynamic(self, dot, atom_a, atom_b, side=+1, sep=0.18, amplitude=0.03, omega=3.0, phase_shift=0.0, noise_amp=0.004):
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
        e = Dot(radius=0.085, color=RED).move_to(pos)
        glow1 = Circle(radius=0.18, stroke_width=2, color=BLUE_C, stroke_opacity=0.8).move_to(pos)
        glow2 = Circle(radius=0.28, stroke_width=2, color=BLUE_C, stroke_opacity=0.5).move_to(pos)

        def pulse(m, dt):
            if not hasattr(m, "t"): m.t = 0
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
        h = Circle(radius=0.11, stroke_width=4, color=YELLOW_C).move_to(pos)
        label = MathTex("h^+", font_size=26, color=YELLOW_C)
        label.add_updater(lambda m: m.next_to(h, UP, buff=0.12))
        return VGroup(h, label)