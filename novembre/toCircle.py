from manim import *
import numpy as np

class SiliconLatticeSpin(Scene):
    def construct(self):
        self.camera.background_color = "#0b0f17"

        # -----------------------------
        # 1) ATOMO SINGOLO CHE RUOTA
        # -----------------------------
        center = ORIGIN
        central_atom = self.make_valence_atom()
        #central_atom.move_to(center)

        self.start_spin(central_atom, omega_spin=2.5)

        self.play(FadeIn(central_atom, scale=1), run_time=1.5)
        #self.play(central_atom.animate.scale(1), run_time=1.0)
        self.wait(6.5)

        # -----------------------------
        # 2) RETICOLO 2D PIÙ GRANDE
        # -----------------------------
        rows = 5
        cols = 7
        spacing = 2.0

        atom_grid = [[None for _ in range(cols)] for _ in range(rows)]
        new_atoms_group = VGroup()

        for i in range(rows):
            for j in range(cols):
                x = (j - (cols - 1)/2) * spacing
                y = (i - (rows - 1)/2) * spacing
                pos = np.array([x, y, 0.0])

                if np.allclose(pos, center):
                    atom = central_atom
                else:
                    atom = self.make_valence_atom()
                    atom.move_to(pos)
                    self.start_spin(atom, omega_spin=2.5)
                    new_atoms_group.add(atom)

                atom_grid[i][j] = atom

        # comparsa graduale degli altri atomi
        self.play(
            LaggedStart(
                *[FadeIn(a, scale=0.9) for a in new_atoms_group],
                lag_ratio=0.08,
                run_time=2.5
            )
        )
        self.wait(5.5)

        # -----------------------------
        # 3) FORMAZIONE DEI LEGAMI
        #    (blob giallini + rimozione
        #     degli elettroni usati)
        # -----------------------------
        bond_electrons = VGroup()

        def make_bond_pair(p1, p2):
            """
            Crea due elettroni di legame attorno al punto medio tra p1 e p2.
            Restano vicini (separazione fissa), ognuno con la sua oscillazione armonica.
            """
            mid = 0.5 * (p1 + p2)

            # direzione lungo il legame
            dir_vec = p2 - p1
            dir_vec[2] = 0
            norm = np.linalg.norm(dir_vec[:2]) or 1.0
            dir_vec /= norm

            # direzione perpendicolare (dove mettiamo i due elettroni)
            perp = np.array([-dir_vec[1], dir_vec[0], 0.0])

            sep = 0.18  # distanza dal centro
            pos1 = mid + sep * perp
            pos2 = mid - sep * perp

            e1 = Dot(pos1, radius=0.07, color="#6EA8FE")
            e2 = Dot(pos2, radius=0.07, color="#6EA8FE")

            # oscillazioni armoniche attorno alle posizioni di equilibrio
            self.add_oscillation_to_bond_electron(e1, base_center=mid, base_offset=sep*perp,
                                                amplitude=0.03, omega=3.0, phase_shift=0.0)
            self.add_oscillation_to_bond_electron(e2, base_center=mid, base_offset=-sep*perp,
                                                amplitude=0.03, omega=3.0, phase_shift=PI)

            return VGroup(e1, e2)

        # Legami orizzontali: RIGHT del sinistro, LEFT del destro
        for i in range(rows):
            for j in range(cols - 1):
                a_left = atom_grid[i][j]
                a_right = atom_grid[i][j+1]

                e_left = a_left.e_by_dir["right"]
                e_right = a_right.e_by_dir["left"]

              

                pair = make_bond_pair(a_left.get_center(), a_right.get_center())

                self.play(
                    FadeIn(pair, scale=0.4),
                    run_time=0.2
                )

                bond_electrons.add(pair)

        # Legami verticali: DOWN del sopra, UP del sotto
        for i in range(rows - 1):
            for j in range(cols):
                a_top = atom_grid[i][j]
                a_bottom = atom_grid[i+1][j]

                e_top = a_top.e_by_dir["down"]
                e_bottom = a_bottom.e_by_dir["up"]

               

                pair = make_bond_pair(a_top.get_center(), a_bottom.get_center())

                self.play(
                    FadeIn(pair, scale=0.4),
                    run_time=0.2
                )

                
        # Dopo aver creato TUTTI i legami (orizzontali e verticali)
        # rimuoviamo tutti gli elettroni dal reticolo in un colpo solo
        all_electrons = VGroup()
        for i in range(rows):
            for j in range(cols):
                atom_grid[i][j].electrons.clear_updaters()
                all_electrons.add(atom_grid[i][j].electrons)

        self.play(FadeOut(all_electrons), run_time=1.0)
        self.wait(4.0)

        # -----------------------------
        # 4) HIGHLIGHT FINALE
        #    SOLO NUCLEI+ORBITE + LEGAMI
        # -----------------------------
        cores = VGroup()
        for i in range(rows):
            for j in range(cols):
                cores.add(atom_grid[i][j].core)

        full_lattice = VGroup(cores, bond_electrons)

        self.play(Indicate(full_lattice, scale_factor=1.02), run_time=2.0)
        self.wait(10.0)

    # -----------------------------
    # HELPER: atomo con 4 elettroni di valenza
    # -----------------------------
    def make_valence_atom(self, orbit_radius=0.8) -> VGroup:
        # nucleo
        nucleus = Circle(
            radius=0.22,
            color="#66E0A3",
            fill_opacity=1.0
        ).set_stroke(width=0)

        # orbita giallina
        orbit = Circle(
            radius=orbit_radius,
            stroke_width=2.0,
            stroke_opacity=0.85,
        )
        orbit.set_stroke(color="#ffe680")  # giallino

        # 4 elettroni di valenza
        directions = {
            "up": UP,
            "right": RIGHT,
            "down": DOWN,
            "left": LEFT,
        }
        electrons = VGroup()
        e_by_dir = {}
        for name, d in directions.items():
            pos = orbit_radius * d
            e = Dot(point=pos, radius=0.07, color="#6EA8FE")
            electrons.add(e)
            e_by_dir[name] = e

        atom = VGroup(orbit, nucleus, electrons)
        atom.orbit = orbit
        atom.nucleus = nucleus
        atom.electrons = electrons
        atom.e_by_dir = e_by_dir
        atom.core = VGroup(orbit, nucleus)  # per highlight finale (senza elettroni)
        return atom

    # -----------------------------
    # HELPER: spin continuo degli elettroni
    # -----------------------------
    def start_spin(self, atom: VGroup, omega_spin: float = 1.5,
               osc_amp: float = 0.03, omega_osc: float = 3.0):
        """
        Fa ruotare gli elettroni di valenza intorno all'atomo
        e aggiunge una leggera oscillazione radiale armonica
        centrata nel tempo (niente drift).
        """
        eg = atom.electrons
        base = eg.copy()  # geometria di riferimento
        t = [0.0]
        phase = np.random.uniform(0, TAU)

        def updater(mobj, dt,
                    a=atom,
                    base_geom=base,
                    w_spin=omega_spin,
                    w_osc=omega_osc,
                    amp=osc_amp,
                    ph=phase):
            t[0] += dt
            # fattore di scala radiale attorno a 1
            scale = 1.0 + amp * np.sin(w_osc * t[0] + ph)
            angle = w_spin * t[0]

            # ricostruiamo ogni volta dalla geometria base
            new = base_geom.copy()
            new.scale(scale, about_point=ORIGIN)
            new.rotate(angle, about_point=ORIGIN)
            new.move_to(a.get_center())
            mobj.become(new)

        eg.add_updater(updater)

    def stop_spin(self, atom: VGroup):
        """Se ti serve, ferma la rotazione."""
        atom.electrons.clear_updaters()

    def add_oscillation_to_blob(self, blob: Mobject, amplitude=0.02, omega=2.0):
        """
        Oscillazione armonica 2D attorno al centro originale del blob.
        Resta sempre centrato, niente drift.
        """
        base_center = blob.get_center().copy()
        phase_x = np.random.uniform(0, TAU)
        phase_y = np.random.uniform(0, TAU)
        t = [0.0]

        def updater(mobj, dt,
                    amp=amplitude,
                    w=omega,
                    px=phase_x,
                    py=phase_y,
                    base=base_center):
            t[0] += dt
            offset = np.array([
                amp * np.sin(w * t[0] + px),
                amp * np.sin(w * t[0] + py),
                0.0
            ])
            mobj.move_to(base + offset)

        blob.add_updater(updater)
        
    def add_oscillation_to_bond_electron(
        self,
        dot: Mobject,
        base_center: np.ndarray,
        base_offset: np.ndarray,
        amplitude: float = 0.03,
        omega: float = 3.0,
        phase_shift: float = 0.0,
        noise_amp: float = 0.01,
    ):
        """
        Oscillazione armonica per un singolo elettrone di legame.
        Rimane in media vicino a base_center + base_offset (niente drift),
        con piccola oscillazione lungo la direzione perpendicolare al legame.
        """
        base_center = np.array(base_center)
        base_offset = np.array(base_offset)
        t = [0.0]

        # direzione di oscillazione = stessa direzione di base_offset, normalizzata
        perp = base_offset.copy()
        if np.linalg.norm(perp[:2]) == 0:
            perp = np.array([1.0, 0.0, 0.0])
        perp = perp / np.linalg.norm(perp[:2])

        def updater(mobj, dt,
                    amp=amplitude,
                    w=omega,
                    ph=phase_shift,
                    center=base_center,
                    off=base_offset,
                    direction=perp,
                    noise=noise_amp):
            t[0] += dt
            osc = amp * np.sin(w * t[0] + ph)
            jitter = noise * np.random.normal(size=3)
            jitter[2] = 0.0
            new_pos = center + off + osc * direction + jitter
            mobj.move_to(new_pos)

        dot.add_updater(updater)