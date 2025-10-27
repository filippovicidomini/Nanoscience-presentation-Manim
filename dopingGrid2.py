from manim import *
import numpy as np
import random

# Palette
SI_COLOR = GRAY_B
P_COLOR  = BLUE_B
E_BOND_COLOR  = BLUE_C   # elettroni di legame
E_FREE_COLOR  = YELLOW   # elettrone donato (quasi libero)

def make_atom(center, label="Si", color=SI_COLOR, r=0.35):
    ring = Circle(radius=r, color=color, stroke_width=4).move_to(center)
    fill = Circle(radius=r, color=color, stroke_width=0).move_to(center)
    fill.set_opacity(0.15)
    txt  = Text(label, font_size=22).move_to(center)
    return VGroup(fill, ring, txt)

def bond(p, q):
    seg = Line(p, q, color=WHITE, stroke_width=5)
    seg.set_opacity(0.25)
    return seg

def electron_on_bond(p, q, speed=2, radius=0.06, color=E_BOND_COLOR):
    """Elettrone che rimbalza avanti/indietro sul legame p<->q (moto di legame)."""
    d = Dot(radius=radius, color=color)
    t = ValueTracker(np.random.rand())
    direction = {"d": 1}
    def upd(m, dt):
        new_t = t.get_value() + direction["d"] * speed * dt
        if new_t > 1.0:
            new_t = 2.0 - new_t
            direction["d"] *= -1
        if new_t < 0.0:
            new_t = -new_t
            direction["d"] *= -1
        t.set_value(new_t)
        pos = (1-new_t) * p + new_t * q
        m.move_to(pos)
    d.add_updater(upd)
    d.move_to((1-t.get_value())*p + t.get_value()*q)
    return d

class NTypeNoFieldAttracted2D(Scene):
    def construct(self):
        self.camera.background_color = "#111111"

        # ----- Layout ---------------------------------------------------------
        TITLE_Y = 3.5
        GRID_CENTER_Y = -0.3

        title = Text("n-type (Fosforo): moto termico con attrazione verso gli atomi (nessun campo)", font_size=30)\
                .move_to(np.array([0, TITLE_Y, 0]))

        # Griglia 3x3
        rows, cols = 5, 5
        dx, dy = 1.8, 1.4
        grid_width = (cols - 1) * dx
        grid_height = (rows - 1) * dy
        x0 = -0.5 * grid_width
        y0 = GRID_CENTER_Y + 0.5 * grid_height
        # Limiti: rettangolo esattamente attorno alla griglia (con un piccolo margine)
        GRID_PAD = 0.15
        grid_min_x = x0 - GRID_PAD
        grid_max_x = x0 + grid_width + GRID_PAD
        grid_max_y = y0 + GRID_PAD
        grid_min_y = y0 - grid_height - GRID_PAD
        
        BOUNDS = {
            "x_min": grid_min_x,
            "x_max": grid_max_x,
            "y_min": grid_min_y,
            "y_max": grid_max_y,
            }       


        centers = [[np.array([x0 + c*dx, y0 - r*dy, 0]) for c in range(cols)] for r in range(rows)]
        centers_flat = [centers[r][c] for r in range(rows) for c in range(cols)]

        atoms = VGroup()
        bonds = VGroup()
        edges = []
        for r in range(rows):
            for c in range(cols):
                atoms.add(make_atom(centers[r][c]))
        for r in range(rows):
            for c in range(cols):
                p = centers[r][c]
                if c+1 < cols:
                    q = centers[r][c+1]
                    bonds.add(bond(p, q)); edges.append((p, q))
                if r+1 < rows:
                    q = centers[r+1][c]
                    bonds.add(bond(p, q)); edges.append((p, q))

        # Elettroni sui legami
        bond_electrons = VGroup()
        for (p, q) in edges:
            bond_electrons.add(
                electron_on_bond(p, q, speed=0.9, radius=0.06, color=E_BOND_COLOR),
                electron_on_bond(q, p, speed=0.9, radius=0.06, color=E_BOND_COLOR)
            )

        # Intro
        self.play(FadeIn(bonds), *[FadeIn(a) for a in atoms], FadeIn(title))
        self.play(FadeIn(bond_electrons))
        self.wait(0.6)

        # Inseriamo Fosforo (donor) in alto-centro (0,1)
        rP, cP = 2, 2
        phosphorus = make_atom(centers[rP][cP], label="P", color=P_COLOR)
        self.play(Transform(atoms[rP*cols + cP], phosphorus))

        # Elettrone extra vicino a P
        extra_e = Dot(radius=0.085, color=E_FREE_COLOR).move_to(centers[rP][cP] + 0.45*UP)
        self.play(FadeIn(extra_e, scale=0.9))

        # -------- Moto "quasi libero" con attrazione verso gli atomi ----------
        # Modello: v' = -gamma*v + k*(x_near - x) + sigma*noise  (con anti-collisione a corto raggio)
        v = np.array([-1.2, 1.0, 0.0])  # velocità iniziale
        gamma = 0.1     # smorzamento
        k_attr = 2.0    # forza di attrazione verso l'atomo più vicino
        sigma = 10.2     # ampiezza rumore (termico)
        r_rep = 1    # raggio di repulsione (evita il nucleo)
        k_rep = 10.0     # intensità repulsione a corto raggio
        v_max = 3.0     # limite superiore velocità

        def nearest_center(x, centers_list):
            d2 = [(np.linalg.norm(x - c), c) for c in centers_list]
            d2.sort(key=lambda t: t[0])
            return d2[0]  # (dist, center)

        def update_free(m: Mobject, dt):
            if dt <= 0: 
                return
            pos = m.get_center()
            dist, c_near = nearest_center(pos, centers_flat)
            # Attrazione verso il centro più vicino
            a_attr = k_attr * (c_near - pos)

            # Repulsione se troppo vicino
            a_rep = np.array([0.0, 0.0, 0.0])
            if dist < r_rep:
                dir_out = (pos - c_near)
                n = np.linalg.norm(dir_out)
                if n > 1e-6:
                    dir_out /= n
                    # più vicino → più repulsione
                    a_rep = k_rep * (r_rep - dist) * dir_out

            # Rumore termico (gaussiano)
            noise = np.array([np.random.randn(), np.random.randn(), 0.0]) * sigma

            # Aggiorna velocità con smorzamento
            nonlocal v
            v += (-gamma * v + a_attr + a_rep + noise) * dt

            # Limita la velocità
            speed = np.linalg.norm(v)
            if speed > v_max:
                v *= (v_max / max(speed, 1e-9))

            # Nuova posizione
            new_pos = pos + v * dt

            # Bordi del frame (rimbalzo morbido)
            # Bordi della GRIGLIA: riflessione elastica e rientro leggero
            if new_pos[0] < BOUNDS["x_min"]:
                v[0] = abs(v[0]) * 0.8
                new_pos[0] = BOUNDS["x_min"] + 1e-3
            elif new_pos[0] > BOUNDS["x_max"]:
                v[0] = -abs(v[0]) * 0.8
                new_pos[0] = BOUNDS["x_max"] - 1e-3

            if new_pos[1] < BOUNDS["y_min"]:
                v[1] = abs(v[1]) * 0.8
                new_pos[1] = BOUNDS["y_min"] + 1e-3
            elif new_pos[1] > BOUNDS["y_max"]:
                v[1] = -abs(v[1]) * 0.8
                new_pos[1] = BOUNDS["y_max"] - 1e-3

            # Evita di salire troppo vicino al titolo
            if new_pos[1] > TITLE_Y - 0.6:
                new_pos[1] = TITLE_Y - 0.6
                v[1] = -abs(v[1])  # rimbalzo verso il basso

            m.move_to(new_pos)

        # Sposta l'elettrone un po' sopra P e attiva updater
        extra_e.move_to(centers[rP][cP] + np.array([0.18, 0.55, 0]))
        extra_e.add_updater(update_free)

        # Lascia evolvere il moto per un po'
        self.wait(10.0)

        # Nota finale
        final = Text("Senza campo: l’elettrone donato resta vicino al reticolo\n(moto termico con lieve attrazione ai siti) → corrente netta ≈ 0", font_size=26)\
                .to_edge(DOWN)
        self.play(FadeIn(final))
        self.wait(1.6)
