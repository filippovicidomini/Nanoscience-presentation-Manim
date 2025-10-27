from manim import *
import numpy as np
import random

# set seed for reproducibility
random.seed(42)
np.random.seed(42)

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
    seg.set_opacity(0.05)
    return seg

def shrink_segment(p, q, gap):
    """
    Restituisce due punti p', q' sul segmento p->q, accorciato di 'gap' a ciascuna estremità.
    Se il segmento è troppo corto, ritorna p, q senza modifiche.
    """
    v = q - p
    L = np.linalg.norm(v)
    if L < 2*gap + 1e-2:
        return p, q
    u = v / L
    return p + gap*u, q - gap*u

def electron_on_bond(p, q, speed=0.9, radius=0.06, color=E_BOND_COLOR,
                     jitter_amp=0.12, jitter_freq=1.5,
                     core_gap=0.42):
    """
    Elettrone di legame che si muove lungo p<->q con piccola ondulazione trasversa,
    SENZA mai passare sui centri atomici grazie al 'core_gap' che accorcia il segmento.
    """
    # Accorcia il segmento per evitare i dischi atomici
    p_s, q_s = shrink_segment(p, q, core_gap)

    d = Dot(radius=radius, color=color)
    t = ValueTracker(np.random.rand())   # 0..1 lungo il segmento accorciato
    direction = {"d": 1}
    phase = random.random() * TAU
    osc_time = [0.0]

    v = q_s - p_s
    nv = np.linalg.norm(v)
    if nv < 1e-9:
        d.move_to(p_s)
        return d
    u = v / nv
    perp = np.array([-u[1], u[0], 0.0])  # perpendicolare nel piano 2D

    def upd(m, dt):
        # avanzamento con rimbalzo
        new_t = t.get_value() + direction["d"] * speed * dt
        if new_t > 1.0:
            new_t = 2.0 - new_t
            direction["d"] *= -1
        if new_t < 0.0:
            new_t = -new_t
            direction["d"] *= -1
        t.set_value(new_t)

        # tempo per l’ondulazione
        osc_time[0] += dt

        # smorzamento dell’ampiezza vicino alle estremità (evita “tocchi” ai cerchi)
        # s = 0 agli estremi, 1 al centro
        s = np.sin(np.pi * new_t)**1.2
        wiggle = (jitter_amp * s) * np.sin(TAU * jitter_freq * osc_time[0] + phase)

        pos = p_s + new_t * v + wiggle * perp
        m.move_to(pos)

    d.add_updater(upd)
    d.move_to(p_s + t.get_value() * (q_s - p_s))
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
        # Elettroni sui legami (due per legame, con parametri leggermente diversi)
        bond_electrons = VGroup()
        for (p, q) in edges:
        # uno in un verso, uno nell’altro; ampiezza/frequenza un po’ diverse
            e1 = electron_on_bond(p, q,
                          speed=0.9 + 0.15*random.random(),
                          radius=0.06,
                          color=E_BOND_COLOR,
                          jitter_amp=0.15 + 0.1*random.random(),
                          jitter_freq=1.2 + 0.8*random.random())
            e2 = electron_on_bond(q, p,
                          speed=0.9 + 0.15*random.random(),
                          radius=0.06,
                          color=E_BOND_COLOR,
                          jitter_amp=0.10 + 0.1*random.random(),
                          jitter_freq=1.2 + 0.8*random.random())
            bond_electrons.add(e1, e2)

        # Intro
        self.play(FadeIn(bonds), *[FadeIn(a) for a in atoms], FadeIn(title))
        self.play(FadeIn(bond_electrons))
        self.wait(0.6)

        # Inseriamo Fosforo (donor) in alto-centro (0,1)
        rP, cP = 2, 2
        phosphorus = make_atom(centers[rP][cP], label="P", color=PURE_BLUE)
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

        


        
        # fade out
        self.play(FadeOut(VGroup(atoms, bonds, bond_electrons, extra_e, title)))
        self.wait(0.5)
