from manim import *
import numpy as np

SI_COLOR = GRAY_B
B_COLOR = RED
E_COLOR = BLUE_B
HOLE_COLOR = PINK

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

def electron_on_bond(p, q, speed=0.8, radius=0.06, color=E_COLOR):
    """ Elettrone che si muove avanti/indietro lungo il legame p<->q. """
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
    # posiziona iniziale
    d.move_to((1-t.get_value())*p + t.get_value()*q)
    return d

class SiBDopingHopping2D(Scene):
    def construct(self):
        self.camera.background_color = "#111111"

        # --- piccola griglia 2x3 ---
        rows, cols = 2, 3
        dx, dy = 2.6, 2.0
        x0, y0 = -dx, +0.6
        centers = [[np.array([x0 + c*dx, y0 - r*dy, 0]) for c in range(cols)] for r in range(rows)]

        title = Text("Silicio: moto degli elettroni lungo i legami • poi doping con B (p-type)", font_size=32).to_edge(UP)

        # Atomi Si e legami (solo primi vicini orizz/vert)
        atoms = VGroup()
        bonds = VGroup()
        for r in range(rows):
            for c in range(cols):
                atoms.add(make_atom(centers[r][c]))

        # legami
        edges = []  # lista di (p,q) per mettere elettroni sopra
        for r in range(rows):
            for c in range(cols):
                p = centers[r][c]
                if c+1 < cols:
                    q = centers[r][c+1]
                    bonds.add(bond(p, q)); edges.append((p, q))
                if r+1 < rows:
                    q = centers[r+1][c]
                    bonds.add(bond(p, q)); edges.append((p, q))

        # Elettroni che “corrono” lungo i legami (visualizzazione del legame covalente)
        electrons = VGroup()
        for (p, q) in edges:
            # due elettroni condivisi per legame → mettiamo 2 dot sfalsati
            e1 = electron_on_bond(p, q, speed=0.9, radius=0.06)
            e2 = electron_on_bond(q, p, speed=0.9, radius=0.06)  # direzione opposta per sfalsare
            electrons.add(e1, e2)

        # Legenda (griglia) in pannello
        legend = VGroup(
            VGroup(Dot(radius=0.09, color=E_COLOR), Text("e⁻ (elettrone)", font_size=22, color=E_COLOR)).arrange(RIGHT, buff=0.2),
            VGroup(Dot(radius=0.09, color=HOLE_COLOR), Text("lacuna (buco)", font_size=22, color=HOLE_COLOR)).arrange(RIGHT, buff=0.2),
        ).arrange_in_grid(rows=2, cols=1, buff=0.2, col_alignments=["l"]).to_edge(DOWN).shift(0.15*DOWN)
        panel = SurroundingRectangle(legend, color=WHITE, buff=0.25, corner_radius=0.15)
        panel.set_opacity(0.15).set_fill(WHITE, opacity=0.06)
        legend_block = VGroup(panel, legend)

        # Intro
        self.play(FadeIn(bonds), *[FadeIn(a) for a in atoms], FadeIn(title))
        self.play(FadeIn(electrons), FadeIn(legend_block))
        self.wait(1.0)

        # --- Inseriamo Boro al centro della griglia (r=0, c=1) per renderlo visibile in mezzo ---
        rB, cB = 0, 1
        boron = make_atom(centers[rB][cB], label="B", color=B_COLOR)
        explain1 = Text("Dopiamo con B: 3 elettroni di valenza → un legame rimane “in deficit”", font_size=26)\
                   .next_to(title, DOWN, buff=0.35)
        self.play(Transform(atoms[rB*cols + cB], boron), FadeIn(explain1))

        # Scegli il legame B–vicino sinistro come “in deficit”
        pB = centers[rB][cB]
        pL = centers[rB][cB-1] if cB-1 >= 0 else pB + LEFT*1.3
        deficit_edge = (pB, pL)

        # Trova un paio di elettroni che stanno su quel legame e rimuovine UNO (visualizza la carenza)
        # Cerca nel gruppo "electrons" quelli più vicini al segmento (pB,pL)
        def seg_dist(m, A, B):
            P = m.get_center()
            AB = B - A
            t = np.clip(np.dot(P-A, AB)/np.dot(AB, AB), 0, 1)
            Q = A + t*AB
            return np.linalg.norm(P - Q)

        # pick: l'elettrone più vicino al legame B-L
        e_candidates = sorted([e for e in electrons], key=lambda e: seg_dist(e, pB, pL))
        if e_candidates:
            e_remove = e_candidates[0]
            # fermiamo il suo updater e dissolviamo
            for up in list(e_remove.updaters):
                e_remove.remove_updater(up)
            self.play(e_remove.animate.set_opacity(0.0), run_time=0.5)
            self.remove(e_remove)

        # Disegna la lacuna su quel legame
        hole_pos = pB + 0.48*(pL - pB)
        hole = Circle(radius=0.12, color=HOLE_COLOR, stroke_width=6).move_to(hole_pos)
        hole.set_opacity(0.95)
        glow = Circle(radius=0.18, color=HOLE_COLOR, stroke_width=0).move_to(hole_pos)
        glow.set_opacity(0.28)
        hole_lbl = Text("buco", font_size=20, color=HOLE_COLOR).next_to(hole, DOWN, buff=0.08)
        self.play(FadeIn(glow), FadeIn(hole), FadeIn(hole_lbl))
        self.wait(0.4)

        # --- Un elettrone dal legame a sinistra riempie il deficit → il buco si sposta ---
        explain2 = Text("Un e⁻ vicino riempie il legame con B → il buco “si sposta” indietro", font_size=26)\
                   .next_to(explain1, DOWN, buff=0.2)
        self.play(FadeIn(explain2))

        # prendi un elettrone da (L–LL) se esiste, altrimenti da (L–sopra) o (L–sotto)
        candidate_neighbors = []
        if cB-2 >= 0:
            candidate_neighbors.append((centers[rB][cB-1], centers[rB][cB-2]))
        if rB-1 >= 0:
            candidate_neighbors.append((centers[rB][cB-1], centers[rB-1][cB-1]))
        if rB+1 < rows:
            candidate_neighbors.append((centers[rB][cB-1], centers[rB+1][cB-1]))

        # trova un elettrone che stia su uno di questi legami
        donor_e = None
        donor_edge = None
        for (P, Q) in candidate_neighbors:
            cand = sorted([e for e in electrons], key=lambda e: seg_dist(e, P, Q))
            for e in cand:
                # se sufficent. vicino al segmento, prendilo
                if seg_dist(e, P, Q) < 0.25:
                    donor_e = e; donor_edge = (P, Q); break
            if donor_e is not None:
                break

        if donor_e is not None:
            # stacca updater, muovi l'elettrone verso la lacuna (riempie)
            for up in list(donor_e.updaters):
                donor_e.remove_updater(up)
            self.play(donor_e.animate.move_to(hole.get_center()), run_time=0.6)
            # la lacuna si sposta lungo la catena (verso il bordo)
            # nuovo target del buco: più indietro sul legame del donatore
            P, Q = donor_edge
            new_hole = P + 0.48*(Q - P)
            self.play(hole.animate.move_to(new_hole), glow.animate.move_to(new_hole), run_time=0.6)
            # l'elettrone che ha riempito scompare (ora legame B-L è completo)
            self.play(donor_e.animate.set_opacity(0.0), run_time=0.25)
            self.remove(donor_e)

        # Campo elettrico + piccolo drift coerente verso sinistra
        efield = Arrow(LEFT*5.8 + UP*3.1, RIGHT*5.8 + UP*3.1, buff=0).set_stroke(width=6)
        efield_lbl = Text("Campo elettrico", font_size=22).next_to(efield, DOWN, buff=0.12)
        self.play(FadeIn(efield), FadeIn(efield_lbl))

        # drift collettivo: sposta leggermente gli elettroni lungo i legami verso sinistra
        drift = 0.25*LEFT
        movers = [e for e in electrons if e in self.mobjects]  # quelli non rimossi
        self.play(*[e.animate.shift(drift) for e in movers],
                  hole.animate.shift(0.4*LEFT),
                  glow.animate.shift(0.4*LEFT),
                  run_time=1.1)

        final = Text("B crea lacune → la loro propagazione genera corrente p-type", font_size=26).to_edge(DOWN)
        self.play(FadeOut(explain1), FadeOut(explain2))
        self.play(FadeIn(final))
        self.wait(1.0)