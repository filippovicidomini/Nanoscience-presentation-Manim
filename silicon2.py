from manim import *
import numpy as np
import random, math

#quality 

# set seed for reproducibility
random.seed(42)
np.random.seed(42)

NUCLEUS_COLOR = GRAY_B
SHELL_COLOR   = GRAY_D
E_CORE_COLOR  = BLUE_C    # elettroni K e L
E_VAL_COLOR   = YELLOW    # elettroni M (valenza)
BOND_COLOR    = WHITE

def electron_on_shell(center, R, speed=1.0, color=E_CORE_COLOR, radius=0.055, phase=0.0, start_active=False):
    dot = Dot(radius=radius, color=color)
    theta = ValueTracker(phase)
    running = [start_active]
    def upd(m, dt):
        if not running[0]:
            return
        dt = min(dt, 1/60)
        theta.increment_value(speed * dt * TAU)   # speed ~ giri/s
        a = theta.get_value()
        m.move_to(center + R*np.array([np.cos(a), np.sin(a), 0]))
    dot.add_updater(upd)
    dot.move_to(center + R*np.array([np.cos(phase), np.sin(phase), 0]))
    dot.running = running
    return dot

class MultiSiBonding2D(Scene):
    def construct(self):
        self.camera.background_color = "#111111"

        # ------- Parametri atomo -------
        R_K, R_L, R_M = 0.42, 0.80, 1.18
        r_nucleus = 0.17
        gap_to_bond = 0.15  # quanto fuori dal guscio M mettiamo la coppia condivisa

        # ------- Griglia di atomi 3x3 -------
        rows, cols = 5, 7
        dx, dy = 2.95, 2.95
        origin = np.array([-(cols-1)*dx/2, +(rows-1)*dy/2, 0])

        # Per ogni cella memorizziamo: nucleus, shells (K,L,M), e- K, e- L, e- M
        atoms = []
        for r in range(rows):
            row_atoms = []
            for c in range(cols):
                C = origin + np.array([c*dx, -r*dy, 0])

                # Nucleo
                nuc_fill = Circle(radius=r_nucleus, color=NUCLEUS_COLOR, stroke_width=0).move_to(C).set_opacity(0.25)
                nuc_ring = Circle(radius=r_nucleus, color=NUCLEUS_COLOR, stroke_width=2).move_to(C)
                nucleus  = VGroup(nuc_fill, nuc_ring)

                # Gusci
                sK = Circle(radius=R_K, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.5)
                sL = Circle(radius=R_L, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.5)
                sM = Circle(radius=R_M, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.6)

                # Elettroni K (2), L (8), M (4)
                eK = VGroup(*[
                    electron_on_shell(C, R_K, speed=0.8, color=E_CORE_COLOR, phase=ph, start_active=False)
                    for ph in [0, np.pi]
                ])
                eL = VGroup(*[
                    electron_on_shell(C, R_L, speed=0.6+0.1*random.random(), color=E_CORE_COLOR,
                                      phase=k*np.pi/4 + 0.1*random.random(), start_active=False)
                    for k in range(8)
                ])
                # quattro direzioni cardinali per la valenza (coerente con legami su griglia 2D)
                val_phases = [0, np.pi/2, np.pi, 3*np.pi/2]
                eM = VGroup(*[
                    electron_on_shell(C, R_M, speed=0.5, color=E_VAL_COLOR, phase=ph, start_active=False)
                    for ph in val_phases
                ])

                group = {
                    "center": C,
                    "nucleus": nucleus,
                    "shells": (sK, sL, sM),
                    "eK": eK,
                    "eL": eL,
                    "eM": eM,
                    "val_phases": val_phases
                }
                row_atoms.append(group)
            atoms.append(row_atoms)

        # ------ Disegna tutto (fase "tutti i gusci") ------
        # Anelli dei legami (visual), poi atomi
        shellsK = VGroup(*[atoms[r][c]["shells"][0] for r in range(rows) for c in range(cols)])
        shellsL = VGroup(*[atoms[r][c]["shells"][1] for r in range(rows) for c in range(cols)])
        shellsM = VGroup(*[atoms[r][c]["shells"][2] for r in range(rows) for c in range(cols)])
        nuclei  = VGroup(*[atoms[r][c]["nucleus"] for r in range(rows) for c in range(cols)])
        eK_all  = VGroup(*[atoms[r][c]["eK"] for r in range(rows) for c in range(cols)])
        eL_all  = VGroup(*[atoms[r][c]["eL"] for r in range(rows) for c in range(cols)])
        eM_all  = VGroup(*[atoms[r][c]["eM"] for r in range(rows) for c in range(cols)])

        self.play(FadeIn(nuclei), FadeIn(shellsK), FadeIn(shellsL), FadeIn(shellsM))
        self.play(FadeIn(eK_all), FadeIn(eL_all), FadeIn(eM_all))
        # attiva tutti gli updaters
        for r in range(rows):
            for c in range(cols):
                for d in list(atoms[r][c]["eK"]) + list(atoms[r][c]["eL"]) + list(atoms[r][c]["eM"]):
                    d.running[0] = True
        self.wait(5)

        # ------ Focus valenza: togli K e L su tutti gli atomi ------
        for r in range(rows):
            for c in range(cols):
                for d in list(atoms[r][c]["eK"]) + list(atoms[r][c]["eL"]):
                    d.running[0] = False
        self.play(FadeOut(shellsK, rate_func=linear), FadeOut(shellsL, rate_func=linear),
                  FadeOut(eK_all, rate_func=linear), FadeOut(eL_all, rate_func=linear),
                  run_time=2)
        self.wait(0.5)

        # ------ Forma i legami covalenti fra vicini (destra e giù per non duplicare) ------
        bonds = VGroup()
        used_valence = {}  # (r,c) -> set di indici valenza già usati (0:0°,1:90°,2:180°,3:270°)

        def pick_valence_index(direction_angle, phases):
            # scegli la fase più vicina alla direzione richiesta
            best_idx, best_score = None, 1e9
            for idx, ph in enumerate(phases):
                d = abs(((direction_angle - ph + np.pi) % (2*np.pi)) - np.pi)
                if d < best_score:
                    best_score, best_idx = d, idx
            return best_idx

        def reserve_valence(r, c, idx):
            used_valence.setdefault((r,c), set()).add(idx)

        def is_used(r, c, idx):
            return idx in used_valence.get((r,c), set())

        bond_pairs = []  # [(e1, e2, mid), ...] per micro-oscillazioni finali

        for r in range(rows):
            for c in range(cols):
                A = atoms[r][c]
                C_A = A["center"]
                phases_A = A["val_phases"]
                eM_A = list(A["eM"])

                # destra
                if c+1 < cols:
                    B = atoms[r][c+1]
                    C_B = B["center"]
                    v = C_B - C_A
                    ang = math.atan2(v[1], v[0])

                    idxA = pick_valence_index(ang, phases_A)
                    idxB = pick_valence_index(ang + np.pi, B["val_phases"])

                    # evita di riusare lo stesso e- per più legami
                    if not is_used(r, c, idxA) and not is_used(r, c+1, idxB):
                        reserve_valence(r, c, idxA)
                        reserve_valence(r, c+1, idxB)

                        eA = eM_A[idxA]; eB = list(B["eM"])[idxB]
                        # spegni updaters per animare la coppia
                        eA.running[0] = False; eB.running[0] = False

                        u = v / (np.linalg.norm(v) + 1e-9)
                        pA = C_A + (R_M + gap_to_bond) * u
                        pB = C_B - (R_M + gap_to_bond) * u
                        mid = (pA + pB) / 2

                        self.play(eA.animate.move_to(pA), eB.animate.move_to(pB), run_time=0.2)
                        seg = Line(pA, pB, color=BOND_COLOR, stroke_width=6).set_opacity(0.25)
                        self.play(Create(seg), run_time=0.05)
                        bonds.add(seg)
                        bond_pairs.append((eA, eB, mid))

                # giù
                if r+1 < rows:
                    B = atoms[r+1][c]
                    C_B = B["center"]
                    v = C_B - C_A
                    ang = math.atan2(v[1], v[0])

                    idxA = pick_valence_index(ang, phases_A)
                    idxB = pick_valence_index(ang + np.pi, B["val_phases"])

                    if not is_used(r, c, idxA) and not is_used(r+1, c, idxB):
                        reserve_valence(r, c, idxA)
                        reserve_valence(r+1, c, idxB)

                        eA = eM_A[idxA]; eB = list(B["eM"])[idxB]
                        eA.running[0] = False; eB.running[0] = False

                        u = v / (np.linalg.norm(v) + 1e-9)
                        pA = C_A + (R_M + gap_to_bond) * u
                        pB = C_B - (R_M + gap_to_bond) * u
                        mid = (pA + pB) / 2

                        self.play(eA.animate.move_to(pA), eB.animate.move_to(pB), run_time=0.2)
                        seg = Line(pA, pB, color=BOND_COLOR, stroke_width=6).set_opacity(0.25)
                        self.play(Create(seg), run_time=0.05)
                        bonds.add(seg)
                        bond_pairs.append((eA, eB, mid))

        # ------ Piccola vibrazione dei legami (qualitativa) ------
        for (e1, e2, mid) in bond_pairs:
            t = ValueTracker(random.random())
            amp = 0.1 + 0.05*random.random()
            freq = 1.0 + 0.9*random.random()
            vdir = e2.get_center() - e1.get_center()
            U = vdir / (np.linalg.norm(vdir) + 1e-9)
            N = np.array([-U[1], U[0], 0.0])
            def make_updater(dot, sign=+1, tt=t, A=amp, f=freq, u=U, n=N, MID=mid):
                def _upd(m, dt):
                    dt = min(dt, 1/60)
                    tt.increment_value(dt)
                    off = A * np.sin(TAU * f * tt.get_value())
                    m.move_to(MID + sign*(0.25*off)*u + (0.5*off)*n)
                return _upd
            e1.add_updater(make_updater(e1, +1))
            e2.add_updater(make_updater(e2, -1))

        self.wait(5)