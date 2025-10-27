from manim import *
import numpy as np
import math, random

# Colori
NUCLEUS_COLOR = GRAY_B
SHELL_COLOR   = GRAY_D
E_CORE_COLOR  = BLUE_C    # elettroni nei gusci interni (K,L)
E_VAL_COLOR   = YELLOW    # elettroni di valenza (M)
BOND_COLOR    = WHITE

# --- Helper: elettrone che orbita su un guscio circolare, con interruttore start/stop ---
def electron_on_shell(center, R, speed=1.0, color=E_CORE_COLOR, radius=0.06, phase=0.0, start_active=False):
    dot = Dot(radius=radius, color=color)
    theta = ValueTracker(phase)     # angolo attuale
    running = [start_active]        # interruttore ON/OFF

    def upd(m, dt):
        if not running[0]:
            return
        dt = min(dt, 1/60)         # evita primo passo enorme
        theta.increment_value(speed * dt * TAU)  # speed ≈ giri/sec
        ang = theta.get_value()
        pos = center + R * np.array([np.cos(ang), np.sin(ang), 0])
        m.move_to(pos)

    dot.add_updater(upd)
    dot.move_to(center + R * np.array([np.cos(phase), np.sin(phase), 0]))
    dot.running = running
    return dot

class SiliconValenceAndBonding2D(Scene):
    def construct(self):
        self.camera.background_color = "#111111"

        title = Text("Silicio: struttura elettronica e legami covalenti", font_size=34).to_edge(UP)

        # Centro atomo principale
        C = np.array([0.0, 0.2, 0.0])

        # Raggi gusci (K, L, M)
        R_K, R_L, R_M = 0.55, 1.05, 1.55

        # Nucleo
        nucleus_fill = Circle(radius=0.22, color=NUCLEUS_COLOR, stroke_width=0).move_to(C).set_opacity(0.25)
        nucleus_ring = Circle(radius=0.22, color=NUCLEUS_COLOR, stroke_width=3).move_to(C)
        nucleus_label = Text("Si", font_size=26).move_to(C)
        nucleus = VGroup(nucleus_fill, nucleus_ring, nucleus_label)

        # Gusci (anelli)
        shellK = Circle(radius=R_K, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.6)
        shellL = Circle(radius=R_L, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.6)
        shellM = Circle(radius=R_M, color=SHELL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.6)

        # Etichette
        labelK = Text("K (2)", font_size=20).next_to(shellK, LEFT, buff=0.15)
        labelL = Text("L (8)", font_size=20).next_to(shellL, LEFT, buff=0.15)
        labelM = Text("M (4 valenza)", font_size=20, color=E_VAL_COLOR).next_to(shellM, LEFT, buff=0.15)

        # Elettroni del Si centrale
        eK = VGroup(*[electron_on_shell(C, R_K, speed=0.7, color=E_CORE_COLOR, phase=ph, start_active=False)
                      for ph in [0, np.pi]])
        eL = VGroup(*[electron_on_shell(C, R_L, speed=0.6, color=E_CORE_COLOR, phase=k*np.pi/4, start_active=False)
                      for k in range(8)])
        eM_phases = [0, np.pi/2, np.pi, 3*np.pi/2]
        eM = VGroup(*[electron_on_shell(C, R_M, speed=0.5, color=E_VAL_COLOR, phase=ph, start_active=False)
                      for ph in eM_phases])

        # Pannello didascalie (destra)
        panel_w, panel_h = 4.2, 3.2
        panel_center = np.array([5.2, 0.6, 0])
        caption_panel = RoundedRectangle(corner_radius=0.15, width=panel_w, height=panel_h, color=WHITE)\
                           .set_fill(WHITE, opacity=0.06).set_stroke(opacity=0.6).move_to(panel_center)
        caption = Text("", font_size=24, line_spacing=0.9).move_to(panel_center)
        def set_caption(t):
            caption.become(Text(t, font_size=24, line_spacing=0.9).move_to(panel_center))

        # Fase 1: atomo singolo
        self.play(FadeIn(title))
        self.play(FadeIn(nucleus), FadeIn(shellK), FadeIn(shellL), FadeIn(shellM),
                  FadeIn(labelK), FadeIn(labelL), FadeIn(labelM))
        set_caption("Atomo di Silicio (Z=14)\nConfigurazione: [Ne] 3s² 3p²\nGusci: K=2, L=8, M=4 (valenza)")
        self.play(FadeIn(caption_panel), FadeIn(caption))

        self.play(FadeIn(eK), FadeIn(eL), FadeIn(eM))
        for d in list(eK) + list(eL) + list(eM):
            d.running[0] = True
        self.wait(1.0)

        # Fase 2: evidenzia valenza
        ring_val = Circle(radius=R_M+0.08, color=E_VAL_COLOR, stroke_width=2).move_to(C).set_stroke(opacity=0.6)
        tip_val = Text("Elettroni di valenza (4)", font_size=26, color=E_VAL_COLOR).next_to(shellM, DOWN, buff=0.25)
        self.play(Create(ring_val), FadeIn(tip_val))
        set_caption("I 4 elettroni di valenza (guscio M)\nsono responsabili dei legami covalenti.")
        self.wait(0.8)

        # Fase 3: porta 4 atomi vicini — ORA COMPLETI (K, L, M + elettroni)
        R_NEI = 2.6
        dirs = [45, 135, -135, -45]
        neighbors_pos = [C + R_NEI * np.array([np.cos(np.deg2rad(a)), np.sin(np.deg2rad(a)), 0]) for a in dirs]

        neighbors = VGroup()
        neighbors_shellsK = VGroup()
        neighbors_shellsL = VGroup()
        neighbors_shellsM = VGroup()
        neighbors_eK = VGroup()
        neighbors_eL = VGroup()
        neighbors_eM = VGroup()

        for P in neighbors_pos:
            # nucleo
            n_f = Circle(radius=0.18, color=NUCLEUS_COLOR, stroke_width=0).move_to(P).set_opacity(0.2)
            n_r = Circle(radius=0.18, color=NUCLEUS_COLOR, stroke_width=2).move_to(P)
            lab = Text("Si", font_size=20).move_to(P)
            # gusci
            sK = Circle(radius=R_K, color=SHELL_COLOR, stroke_width=2).move_to(P).set_stroke(opacity=0.5)
            sL = Circle(radius=R_L, color=SHELL_COLOR, stroke_width=2).move_to(P).set_stroke(opacity=0.5)
            sM = Circle(radius=R_M, color=SHELL_COLOR, stroke_width=2).move_to(P).set_stroke(opacity=0.5)
            neighbors_shellsK.add(sK)
            neighbors_shellsL.add(sL)
            neighbors_shellsM.add(sM)
            # elettroni K (2), L(8), M(4)
            eK_n = VGroup(*[electron_on_shell(P, R_K, speed=0.7, color=E_CORE_COLOR, phase=ph, start_active=False)
                            for ph in [0, np.pi]])
            eL_n = VGroup(*[electron_on_shell(P, R_L, speed=0.6, color=E_CORE_COLOR, phase=k*np.pi/4, start_active=False)
                            for k in range(8)])
            eM_n = VGroup(*[electron_on_shell(P, R_M, speed=0.5, color=E_VAL_COLOR, phase=ph, start_active=False)
                            for ph in [0, np.pi/2, np.pi, 3*np.pi/2]])
            neighbors_eK.add(eK_n)
            neighbors_eL.add(eL_n)
            neighbors_eM.add(eM_n)
            neighbors.add(VGroup(n_f, n_r, lab, sK, sL, sM, eK_n, eL_n, eM_n))

        self.play(*[FadeIn(n) for n in neighbors])
        # attiva moto di TUTTI gli elettroni dei vicini
        for grp in list(neighbors_eK) + list(neighbors_eL) + list(neighbors_eM):
            for d in grp:
                d.running[0] = True

        set_caption("Avvicinando atomi di Silicio completi (K,L,M),\nle orbite esterne si avvicinano: la valenza può condividere coppie.")
        self.wait(1.0)

        # Fase 4: FOCUS sulla valenza — TOGLIAMO K e L (centrale + vicini)
        # Spegni updaters K,L e fai fade out pulito (lascia M)
        for d in list(eK) + list(eL):
            d.running[0] = False
        for grp in list(neighbors_eK) + list(neighbors_eL):
            for d in grp:
                d.running[0] = False

        self.play(
            *[FadeOut(m, rate_func=linear) for m in [shellK, shellL, neighbors_shellsK, neighbors_shellsL, eK, eL, neighbors_eK, neighbors_eL]],
            run_time=0.8
        )
        set_caption("Focus sulla valenza: lasciamo SOLO i 4 elettroni M\n(per il centrale e per i vicini).")
        self.wait(0.6)

        # Fase 5: forma i 4 legami covalenti (come prima)
        R_NEI = 2.6
        neighbors_pos = [C + R_NEI * np.array([np.cos(np.deg2rad(a)), np.sin(np.deg2rad(a)), 0]) for a in dirs]

        def valence_e_pos_list(center, phases):
            return [center + R_M*np.array([np.cos(ph), np.sin(ph), 0]) for ph in phases]

        central_valence_electrons = list(eM)
        bonds_drawn = VGroup()
        bond_pairs  = []

        for k, P in enumerate(neighbors_pos):
            v = P - C
            ang_dir = math.atan2(v[1], v[0])
            # centrale: scegli e- M più allineato
            best_idx, best_score = None, 1e9
            for idx, ph in enumerate(eM_phases):
                d = abs(((ang_dir - ph + np.pi) % (2*np.pi)) - np.pi)
                if d < best_score:
                    best_score, best_idx = d, idx
            e_central = central_valence_electrons[best_idx]

            # vicino: prendi l'e- M più allineato verso il centro
            ang_to_center = math.atan2((C - P)[1], (C - P)[0])
            best_n_idx, best_n_score = None, 1e9
            phases_nei = [0, np.pi/2, np.pi, 3*np.pi/2]
            nei_group = neighbors_eM[k]
            for n_idx, ph in enumerate(phases_nei):
                d = abs(((ang_to_center - ph + np.pi) % (2*np.pi)) - np.pi)
                if d < best_n_score:
                    best_n_score, best_n_idx = d, n_idx
            e_neighbor = nei_group[best_n_idx]

            # Spegni updater prima di animare la coppia verso il legame
            e_central.running[0] = False
            e_neighbor.running[0] = False

            gap = 0.18
            u = v / np.linalg.norm(v)
            p_c_target = C + (R_M + gap) * u
            p_n_target = P - (R_M + gap) * u
            bond_mid   = (p_c_target + p_n_target) / 2

            bond_line = Line(p_c_target, p_n_target, color=BOND_COLOR, stroke_width=6).set_opacity(0.25)

            self.play(
                e_central.animate.move_to(p_c_target),
                e_neighbor.animate.move_to(p_n_target),
                run_time=0.8
            )
            self.play(Create(bond_line), run_time=0.5)

            bonds_drawn.add(bond_line)
            bond_pairs.append((e_central, e_neighbor, bond_mid))

        set_caption("Condivisione coppie di valenza → 4 legami covalenti Si–Si.")
        self.wait(1.0)

        # Fase 6: vibrazione qualitativa dei legami (fononi)
        for (e1, e2, mid) in bond_pairs:
            t = ValueTracker(random.random())
            amp = 0.05 + 0.03*random.random()
            freq = 1.2 + 0.6*random.random()
            vdir = e2.get_center() - e1.get_center()
            u = vdir / (np.linalg.norm(vdir) + 1e-9)
            n = np.array([-u[1], u[0], 0.0])

            def make_updater(dot, sign=+1, tt=t, A=amp, f=freq, U=u, N=n, MID=mid):
                def _upd(m, dt):
                    dt = min(dt, 1/60)
                    tt.increment_value(dt)
                    off = A * np.sin(TAU * f * tt.get_value())
                    pos = MID + sign * (0.25*off)*U + (0.5*off)*N
                    m.move_to(pos)
                return _upd

            e1.add_updater(make_updater(e1, +1))
            e2.add_updater(make_updater(e2, -1))

        set_caption("I legami vibrano (fononi): qui li rendiamo visibili in modo qualitativo.")
        final = Text("Si: 2, 8, 4 → i 4 di valenza formano 4 legami covalenti", font_size=26).to_edge(DOWN)
        self.play(FadeIn(final))
        self.wait(1.2)