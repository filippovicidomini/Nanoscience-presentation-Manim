from manim import *
import numpy as np

# Palette
SI_COLOR = GRAY_B
B_COLOR = RED
E_COLOR = BLUE_B
HOLE_COLOR = PINK

def make_atom(center, radius=0.45, color=SI_COLOR, label_text="Si"):
    atom = Circle(radius=radius, color=color, stroke_width=4).move_to(center)
    fill = Circle(radius=radius, color=color, stroke_width=0).move_to(center)
    fill.set_opacity(0.15)
    label = Text(label_text, font_size=28, weight=BOLD).move_to(center)
    return VGroup(fill, atom, label)

def make_orbiters(center, n=3, R=0.75, speed=1.0, color=E_COLOR, size=0.08):
    """
    Crea n elettroni che orbitano attorno al centro su cerchi di raggio R.
    Usa updaters per rotazione continua.
    """
    dots = VGroup()
    for k in range(n):
        theta0 = 2*np.pi * k / n
        tracker = ValueTracker(theta0)

        def make_dot(tr=tracker):
            return always_redraw(lambda: Dot(
                point = center + R*np.array([np.cos(tr.get_value()), np.sin(tr.get_value()), 0]),
                radius = size,
                color = color
            ))

        d = make_dot(tracker)
        # updater di rotazione costante
        d.add_updater(lambda m, dt, tr=tracker: tr.increment_value(speed*dt))
        dots.add(d)
    return dots

def bond(p, q):
    seg = Line(p, q, stroke_width=6, color=WHITE)
    seg.set_opacity(0.25)
    return seg

class PTypeDoping2D(Scene):
    def construct(self):
        self.camera.background_color = "#111111"

        # --- Layout: 4 atomi in riga
        xs = [-4, -1.3, 1.4, 4.1]
        y = 0
        centers = [np.array([x, y, 0]) for x in xs]

        # --- Stato iniziale: Silicio intrinseco (tutti "Si")
        atoms = VGroup(*[make_atom(c) for c in centers])
        bonds = VGroup(*[
            bond(centers[i], centers[i+1]) for i in range(len(centers)-1)
        ])

        # Elettroni: 3 orbitali per atomo (astrazione didattica)
        orbiters = VGroup(*[
            make_orbiters(c, n=3, R=0.75, speed=1.2) for c in centers
        ])

        title = Text("Silicio (intrinseco): elettroni in equilibrio dinamico", font_size=32)\
                .to_edge(UP)

        legend = VGroup(
    VGroup(
        Dot(radius=0.09, color=E_COLOR),
        Text("e⁻ (elettrone)", font_size=22, color=E_COLOR),
    ).arrange(RIGHT, buff=0.2),
    VGroup(
        Dot(radius=0.09, color=HOLE_COLOR),
        Text("lacuna (buco)", font_size=22, color=HOLE_COLOR),
    ).arrange(RIGHT, buff=0.2),
).arrange(DOWN, aligned_edge=LEFT, buff=0.15).to_edge(DOWN).shift(0.2*DOWN)

        self.play(FadeIn(bonds), *[FadeIn(a) for a in atoms], FadeIn(title))
        for grp in orbiters:
            for d in grp:
                d.set_z_index(3)
        self.play(*[FadeIn(grp) for grp in orbiters], FadeIn(legend))
        self.wait(1.0)

        # --- Introduzione del boro (p-type): sostituisco l'ultimo atomo
        boron = make_atom(centers[-1], color=B_COLOR, label_text="B")
        explain1 = Text("Dopiamo con Boro: B ha 3 e⁻ di valenza → manca un legame", font_size=26)\
                   .next_to(title, DOWN, buff=0.4)
        self.play(Transform(atoms[-1], boron), FadeIn(explain1))
        self.wait(0.6)

        # Rimuovo 1 elettrone dall'orbita dell'atomo B per visualizzare la "carenza"
        b_orbiters = orbiters[-1]
        if len(b_orbiters) > 0:
            removed_e = b_orbiters[-1]
            orbiters[-1].remove(removed_e)
            self.play(removed_e.animate.scale(0.6).set_opacity(0.0), run_time=0.6)

        # Creo una "lacuna" vicina a B (vuoto di elettrone su legame verso sinistra)
        hole_pos = centers[-1] + 0.45*(centers[-2]-centers[-1])  # verso il vicino a sinistra
        hole = Circle(radius=0.12, color=HOLE_COLOR, stroke_width=6).move_to(hole_pos)
        hole.set_opacity(0.9)
        hole_glow = Circle(radius=0.18, stroke_width=0, color=HOLE_COLOR).move_to(hole_pos)
        hole_glow.set_opacity(0.25)
        hole_lbl = Text("buco", font_size=20, color=HOLE_COLOR).next_to(hole, DOWN, buff=0.1)
        self.play(FadeIn(hole_glow), FadeIn(hole), FadeIn(hole_lbl))
        self.wait(0.6)

        # --- Ridistribuzione: un e⁻ dal vicino "riempie" B → il buco si sposta a sinistra
        explain2 = Text("Un e⁻ dal vicino “riempie” il legame con B → il buco si sposta", font_size=26)\
                   .next_to(explain1, DOWN, buff=0.25)
        self.play(FadeIn(explain2))

        # Prendo un elettrone dall'atomo N-1 e lo faccio saltare verso B
        donor_atom_idx = -2  # vicino a sinistra di B
        donor_orb = orbiters[donor_atom_idx][0]
        jump_to_b = centers[-1] + 0.35*(centers[-2]-centers[-1])  # punto sul legame verso B

        # Per "staccare" l'elettrone dal moto orbitale, rimuovo updater e lo animo
        for up in donor_orb.updaters:
            donor_orb.remove_updater(up)
        self.play(donor_orb.animate.move_to(jump_to_b), run_time=0.6)
        self.play(donor_orb.animate.move_to(centers[-1] + 0.5*(centers[-2]-centers[-1])), run_time=0.5)

        # Il "buco" salta al posto lasciato dall'elettrone
        new_hole_pos = centers[-2] + 0.45*(centers[-1]-centers[-2])
        self.play(hole.animate.move_to(new_hole_pos), hole_glow.animate.move_to(new_hole_pos), run_time=0.6)
        self.wait(0.3)

        # Rimetto il donatore su un'orbita attorno a B (ora legame con B è "completo")
        # Creo un nuovo orbiter ancorato a B (così non devo reinserire updater a mano)
        new_e_around_B = make_orbiters(centers[-1], n=1, R=0.55, speed=1.4, color=E_COLOR, size=0.085)
        self.add(*new_e_around_B)
        self.play(donor_orb.animate.set_opacity(0.0), run_time=0.4)
        self.remove(donor_orb)

        # --- Propagazione del buco: ripeti il "salto" a cascata verso sinistra
        chain_text = Text("La lacuna si propaga a sinistra: conduzione p-type", font_size=26)\
                     .next_to(explain2, DOWN, buff=0.25)
        self.play(FadeIn(chain_text))

        # Lista di target per il buco: legami successivi verso sinistra
        hole_targets = [
            centers[-3] + 0.45*(centers[-2]-centers[-3]),
            centers[-4] + 0.45*(centers[-3]-centers[-4]),
        ]

        # Per ogni salto, prendo un orbiter dal vicino di sinistra e lo faccio "riempire" il buco
        donor_indices = [-3, -4]
        for idx, hpos in zip(donor_indices, hole_targets):
            # prendo un elettrone dall'atomo idx
            if len(orbiters[idx]) == 0:
                continue
            e = orbiters[idx][0]
            # stacco l'updater
            for up in list(e.updaters):
                e.remove_updater(up)
            # animo: e⁻ va verso la vecchia posizione del buco
            self.play(e.animate.move_to(hole.get_center()), run_time=0.6)
            # il buco si sposta alla posizione successiva
            self.play(hole.animate.move_to(hpos), hole_glow.animate.move_to(hpos), run_time=0.6)
            # l'elettrone che ho mosso "scompare" (rappresenta il fatto che ha riempito quel legame)
            self.play(e.animate.set_opacity(0.0), run_time=0.3)
            self.remove(e)

        self.wait(0.4)

        # --- Campo elettrico e direzioni apparenti di moto
        efield = Arrow(LEFT*5.5 + UP*2.9, RIGHT*5.5 + UP*2.9, buff=0).set_stroke(width=6)
        efield_lbl = Text("Campo elettrico", font_size=24).next_to(efield, DOWN, buff=0.15)
        current_lbl = Text("Corrente p-type: le lacune avanzano", font_size=26, color=HOLE_COLOR)\
                      .to_edge(DOWN).shift(0.6*UP)
        self.play(FadeIn(efield), FadeIn(efield_lbl), FadeIn(current_lbl))

        # Piccolo drift: sposta lievemente gli elettroni verso sinistra mentre il buco va a sinistra
        drift = 0.25*LEFT
        drift_group = VGroup(*orbiters[:-1], *new_e_around_B)  # tutti tranne B originali già usati
        self.play(
            AnimationGroup(
                *[AnimationGroup(*[d.animate.shift(drift) for d in grp], lag_ratio=0.1) for grp in drift_group],
                hole.animate.shift(0.5*LEFT),
                hole_glow.animate.shift(0.5*LEFT),
                lag_ratio=0.1
            ),
            run_time=1.2
        )

        # --- Chiusura
        final_note = Text("Doping con Boro → crea lacune → conduzione per buco (p-type)", font_size=28)\
                     .to_edge(DOWN)
        self.play(FadeOut(explain1), FadeOut(explain2), FadeOut(chain_text), FadeOut(current_lbl))
        self.play(FadeIn(final_note))
        self.wait(1.0)