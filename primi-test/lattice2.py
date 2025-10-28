from manim import *

# --- Diagramma a bande 2D ----------------------------------------------------
class BandDiagram(VGroup):
    def __init__(self, width=4, height=3, label="n-type (donor)", **kwargs):
        super().__init__(**kwargs)
        self.width = width
        self.height = height

        y_cb = +height/2 * 0.85
        y_vb = -height/2 * 0.85

        cb = Line(LEFT*width/2 + UP*y_cb, RIGHT*width/2 + UP*y_cb)
        vb = Line(LEFT*width/2 + UP*y_vb, RIGHT*width/2 + UP*y_vb)

        cb_lbl = Text("Conduction band (CB)", font_size=24).next_to(cb, UP, buff=0.1)
        vb_lbl = Text("Valence band (VB)", font_size=24).next_to(vb, DOWN, buff=0.1)

        axis = Line(UP*height/2, DOWN*height/2).set_stroke(opacity=0.3)
        e_lbl = Text("Energy", font_size=20).rotate(PI/2).next_to(axis, LEFT, buff=0.2)

        self.level = Line(ORIGIN, ORIGIN)
        self.level_lbl = Text("", font_size=24)

        self.add(axis, e_lbl, cb, vb, cb_lbl, vb_lbl, self.level, self.level_lbl)
        self.cb, self.vb = cb, vb
        self.y_cb, self.y_vb = y_cb, y_vb

    def set_donor_level(self, gap_fraction=0.12):
        y = self.y_cb - (self.y_cb - self.y_vb) * gap_fraction
        lvl = Line(LEFT*self.width/2, RIGHT*self.width/2).shift(UP*y)
        lbl = Text("Donor level (P)", font_size=24).next_to(lvl, DOWN, buff=0.1)
        self.level.become(lvl)
        self.level_lbl.become(lbl)
        return self

    def set_acceptor_level(self, gap_fraction=0.12):
        y = self.y_vb + (self.y_cb - self.y_vb) * gap_fraction
        lvl = Line(LEFT*self.width/2, RIGHT*self.width/2).shift(UP*y)
        lbl = Text("Acceptor level (B)", font_size=24).next_to(lvl, UP, buff=0.1)
        self.level.become(lvl)
        self.level_lbl.become(lbl)
        return self

# --- Reticolo 3D semplificato ------------------------------------------------
def make_silicon_lattice(n=5, spacing=0.8, dopant_coord=(3,2,2), dopant_color=BLUE):
    atoms = VGroup()
    bonds = VGroup()

    def idx_to_point(i,j,k):
        return np.array([
            (i-(n-1)/2)*spacing,
            (j-(n-1)/2)*spacing,
            (k-(n-1)/2)*spacing
        ])

    for i in range(n):
        for j in range(n):
            for k in range(n):
                p = idx_to_point(i,j,k)
                is_dop = (i,j,k) == dopant_coord
                dot = Dot3D(p, radius=0.06, color=(dopant_color if is_dop else DARK_GRAY))
                atoms.add(dot)

    # legami verso i vicini (x,y,z) con stroke sottile e opacità bassa
    for i in range(n):
        for j in range(n):
            for k in range(n):
                p = idx_to_point(i,j,k)
                for di,dj,dk in [(1,0,0),(0,1,0),(0,0,1)]:
                    ni, nj, nk = i+di, j+dj, k+dk
                    if 0 <= ni < n and 0 <= nj < n and 0 <= nk < n:
                        q = idx_to_point(ni,nj,nk)
                        seg = Line3D(p, q, color=GRAY)
                        seg.set_stroke(width=1).set_opacity(0.35)
                        bonds.add(seg)

    return VGroup(bonds, atoms)

# --- Scena -------------------------------------------------------------------
class DopingInSilicon(ThreeDScene):
    def construct(self):
        self.camera.background_color = "#111111"
        self.set_camera_orientation(phi=60*DEGREES, theta=-45*DEGREES, zoom=1.0)

        # Reticolo n-type
        lattice_n = make_silicon_lattice(n=5, spacing=0.8, dopant_coord=(3,2,2), dopant_color=BLUE)
        title3d = Text("Silicon lattice (n-type: Phosphorus donor)", font_size=28).to_edge(UP).set_z_index(5)
        self.play(FadeIn(lattice_n, shift=0.5*OUT), FadeIn(title3d))

        # Evidenzia dopante (rettangolo 2D fissato alla camera)
        dopant_atom = [m for m in lattice_n.submobjects[1] if isinstance(m, Dot3D) and m.get_color()==BLUE][0]
        dopant_glow = SurroundingRectangle(Dot(dopant_atom.get_center()), color=BLUE, buff=0.12).set_z_index(4)
        dopant_lbl = Text("P (donor)", font_size=24, color=BLUE).next_to(dopant_glow, UP).set_z_index(4)
        self.add_fixed_in_frame_mobjects(dopant_glow, dopant_lbl)
        self.play(Create(dopant_glow), FadeIn(dopant_lbl))

        # Elettrone legato
        e_bound = Dot3D(dopant_atom.get_center()+0.1*UP, radius=0.05, color=YELLOW)
        self.play(FadeIn(e_bound, scale=0.5))

        # Band diagram 2D
        band = BandDiagram(width=4.2, height=3.0).to_corner(UR).set_z_index(10)
        band.set_donor_level(gap_fraction=0.12)
        self.add_fixed_in_frame_mobjects(band)
        self.play(FadeIn(band))

        # Elettrone dal donor alla CB (overlay 2D)
        donor_point = band.level.get_center()
        cb_point = band.cb.get_center()
        e2d = Dot(donor_point, radius=0.06, color=YELLOW).set_z_index(11)
        self.add_fixed_in_frame_mobjects(e2d)
        self.play(FadeIn(e2d, scale=0.5))
        self.play(e2d.animate.move_to(cb_point), run_time=1.2)

        # Campo elettrico: freccia 2D fissa
        field_arrow = Arrow(LEFT*2.5+DOWN*2.8, RIGHT*2.5+DOWN*2.8).set_z_index(6)
        field_lbl = Text("E-field", font_size=22).next_to(field_arrow, DOWN).set_z_index(6)
        self.add_fixed_in_frame_mobjects(field_arrow, field_lbl)
        self.play(FadeIn(field_arrow), FadeIn(field_lbl))

        # Elettrone libero che si muove nel reticolo
        path_points = [
            dopant_atom.get_center()+0.2*RIGHT,
            dopant_atom.get_center()+RIGHT*1.2+0.1*UP,
            dopant_atom.get_center()+RIGHT*2.0+0.2*DOWN,
            dopant_atom.get_center()+RIGHT*3.0+0.2*UP,
        ]
        e_free = Dot3D(path_points[0], radius=0.05, color=YELLOW)
        self.play(FadeIn(e_free))
        for a,b in zip(path_points, path_points[1:]):
            self.play(e_free.animate.move_to(b), run_time=0.8, rate_func=linear)

        # Passa a p-type
        self.play(FadeOut(e_free), FadeOut(e_bound))
        title_p = Text("Silicon lattice (p-type: Boron acceptor)", font_size=28).to_edge(UP).set_z_index(5)
        self.play(Transform(title3d, title_p))

        lattice_p = make_silicon_lattice(n=5, spacing=0.8, dopant_coord=(1,2,2), dopant_color=RED)
        self.play(ReplacementTransform(lattice_n, lattice_p))

        b_dopant_atom = [m for m in lattice_p.submobjects[1] if isinstance(m, Dot3D) and m.get_color()==RED][0]
        b_glow = SurroundingRectangle(Dot(b_dopant_atom.get_center()), color=RED, buff=0.12).set_z_index(4)
        b_lbl = Text("B (acceptor)", font_size=24, color=RED).next_to(b_glow, UP).set_z_index(4)
        self.add_fixed_in_frame_mobjects(b_glow, b_lbl)
        self.play(Create(b_glow), FadeIn(b_lbl))

        # Aggiorna band diagram con livello acceptor
        self.play(FadeOut(band.level), FadeOut(band.level_lbl))
        band.set_acceptor_level(gap_fraction=0.12)
        self.play(FadeIn(band.level), FadeIn(band.level_lbl))

        # Elettrone dalla VB al livello acceptor (2D)
        vb_point = band.vb.get_center()
        acc_point = band.level.get_center()
        e_from_vb = Dot(vb_point + 1.0*LEFT, radius=0.06, color=YELLOW).set_z_index(11)
        hole = Dot(vb_point + 1.0*LEFT, radius=0.06, color=PURE_RED).set_z_index(11)
        hole.set_opacity(0.0)
        self.add_fixed_in_frame_mobjects(e_from_vb, hole)
        self.play(FadeIn(e_from_vb, scale=0.5))
        self.play(e_from_vb.animate.move_to(acc_point), run_time=1.2)
        hole.set_opacity(1.0)
        hole_lbl = Text("hole", font_size=20, color=PURE_RED).next_to(hole, DOWN, buff=0.05)
        self.add_fixed_in_frame_mobjects(hole_lbl)
        self.play(FadeIn(hole, hole_lbl))

        # Moto del buco lungo la VB (2D)
        vb_left = band.vb.get_start() + 0.3*RIGHT
        vb_right = band.vb.get_end() - 0.3*RIGHT
        for target in [vb_left, vb_right, vb_left + 0.8*RIGHT]:
            self.play(hole.animate.move_to(np.array([target[0], vb_point[1], 0])), run_time=0.7, rate_func=linear)

        # “Buco” nel reticolo (3D) semitrasparente
        hole3d = Dot3D(b_dopant_atom.get_center()+0.2*LEFT, radius=0.06, color=PURE_RED)
        hole3d.set_opacity(0.7)
        path2 = [
            b_dopant_atom.get_center()+0.2*LEFT,
            b_dopant_atom.get_center()+LEFT*1.0+0.1*UP,
            b_dopant_atom.get_center()+LEFT*2.0+0.2*DOWN,
            b_dopant_atom.get_center()+LEFT*3.0+0.1*UP,
        ]
        self.play(FadeIn(hole3d))
        for a,b in zip(path2, path2[1:]):
            self.play(hole3d.animate.move_to(b), run_time=0.8, rate_func=linear)

        # Riepilogo
        recap = VGroup(
            Text("n-type: donor → elettroni liberi in CB", font_size=26, color=YELLOW),
            Text("p-type: acceptor → buchi mobili in VB", font_size=26, color=PURE_RED)
        ).arrange(DOWN, aligned_edge=LEFT).to_corner(DL).set_z_index(6)
        self.add_fixed_in_frame_mobjects(recap)
        self.play(FadeIn(recap))
        self.wait(1.0)

        self.play(*[FadeOut(m) for m in [hole3d, e_from_vb, hole, hole_lbl, band, title3d, dopant_glow, dopant_lbl, b_glow, b_lbl]])
        self.play(FadeOut(lattice_p))
        self.wait(0.2)