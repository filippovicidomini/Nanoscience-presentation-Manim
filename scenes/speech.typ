
#set par(justify: true)
#set par(first-line-indent: 1.5em)
#set text(size: 14pt, font: "Libertinus Serif")
#set page(numbering: "1")
#set document(
  title: [From Silicon to Light ],
  author: "Filippo Vicidomini",
)



#align(center)[
  #title()
  #text(font: "Libertinus Serif", style: "italic")[
    Filippo Vicidomini\
    University Ca' Foscari of Venice\
    Nanomaterials and Nanotechnology\

  ]
]

#outline(depth: 1)

⸻
#set heading(numbering: "1.1")


= Construction of the Silicon Lattice

In this first part, we introduce the atomic structure of crystalline silicon.
Starting from a single atom, we build a simplified two-dimensional representation
of the silicon lattice and show how covalent bonds arise from the sharing of valence electrons.
Although real silicon is a three-dimensional crystal,
this model captures the essential bonding mechanism
that holds the solid together.

== on video
[0–6 s] – Single atom

“Let’s start with a single silicon atom.
At the center is a positively charged nucleus,
surrounded by four valence electrons.”

[6–14 s] – Building the lattice

“Now imagine many silicon atoms coming together.
We build a simplified two-dimensional lattice,
where each atom is connected to four neighbors.
Although real silicon is three-dimensional,
this view helps us clearly see how bonds form.”

[14–24 s] – Covalent bonds

“When two atoms bond, they share pairs of valence electrons.
Each silicon atom forms four covalent bonds,
sharing its electrons with its neighbors.
These shared electrons hold the crystal together.”

[24–30 s] – Final idea

“The result is an ordered silicon lattice.
In reality, these electrons are not localized,
but spread through the crystal forming energy bands.
Still, this model captures the essential physics of bonding in silicon.”

⸻

= n-Type Doping: Phosphorus in Silicon

When silicon is doped with phosphorus, a group-V element, one silicon atom in the lattice is replaced by a phosphorus atom. While silicon has four valence electrons, phosphorus has five. Four of these electrons participate in covalent bonding with neighboring silicon atoms, preserving the local bonding structure of the lattice. The fifth electron, however, is not required to complete the bonds and remains only weakly bound to the donor site.

This situation is commonly described by the donor ionization process $P_(S i) = P_(S i)^+ + e^-$,
where the phosphorus atom becomes positively charged after releasing an electron into the crystal. The released electron can move through the lattice and contributes to electrical conduction.

From an energy-band perspective, this extra electron originates from a donor energy level $E_D$ located slightly below the conduction band edge $E_c$. The energy required to free the electron is therefore small and given by $E_(i o n)= E_c - E_D$.
At room temperature, this energy is typically low enough that most donor atoms are ionized, providing a significant number of free electrons.









== on video

[0–5 s] – Starting point

“We start from the silicon lattice we have just built,
where atoms are connected by shared electron pairs
forming a stable covalent network.”

[5–12 s] – Introducing the dopant

“Now we replace one silicon atom with a phosphorus atom.
Phosphorus has five valence electrons,
one more than silicon.”

[12–20 s] – Donor electron

“Four of these electrons form covalent bonds with the lattice.
The fifth electron is not needed for bonding,
and is only weakly bound to the crystal.”

[20–26 s] – Free electron

“This extra electron can move through the lattice,
acting as a mobile charge carrier.”

[26–30 s] – Meaning

“This is an n-type semiconductor:
electrons donated by impurity atoms
make electrical conduction much easier.”

⸻

= p-Type Doping: Boron in Silicon

In p-type doping, silicon is doped with acceptor impurities, such as boron, which belong to group III of the periodic table. These atoms have three valence electrons, one fewer than silicon. When an acceptor atom replaces a silicon atom in the crystal lattice, it can form only three covalent bonds with neighboring silicon atoms. As a result, one bond remains incomplete, corresponding to the absence of an electron.

This process is described by the acceptor ionization reaction
$A_(S i) = A_(S i)^- + h^+$,
where the acceptor atom captures an electron from the valence band and becomes negatively charged, leaving behind a hole. The hole is not a physical particle, but it behaves as a mobile positive charge carrier within the lattice.

From the energy-band perspective, acceptor impurities introduce an energy level E_A located slightly above the valence band edge E_v. The energy required to promote an electron from the valence band into the acceptor level is given by
$E_(i o n) = E_A - E_v$.
At typical operating temperatures, this ionization energy is small compared to the thermal energy $k_B T$, so most acceptor atoms are ionized, generating a high concentration of holes.

== on video
[0–5 s] – Starting point

“Let us again start from a perfect silicon lattice,
where all valence electrons are paired in covalent bonds.”

[5–12 s] – Introducing the acceptor

“Now one silicon atom is replaced by a boron atom.
Boron has three valence electrons,
one fewer than silicon.”

[12–20 s] – Hole creation

“When boron forms bonds with its neighbors,
one electron is missing to complete the four bonds.
This absence is called a hole.”

[20–26 s] – Hole motion

“The hole behaves like a positive charge.
As nearby electrons move to fill the missing bond,
the hole effectively moves through the lattice.”

[26–30 s] – Meaning

“This defines a p-type semiconductor:
charge transport is dominated by mobile holes,
rather than free electrons.”


= Formation of the PN Junction
We start with two separate silicon lattices: an n-type region on the left and a p-type region on the right.
As the two regions are brought closer together, their crystal lattices connect seamlessly at the interface. New covalent bonds form between atoms at the boundary, creating a continuous silicon structure.
Once the junction is formed, carriers begin to diffuse due to concentration gradients. Electrons diffuse from the n-type region into the p-type region, while holes diffuse in the opposite direction, from p-type to n-type.
Near the interface, electrons and holes recombine. This recombination removes mobile charge carriers from the junction region, leaving behind fixed ionized donors on the n-side and fixed ionized acceptors on the p-side.
The result is the formation of a depletion region, a zone depleted of free carriers. The exposed fixed charges generate an internal electric field pointing from the n-type side to the p-type side.

== on video
[0–7 s] – Separate regions

“We begin with two separate silicon regions:
an n-type region on the left, rich in electrons,
and a p-type region on the right, rich in holes.”

[7–15 s] – Bringing regions together

“When the two regions are brought into contact,
their crystal lattices connect at the interface.
New covalent bonds form,
creating a continuous silicon crystal.”

[15–26 s] – Carrier diffusion

“After the junction forms, charge carriers begin to diffuse
due to concentration differences.
Electrons move from the n-type side toward the p-type side,
while holes diffuse in the opposite direction.”

[26–34 s] – Recombination

“Near the interface, electrons and holes recombine.
This removes mobile carriers from the junction region,
leaving behind fixed ionized donors on the n-side
and fixed ionized acceptors on the p-side.”

[34–40 s] – Depletion region

“This creates the depletion region,
a zone with almost no free carriers.
The fixed charges generate an internal electric field,
directed from the n-type side toward the p-type side.”


= Forward Bias
Now we apply forward bias: the p-side is connected to the positive terminal and the n-side to the negative terminal. This external voltage partially cancels the built-in potential.
Visually, the depletion region shrinks because the barrier is reduced. With a lower barrier, majority carriers can cross the junction more easily: electrons are injected from the n-side into the p-side, and holes are injected from the p-side into the n-side.
Once injected, these carriers move across the junction and recombine, and the result is a large current. This is why a diode conducts strongly in forward bias: the junction is no longer “blocking” carriers, it’s allowing injection and recombination.

== on video
[0–8 s] – Applying forward bias
“Now we apply forward bias: the p-side is connected to the positive terminal
and the n-side to the negative terminal.
This external voltage partially cancels the built-in potential.”

[8–16 s] – Depletion region shrinks
“Visually, the depletion region shrinks because the barrier is reduced.
With a lower barrier, majority carriers can cross the junction more easily:
electrons are injected from the n-side into the p-side,
and holes are injected from the p-side into the n-side.”

[16–28 s] – Carrier injection
“Once injected, these carriers move across the junction and recombine.
The result is a large current.
This is why a diode conducts strongly in forward bias:
the junction is no longer blocking carriers, it’s allowing injection and recombination.”


= Reverse Bias
Now we apply reverse bias: the p-side is connected to the negative terminal and the n-side to the positive terminal. This external voltage adds to the built-in potential.
Visually, the depletion region widens because the barrier is increased. With a higher barrier, majority carriers find it much harder to cross the junction: electrons on the n-side are pulled away from the junction, and holes on the p-side are also pulled away.
As a result, very few carriers can cross the junction, leading to a very small current. This is why a diode blocks current in reverse bias: the junction acts as a strong barrier preventing carrier flow.
== on video
[0–8 s] – Applying reverse bias
“Now we apply reverse bias: the p-side is connected to the negative terminal
and the n-side to the positive terminal.
This external voltage adds to the built-in potential.”

[8–16 s] – Depletion region widens
“Visually, the depletion region widens because the barrier is increased.
With a higher barrier, majority carriers find it much harder to cross the junction:
electrons on the n-side are pulled away from the junction,
and holes on the p-side are also pulled away.”

[16–28 s] – Carrier blocking
“As a result, very few carriers can cross the junction, leading to a very small current.
This is why a diode blocks current in reverse bias:
the junction acts as a strong barrier preventing carrier flow.”




= Energy Band Diagram of the PN Junction
This scene translates the PN junction into the energy-band diagram, where the vertical axis is energy and the horizontal axis is position across the junction.

Before contact, the n-type and p-type regions are still separate. Their conduction and valence bands are flat in each region, but their Fermi levels are different: the n-type side has a higher Fermi level (more electrons available), while the p-type side has a lower Fermi level (more holes).

After contact, at equilibrium, carriers diffuse and leave behind fixed charges in the depletion region. The electrostatic potential that forms inside the junction appears here as band bending: both the conduction band edge and the valence band edge bend across the junction. At equilibrium the system settles into a single, flat Fermi level across the entire device, meaning there is no net current. The vertical brace highlights the potential barrier that carriers must overcome.

Under forward bias, the external voltage reduces the built-in potential. In the band diagram, that means the bending is weaker and the barrier height decreases. With a lower barrier, carriers can cross the junction more easily: electrons from the n-side and holes from the p-side are injected across, producing a large current.
Under reverse bias, the external voltage strengthens the built-in potential. The band bending becomes stronger and the barrier increases, widening the depletion region and making carrier injection extremely unlikely. In this condition, the junction blocks current (aside from a tiny leakage current not shown here).

= on video
[0–7 s] – Before contact

“We now describe the PN junction using the energy-band diagram,
where energy is plotted versus position.
Before contact, the n-type and p-type regions are separate.
Their bands are flat,
but their Fermi levels are different:
higher on the n-side, lower on the p-side.”

[7–16 s] – After contact (equilibrium)

“After contact, carriers diffuse and leave fixed charges in the depletion region.
This creates an internal electric potential,
which appears here as band bending.
At equilibrium, the Fermi level becomes flat across the junction,
meaning there is no net current.
The band bending forms a potential barrier.”

[16–25 s] – Forward bias

“Under forward bias, the applied voltage reduces this barrier.
The band bending weakens,
allowing electrons and holes to cross the junction easily
and produce a large current.”

[25–33 s] – Reverse bias

“Under reverse bias, the applied voltage increases the barrier.
The bands bend more strongly,
the depletion region widens,
and current is effectively blocked.”




= Light Emitting Diode (LED)
In this scene we use the energy-band picture to explain how a PN junction becomes a light source. The conduction band E_c and the valence band E_v are shown versus position, with a small band bending consistent with forward bias.

Under forward bias, electrons are injected into the p-side and holes are injected into the n-side. In the active region near the junction, an electron can drop from the conduction band down toward the valence band and recombine with a hole.
That transition releases energy. In a direct bandgap semiconductor, the released energy is emitted efficiently as a photon. The expanding circular wave represents that light emission event.
The crucial idea is that the photon energy is approximately the bandgap energy:
$E_gamma approx E_g = E_c - E_v$
So when the bandgap is small, the photon has lower energy and longer wavelength (red). When the bandgap increases, the emitted photons become higher energy and shift to shorter wavelengths (green, then blue).
The last part of the animation makes this explicit by increasing the $E_c–E_v$ separation and showing the corresponding color change in the emitted light.
== on video
[0–10 s] – Forward bias and carrier injection
“In this scene we use the energy-band picture to explain how a PN junction becomes a light source.
The conduction band E_c and the valence band E_v are shown versus position, with a small band bending consistent with forward bias.
Under forward bias, electrons are injected into the p-side and holes are injected into the n-side.
In the active region near the junction, an electron can drop from the conduction band down toward the valence band and recombine with a hole.”

[10–22 s] – Photon emission
“That transition releases energy.
In a direct bandgap semiconductor, the released energy is emitted efficiently as a photon.
The expanding circular wave represents that light emission event.”

[22–34 s] – Photon energy and bandgap
“The crucial idea is that the photon energy is approximately the bandgap energy:
$E_gamma approx E_g = E_c - E_v$
So when the bandgap is small, the photon has lower energy and longer wavelength (red).
When the bandgap increases, the emitted photons become higher energy and shift to shorter wavelengths (green, then blue).”

[34–48 s] – Color change with bandgap
“The last part of the animation makes this explicit by increasing the $E_c–E_v$ separation
and showing the corresponding color change in the emitted light.”


= I-V Characteristic of the PN Junction
This scene shows the macroscopic signature of everything we saw at the junction: diffusion, drift, the depletion region, and the potential barrier.

In reverse bias (negative voltage), the depletion region widens and the barrier increases. Majority carriers are pulled away from the junction, so almost no current flows — only a tiny leakage current. At sufficiently large reverse voltage, the diode can enter reverse breakdown, where the reverse current increases rapidly (shown here as a conceptual marker).

In forward bias (positive voltage), the external voltage lowers the barrier and narrows the depletion region. Once the applied voltage approaches the knee voltage (around $0.7 V$ for silicon), carrier injection becomes efficient and the current rises exponentially.
The zoom around $0 V$ highlights that the current does not “jump” abruptly: it increases smoothly, but the exponential growth makes it look like a sharp turn on a large scale.

This is why the I–V characteristic is a compact, measurable summary of the junction’s internal potential barrier.
== on video
[0–10 s] – Reverse bias
“In reverse bias (negative voltage), the depletion region widens and the barrier increases.
Majority carriers are pulled away from the junction, so almost no current flows — only a tiny leakage current.
At sufficiently large reverse voltage, the diode can enter reverse breakdown, where the reverse current increases rapidly (shown here as a conceptual marker).”

[10–22 s] – Forward bias
“In forward bias (positive voltage), the external voltage lowers the barrier and narrows the depletion region.
Once the applied voltage approaches the knee voltage (around $0.7 V$ for silicon), carrier injection becomes efficient and the current rises exponentially.”

[22–34 s] – Zoom around 0 V
“The zoom around $0 V$ highlights that the current does not jump abruptly: it increases smoothly,
but the exponential growth makes it look like a sharp turn on a large scale.”

[34–48 s] – Summary
“This is why the I–V characteristic is a compact, measurable summary of the junction’s internal potential barrier.”










#pagebreak()
#set bibliography(full: true, title: "References")
#bibliography("../Nanomaterials-Presentation.bib")
