
#set par(justify: true)
#set par(first-line-indent: 1.5em)
#set text(size: 14pt, font: "Libertinus Serif")
#set page(numbering: "1")
#set document(
  title: [From Doping Semiconductors \n to PN junctions],
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
[0–8 s] – Single atom

“Let’s begin with a single silicon atom.
At the center we have the positively charged nucleus,
and around it four valence electrons.”

[8–16 s] – Building the lattice

“Now, instead of a single atom, imagine an entire silicon crystal.
Here we build a simplified two-dimensional lattice,
where each silicon atom is surrounded by four neighbors.
In reality the structure is three-dimensional,
but this 2D view lets us clearly see how bonds form.”

[16–29 s] – Covalent bonds

“When two silicon atoms are close,
each provides one valence electron to the shared bond.
Each bond is made of a pair of electrons,
As the lattice forms, every atom shares four electron pairs
with its neighbors.
The individual valence electrons that originally orbited each atom
are now absorbed into these covalent bonds,
holding the entire lattice together.”

[29–37 s] – Final highlight

“The final result is an ordered silicon lattice.
Each atom is connected to its neighbors through four shared electron pairs.
In this classical picture the electrons appear localized,
but in a real crystal their wavefunctions spread over many atoms
and form the electronic bands of the solid.
Still, this model captures the essential idea:
the silicon lattice is held together by a network of shared valence electrons.”

⸻

= n-Type Doping: Phosphorus in Silicon

When silicon is doped with phosphorus, a group-V element, one silicon atom in the lattice is replaced by a phosphorus atom. While silicon has four valence electrons, phosphorus has five. Four of these electrons participate in covalent bonding with neighboring silicon atoms, preserving the local bonding structure of the lattice. The fifth electron, however, is not required to complete the bonds and remains only weakly bound to the donor site.

This situation is commonly described by the donor ionization process $P_(S i) = P_(S i)^+ + e^-$,
where the phosphorus atom becomes positively charged after releasing an electron into the crystal. The released electron can move through the lattice and contributes to electrical conduction.

From an energy-band perspective, this extra electron originates from a donor energy level $E_D$ located slightly below the conduction band edge $E_c$. The energy required to free the electron is therefore small and given by $E_(i o n)= E_c - E_D$.
At room temperature, this energy is typically low enough that most donor atoms are ionized, providing a significant number of free electrons.









== on video

[0–6 s] – Starting point

“We start from the same silicon lattice we have just built.
Each atom is bonded to its neighbors through shared electron pairs,
forming a stable covalent network.”

[6–14 s] – Introducing the dopant

“Now, we replace one silicon atom with a phosphorus atom.
Phosphorus belongs to group five of the periodic table,
so it has one extra valence electron compared to silicon.”

[14–24 s] – Donor electron

“Four of phosphorus’ valence electrons participate in covalent bonds,
just like silicon.
The fifth electron, however, is not needed to complete the bonds.
It is only weakly bound to the lattice
and can move freely through the crystal.”

[24–38 s] – Free electron motion

“This extra electron becomes a mobile charge carrier.
It can wander through the lattice,
being slightly influenced by the surrounding atoms
but not tied to any specific bond.
This is what defines an n-type semiconductor:
the presence of free electrons donated by impurity atoms.”

[38–48 s] – Physical interpretation

“In terms of energy bands,
this electron occupies states close to the conduction band,
making electrical conduction much easier.
Even a small concentration of donor atoms
can significantly increase the conductivity of silicon.”

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
[0–6 s] – Starting point

“Let us again start from a perfect silicon lattice,
where all valence electrons are paired in covalent bonds.”

[6–14 s] – Introducing the acceptor

“Now we replace one silicon atom with a boron atom.
Boron belongs to group three,
so it has one valence electron less than silicon.”

[14–24 s] – Creation of a hole

“When boron forms bonds with neighboring silicon atoms,
one electron is missing to complete the four covalent bonds.
This absence of an electron creates what we call a hole.”

[24–40 s] – Hole motion

“The hole is not a physical particle,
but it behaves like a positive charge.
When a nearby electron moves to fill the missing bond,
it leaves behind another hole.
In this way, the hole effectively moves through the lattice,
even though only electrons are actually moving.”

[40–52 s] – Physical interpretation

“This mechanism defines a p-type semiconductor.
Charge transport is dominated by holes,
which act as mobile positive carriers.
In the band picture,
this corresponds to empty states near the top of the valence band,
where electrons can easily move into.”


#pagebreak()
#set bibliography(full: true, title: "References")
#bibliography("../Nanomaterials-Presentation.bib")
