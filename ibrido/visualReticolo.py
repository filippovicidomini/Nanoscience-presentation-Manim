#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# === PERCORSI GIÀ INSERITI (controlla che esistano) ===
# --- PATCH PERCORSI (usa QUESTI) ---
from pathlib import Path

# 1) Assoluto (come nello screenshot: file in radice del repo)
SRC_OBJ = "/Users/filippovicidomini/Library/CloudStorage/GoogleDrive-filippo.vicidomini@alumni.unive.it/My Drive/University/Nanoscience-presentation-Manim/tetrahedral_touch_calibrated.obj"

# 2) Output: mettilo dove vuoi (qui lo scrivo in ibrido/)
OUT_OBJ = str(Path(SRC_OBJ).with_name("edge_dimer_lattice.obj"))
# --- fallback automatico se il path assoluto non esiste ---
if not Path(SRC_OBJ).exists():
    # se stai eseguendo lo script da ibrido/, risale di una cartella e prende il file in radice
    SRC_OBJ = str(Path(__file__).resolve().parent.parent / "tetrahedral_touch_calibrated.obj")

print("USO SRC_OBJ =", SRC_OBJ)
print("USO OUT_OBJ =", OUT_OBJ)


# === PARAMETRI PRINCIPALI (modifica qui se vuoi) ===
TARGET_R   = 3     # raggio "visivo" del tetraedro normalizzato
GAP        = -0.01   # <0: leggera compenetrazione sullo spigolo; >0: piccola fessura
MARGIN     = 1    # separazione tra dimeri (>=1.4 per evitare contatti indesiderati)
N_X, N_Y, N_Z = 3, 3, 3   # reticolo di dimeri (piano per avere 1 solo contatto per tetra)
SIMPLIFY_TO = 4000   # 0 = nessuna decimazione; altrimenti facce target globali

# ======================================================

import sys
from pathlib import Path
import numpy as np
import trimesh
from trimesh.transformations import rotation_matrix, translation_matrix
from trimesh.geometry import align_vectors

def load_mesh_center_scale(path: Path, target_r: float) -> trimesh.Trimesh:
    if not path.exists():
        print(f"[ERRORE] File sorgente non trovato:\n{path}", file=sys.stderr)
        sys.exit(1)
    tm = trimesh.load(path, process=False)
    if isinstance(tm, trimesh.Scene):
        if len(tm.geometry) == 0:
            print("[ERRORE] La scena sorgente è vuota.", file=sys.stderr)
            sys.exit(1)
        tm = trimesh.util.concatenate([g for g in tm.geometry.values()])
    V = np.asarray(tm.vertices, dtype=np.float64).copy()
    F = np.asarray(tm.faces, dtype=np.int64).copy()
    if V.size == 0 or F.size == 0:
        print("[ERRORE] La mesh sorgente non ha vertici/facce.", file=sys.stderr)
        sys.exit(1)
    # centra e scala
    V -= V.mean(axis=0)
    r = np.linalg.norm(V, axis=1).max()
    if r > 1e-12:
        V *= (target_r / r)
    tm = trimesh.Trimesh(vertices=V, faces=F, process=False)
    # triangola se necessario
    if tm.faces.shape[1] != 3:
        try:
            tm = tm.triangulate()
        except Exception:
            tm = tm.subdivide()
    return tm

def longest_edge_indices(mesh: trimesh.Trimesh) -> tuple[int,int]:
    edges = mesh.edges_unique
    ev = mesh.vertices[edges]
    lengths = np.linalg.norm(ev[:,0,:] - ev[:,1,:], axis=1)
    k = int(np.argmax(lengths))
    i, j = int(edges[k,0]), int(edges[k,1])
    return (i, j)

def build_edge_dimer(base: trimesh.Trimesh, gap: float = 0.0) -> trimesh.Trimesh:
    # scegli lo spigolo di contatto (quello più lungo per robustezza)
    i, j = longest_edge_indices(base)
    Va = base.vertices
    ea0, ea1 = Va[i], Va[j]
    mid = 0.5*(ea0 + ea1)
    u = ea1 - ea0
    nu = np.linalg.norm(u)
    if nu < 1e-12:
        raise ValueError("Spigolo nullo.")
    u /= nu

    # B = copia del tetra; allinea il SUO spigolo i-j allo stesso asse u
    B = base.copy()
    Vb = B.vertices
    eb0, eb1 = Vb[i], Vb[j]
    v = eb1 - eb0
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        raise ValueError("Spigolo B nullo.")
    v /= nv

    # 1) allinea direzione spigolo (v -> u)
    R_align = align_vectors(v, u)  # 4x4
    B.apply_transform(R_align)

    # 2) porta il midpoint dello spigolo di B su quello di A, con piccolo offset lungo u (gap)
    Vb = B.vertices
    eb0p, eb1p = Vb[i], Vb[j]
    mid_b = 0.5*(eb0p + eb1p)
    T_match = translation_matrix((mid - mid_b) + (gap * u))
    B.apply_transform(T_match)

    # 3) ruota di 180° attorno allo spigolo (asse u, punto mid) per metterlo “dall'altro lato”
    R_flip = rotation_matrix(np.pi, direction=u, point=mid)
    B.apply_transform(R_flip)

    A = base.copy()
    return trimesh.util.concatenate([A, B])

def tile_dimer(dimer: trimesh.Trimesh, nx: int, ny: int, nz: int, margin: float = 1.5) -> trimesh.Trimesh:
    mins, maxs = dimer.bounds
    extent = (maxs - mins)
    step = float(np.max(extent)) * float(margin)
    offx, offy, offz = (nx - 1) / 2.0, (ny - 1) / 2.0, (nz - 1) / 2.0
    meshes = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                g = dimer.copy()
                T = translation_matrix([
                    (ix - offx) * step,
                    (iy - offy) * step,
                    (iz - offz) * step
                ])
                g.apply_transform(T)
                meshes.append(g)
    return trimesh.util.concatenate(meshes)

def simplify_mesh(mesh: trimesh.Trimesh, target_faces: int | None) -> trimesh.Trimesh:
    if not target_faces or target_faces <= 0:
        return mesh
    before = len(mesh.faces)
    try:
        m2 = mesh.simplify_quadratic_decimation(int(target_faces))
        if len(m2.faces) < before:
            return m2
    except Exception:
        pass
    return mesh  # se non decima, proseguiamo comunque

def main():
    src = Path(SRC_OBJ)
    outp = Path(OUT_OBJ)

    print("== EDGE-TO-EDGE LATTICE ==")
    print("Sorgente :", src)
    print("Output   :", outp)
    print(f"N        : {N_X} x {N_Y} x {N_Z}")
    print(f"targetR  : {TARGET_R}   | gap: {GAP}   | margin: {MARGIN}")
    print(f"simplify : {SIMPLIFY_TO}")

    base = load_mesh_center_scale(src, target_r=TARGET_R)
    dimer = build_edge_dimer(base, gap=GAP)
    lattice = tile_dimer(dimer, N_X, N_Y, N_Z, margin=MARGIN)
    print("Pre-decimazione:", len(lattice.vertices), "vertici,", len(lattice.faces), "facce")

    lattice = simplify_mesh(lattice, target_faces=SIMPLIFY_TO)
    print("Post-decimazione:", len(lattice.vertices), "vertici,", len(lattice.faces), "facce")

    outp.parent.mkdir(parents=True, exist_ok=True)
    lattice.export(outp)
    print(">> Scritto:", outp.resolve())
    # salvalo anche in formato .ply per visualizzazione rapida
    ply_out = outp.with_suffix('.ply')
    lattice.export(ply_out)
    print(">> Scritto:", ply_out.resolve())

if __name__ == "__main__":
    main()