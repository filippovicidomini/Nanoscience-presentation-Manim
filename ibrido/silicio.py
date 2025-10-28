#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Orbitali atomici reali (s, p) con Z_eff (Slater), somma densità e isosuperficie unica.
- Supporta: modalità 'atomic' (configurazione elettronica) e 'sp3' (ibridi tetraedrici).
- Export: .obj di una singola mesh triangolare (isosuperficie della densità totale).
- Vista interattiva opzionale con Plotly (se disponibile).

Dipendenze obbligatorie: numpy, scipy, scikit-image
Opzionali: plotly (solo per preview interattivo)
"""

import re
import argparse
import numpy as np
from math import sqrt, pi
from scipy.special import genlaguerre, gamma
from skimage.measure import marching_cubes
from scipy.ndimage import gaussian_filter

# --------- Utilità base ---------

A0 = 1.0  # unità atomiche: r in multipli di a0

def hydrogenic_R_nl(r, n, l, Z_eff):
    """
    Radiale idrogeno-like (unità atomiche, a0=1), con Z_eff (scalatura di Bohr).
    R_{nl}(r) = (2Z/(na0))^{3/2} * sqrt((n-l-1)! / (2n (n+l)!)) * e^{-ρ/2} * ρ^l * L_{n-l-1}^{2l+1}(ρ)
    ρ = 2Zr/(na0)
    """
    r = np.asarray(r)
    rho = 2.0 * Z_eff * r / float(n)
    k = n - l - 1
    if k < 0:
        # n < l+1 non fisico: restituisco zeri per sicurezza
        return np.zeros_like(r)
    # Normalizzazione
    norm = (2.0 * Z_eff / float(n))**1.5 * np.sqrt(gamma(k+1) / (2.0*n*gamma(n+l+1)))
    L = genlaguerre(k, 2*l+1)(rho)  # polinomio associato di Laguerre
    R = norm * np.exp(-rho/2.0) * (rho**l) * L
    return R

# Sferiche armoniche REALI per l=0,1 (basi s, p_x, p_y, p_z)
# Y_00 = 1/sqrt(4π)
# p_x ∝ sinθ cosφ = x/r, p_y ∝ y/r, p_z ∝ cosθ = z/r, tutti * sqrt(3/4π)

C_S = 1.0/np.sqrt(4.0*pi)
C_P = np.sqrt(3.0/(4.0*pi))

def ang_s(x, y, z, r):
    return C_S * np.ones_like(r)

def ang_px(x, y, z, r):
    return C_P * np.where(r>0, x/r, 0.0)

def ang_py(x, y, z, r):
    return C_P * np.where(r>0, y/r, 0.0)

def ang_pz(x, y, z, r):
    return C_P * np.where(r>0, z/r, 0.0)

# --------- Regole di Slater (solo ns/np, sufficienti per valenza di molti elementi) ---------

AUFBAU_ORDER = [
    (1,'s',2), (2,'s',2), (2,'p',6), (3,'s',2), (3,'p',6),
    (4,'s',2), (3,'d',10), (4,'p',6), (5,'s',2), (4,'d',10),
    (5,'p',6), (6,'s',2), (4,'f',14), (5,'d',10), (6,'p',6),
    (7,'s',2), (5,'f',14), (6,'d',10), (7,'p',6)
]
# Mappa simbolo -> Z e configurazione standard (fino ad Ar bastano s/p per molti usi)
ELEMENTS_SP = {
    'H':  (1,  "1s1"),
    'He': (2,  "1s2"),
    'Li': (3,  "1s2 2s1"),
    'Be': (4,  "1s2 2s2"),
    'B':  (5,  "1s2 2s2 2p1"),
    'C':  (6,  "1s2 2s2 2p2"),
    'N':  (7,  "1s2 2s2 2p3"),
    'O':  (8,  "1s2 2s2 2p4"),
    'F':  (9,  "1s2 2s2 2p5"),
    'Ne': (10, "1s2 2s2 2p6"),
    'Na': (11, "1s2 2s2 2p6 3s1"),
    'Mg': (12, "1s2 2s2 2p6 3s2"),
    'Al': (13, "1s2 2s2 2p6 3s2 3p1"),
    'Si': (14, "1s2 2s2 2p6 3s2 3p2"),
    'P':  (15, "1s2 2s2 2p6 3s2 3p3"),
    'S':  (16, "1s2 2s2 2p6 3s2 3p4"),
    'Cl': (17, "1s2 2s2 2p6 3s2 3p5"),
    'Ar': (18, "1s2 2s2 2p6 3s2 3p6"),
}

def parse_config(config_str):
    """
    Parse stringa tipo '1s2 2s2 2p6 3s2 3p2' -> lista di dict {'n':3, 'l':0/1/2, 'sub':'s','occ':2}
    """
    out = []
    tokens = re.findall(r'(\d)([spdf])(\d+)', config_str.replace(',', ' '))
    lmap = {'s':0, 'p':1, 'd':2, 'f':3}
    for n, sub, occ in tokens:
        out.append({'n': int(n), 'sub': sub, 'l': lmap[sub], 'occ': int(occ)})
    # ordina per n crescente, poi s,p,d,f
    order = {'s':0,'p':1,'d':2,'f':3}
    out.sort(key=lambda t: (t['n'], order[t['sub']]))
    return out

def slater_Zeff_for_ns_np(target_n, target_sub, config, Z):
    """
    Z_eff per un elettrone in ns/np (regole di Slater):
    - Stessa shell (n): 0.35 ciascuno (0.30 per 1s); escludi l'elettrone bersaglio
    - n-1: 0.85 ciascuno
    - n-2 o minori: 1.00 ciascuno
    Ignora d/f nel dettaglio (ok per valenza s/p).
    """
    S = 0.0
    # conta elettroni per shell
    by_shell = {}
    for sh in config:
        by_shell.setdefault(sh['n'], 0)
        by_shell[sh['n']] += sh['occ']

    # stessa shell
    same_n_e = by_shell.get(target_n, 0)
    same_coeff = 0.30 if target_n == 1 else 0.35
    S += same_coeff * max(same_n_e - 1, 0)  # escludi l'elettrone stesso

    # n-1
    S += 0.85 * by_shell.get(target_n - 1, 0)

    # n-2 o meno
    lower = [n for n in by_shell.keys() if n <= target_n - 2]
    S += sum(by_shell[n] for n in lower) * 1.0

    Zeff = max(Z - S, 0.5)  # evita degenerazioni
    return Zeff

# --------- Costruzione degli orbitali reali richiesti ---------

def build_orbitals_atomic(config, Z, only_valence=True, hund_for_p=True):
    """
    Dalla configurazione elettronica costruisce la lista degli orbitali reali da sommare.
    Ogni elemento: dict con chiavi:
        {'label','n','l','type' in {'s','px','py','pz'}, 'occ', 'Zeff'}
    Se only_valence=True, prende solo shell con n massimo.
    Se hund_for_p=True, riempie i p come px,py,pz in ordine ciclico (1 elettrone per orbitale).
    """
    orbitals = []
    if only_valence:
        max_n = max(sh['n'] for sh in config)
        shells = [sh for sh in config if sh['n'] == max_n]
    else:
        shells = config

    for sh in shells:
        n, sub, l, occ = sh['n'], sh['sub'], sh['l'], sh['occ']
        if sub not in ('s','p'):  # qui limitiamo a s/p per visual valenza
            continue
        Zeff = slater_Zeff_for_ns_np(n, sub, config, Z)

        if sub == 's':
            if occ > 0:
                orbitals.append({'label': f'{n}s', 'n': n, 'l': 0, 'type': 's', 'occ': occ, 'Zeff': Zeff})
        elif sub == 'p':
            # distribuisci px, py, pz secondo Hund
            seq = ['px','py','pz','px','py','pz']  # fino a 6
            for i in range(min(occ, 6)):
                orbitals.append({'label': f'{n}p-{seq[i]}', 'n': n, 'l': 1, 'type': seq[i], 'occ': 1, 'Zeff': Zeff})
    return orbitals

def build_orbitals_sp3(n_s, n_p, Z, config_all):
    """
    Quattro ibridi sp3 orientati verso i vertici del tetraedro: (±1,±1,±1)/√3 con numero pari di -.
    Ogni ibrido: ψ = (1/2) s_n + (√3/2)(n_x p_x + n_y p_y + n_z p_z), con R_s(n_s), R_p(n_p).
    Per il raggio usa Z_eff calcolati dalle Slater per ns e np (dalla config completa).
    """
    dirs = np.array([
        [ 1,  1,  1],
        [-1, -1,  1],
        [-1,  1, -1],
        [ 1, -1, -1],
    ], dtype=float)
    dirs = dirs / np.linalg.norm(dirs, axis=1, keepdims=True)

    Zeff_s = slater_Zeff_for_ns_np(n_s, 's', config_all, Z)
    Zeff_p = slater_Zeff_for_ns_np(n_p, 'p', config_all, Z)

    hybrids = []
    for i, v in enumerate(dirs):
        hybrids.append({
            'label': f'sp3_{i+1}',
            'type': 'sp3',
            'n_s': n_s, 'n_p': n_p,
            'Zeff_s': Zeff_s, 'Zeff_p': Zeff_p,
            'dir': v, 'occ': 1  # un elettrone per ibrido (4 in totale)
        })
    return hybrids

# --------- Valutazione campo d'onda reale e densità ---------

def psi_orbital_real(x, y, z, n, l, typ, Zeff):
    r = np.sqrt(x*x + y*y + z*z)
    if typ == 's':
        R = hydrogenic_R_nl(r, n, 0, Zeff)
        ang = ang_s(x,y,z,r)
        return R * ang
    elif typ == 'px':
        R = hydrogenic_R_nl(r, n, 1, Zeff)
        return R * ang_px(x,y,z,r)
    elif typ == 'py':
        R = hydrogenic_R_nl(r, n, 1, Zeff)
        return R * ang_py(x,y,z,r)
    elif typ == 'pz':
        R = hydrogenic_R_nl(r, n, 1, Zeff)
        return R * ang_pz(x,y,z,r)
    else:
        raise ValueError(f"Tipo orbitale non supportato: {typ}")

def psi_sp3(x, y, z, n_s, n_p, Zeff_s, Zeff_p, vdir):
    # ψ = (1/2) s + (√3/2)(v·p), con p_x = R_{n_p,1} * ang_px, etc.
    r = np.sqrt(x*x + y*y + z*z)
    s_part = 0.5 * hydrogenic_R_nl(r, n_s, 0, Zeff_s) * ang_s(x,y,z,r)
    Rp = hydrogenic_R_nl(r, n_p, 1, Zeff_p)
    px = Rp * ang_px(x,y,z,r)
    py = Rp * ang_py(x,y,z,r)
    pz = Rp * ang_pz(x,y,z,r)
    p_dot = vdir[0]*px + vdir[1]*py + vdir[2]*pz
    return s_part + (sqrt(3)/2.0)*p_dot

def total_density(grid, orbitals=None, hybrids=None, smooth_sigma=0.0):
    """
    Somma densità: Σ occ * |ψ|^2 su griglia.
    orbitals: lista di orbitali reali (s, px, py, pz)
    hybrids: lista di ibridi sp3
    """
    x, y, z = grid
    rho = np.zeros_like(x, dtype=np.float64)

    if orbitals:
        for ob in orbitals:
            psi = psi_orbital_real(x, y, z, ob['n'], ob['l'], ob['type'], ob['Zeff'])
            rho += ob['occ'] * (psi*psi)  # |ψ|^2 (reale)

    if hybrids:
        for hb in hybrids:
            psi = psi_sp3(x, y, z, hb['n_s'], hb['n_p'], hb['Zeff_s'], hb['Zeff_p'], hb['dir'])
            rho += hb['occ'] * (psi*psi)

    if smooth_sigma > 0:
        rho = gaussian_filter(rho, sigma=smooth_sigma, mode='nearest')

    return rho

# --------- Mesh export .OBJ (senza dipendenze esterne) ---------

def save_obj(path, verts, faces):
    with open(path, 'w') as f:
        for v in verts:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for tri in faces:
            # OBJ è 1-based
            f.write(f"f {tri[0]+1} {tri[1]+1} {tri[2]+1}\n")

# --------- Plotly opzionale ---------

def try_show_plotly(x, y, z, rho, iso):
    try:
        import plotly.graph_objects as go
    except Exception:
        print("Plotly non installato: salto la preview interattiva.")
        return
    # Plotly Isosurface vuole vettori 1D per assi e un volume 3D
    fig = go.Figure(data=go.Isosurface(
        x=np.repeat(x[:,0,0], y.shape[1]*z.shape[2]),
        y=np.tile(np.repeat(y[0,:,0], z.shape[2]), x.shape[0]),
        z=np.tile(z[0,0,:], x.shape[0]*y.shape[1]),
        value=rho.flatten(order='C'),
        isomin=iso, isomax=rho.max(),
        caps=dict(x_show=False, y_show=False, z_show=False),
        surface_count=1
    ))
    fig.update_layout(scene_aspectmode='data', title=f'Isosuperficie densità (iso={iso})')
    fig.show()

# --------- CLI ---------

def main():
    parser = argparse.ArgumentParser(description="Somma densità orbitali e isosuperficie unica")
    parser.add_argument('--element', type=str, default='Si', help="Simbolo elemento (es. Si) oppure 'Z=14'")
    parser.add_argument('--config', type=str, default=None, help="Config esplicita, es: '1s2 2s2 2p6 3s2 3p2'")
    parser.add_argument('--mode', type=str, default='atomic', choices=['atomic','sp3'],
                        help="atomic = usa configurazione; sp3 = 4 ibridi tetraedrici")
    parser.add_argument('--only-valence', action='store_true', help="Somma solo gli orbitali di valenza")
    parser.add_argument('--grid', type=int, default=128, help="Risoluzione per lato (NxNxN)")
    parser.add_argument('--rmax', type=float, default=8.0, help="Raggio massimo in a0 per coprire la densità")
    parser.add_argument('--iso', type=float, default=5e-4, help="Valore di iso per marching cubes (densità)")
    parser.add_argument('--smooth', type=float, default=0.0, help="Sigma smoothing Gaussian per volume")
    parser.add_argument('--export', type=str, default=None, help="Percorso file .obj da esportare")
    parser.add_argument('--show', action='store_true', help="Preview interattiva (richiede Plotly)")
    args = parser.parse_args()

    # Z & configurazione
    if args.config:
        config = parse_config(args.config)
        if args.element.upper().startswith('Z='):
            Z = int(args.element.split('=')[1])
        else:
            Z = ELEMENTS_SP.get(args.element, (None, None))[0]
            if Z is None:
                raise ValueError("Specifica anche Z=... se usi --config con elemento non in ELEMENTS_SP.")
    else:
        if args.element.upper().startswith('Z='):
            Z = int(args.element.split('=')[1])
            raise ValueError("Per Z puro, serve anche --config con la configurazione elettronica.")
        else:
            Z, cfg = ELEMENTS_SP.get(args.element, (None, None))
            if Z is None:
                raise ValueError(f"Elemento non supportato nella tabella rapida: {args.element}")
            config = parse_config(cfg)

    # Griglia cartesiana
    N = args.grid
    lin = np.linspace(-args.rmax, args.rmax, N)
    x, y, z = np.meshgrid(lin, lin, lin, indexing='ij')

    # Costruisci orbitali/hybrids
    orbitals = None
    hybrids = None
    if args.mode == 'atomic':
        orbitals = build_orbitals_atomic(config, Z, only_valence=args.only_valence, hund_for_p=True)
        labels = ", ".join(o['label'] for o in orbitals)
        print(f"[atomic] Z={Z}  orbitals: {labels}")
    else:
        # tipicamente per Si: n_s=3, n_p=3
        max_n = max(sh['n'] for sh in config)
        # se in valenza abbiamo sia s che p con lo stesso n (es. 3s/3p), usiamo quelli
        n_s = next((sh['n'] for sh in reversed(config) if sh['sub']=='s'), max_n)
        n_p = next((sh['n'] for sh in reversed(config) if sh['sub']=='p'), max_n)
        hybrids = build_orbitals_sp3(n_s=n_s, n_p=n_p, Z=Z, config_all=config)
        print(f"[sp3] Z={Z}  n_s={n_s}  n_p={n_p}")

    # Densità totale
    rho = total_density((x,y,z), orbitals=orbitals, hybrids=hybrids, smooth_sigma=args.smooth)

    # Marching cubes su livello 'iso'
    # spacing per scalare correttamente in coordinate fisiche (qui già in a0)
    spacing = (lin[1]-lin[0], lin[1]-lin[0], lin[1]-lin[0])
    verts, faces, _, _ = marching_cubes(rho, level=args.iso, spacing=spacing)
    print(f"Mesh: {len(verts)} vertici, {len(faces)} triangoli (iso={args.iso})")

    # Esporta .obj
    if args.export:
        save_obj(args.export, verts, faces)
        print(f"Salvato OBJ -> {args.export}")

    # Preview interattiva opzionale 
    if args.show:
        try_show_plotly(x, y, z, rho, args.iso)

if __name__ == '__main__':
    main() 