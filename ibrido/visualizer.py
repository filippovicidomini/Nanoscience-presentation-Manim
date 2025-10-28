# ibrido/ibrido.py

# ... [previous imports and code] ...

# Unified wrapper for sph_harm (new/old API)
def _sphY(l, m, theta, phi):
    if HAS_SPH_Y:
        return sphY_new(l, m, theta, phi)  # (l,m,theta,phi)
    else:
        return sphY_old(m, l, phi, theta)  # (m,l,phi,theta)

def Y_real(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    if not SCIPY_OK:
        raise RuntimeError("SciPy è richiesto per sferiche e Laguerre.")
    Y = _sphY(l, m, theta, phi)
    return np.real(Y)

def psi_field(grid: GridSpec, which: str, use_slater=True):
    assert which in ORBITALS
    n, l, m = ORBITALS[which]
    X, Y, Z, r, th, ph, spacing = make_grid(grid.N, grid.Rmax)

    Z_eff = Z_EFF_SI_3SP if (use_slater and n == 3 and l in (0,1)) else Z_SI
    Rnl = R_nl_hydrogenic(n, l, r, Z_eff)

    if l == 1:  # p-orbitals with explicit real combinations
        Yp1 = _sphY(1, +1, th, ph)
        Y10 = _sphY(1, 0, th, ph)
        if which == '3px':
            Ylm = np.sqrt(2.0) * np.real(Yp1)
        elif which == '3py':
            Ylm = np.sqrt(2.0) * np.imag(Yp1)
        else:  # '3pz'
            Ylm = np.real(Y10)
    else:
        Ylm = np.real(_sphY(l, m, th, ph))

    psi = Rnl * Ylm
    return psi, spacing

# ---- sp3 hybrids ----
TETRA_DIRS = [
    np.array([+1, +1, +1], dtype=float),
    np.array([+1, -1, -1], dtype=float),
    np.array([-1, +1, -1], dtype=float),
    np.array([-1, -1, +1], dtype=float),
]

def psi_sp3(direction: np.ndarray, grid: GridSpec, s_weight=0.5, use_slater=True):
    # Normalize direction
    u = np.array(direction, dtype=float)
    u /= np.linalg.norm(u)

    X, Y, Z, r, th, ph, spacing = make_grid(grid.N, grid.Rmax)

    # Radiali per 3s e 3p con lo stesso Z_eff
    Z_eff = Z_EFF_SI_3SP if use_slater else Z_SI
    R3s = R_nl_hydrogenic(3, 0, r, Z_eff)
    R3p = R_nl_hydrogenic(3, 1, r, Z_eff)

    # Armoniche reali base
    Y00 = _sphY(0, 0, th, ph)  # reale
    Yp1 = _sphY(1, +1, th, ph)
    Y10 = _sphY(1, 0, th, ph)

    # Costruisci componente p direzionale: u_x px + u_y py + u_z pz
    Ypx = np.sqrt(2.0) * np.real(Yp1)   # ∝ sinθ cosφ
    Ypy = np.sqrt(2.0) * np.imag(Yp1)   # ∝ sinθ sinφ
    Ypz = np.real(Y10)                  # ∝ cosθ
    Y_pu = u[0]*Ypx + u[1]*Ypy + u[2]*Ypz

    # Mixing coefficienti: s_weight per s, p_weight dedotto per normalizzazione
    p_weight = math.sqrt(1.0 - s_weight**2)  # mantiene ||ψ||≈1 se base ortonormale

    psi = s_weight * R3s * np.real(Y00) + p_weight * R3p * Y_pu
    return psi, spacing

def export_sp3_tetra(grid: GridSpec, iso: float | None, auto_iso_q: float | None,
                     use_slater=True, do_simplify=True, s_weight=0.5):
    outputs = []
    for idx, u in enumerate(TETRA_DIRS, start=1):
        psi, spacing = psi_sp3(u, grid, s_weight=s_weight, use_slater=use_slater)
        # auto-iso
        if iso is None:
            q = 0.995 if auto_iso_q is None else float(auto_iso_q)
            level = float(np.quantile(np.abs(psi), q)) * 0.75
        else:
            level = float(iso)
        mP, mN = mc_signed(psi, level, spacing)
        if do_simplify:
            mP = simplify(mP)
            mN = simplify(mN)
        tag = ["111","1-1-1","-11-1","-1-1 1"][idx-1]
        base = f"Si_sp3_{tag.replace(' ', '')}"
        outP = f"{base}_pos.obj"; mP.export(outP)
        outN = f"{base}_neg.obj"; mN.export(outN)
        print(f"[OK] sp3 {tag}: +{len(mP.vertices)}/{len(mP.faces)}  -{len(mN.vertices)}/{len(mN.faces)}  |ψ|={level:.3e}")
        outputs.append((outP, outN))
    return outputs

def parse_args():
    p = argparse.ArgumentParser(description="Calcolo orbitali silicio")
    p.add_argument('which', choices=['3s','3px','3py','3pz','all','sp3tetra'])
    # ... rest of parse_args unchanged ...

def main():
    args = parse_args()
    grid = GridSpec(N=args.N, Rmax=args.Rmax)
    if args.which == 'sp3tetra':
        export_sp3_tetra(
            grid,
            iso=args.iso,
            auto_iso_q=args.auto_iso,
            use_slater=(not args.no_slater),
            do_simplify=(not args.no_simplify),
            s_weight=0.5,
        )
        return
    # ... rest of main unchanged ...