# precompute_orbital_mesh.py
import numpy as np
from skimage import measure
import trimesh
import sys

# SciPy: usa nuova API se c'è, altrimenti fallback
try:
    from scipy.special import sph_harm_y as sphY    # SciPy >= 1.15
    def Y_real(l, m, theta, phi):
        Y = sphY(l, m, theta, phi)
        if l == 1 and m == +1:   # p_x
            return np.sqrt(2)*np.real(Y)
        if l == 1 and m == -1:   # p_y
            return np.sqrt(2)*np.imag(Y)
        return np.real(Y)        # s, p_z
except ImportError:
    from scipy.special import sph_harm as sphY      # SciPy < 1.15
    def Y_real(l, m, theta, phi):
        Y = sphY(m, l, phi, theta)                  # firma vecchia
        if l == 1 and m == +1:
            return np.sqrt(2)*np.real(Y)
        if l == 1 and m == -1:
            return np.sqrt(2)*np.imag(Y)
        return np.real(Y)

# -------- Parametri principali --------
n, l, m = 5, 2, 0        # 3p_z (cambia qui: p_x -> (3,1,+1), p_y -> (3,1,-1), 3s -> (3,0,0))
a0      = 1.0
grid_N  = 128            # 112–160: compromesso qualità/tempo
Rmax    = 10.0           # estensione spaziale
iso     = 3e-4           # livello isosuperficie (aumenta se la mesh è troppo grossa)

# -------- Griglia cartesiana --------
xs = np.linspace(-Rmax, Rmax, grid_N)
ys = np.linspace(-Rmax, Rmax, grid_N)
zs = np.linspace(-Rmax, Rmax, grid_N)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
r = np.sqrt(X*X + Y*Y + Z*Z) + 1e-12
theta = np.arccos(np.clip(Z/r, -1.0, 1.0))
phi   = (np.arctan2(Y, X) + 2*np.pi) % (2*np.pi)

# -------- Radiali idrogeno-like (qualitativi ma corretti nei nodi) --------
x = r/a0
if (n, l) == (3, 0):      # 3s: due nodi radiali
    Rnl = (1 - (2*x)/3 + (2*x**2)/27) * np.exp(-x/3)
elif (n, l) == (3, 1):    # 3p: un nodo radiale
    Rnl = (x) * (1 - x/6) * np.exp(-x/3)
else:                     # fallback semplice
    Rnl = np.exp(-r/(n*a0))

Ylm_real = Y_real(l, m, theta, phi)
psi2 = (Rnl * Ylm_real)**2

# -------- Marching Cubes --------
spacing = (xs[1]-xs[0], ys[1]-ys[0], zs[1]-zs[0])
verts, faces, normals, _ = measure.marching_cubes(psi2, level=iso, spacing=spacing)

mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)

# -------- Semplificazione (Open3D) --------
import open3d as o3d
o3 = o3d.geometry.TriangleMesh(
    o3d.utility.Vector3dVector(np.asarray(mesh.vertices)),
    o3d.utility.Vector3iVector(np.asarray(mesh.faces))
)
o3.compute_vertex_normals()
target = max(3000, int(len(mesh.faces) * 0.25))   # ~25% o minimo 3k triangoli
o3 = o3.simplify_quadric_decimation(target_number_of_triangles=target)

V = np.asarray(o3.vertices)
F = np.asarray(o3.triangles)

out = f"orbital_n{n}_l{l}_m{m}.obj"
trimesh.Trimesh(vertices=V, faces=F, process=False).export(out)
print(f"Saved {out} — verts: {len(V)}, faces: {len(F)}")