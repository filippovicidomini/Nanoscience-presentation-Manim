import numpy as np
from skimage import measure
import trimesh
import open3d as o3d
from scipy.special import sph_harm

# --- Funzione radiale approssimata (4d) ---
def R_nl(n, l, r, a0=1.0):
    x = r / a0
    if (n, l) == (4, 2):
        return (x**2) * np.exp(-x/3.0)
    return np.exp(-x/n)

# --- Armoniche sferiche reali ---
def Y_real(l, m, theta, phi):
    Y = sph_harm(m, l, phi, theta)
    if m > 0:
        return np.sqrt(2) * np.real(Y)
    elif m < 0:
        return np.sqrt(2) * np.imag(Y)
    else:
        return np.real(Y)

# --- Parametri principali ---
n, l, m = 4, 2, 0
a0 = 1.00
grid_N = 160
Rmax = 20.0
iso_fraction = 0.32      # frazione del valore massimo per la superficie
asym = 0.3              # asimmetria dei lobi (1 = simmetrico)

# --- Griglia cartesiana ---
xs = np.linspace(-Rmax, Rmax, grid_N)
ys = np.linspace(-Rmax, Rmax, grid_N)
zs = np.linspace(-Rmax, Rmax, grid_N)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")

# --- Coordinate sferiche ---
r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12
theta = np.arccos(np.clip(Z / r, -1.0, 1.0))
phi = (np.arctan2(Y, X) + 2*np.pi) % (2*np.pi)

# --- Funzione d'onda reale ψ(r,θ,φ) ---
Rpart = R_nl(n, l, r, a0)
Ypart = Y_real(l, m, theta, phi)
direction_factor = 1 + (1 - asym) * np.cos(theta)
psi = Rpart * Ypart * direction_factor

# --- Densità di probabilità ---
rho = psi**2
rho -= rho.min()
rho /= rho.max()

stretch_z = 2    # >1 allunga lungo z
stretch_xy = 1   # <1 stringe sui lati

X_scaled = X * stretch_xy
Y_scaled = Y * stretch_xy
Z_scaled = Z * stretch_z
# --- Marching Cubes ---
iso = iso_fraction
spacing =  (
    (xs[1] - xs[0]) * stretch_xy,
    (ys[1] - ys[0]) * stretch_xy,
    (zs[1] - zs[0]) * stretch_z
)
verts, faces, _, _ = measure.marching_cubes(rho, level=iso, spacing=spacing)
mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)

# --- Semplificazione della mesh ---
o3mesh = o3d.geometry.TriangleMesh(
    o3d.utility.Vector3dVector(np.asarray(mesh.vertices)),
    o3d.utility.Vector3iVector(np.asarray(mesh.faces))
)
o3mesh.compute_vertex_normals()
target = max(4000, int(len(mesh.faces) * 0.75))
o3mesh = o3mesh.simplify_quadric_decimation(target_number_of_triangles=target)

V = np.asarray(o3mesh.vertices)
F = np.asarray(o3mesh.triangles)

# --- Esporta il modello 3D ---
out = f"orbital_n{n}_l{l}_m{m}_asym{asym:.2f}.obj"
trimesh.Trimesh(vertices=V, faces=F, process=False).export(out)
print(f"Saved {out} — verts: {len(V)}, faces: {len(F)}")