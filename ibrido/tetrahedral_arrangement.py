import numpy as np
import trimesh
from scipy.spatial.transform import Rotation as R

# === File di input ===
input_file = "orbital_n4_l2_m0_asym0.30.obj"
output_file = "tetrahedral_touch_calibrated.obj"

# === Direzioni tetraedriche ===
directions = np.array([
    [ 1,  1,  1],
    [ 1, -1, -1],
    [-1,  1, -1],
    [-1, -1,  1]
], dtype=float)
directions /= np.linalg.norm(directions, axis=1)[:, None]

# === Carica il modello base ===
base = trimesh.load(input_file, process=False)

# --- centra il modello sull’origine ---
base.vertices -= base.center_mass

# --- misura le dimensioni ---
bbox = base.bounds
size = bbox[1] - bbox[0]
extent = np.linalg.norm(size) / 2        # raggio caratteristico
print(f"Bounding box size: {size},  extent ~ {extent:.3f}")

# === calcola la traslazione ideale ===
# 0.25–0.35 dell'estensione dà contatto dei vertici al centro
touch_factor = 0.8
t = extent * touch_factor
print(f"Using translation distance t = {t:.3f}")

# === Rotazione e traslazione ===
def orient_and_translate(mesh, direction, t):
    z = np.array([0, 0, 1])
    target = direction / np.linalg.norm(direction)
    if np.allclose(z, target):
        rotated = mesh.vertices.copy()
    else:
        axis = np.cross(z, target)
        axis /= np.linalg.norm(axis)
        angle = np.arccos(np.dot(z, target))
        rot = R.from_rotvec(axis * angle)
        rotated = rot.apply(mesh.vertices)
    translated = rotated + target * t
    return trimesh.Trimesh(vertices=translated, faces=mesh.faces, process=False)

# === Crea le quattro copie ===
meshes = []
for d in directions:
    orb = orient_and_translate(base, d, t)
    meshes.append(orb)

# --- (opzionale) nucleo centrale ---
core = trimesh.creation.icosphere(subdivisions=3, radius=extent * 0.1)
meshes.append(core)

# === Unisci ed esporta ===
combined = trimesh.util.concatenate(meshes)
combined.export(output_file)
print(f"Saved {output_file}")