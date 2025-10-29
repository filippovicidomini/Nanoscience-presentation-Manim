# ====================================================
from pathlib import Path
import sys, numpy as np, trimesh
from trimesh.transformations import rotation_matrix, translation_matrix
from trimesh.geometry import align_vectors