"""gridio — minimal reader for gRASPA 3D density grids (.cube and .vtk).

Dependency-light: numpy only.

Both writers describe the SAME scalar field (adsorbate-atom number density,
molecules-of-atoms / Angstrom^3) but with different geometry conventions:

  * .vtk  DATASET STRUCTURED_POINTS : SPACING is in Angstrom, values are written
          x-fastest (VTK point order: idx = x + nx*(y + ny*z)). Assumes an
          orthorhombic cell.
  * .cube Gaussian cube             : voxel vectors + origin are in BOHR, values
          are written z-fastest (cube order: idx = z + nz*(y + ny*x)). Triclinic
          correct (full voxel vectors).

load_grid() returns everything in a single canonical convention so the three CLI
tools can treat the two formats interchangeably:

  density : float64 ndarray, shape (nx, ny, nz), axes = (x, y, z)
  voxvol  : scalar voxel volume in Angstrom^3
  meta    : dict with 'format', 'dims', 'path', and (vtk) 'spacing'
"""
from __future__ import annotations
import os
import numpy as np

BOHR_TO_ANG = 0.52917721067  # CODATA-ish; only used for .cube geometry


def _read_vtk(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    dims = spacing = None
    data_start = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("DIMENSIONS"):
            dims = tuple(int(v) for v in s.split()[1:4])
        elif s.startswith("SPACING"):
            spacing = tuple(float(v) for v in s.split()[1:4])
        elif s.startswith("LOOKUP_TABLE"):
            data_start = i + 1
            break
    if dims is None or spacing is None or data_start is None:
        raise ValueError(f"{path}: not a STRUCTURED_POINTS VTK with SPACING/LOOKUP_TABLE")
    nx, ny, nz = dims
    vals = np.fromstring(" ".join(lines[data_start:]), sep=" ")
    if vals.size != nx * ny * nz:
        raise ValueError(f"{path}: expected {nx*ny*nz} values, got {vals.size}")
    # VTK is x-fastest -> reshape (z,y,x) then transpose to canonical (x,y,z)
    dens = vals.reshape((nz, ny, nx)).transpose(2, 1, 0)
    voxvol = float(spacing[0] * spacing[1] * spacing[2])
    return np.ascontiguousarray(dens, dtype=np.float64), voxvol, {
        "format": "vtk", "dims": (nx, ny, nz), "spacing": spacing, "path": path}


def _read_cube(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    natoms = int(float(lines[2].split()[0]))
    nabs = abs(natoms)
    def _n_and_vec(line):
        p = line.split()
        return int(p[0]), np.array(list(map(float, p[1:4])))
    nx, v1 = _n_and_vec(lines[3])
    ny, v2 = _n_and_vec(lines[4])
    nz, v3 = _n_and_vec(lines[5])
    voxvol_bohr = abs(np.linalg.det(np.stack([v1, v2, v3])))
    voxvol = float(voxvol_bohr * BOHR_TO_ANG ** 3)
    body = " ".join(lines[6 + nabs:])
    vals = np.fromstring(body, sep=" ")
    # a cube may carry NVAL>1 per point; gRASPA writes 1, but be defensive
    if vals.size != nx * ny * nz:
        # trim any trailing DSET_IDS / extras if size is a clean multiple
        if vals.size > nx * ny * nz:
            vals = vals[: nx * ny * nz]
        else:
            raise ValueError(f"{path}: expected {nx*ny*nz} values, got {vals.size}")
    # cube is z-fastest -> reshape (x,y,z) is already canonical
    dens = vals.reshape((nx, ny, nz))
    return np.ascontiguousarray(dens, dtype=np.float64), voxvol, {
        "format": "cube", "dims": (nx, ny, nz), "path": path}


def load_grid(path):
    """Load a .cube or .vtk gRASPA density grid -> (density[x,y,z], voxvol_A3, meta)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".vtk":
        return _read_vtk(path)
    if ext == ".cube":
        return _read_cube(path)
    # sniff
    head = open(path).read(256)
    if "STRUCTURED_POINTS" in head or head.startswith("# vtk"):
        return _read_vtk(path)
    return _read_cube(path)


def find_grid(run_dir, prefer="cube"):
    """Return the density grid path under a run dir (DensityGrid/System_0/*).

    Defaults to the .cube: its voxel vectors are triclinic-correct, whereas the
    .vtk STRUCTURED_POINTS geometry assumes an orthorhombic cell (its declared
    SPACING over-states the voxel volume for non-orthorhombic frameworks, so
    re-integrating a .vtk from a triclinic MOF overcounts by |a||b||c|/V_cell).
    The scalar density VALUES are identical in both files.
    """
    exts = (".vtk", ".cube") if prefer == "vtk" else (".cube", ".vtk")
    for root, _dirs, files in os.walk(run_dir):
        if os.path.basename(root) != "System_0" and "DensityGrid" not in root:
            continue
        for ext in exts:
            hits = sorted(f for f in files if f.lower().endswith(ext) and "Density" in f)
            if hits:
                return os.path.join(root, hits[0])
    # fallback: any grid anywhere under run_dir
    for ext in exts:
        for root, _dirs, files in os.walk(run_dir):
            hits = sorted(f for f in files if f.lower().endswith(ext))
            if hits:
                return os.path.join(root, hits[0])
    return None
