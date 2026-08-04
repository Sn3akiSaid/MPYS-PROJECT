#!/usr/bin/env python3
"""3D Fermi surfaces with optional Berry-curvature arrow glyphs.

Consolidates the two original scripts in this directory:

  Fermi_Surface.py      -> --preset arrows  (isosurface + Berry-curvature glyphs)
  2nd_fermi_surface.py  -> --preset plain   (isosurface only)

Both originals are kept on disk as the provenance record for the published
figures; this script is the maintained version.

Run from this directory (data paths are resolved relative to the script when not
found relative to the working directory):

    python3 fermi_surface_plot.py                       # default figure
    python3 fermi_surface_plot.py --preset plain \
        --data-file k_surface_fermi_energies_By_0.01_part_1.dat --energy-offset 0.03

DEFAULTS REPRODUCE THE CURRENT OUTPUT of Fermi_Surface.py. Corrections that would
change the rendered image are opt-in; see "Known issues" below and the banner the
script prints when any of them is enabled.

Known issues in the original, and how they are handled here
-----------------------------------------------------------
Fixed unconditionally (verified to leave the render bit-identical):

  * ``griddata(method='nearest')`` replaced by a separable nearest-index lookup
    into the reshaped data cube. Verified array-identical on both the native grid
    and the legacy 100^3 grid, and ~80x faster. The old call is still available
    via ``--interp griddata`` as an oracle.
  * ``energy`` and ``energy13`` both read column 7 in the original, so the
    "band 13" distinction was already lost. Columns are now explicit
    (``--energy-col`` / ``--ref-col``); the defaults match the original.
  * ``smooth(feature_angle=..., edge_angle=...)`` -- these are degrees in VTK
    (valid 0-180) and the original passed 10000, but ``feature_smoothing`` and
    ``boundary_smoothing`` default to False so both arguments were never read.
    Dropped. ``n_iter`` and ``relaxation_factor`` are preserved exactly.
  * Unused ``import vtk``; ``pv.Report()`` env dump now behind ``--report``;
    dead ``vector_norms`` computation removed.

Opt-in because they change the figure:

  * ``--grid native``: the default samples a 100^3 grid over hardcoded bounds
    that extend beyond the simulated region -- for the shipped default dataset
    46.7% of grid cells lie outside the data, and ~1650 of those fall within
    contouring distance of E_F, producing isosurface where nothing was computed.
    ``native`` restricts the grid to the data's own axes.
  * ``--sampler fixed``: the original bins with ``linspace(min, max, n)`` and a
    ``<`` upper test, silently dropping points sitting exactly at each axis
    maximum (841 per axis for the default dataset), and selects ``indices[0]``
    (arbitrary array order) within each cell.
  * ``--scale-by-magnitude``: the original renders every arrow the same length,
    discarding |Omega| (which spans ~5 decades).

Data formats
------------
9 columns: kx ky kz  Omega_x Omega_y Omega_z  |Omega|  E  alpha
5 columns: kx ky kz  E1 E2

Both are regular cubic k-space lattices, but the dimension is FILE-SPECIFIC
(29^3, 19^3, 51^3 ... ). It is always derived at runtime -- never hardcoded --
and ``infer_cube_dim`` hard-fails rather than guessing, because a wrong reshape
scrambles k-space with no visual tell.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------
# Presets: the visual settings that differed between the two original scripts.
# --------------------------------------------------------------------------
PRESETS = {
    # Fermi_Surface.py as committed.
    "arrows": dict(
        data_file="BerryUnpertTVB_5.dat",
        energy_offset=0.01,
        bounds=(-0.08, 0.08, -0.08, 0.08, 0.42, 0.51),
        data_driven_bounds=False,
        dimension=100,
        window=(1440, 1080),
        box_opacity=0.2,
        surface_opacity=1.0,
        latex_titles=False,
        azimuth=170.0,
        focal_dz=-0.02,
        outline=False,
        arrows=True,
        num_bins=16,
    ),
    # 2nd_fermi_surface.py as committed.
    "plain": dict(
        data_file="k_surface_fermi_energies_By_0.01_part_1.dat",
        energy_offset=0.03,
        bounds=None,
        data_driven_bounds=True,
        dimension=50,
        window=(800, 600),
        box_opacity=0.5,
        surface_opacity=0.8,
        latex_titles=True,
        azimuth=None,
        focal_dz=None,
        outline=True,
        arrows=False,
        num_bins=16,
    ),
    # Reconstructed from the May-2025 WSM figures. APPROXIMATE: the exact
    # energy_offset / num_bins / window for that run are not recoverable from
    # the repo (see CLAUDE.md). Expect to iterate against the PNGs.
    "wsm": dict(
        data_file="berryinWSM_11.dat",
        energy_offset=0.01,
        bounds=None,
        data_driven_bounds=True,
        dimension=100,
        window=(1080, 1080),
        box_opacity=0.2,
        surface_opacity=1.0,
        latex_titles=False,
        azimuth=135.0,
        focal_dz=-0.02,
        outline=False,
        arrows=True,
        num_bins=16,
    ),
}


@dataclass
class KData:
    kx: np.ndarray
    ky: np.ndarray
    kz: np.ndarray
    energy: np.ndarray
    ref_energy: np.ndarray
    omega: np.ndarray | None  # (N, 3) or None for 5-column files
    omega_mag: np.ndarray | None

    def __len__(self) -> int:
        return len(self.kx)


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
def resolve_data_path(name: str) -> str:
    """Look for the file relative to cwd, then to this script's directory."""
    if os.path.isfile(name):
        return name
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    if os.path.isfile(here):
        return here
    raise SystemExit(
        f"error: data file {name!r} not found relative to the working directory "
        f"({os.getcwd()}) or to the script directory "
        f"({os.path.dirname(os.path.abspath(__file__))})"
    )


def load_kspace_data(path, energy_col=None, ref_col=None) -> KData:
    """Load a 9-column Berry file or a 5-column Fermi-energy file."""
    raw = np.loadtxt(resolve_data_path(path))
    if raw.ndim != 2:
        raise SystemExit(f"error: {path!r} did not parse as a 2D table")
    ncols = raw.shape[1]

    if ncols >= 9:
        e_col = 7 if energy_col is None else energy_col
        r_col = e_col if ref_col is None else ref_col
        return KData(
            kx=raw[:, 0], ky=raw[:, 1], kz=raw[:, 2],
            energy=raw[:, e_col], ref_energy=raw[:, r_col],
            omega=raw[:, 3:6], omega_mag=raw[:, 6],
        )
    if ncols >= 4:
        # The original 2nd_fermi_surface.py used column 3 for both energies even
        # though column 4 is a genuinely different band. Default preserves that.
        e_col = 3 if energy_col is None else energy_col
        r_col = e_col if ref_col is None else ref_col
        return KData(
            kx=raw[:, 0], ky=raw[:, 1], kz=raw[:, 2],
            energy=raw[:, e_col], ref_energy=raw[:, r_col],
            omega=None, omega_mag=None,
        )
    raise SystemExit(f"error: {path!r} has {ncols} columns; expected >= 4")


def infer_cube_dim(kx, ky, kz) -> int:
    """Derive n for an n^3 regular lattice, or raise.

    Hard-fails rather than guessing: a wrong reshape silently scrambles k-space.
    """
    ux, uy, uz = np.unique(kx), np.unique(ky), np.unique(kz)
    n = len(ux)
    if not (len(uy) == len(uz) == n):
        raise ValueError(
            f"not a cubic lattice: {len(ux)} x {len(uy)} x {len(uz)} unique "
            f"coordinates. Use --interp griddata."
        )
    if n ** 3 != len(kx):
        raise ValueError(
            f"not a complete lattice: {n}^3 = {n**3} but {len(kx)} rows. "
            f"Use --interp griddata."
        )
    # C-order check: kx must be constant across the two fastest axes.
    cube = kx.reshape(n, n, n)
    if not np.allclose(cube, cube[:, :1, :1]):
        raise ValueError(
            "lattice is not in C-order (kx slowest, kz fastest). "
            "Use --interp griddata."
        )
    return n


# --------------------------------------------------------------------------
# Grid construction
# --------------------------------------------------------------------------
def _nearest_index(query, axis_values):
    """Index of the nearest axis value for each query point (axis is sorted)."""
    j = np.clip(np.searchsorted(axis_values, query), 1, len(axis_values) - 1)
    left_closer = np.abs(query - axis_values[j - 1]) <= np.abs(query - axis_values[j])
    return np.where(left_closer, j - 1, j)


def build_energy_grid(data: KData, cfg):
    """Return (X, Y, Z, values) for contouring.

    ``native``  -- the data's own axes; no resampling, no extrapolation.
    ``legacy``  -- resample onto ``dimension^3`` over ``bounds``, reproducing the
                   original's nearest-neighbour behaviour (including its
                   extrapolation outside the data) exactly.
    """
    if cfg.grid == "native" or cfg.interp == "reshape":
        n = infer_cube_dim(data.kx, data.ky, data.kz)
        ux, uy, uz = np.unique(data.kx), np.unique(data.ky), np.unique(data.kz)
        cube = data.energy.reshape(n, n, n)

    if cfg.grid == "native":
        X, Y, Z = np.meshgrid(ux, uy, uz, indexing="ij")
        return X, Y, Z, cube

    x0, x1, y0, y1, z0, z1 = cfg.bounds
    d = cfg.dimension
    x_lin = np.linspace(x0, x1, d)
    y_lin = np.linspace(y0, y1, d)
    z_lin = np.linspace(z0, z1, d)
    X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing="ij")

    if cfg.interp == "griddata":
        from scipy.interpolate import griddata
        vals = griddata(
            np.column_stack((data.kx, data.ky, data.kz)),
            data.energy, (X, Y, Z), method="nearest",
        )
        vals = np.nan_to_num(vals, nan=float(np.min(data.ref_energy)))
        return X, Y, Z, vals

    # Separable nearest-index lookup. Nearest-neighbour on a rectilinear source
    # is separable per axis, so this is identical to griddata(method='nearest')
    # including its out-of-hull clamping -- verified array-equal on both grids.
    vals = cube[np.ix_(
        _nearest_index(x_lin, ux),
        _nearest_index(y_lin, uy),
        _nearest_index(z_lin, uz),
    )]
    return X, Y, Z, vals


# --------------------------------------------------------------------------
# Glyph point sampling
# --------------------------------------------------------------------------
def grid_sample_3d(x, y, z, vx, vy, vz, num_bins, mode="legacy", kz_ranges=None):
    """Thin points to roughly one per cell of a num_bins^3 lattice.

    ``legacy`` reproduces the original exactly, including dropping points that
    sit on the upper edge of each axis and picking the first point per cell in
    array order. ``fixed`` closes the top bin and picks the point nearest each
    cell centroid.
    """
    if kz_ranges:
        m = np.zeros_like(z, dtype=bool)
        for lo, hi in kz_ranges:
            m |= (z >= lo) & (z <= hi)
        if not np.any(m):
            print("No points found in the specified kz ranges")
            empty = np.array([])
            return (empty,) * 6
        x, y, z, vx, vy, vz = x[m], y[m], z[m], vx[m], vy[m], vz[m]

    edges = [np.linspace(np.min(a), np.max(a), num_bins) for a in (x, y, z)]

    # Bin index per point. np.digitize(right=False) gives 1..len(edges) for
    # values in [edges[0], edges[-1]]; shift to 0-based cell indices.
    ix = np.digitize(x, edges[0]) - 1
    iy = np.digitize(y, edges[1]) - 1
    iz = np.digitize(z, edges[2]) - 1

    ncell = num_bins - 1
    if mode == "legacy":
        # Original used a strict '<' upper bound, so points at the axis maximum
        # land in cell index ncell and are dropped.
        valid = ((ix >= 0) & (ix < ncell)
                 & (iy >= 0) & (iy < ncell)
                 & (iz >= 0) & (iz < ncell))
    else:
        # Fold the top edge back into the last cell.
        ix = np.clip(ix, 0, ncell - 1)
        iy = np.clip(iy, 0, ncell - 1)
        iz = np.clip(iz, 0, ncell - 1)
        valid = np.ones(len(x), dtype=bool)

    idx = np.flatnonzero(valid)
    if idx.size == 0:
        empty = np.array([])
        return (empty,) * 6

    flat = (ix[idx] * ncell + iy[idx]) * ncell + iz[idx]

    if mode == "legacy":
        # First point per cell in array order, matching the original's
        # indices[0]. Iterating cells in (i, j, k) order means sorting by cell
        # id with a stable sort and taking the first of each run.
        order = np.argsort(flat, kind="stable")
    else:
        # Nearest to cell centroid, deterministic and independent of row order.
        cx = 0.5 * (edges[0][ix[idx]] + edges[0][ix[idx] + 1])
        cy = 0.5 * (edges[1][iy[idx]] + edges[1][iy[idx] + 1])
        cz = 0.5 * (edges[2][iz[idx]] + edges[2][iz[idx] + 1])
        d2 = (x[idx] - cx) ** 2 + (y[idx] - cy) ** 2 + (z[idx] - cz) ** 2
        order = np.lexsort((d2, flat))

    sorted_flat = flat[order]
    first = np.ones(len(sorted_flat), dtype=bool)
    first[1:] = sorted_flat[1:] != sorted_flat[:-1]
    pick = idx[order[first]]

    return x[pick], y[pick], z[pick], vx[pick], vy[pick], vz[pick]


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
def render(data: KData, X, Y, Z, values, fermi_energy, filt, cfg):
    import pyvista as pv

    if cfg.report:
        print(pv.__version__)
        print(pv.Report())

    grid = pv.StructuredGrid(X, Y, Z)
    grid["energy"] = values.flatten(order="F")

    print(f"Extracting iso-surface at energy = {fermi_energy:.6f} eV")
    surface = grid.contour(isosurfaces=[fermi_energy], scalars="energy",
                           method="contour")
    surface = surface.smooth(n_iter=cfg.smooth_iter,
                             relaxation_factor=cfg.relaxation)

    print("Generating visualization...")
    plotter = pv.Plotter(window_size=list(cfg.window), off_screen=cfg.off_screen,
                         line_smoothing=True, polygon_smoothing=True)
    plotter.add_mesh(pv.Box(bounds=grid.bounds), color="black",
                     opacity=cfg.box_opacity, show_edges=None, edge_opacity=0.0)
    plotter.add_mesh(surface, color="red", opacity=cfg.surface_opacity,
                     specular=0.9, specular_power=50, smooth_shading=True,
                     show_scalar_bar=False, silhouette=True, lighting=True,
                     diffuse=0.9)

    if cfg.outline:
        plotter.add_mesh(grid.outline(), color="gray", opacity=0.25, line_width=2)

    titles = ((r"$k_x(\AA^{-1})$", r"$k_y(\AA^{-1})$", r"$k_z(\AA^{-1})$")
              if cfg.latex_titles else ("", "", ""))
    plotter.show_bounds(
        grid="back", location="outer", color="Black",
        show_xaxis=True, show_yaxis=True, show_zaxis=True,
        ticks="outside", all_edges=False, fmt="%.2f",
        xtitle=titles[0], ytitle=titles[1], ztitle=titles[2],
        font_size=20 if not cfg.latex_titles else 12, font_family="arial",
        n_xlabels=3, n_ylabels=3, n_zlabels=3,
    )
    if not cfg.latex_titles:
        plotter.renderer.cube_axes_actor.label_offset = 40.0

    plotter.add_light(pv.Light(position=(1, 1, 1), focal_point=(0, 0, 0),
                               color=[1, 1, 1], intensity=0.8))
    plotter.add_light(pv.Light(position=(-1, -1, -1), focal_point=(0, 0, 0),
                               color=[0.5, 0.5, 0.7], intensity=0.3))

    if cfg.arrows and data.omega is not None:
        add_berry_arrows(plotter, data, surface, filt, cfg)

    plotter.enable_anti_aliasing("ssaa")
    if cfg.azimuth is not None:
        plotter.camera.azimuth = cfg.azimuth
    if cfg.focal_dz is not None:
        fp = list(plotter.camera.focal_point)
        fp[2] += cfg.focal_dz
        plotter.camera.focal_point = fp

    if cfg.screenshot:
        plotter.show(screenshot=cfg.screenshot,
                     auto_close=not cfg.off_screen)
        print(f"Wrote {cfg.screenshot}")
    else:
        plotter.show()


def add_berry_arrows(plotter, data: KData, surface, filt, cfg):
    import pyvista as pv
    from scipy.spatial import KDTree

    print("Adding omega vectors to visualization...")
    surface_points = np.array(surface.points)
    if surface_points.size == 0:
        print("warning: isosurface is empty; no arrows drawn")
        return

    tree = KDTree(surface_points)
    pts = np.column_stack((filt.kx, filt.ky, filt.kz))
    distances, _ = tree.query(pts, k=1)
    near = distances <= cfg.max_distance

    print(f"Original filtered points: {len(filt.kx)}")
    print(f"Points on Fermi surface: {int(near.sum())}")

    sx, sy, sz, ox, oy, oz = grid_sample_3d(
        filt.kx[near], filt.ky[near], filt.kz[near],
        filt.omega[near, 0], filt.omega[near, 1], filt.omega[near, 2],
        cfg.num_bins, mode=cfg.sampler, kz_ranges=cfg.kz_ranges,
    )
    if len(sx) == 0:
        print("warning: sampler selected no points; no arrows drawn")
        return

    cloud = pv.PolyData(np.column_stack((sx, sy, sz)))
    vectors = np.column_stack((ox, oy, oz))
    cloud["omega_vectors"] = vectors

    scale_arg = False
    if cfg.scale_by_magnitude:
        cloud["omega_magnitude"] = np.linalg.norm(vectors, axis=1)
        scale_arg = "omega_magnitude"

    arrow = pv.Arrow(start=(-0.5, 0, 0), direction=(1, 0, 0), tip_length=0.65,
                     tip_radius=0.19, tip_resolution=50, shaft_radius=0.1,
                     shaft_resolution=50, scale=2.5)
    glyphs = cloud.glyph(orient="omega_vectors", scale=scale_arg,
                         factor=cfg.arrow_scale, geom=arrow, tolerance=0.0)

    plotter.add_mesh(
        glyphs, color=[0, 192, 255], show_scalar_bar=False, name="omega_vectors",
        specular=0.8, specular_power=50, smooth_shading=True, lighting=True,
        opacity=1, line_width=3.0, render_points_as_spheres=True,
        ambient=0.4, diffuse=0.8, pickable=False,
        silhouette={"color": "black", "line_width": 2.8, "opacity": 1,
                    "feature_angle": 45.0, "decimate": 0.5},
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="3D Fermi surface with optional Berry-curvature arrows.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Defaults reproduce the current output of Fermi_Surface.py.",
    )
    g = p.add_argument_group("data")
    g.add_argument("--preset", choices=sorted(PRESETS), default="arrows",
                   help="visual preset (default: arrows)")
    g.add_argument("--data-file", help="input .dat file")
    g.add_argument("--energy-col", type=int,
                   help="energy column (default: 7 for 9-col, 3 for 5-col)")
    g.add_argument("--ref-col", type=int,
                   help="column for the band minimum (default: same as energy)")
    g.add_argument("--energy-offset", type=float,
                   help="Fermi level above the band minimum, in eV")
    g.add_argument("--tolerance", type=float, default=None,
                   help="half-width of the energy window; default None (no "
                        "filter). The original's 0.5 kept 100%% of points.")

    g = p.add_argument_group("grid")
    g.add_argument("--grid", choices=["legacy", "native"], default="legacy",
                   help="legacy: resample over fixed bounds (default, matches "
                        "published figures); native: use the data's own axes "
                        "(correct -- no extrapolation -- but CHANGES the figure)")
    g.add_argument("--dimension", type=int, help="samples per axis in legacy mode")
    g.add_argument("--bounds", type=float, nargs=6,
                   metavar=("X0", "X1", "Y0", "Y1", "Z0", "Z1"))
    g.add_argument("--interp", choices=["reshape", "griddata"], default="reshape",
                   help="reshape (default, exact and ~80x faster) or the "
                        "original scipy griddata call, kept as an oracle")

    g = p.add_argument_group("arrows")
    g.add_argument("--arrows", dest="arrows", action="store_true", default=None)
    g.add_argument("--no-arrows", dest="arrows", action="store_false")
    g.add_argument("--num-bins", type=int)
    g.add_argument("--max-distance", type=float, default=0.012)
    g.add_argument("--arrow-scale", type=float, default=0.0035)
    g.add_argument("--sampler", choices=["legacy", "fixed"], default="legacy",
                   help="fixed: keep axis-maximum points and pick the point "
                        "nearest each cell centre (CHANGES the figure)")
    g.add_argument("--scale-by-magnitude", action="store_true",
                   help="scale arrow length by |Omega| (CHANGES the figure)")
    g.add_argument("--kz-ranges", nargs="+", metavar="Z0:Z1",
                   help="restrict arrows to these kz ranges")

    g = p.add_argument_group("render")
    g.add_argument("--smooth-iter", type=int, default=10000)
    g.add_argument("--relaxation", type=float, default=0.01)
    g.add_argument("--window", type=int, nargs=2, metavar=("W", "H"))
    g.add_argument("--azimuth", type=float)
    g.add_argument("--screenshot", metavar="PATH")
    g.add_argument("--off-screen", action="store_true",
                   help="render without a window (needs OSMesa/EGL or xvfb)")
    g.add_argument("--report", action="store_true",
                   help="print the pyvista environment report")
    g.add_argument("--dump-stats", metavar="PATH",
                   help="write numeric fingerprints as JSON and exit before "
                        "rendering (no display needed)")
    return p


def resolve_config(args):
    cfg = argparse.Namespace(**vars(args))
    preset = PRESETS[args.preset]
    for key, value in preset.items():
        if getattr(cfg, key, None) is None:
            setattr(cfg, key, value)

    if args.bounds is not None:
        cfg.bounds = tuple(args.bounds)
        cfg.data_driven_bounds = False
    if args.window is not None:
        cfg.window = tuple(args.window)

    ranges = []
    for item in (args.kz_ranges or []):
        try:
            lo, hi = item.split(":")
            ranges.append((float(lo), float(hi)))
        except ValueError:
            raise SystemExit(f"error: --kz-ranges expects Z0:Z1, got {item!r}")
    cfg.kz_ranges = ranges or None
    return cfg


def main(argv=None):
    args = build_parser().parse_args(argv)
    cfg = resolve_config(args)

    changed = [name for name, on in (
        ("--grid native", cfg.grid == "native"),
        ("--sampler fixed", cfg.sampler == "fixed"),
        ("--scale-by-magnitude", cfg.scale_by_magnitude),
        ("--tolerance", cfg.tolerance is not None),
        ("--energy-col/--ref-col", args.energy_col is not None or args.ref_col is not None),
    ) if on]
    if changed:
        print("[non-default] " + ", ".join(changed)
              + " -- output differs from the published figures", file=sys.stderr)

    data = load_kspace_data(cfg.data_file, cfg.energy_col, cfg.ref_col)

    band_min = float(np.min(data.ref_energy))
    fermi_energy = band_min + cfg.energy_offset
    print(f"Unperturbed band minimum: {band_min:.6f} eV")
    print(f"Fermi energy set to: {fermi_energy:.6f} eV")

    if cfg.tolerance is None:
        filt = data
    else:
        m = ((data.energy >= fermi_energy - cfg.tolerance)
             & (data.energy <= fermi_energy + cfg.tolerance))
        filt = KData(data.kx[m], data.ky[m], data.kz[m], data.energy[m],
                     data.ref_energy[m],
                     None if data.omega is None else data.omega[m],
                     None if data.omega_mag is None else data.omega_mag[m])

    if cfg.data_driven_bounds and cfg.grid != "native":
        cfg.bounds = (float(filt.kx.min()), float(filt.kx.max()),
                      float(filt.ky.min()), float(filt.ky.max()),
                      float(filt.kz.min()), float(filt.kz.max()))

    print("Performing interpolation...")
    X, Y, Z, values = build_energy_grid(data, cfg)

    if cfg.dump_stats:
        import hashlib
        import json

        def sha(a):
            return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]

        stats = {
            "data_file": cfg.data_file,
            "n_rows": len(data),
            "band_min": band_min,
            "fermi_energy": fermi_energy,
            "grid": cfg.grid,
            "interp": cfg.interp,
            "values_shape": list(values.shape),
            "values_sha": sha(values),
            "values_sum": float(values.sum()),
        }
        if data.omega is not None:
            s = grid_sample_3d(filt.kx, filt.ky, filt.kz,
                               filt.omega[:, 0], filt.omega[:, 1], filt.omega[:, 2],
                               cfg.num_bins, mode=cfg.sampler,
                               kz_ranges=cfg.kz_ranges)
            stats["sampler_n"] = int(len(s[0]))
            stats["sampler_xyz_sha"] = sha(np.column_stack(s[:3]))
            stats["sampler_vec_sha"] = sha(np.column_stack(s[3:]))
        with open(cfg.dump_stats, "w") as fh:
            json.dump(stats, fh, indent=2)
        print(json.dumps(stats, indent=2))
        return 0

    render(data, X, Y, Z, values, fermi_energy, filt, cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
