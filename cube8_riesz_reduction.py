#!/usr/bin/env python3
"""Worked 8-vertex reduction example for planar Riesz s=2 subset selection.

This script is a computational companion to the manuscript
"Riesz s-Energy Subset Selection is NP Hard in the Euclidean Plane".

It uses the planar cubic cube graph Q_3 (8 vertices, 12 edges).  A cubic
(simple) graph cannot have 7 vertices because the degree sum 3n must be even.

What the script checks:
  1. Enumerates all 2^8 Ising states and verifies alpha(Q_3)=4 and min H_B=-12.
  2. Builds an explicit orthogonal selector-tree layout with one objective box
     per graph edge.
  3. Computes the exact two-state Riesz interaction of every selector pair
     using high-precision arithmetic for the inverse-square kernel.
  4. Applies the manuscript's first-order local field compensation, rounded to
     finite decimals so all emitted candidate coordinates are rational.
  5. Verifies the two-of-four selector-forcing inequality for the finite demo.
  6. Audits every consistency-tree edge by a worst-case component-flip bound.
  7. Collapses each tree to one source spin, enumerates all 256 normalized
     states, and checks that a threshold between the H_B=-12 and H_B=-8
     levels separates them in the actual finite Riesz geometry.

The constants here are deliberately human-scale demonstration constants, not
those used by the asymptotic theorem.  The general proof uses much larger
feature separation and a much smaller b.  This script is an executable finite
sanity check of the same architecture, not a substitute for the proof.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from pathlib import Path

import mpmath as mp
import numpy as np

try:
    import matplotlib.pyplot as plt
except Exception:  # plotting is optional
    plt = None


# ---------------------------------------------------------------------------
# Source graph: planar cubic cube graph Q_3
# ---------------------------------------------------------------------------

VERTICES = tuple(range(8))
EDGES = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
)
EDGE_SET = {tuple(sorted(e)) for e in EDGES}

# A planar orthogonal drawing.  Outer vertices are 0..3; inner vertices 4..7.
BASE_ROOT = {
    0: (-10, -10), 1: (10, -10), 2: (10, 10), 3: (-10, 10),
    4: (-4, -4), 5: (4, -4), 6: (4, 4), 7: (-4, 4),
}


def independent_set(mask: int) -> bool:
    return all(not (((mask >> u) & 1) and ((mask >> v) & 1)) for u, v in EDGES)


def selected_count(mask: int) -> int:
    return mask.bit_count()


def spin(mask: int, v: int) -> int:
    # +1 = selected, -1 = not selected
    return 1 if ((mask >> v) & 1) else -1


def barahona_energy(mask: int) -> int:
    sig = [spin(mask, v) for v in VERTICES]
    return sum(sig[u] * sig[v] for u, v in EDGES) + sum(sig)


def induced_edge_count(mask: int) -> int:
    return sum(1 for u, v in EDGES if ((mask >> u) & 1) and ((mask >> v) & 1))


# ---------------------------------------------------------------------------
# Orthogonal routing and selector trees
# ---------------------------------------------------------------------------


def scale_point(p: tuple[int, int], scale: int) -> tuple[int, int]:
    return (p[0] * scale, p[1] * scale)


def build_roots(scale: int) -> dict[int, tuple[int, int]]:
    return {v: scale_point(p, scale) for v, p in BASE_ROOT.items()}


def build_routes(scale: int, root: dict[int, tuple[int, int]]):
    s = lambda p: scale_point(p, scale)
    return {
        (0, 1): [root[0], s((-10, -14)), s((10, -14)), root[1]],
        (1, 2): [root[1], s((14, -10)), s((14, 10)), root[2]],
        (2, 3): [root[2], s((10, 14)), s((-10, 14)), root[3]],
        (3, 0): [root[3], s((-14, 10)), s((-14, -10)), root[0]],
        (4, 5): [root[4], root[5]],
        (5, 6): [root[5], root[6]],
        (6, 7): [root[6], root[7]],
        (7, 4): [root[7], root[4]],
        (0, 4): [root[0], s((-4, -10)), root[4]],
        (1, 5): [root[1], s((4, -10)), root[5]],
        (2, 6): [root[2], s((4, 10)), root[6]],
        (3, 7): [root[3], s((-4, 10)), root[7]],
    }


def unit_lattice_points(polyline: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out = [polyline[0]]
    for a, z in zip(polyline, polyline[1:]):
        dx, dy = z[0] - a[0], z[1] - a[1]
        assert dx == 0 or dy == 0
        length = abs(dx) + abs(dy)
        d = (0 if dx == 0 else (1 if dx > 0 else -1),
             0 if dy == 0 else (1 if dy > 0 else -1))
        out.extend((a[0] + d[0] * k, a[1] + d[1] * k)
                   for k in range(1, length + 1))
    return out


def split_route_with_objective(polyline: list[tuple[int, int]], D: int):
    """Split one routed graph edge into two trees and a 45-degree objective.

    On the longest straight segment, the u-side ends at p_0.  The v-side ends
    at q_0=p_0+d+n, where d is the corridor direction and n is a unit normal.
    The v-side then contains q_j=q_0+j d for j=0..D before returning by one
    unit to the original corridor.  The u-side contains p_i=p_0-i d.
    """
    lengths = [abs(z[0] - a[0]) + abs(z[1] - a[1])
               for a, z in zip(polyline, polyline[1:])]
    seg = max(range(len(lengths)), key=lengths.__getitem__)
    a, z = polyline[seg], polyline[seg + 1]
    length = lengths[seg]
    d = ((z[0] - a[0]) // length, (z[1] - a[1]) // length)
    n = (-d[1], d[0])

    # Need D points behind p_0 and D points beyond q_0, plus rejoin room.
    t = (length - (D + 1)) // 2
    if not (t >= D + 1 and t + D + 1 <= length - (D + 1)):
        raise ValueError(f"segment length {length} too short for D={D}")

    p0 = (a[0] + d[0] * t, a[1] + d[1] * t)
    rejoin = (p0[0] + d[0] * (D + 1), p0[1] + d[1] * (D + 1))
    q0 = (p0[0] + d[0] + n[0], p0[1] + d[1] + n[1])

    route_points = unit_lattice_points(polyline)
    p_idx = route_points.index(p0)
    r_idx = route_points.index(rejoin)

    u_branch = route_points[:p_idx + 1]
    v_branch = list(reversed(route_points[r_idx:]))
    offset_end = (rejoin[0] + n[0], rejoin[1] + n[1])
    v_branch.append(offset_end)
    for j in range(D - 1, -1, -1):
        v_branch.append((q0[0] + j * d[0], q0[1] + j * d[1]))

    p_prefix = [(p0[0] - i * d[0], p0[1] - i * d[1]) for i in range(D + 1)]
    q_prefix = [(q0[0] + j * d[0], q0[1] + j * d[1]) for j in range(D + 1)]
    return u_branch, v_branch, p_prefix, q_prefix


def build_selector_layout(scale: int = 2, D: int = 2):
    roots = build_roots(scale)
    routes = build_routes(scale, roots)
    cells: dict[int, set[tuple[int, int]]] = defaultdict(set)
    tree_edges: dict[int, set[tuple[tuple[int, int], tuple[int, int]]]] = defaultdict(set)
    objectives = []

    for (u, v), poly in routes.items():
        u_branch, v_branch, p_prefix, q_prefix = split_route_with_objective(poly, D)
        for owner, seq in ((u, u_branch), (v, v_branch)):
            cells[owner].update(seq)
            for a, z in zip(seq, seq[1:]):
                tree_edges[owner].add(tuple(sorted((a, z))))
        objectives.append({"edge": (u, v), "p_prefix": p_prefix, "q_prefix": q_prefix})

    # Each source-vertex structure must be a tree: |E|=|V|-1 and connected.
    for v in VERTICES:
        assert len(tree_edges[v]) == len(cells[v]) - 1

    centers, owners = [], []
    for v in VERTICES:
        for c in sorted(cells[v]):
            owners.append(v)
            centers.append(c)

    # No two different source trees may share a selector center.
    seen = {}
    for i, (v, c) in enumerate(zip(owners, centers)):
        if c in seen and seen[c] != v:
            raise AssertionError(f"selector-center collision at {c}: {seen[c]} vs {v}")
        seen[c] = v

    return roots, routes, cells, tree_edges, objectives, centers, owners


# ---------------------------------------------------------------------------
# Riesz s=2 selector interactions
# ---------------------------------------------------------------------------


def state_points(center, sigma: int, b: mp.mpf, eta: mp.mpf = mp.mpf("0")):
    x, y = mp.mpf(center[0]), mp.mpf(center[1])
    if sigma == 1:
        a = b + eta
        return ((x + a, y + a), (x - a, y - a))
    return ((x + b, y - b), (x - b, y + b))


def pair_fourier(c1, c2, b: mp.mpf,
                 eta1: mp.mpf = mp.mpf("0"), eta2: mp.mpf = mp.mpf("0")):
    vals = {}
    for s in (1, -1):
        p1 = state_points(c1, s, b, eta1)
        for t in (1, -1):
            p2 = state_points(c2, t, b, eta2)
            z = mp.mpf("0")
            for a in p1:
                for q in p2:
                    dx, dy = a[0] - q[0], a[1] - q[1]
                    z += 1 / (dx * dx + dy * dy)
            vals[(s, t)] = z
    pp, pm = vals[(1, 1)], vals[(1, -1)]
    mpv, mm = vals[(-1, 1)], vals[(-1, -1)]
    A = (pp + pm + mpv + mm) / 4
    h1 = (pp + pm - mpv - mm) / 4
    h2 = (pp - pm + mpv - mm) / 4
    J = (pp - pm - mpv + mm) / 4
    return A, h1, h2, J


def objective_lambda(D: int, b: mp.mpf) -> mp.mpf:
    total = mp.mpf("0")
    for i in range(D + 1):
        for j in range(D + 1):
            total += pair_fourier((-i, 0), (1 + j, 1), b)[3]
    return total


def round_mpf_to_decimal_rational(x: mp.mpf, digits: int) -> tuple[Decimal, mp.mpf]:
    """Round to a finite decimal, hence an exact rational number."""
    getcontext().prec = max(80, digits + 40)
    d = Decimal(mp.nstr(x, 90))
    quantum = Decimal(1).scaleb(-digits)
    q = d.quantize(quantum, rounding=ROUND_HALF_EVEN)
    return q, mp.mpf(str(q))


def decimal_string(x: mp.mpf, digits: int = 30) -> str:
    q, _ = round_mpf_to_decimal_rational(x, digits)
    return format(q, "f")


# ---------------------------------------------------------------------------
# Numerical audit of the complete finite instance
# ---------------------------------------------------------------------------


def run_audit(outdir: Path, scale: int, D: int, b_den: int,
              dps: int, eta_digits: int, threshold_digits: int):
    mp.mp.dps = dps
    b = mp.mpf(1) / b_den

    roots, routes, cells, tree_edges, objectives, centers, owners = build_selector_layout(scale, D)
    M = len(centers)
    k = 2 * M
    index = {(owners[i], centers[i]): i for i in range(M)}

    # Source enumeration.
    source_rows = []
    alpha = 0
    min_H = None
    for mask in range(1 << 8):
        t = selected_count(mask)
        q = induced_edge_count(mask)
        H = barahona_energy(mask)
        indep = independent_set(mask)
        if indep:
            alpha = max(alpha, t)
        min_H = H if min_H is None else min(min_H, H)
        source_rows.append((mask, t, q, indep, H))

    assert alpha == 4
    assert min_H == -12
    assert min(H for _, _, _, _, H in source_rows if H > min_H) == -8

    # Standard local objective coefficient.
    lam = objective_lambda(D, b)

    # Unperturbed geometric one-spin fields G_i.
    G = [mp.mpf("0") for _ in range(M)]
    for i in range(M):
        for j in range(i + 1, M):
            _, h1, h2, _ = pair_fourier(centers[i], centers[j], b)
            G[i] += h1
            G[j] += h2

    beta = -1 / (8 * b**3)
    targets = [lam if centers[i] == roots[owners[i]] else mp.mpf("0") for i in range(M)]
    eta_star = [(targets[i] - G[i]) / beta for i in range(M)]

    # Emit finite-decimal eta values.  These are exact rationals and are used
    # for all subsequent checks, not merely for output formatting.
    eta_dec, eta = zip(*(round_mpf_to_decimal_rational(x, eta_digits) for x in eta_star))
    eta = list(eta)

    # Recompute the exact two-state effective Hamiltonian numerically at high
    # precision for the rounded rational coordinates.
    h_cell = [mp.mpf("0") for _ in range(M)]
    J_float = np.zeros((M, M), dtype=float)
    A_total = mp.mpf("0")

    for i in range(M):
        eplus = 1 / (8 * (b + eta[i])**2)
        eminus = 1 / (8 * b**2)
        A_total += (eplus + eminus) / 2
        h_cell[i] += (eplus - eminus) / 2

    # Source-level fields/couplings and baseline after tree normalization.
    source_h = [mp.mpf("0") for _ in VERTICES]
    source_K = [[mp.mpf("0") for _ in VERTICES] for _ in VERTICES]
    E0 = A_total

    # Pair contributions.
    for i in range(M):
        for j in range(i + 1, M):
            A, h1, h2, J = pair_fourier(centers[i], centers[j], b, eta[i], eta[j])
            h_cell[i] += h1
            h_cell[j] += h2
            J_float[i, j] = J_float[j, i] = float(J)
            E0 += A
            u, v = owners[i], owners[j]
            if u == v:
                # In a normalized tree sigma_i sigma_j = 1, so this is baseline.
                E0 += J
            else:
                if u > v:
                    u, v = v, u
                source_K[u][v] += J

    for i, v in enumerate(owners):
        source_h[v] += h_cell[i]

    # Build tree adjacency in cell-index coordinates and verify connectivity.
    adjacency: dict[int, list[int]] = defaultdict(list)
    tree_index_edges = []
    for v in VERTICES:
        for a, z in tree_edges[v]:
            ia, iz = index[(v, a)], index[(v, z)]
            adjacency[ia].append(iz)
            adjacency[iz].append(ia)
            tree_index_edges.append((v, ia, iz))

        root_idx = index[(v, roots[v])]
        seen = {root_idx}
        stack = [root_idx]
        while stack:
            x = stack.pop()
            for y in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        assert len(seen) == len(cells[v])

    # Worst-case component-flip audit for every consistency-tree edge.
    absJ = np.abs(J_float)
    abs_h = np.abs(np.array([float(x) for x in h_cell]))
    min_flip_margin = math.inf
    worst_flip = None

    def component_after_cut(start: int, a: int, z: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            x = stack.pop()
            for y in adjacency[x]:
                if (x == a and y == z) or (x == z and y == a):
                    continue
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen

    for v, a, z in tree_index_edges:
        root_idx = index[(v, roots[v])]
        ca = component_after_cut(a, a, z)
        comp = ca if root_idx not in ca else component_after_cut(z, a, z)
        mask = np.zeros(M, dtype=bool)
        mask[list(comp)] = True
        favorable = abs(J_float[a, z])
        other_cross = absJ[np.ix_(mask, ~mask)].sum() - favorable
        field_load = abs_h[mask].sum()
        margin = favorable - other_cross - field_load
        if margin < min_flip_margin:
            min_flip_margin = margin
            worst_flip = {
                "owner": v,
                "edge": [centers[a], centers[z]],
                "component_size": len(comp),
                "favorable": favorable,
                "other_cross": other_cross,
                "field_load": field_load,
            }

    # Source-level normalized states.
    normalized_rows = []
    max_remainder = mp.mpf("0")
    min_yes = None
    min_no = None
    max_yes = None
    for mask in range(1 << 8):
        sig = [spin(mask, v) for v in VERTICES]
        H = barahona_energy(mask)
        E = E0
        E += sum(source_h[v] * sig[v] for v in VERTICES)
        E += sum(source_K[u][v] * sig[u] * sig[v]
                 for u in VERTICES for v in VERTICES if u < v)
        ideal = E0 + lam * H
        remainder = E - ideal
        max_remainder = max(max_remainder, abs(remainder))
        scaled = (E - E0) / lam
        normalized_rows.append((mask, H, E, scaled, remainder))
        if H == -12:
            min_yes = E if min_yes is None else min(min_yes, E)
            max_yes = E if max_yes is None else max(max_yes, E)
        elif H >= -8:
            min_no = E if min_no is None else min(min_no, E)

    # Midpoint source threshold for K=4: H0+2 = -10.
    tau_ideal = E0 - 10 * lam
    tau_dec, tau = round_mpf_to_decimal_rational(tau_ideal, threshold_digits)
    assert max_yes < tau < min_no

    # Finite-demo selector forcing check from Lemma 3.1.
    selector_lhs = (mp.mpf(1) / 16) * b**(-2)
    selector_rhs = 4 * M
    max_eta_ratio = max(abs(e) for e in eta) / b

    edge_ratios = [source_K[min(u, v)][max(u, v)] / lam for u, v in EDGES]
    nonedge_ratios = [abs(source_K[u][v] / lam)
                      for u in VERTICES for v in VERTICES if u < v and (u, v) not in EDGE_SET]

    report = {
        "source_graph": "cube graph Q3",
        "n": 8,
        "m": 12,
        "alpha": alpha,
        "min_barahona_energy": min_H,
        "next_barahona_level": -8,
        "selector_cells_M": M,
        "candidate_points": 4 * M,
        "subset_size_k": k,
        "scale": scale,
        "objective_prefix_D": D,
        "b": f"1/{b_den}",
        "lambda": mp.nstr(lam, 30),
        "lambda_over_b4": mp.nstr(lam / b**4, 20),
        "selector_force_lhs": mp.nstr(selector_lhs, 20),
        "selector_force_rhs": selector_rhs,
        "selector_force_holds": bool(selector_lhs > selector_rhs),
        "max_eta_over_b": mp.nstr(max_eta_ratio, 20),
        "field_ratio_min": mp.nstr(min(h / lam for h in source_h), 20),
        "field_ratio_max": mp.nstr(max(h / lam for h in source_h), 20),
        "edge_coupling_ratio_min": mp.nstr(min(edge_ratios), 20),
        "edge_coupling_ratio_max": mp.nstr(max(edge_ratios), 20),
        "max_nonedge_coupling_over_lambda": mp.nstr(max(nonedge_ratios), 20),
        "max_remainder_over_lambda": mp.nstr(max_remainder / lam, 20),
        "min_component_flip_margin_over_b4": f"{min_flip_margin / float(b**4):.12g}",
        "threshold_tau_decimal_rational": format(tau_dec, "f"),
        "max_yes_scaled_energy": mp.nstr((max_yes - E0) / lam, 20),
        "min_no_scaled_energy": mp.nstr((min_no - E0) / lam, 20),
        "threshold_scaled_energy": "-10",
        "worst_flip": worst_flip,
        "note": "Finite high-precision audit of a pedagogical-scale instance; the theorem uses more conservative asymptotic scales."
    }

    # ------------------------------------------------------------------
    # Emit files
    # ------------------------------------------------------------------
    outdir.mkdir(parents=True, exist_ok=True)

    with (outdir / "cube8_reduction_report.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    with (outdir / "cube8_reduction_report.txt").open("w", encoding="utf-8") as f:
        f.write("Worked cube-graph reduction to planar Riesz s=2 subset selection\n")
        f.write("=" * 70 + "\n\n")
        f.write("Source: cube graph Q3, 8 vertices, 12 edges, degree 3 at every vertex.\n")
        f.write("Maximum independent-set size alpha = 4.\n")
        f.write("Barahona ground energy = -12; next source level = -8.\n\n")
        f.write(f"Selector cells M = {M}; candidate points = {4*M}; k = {k}.\n")
        f.write(f"Demo b = 1/{b_den}; objective prefix D = {D}; orthogonal scale = {scale}.\n")
        f.write(f"lambda / b^4 = {mp.nstr(lam / b**4, 14)}.\n")
        f.write(f"Selector forcing: (1/16)b^-2 = {mp.nstr(selector_lhs, 12)} > 4M = {selector_rhs}.\n")
        f.write(f"max |eta_i|/b = {mp.nstr(max_eta_ratio, 12)} (< 1e-6).\n")
        f.write(f"Source field coefficients / lambda lie in [{mp.nstr(min(h/lam for h in source_h), 10)}, {mp.nstr(max(h/lam for h in source_h), 10)}].\n")
        f.write(f"Graph-edge couplings / lambda lie in [{mp.nstr(min(edge_ratios), 10)}, {mp.nstr(max(edge_ratios), 10)}].\n")
        f.write(f"Largest nonedge coupling / lambda = {mp.nstr(max(nonedge_ratios), 10)}.\n")
        f.write(f"Worst normalized remainder |R|/lambda = {mp.nstr(max_remainder/lam, 10)}.\n")
        f.write(f"Minimum component-flip margin = {min_flip_margin/float(b**4):.8f} b^4.\n")
        f.write(f"Ground-state Riesz scaled energies are at most {mp.nstr((max_yes-E0)/lam, 10)}.\n")
        f.write(f"All next-level states have scaled energy at least {mp.nstr((min_no-E0)/lam, 10)}.\n")
        f.write("The midpoint threshold is scaled energy -10.\n")
        f.write(f"A finite-decimal rational threshold is\n  tau = {format(tau_dec, 'f')}\n\n")
        f.write("This is a finite numerical audit of the proof architecture.  It is not\n")
        f.write("used as a proof of the general theorem.\n")

    # Candidate point set.  The printed coordinates are finite decimals, hence rational.
    b_dec = Decimal(1) / Decimal(b_den)
    with (outdir / "cube8_riesz_points.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["point_id", "cell_id", "source_vertex", "center_x", "center_y", "role", "x", "y"])
        pid = 0
        for i, ((cx, cy), v) in enumerate(zip(centers, owners)):
            e = eta_dec[i]
            a = b_dec + e
            points = [
                ("plus_a", Decimal(cx) + a, Decimal(cy) + a),
                ("plus_b", Decimal(cx) - a, Decimal(cy) - a),
                ("minus_a", Decimal(cx) + b_dec, Decimal(cy) - b_dec),
                ("minus_b", Decimal(cx) - b_dec, Decimal(cy) + b_dec),
            ]
            for role, x, y in points:
                w.writerow([pid, i, v, cx, cy, role, format(x, "f"), format(y, "f")])
                pid += 1

    with (outdir / "cube8_normalized_states.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["mask", "selected_vertices", "independent", "t", "q", "H_B", "scaled_riesz_energy", "remainder_over_lambda", "below_threshold"])
        norm_by_mask = {row[0]: row for row in normalized_rows}
        for mask, t, q, indep, H in source_rows:
            _, _, E, scaled, rem = norm_by_mask[mask]
            selected = " ".join(str(v) for v in VERTICES if (mask >> v) & 1)
            w.writerow([
                mask, selected, indep, t, q, H,
                mp.nstr(scaled, 25), mp.nstr(rem / lam, 25), bool(E < tau)
            ])

    metadata = {
        "s": 2,
        "k": k,
        "point_file": "cube8_riesz_points.csv",
        "threshold": format(tau_dec, "f"),
        "threshold_is_finite_decimal_rational": True,
        "source_target_K": 4,
        "source_graph_edges": [list(e) for e in EDGES],
        "root_centers": {str(v): list(roots[v]) for v in VERTICES},
        "selector_cells": M,
        "candidate_points": 4 * M,
        "eta_decimal_places": eta_digits,
        "arithmetic_dps": dps,
    }
    with (outdir / "cube8_instance_metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    # Plot selector-center trees and objective terminal pairs.
    if plt is not None:
        fig, ax = plt.subplots(figsize=(8, 8))
        for v in VERTICES:
            for a, z in tree_edges[v]:
                ax.plot([a[0], z[0]], [a[1], z[1]], linewidth=0.8, alpha=0.7)
        for obj in objectives:
            p0, q0 = obj["p_prefix"][0], obj["q_prefix"][0]
            ax.plot([p0[0], q0[0]], [p0[1], q0[1]], linestyle="--", linewidth=1.4)
        for v in VERTICES:
            ax.scatter([roots[v][0]], [roots[v][1]], s=55, zorder=5)
            ax.text(roots[v][0] + 0.6, roots[v][1] + 0.6, str(v), fontsize=10)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("Cube Q3 selector trees and 45-degree objective terminals")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(True, linewidth=0.25, alpha=0.35)
        fig.tight_layout()
        fig.savefig(outdir / "cube8_selector_layout.png", dpi=180)
        plt.close(fig)

    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("cube8_output"))
    ap.add_argument("--scale", type=int, default=2)
    ap.add_argument("--D", type=int, default=2)
    ap.add_argument("--b-den", type=int, default=500, help="use b=1/b_den")
    ap.add_argument("--dps", type=int, default=70)
    ap.add_argument("--eta-digits", type=int, default=30)
    ap.add_argument("--threshold-digits", type=int, default=30)
    args = ap.parse_args()

    report = run_audit(args.outdir, args.scale, args.D, args.b_den,
                       args.dps, args.eta_digits, args.threshold_digits)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
