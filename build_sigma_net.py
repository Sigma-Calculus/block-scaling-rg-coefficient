#!/usr/bin/env python3
"""
build_sigma_net_fast.py  –  Periodisches T^4-Gitter + 2×2×2-Blocksterne
skalieren linear in der Zahl der Feinkanten.

Beispiel:
    python build_sigma_net_fast.py --N0 32 --N1 64 --N2 64 --N3 64 \
           --out_edges fine_graph.edgelist --out_stars block_stars.json
"""
import argparse, json, sys
from collections import defaultdict
from itertools import product

# ----------------------------------------------------------------------
def parse():
    p = argparse.ArgumentParser()
    p.add_argument("--N0", type=int, default=4)
    p.add_argument("--N1", type=int, default=4)
    p.add_argument("--N2", type=int, default=4)
    p.add_argument("--N3", type=int, default=4)
    p.add_argument("--out_edges", default="fine_graph.edgelist",
                   help="Datei: src dst mu")
    p.add_argument("--out_stars", default="block_stars.json",
                   help="JSON: coarse_id → [fine_edge_ids]")
    return p.parse_args()

# ----------------------------------------------------------------------
def vid(coord, N):
    """4-Tupel → eindeutige Ganzzahl-ID."""
    n0, n1, n2, n3 = coord
    return (((n0*N[1] + n1)*N[2] + n2)*N[3] + n3)

def decode(v, N):
    """Ganzzahl-ID → 4-Tupel."""
    n3 = v % N[3]; v //= N[3]
    n2 = v % N[2]; v //= N[2]
    n1 = v % N[1]; v //= N[1]
    n0 = v
    return (n0, n1, n2, n3)

directions = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]

# ----------------------------------------------------------------------
def even_candidates(c, Ndim):
    """
    Liefert die even-Koordinaten im Abstand ≤ 1 (periodisch).
    Maximal 2 Werte pro Dimension.
    """
    base = (c // 2) * 2
    cand = [base]
    for delta in (-2, 2):
        e = (base + delta) % Ndim
        if (c - e) % Ndim <= 1 or (e - c) % Ndim <= 1:
            cand.append(e)
    return cand

# ----------------------------------------------------------------------
def main():
    args = parse()
    N = [args.N0, args.N1, args.N2, args.N3]
    print("[σ-net] lattice size", N, file=sys.stderr)

    # ---------- 1  Feinkanten (O(|E_f|)) ------------------------------
    edge_set = set()                           # vermeidet Duplikate
    for coord in product(range(N[0]), range(N[1]),
                         range(N[2]), range(N[3])):
        src = vid(coord, N)
        for mu, d in enumerate(directions):
            dst = vid(tuple((coord[i]+d[i]) % N[i] for i in range(4)), N)
            edge_set.add((*sorted((src, dst)), mu))   # undirektional

    edges = sorted(edge_set)                  # (src, dst, mu)
    E_f   = len(edges)
    print("[σ-net] fine edges :", E_f, file=sys.stderr)

    # Mapping (src,dst) → index (beide Orientierungen)
    fine_id = {}
    for idx, (u, v, _) in enumerate(edges):
        fine_id[(u, v)] = fine_id[(v, u)] = idx

    # ---------- 2  coarse edges ---------------------------------------
    coarse_edges, key2cid = [], {}
    for (u, v, mu) in edges:
        if all(c % 2 == 0 for c in decode(u, N)):          # Start even×4
            cid = len(coarse_edges)
            coarse_edges.append((u, v, mu))
            key2cid[(decode(u, N), mu)] = cid
    print("[σ-net] coarse edges:", len(coarse_edges), file=sys.stderr)

    # ---------- 3  Block-Sterne (linear) ------------------------------
    stars = defaultdict(list)                 # cid → [fine-IDs]
    for (u, v, mu_edge) in edges:
        cu, cv = decode(u, N), decode(v, N)
        cand_lists = [even_candidates(min(cu[d], cv[d]), N[d]) for d in range(4)]

        # höchstens 2^4 = 16 Anker pro Fine-Kante
        for anchor in product(*cand_lists):
            for mu_star in range(4):
                cid = key2cid.get((anchor, mu_star))
                if cid is not None:
                    stars[cid].append(fine_id[(u, v)])

    # ---------- 4  Dateien schreiben ---------------------------------
    with open(args.out_edges, "w") as f:
        for u, v, mu in edges:
            f.write(f"{u} {v} {mu}\n")

    with open(args.out_stars, "w") as f:
        json.dump({str(k): v for k, v in stars.items()}, f)

    print("[σ-net] written", args.out_edges, args.out_stars, file=sys.stderr)

# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
