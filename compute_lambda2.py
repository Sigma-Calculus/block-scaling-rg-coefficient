#!/usr/bin/env python3
import json, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

# ------------------------------------------------------------------ #
print("• Dateien laden …")
edges = np.loadtxt("fine_graph.edgelist", dtype=np.int64)

with open("block_stars.json") as f:
    stars_raw = json.load(f)
stars_raw = {int(k): v for k, v in stars_raw.items()}

E_f = edges.shape[0]
star_ids = sorted(stars_raw)            # lückenlos dichte Indizes
id2row = {sid: i for i, sid in enumerate(star_ids)}
I = len(star_ids)

print(f"  fine edges : {E_f:,}")
print(f"  block stars: {I:,} (dichte Rows)")

# ------------------------------------------------------------------ #
print("• M-Matrix bauen …")
row, col, data = [], [], []
for old_id, fine_list in stars_raw.items():
    r = id2row[old_id]
    w = 1.0 / len(fine_list)
    row.extend([r] * len(fine_list))
    col.extend(fine_list)
    data.extend([w] * len(fine_list))

M = sp.coo_matrix((data, (row, col)), shape=(I, E_f), dtype=np.float64).tocsr()

# ------------------------------------------------------------------ #
print("• Laplacian L_sym = I – D^{-1/2} S D^{-1/2} …")
S = M @ M.T
d = np.array(S.sum(axis=1)).ravel()
Dinvsqrt = sp.diags(1.0 / np.sqrt(d))
L = sp.eye(I, format="csr") - Dinvsqrt @ S @ Dinvsqrt
del M, S, d

# ------------------------------------------------------------------ #
print("• zweitkleinster Eigenwert λ₂ …")
lam, _ = eigsh(L, k=2, which="SM", tol=1e-3, maxiter=300)
theta_new = lam[1]

print(f"\nϑ_new  =  {theta_new:.6f}")
print("fertig.")
