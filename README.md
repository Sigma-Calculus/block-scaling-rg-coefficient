# block-scaling-rg-coefficient
Scripts and numerical tests for validating the universal RG coefficient  C∼1/L^2   in the block Laplacian framework.

# Block Scaling and the Universal RG Coefficient \( C \)

This repository contains Python scripts and data to validate the scaling law  
λ₂(L) ~ C / L²
for the spectral gap of the block-averaged Laplacian in a coarse-grained, periodic 4D lattice.

The coefficient \( C \) is analytically fixed by the topology and physical units of the model.  
There is no tunable parameter — only numerical confirmation.

---

## 📁 Contents

- `build_sigma_net_fast.py`  
  Constructs a 4D periodic lattice and the associated block-star mapping for coarse-graining.  
  Produces:
  - `fine_graph.edgelist`: full edge list of the fine lattice  
  - `block_stars.json`: coarse block mapping for Laplacian construction

- `compute_lambda2.py`  
  Computes the second eigenvalue \( \lambda_2 \) of the block-averaged normalized Laplacian using sparse eigensolvers.

- `plot_scaling.py`  
  Generates log–log plots comparing \( \lambda_2 \) to the theoretical form \( C/L^2 \).

---

## ⚙️ Requirements

- Python 3.8+
- NumPy
- SciPy
- Matplotlib

Install dependencies (e.g. via pip):

```bash
pip install numpy scipy matplotlib

# Step 1: Generate fine lattice and block mapping
python build_sigma_net_fast.py --N0 32 --N1 64 --N2 64 --N3 64 \
    --out_edges fine_graph.edgelist --out_stars block_stars.json

# Step 2: Compute spectral gap of coarse Laplacian
python compute_lambda2.py

# Step 3: Plot numerical λ₂(L) against analytic C / L²
python plot_scaling.py

### 📄 Citation

If you use this code, please cite the Zenodo archive:  
https://zenodo.org/records/15694883

