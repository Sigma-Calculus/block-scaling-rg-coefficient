#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 22:30:32 2025

@author: oliversievers
"""

import matplotlib.pyplot as plt
import numpy as np

# Data
L = np.array([8, 16, 32, 64])
lam = np.array([0.351938, 0.110279, 0.029751, 0.007611])

plt.figure()
plt.loglog(L, lam, marker='o', linestyle='-', label=r'$\lambda_2(L)$ (numeric)')

# Fit line for slope -2
C = (lam * L**2).mean()
fit = C / L**2
plt.loglog(L, fit, linestyle='--', label=r'$C/L^2$ fit')

plt.xlabel(r'Box size $L$')
plt.ylabel(r'Spectral gap $\lambda_2$')
plt.title(r'Log–log plot of $\lambda_2$ vs $L$')
plt.legend()
plt.grid(True, which='both', ls=':')
plt.tight_layout()

