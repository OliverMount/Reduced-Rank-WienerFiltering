"""
Adaptive Noise Reduction - Python translation of MATLAB script
================================================================
Sections:
  1. Full-rank (Conventional) Wiener filter
  2. Principal Components / Eigen Decomposition rank reduction
  3. Cross-Spectral (CS) metric rank reduction
  4. Lanczos-based MSNWF (Multi-Stage Nested Wiener Filter)
  5. Filter order selection via Levinson recursion

All helper functions are direct translations of the supplied .m files.
Random seed is NOT fixed so every run gives different realisations,
exactly like MATLAB with no seed set.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz
from itertools import combinations


# =============================================================================
# Helper functions  (direct translations of the .m files)
# =============================================================================

from scipy.linalg import convolution_matrix

def convolmtx(x, M): 
    return convolution_matrix(x, M, mode='full')

 
def convol(x, h):
    """ 
    Convolution via convolution matrix; result is 1-D.
    """ 
    return convolmtx(x, len(h)) @ h  # performs matrix multiplications


def TimeAC(x, lag='half'):
    M = len(x)
    full = np.correlate(x, x, mode='full') / M
    if lag == 'full':
        return full                    # negative and positive lags: -(M-1)..0..+(M-1)
    else:
        return full[M-1:]             # positive lags including zero: 0..+(M-1)

def TimeCC(x, y, lag='half'):
    M = len(x)
    full = np.correlate(x, y, mode='full') / M
    if lag == 'half':
        return full[M-1:]             # positive lags including zero: 0..+(M-1)
    else:
        return full  

def AutoRnoise(a, noise, length):
    """
    Generates AR(p) process
    """  
    p     = len(a) - 1
    a1    = -np.array(a[1:])
    x     = np.zeros(length)
    state = np.zeros(p)

    for i in range(length):
        x[i]      = (state @ a1) + noise[i]
        state[1:] = state[:p-1]
        state[0]  = x[i]
    return x


def Levinson(r, b):
    """
    MATLAB Levinson(r, b)
    Levinson-Durbin: solve Toeplitz(r)*x = b.
    """
    r = np.asarray(r, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    p = len(b)
    a    = np.array([1.0])
    x    = np.array([b[0] / r[0]])
    epsi = np.zeros(p);  epsi[0] = r[0]
    for j in range(1, p):
        g       = np.dot(r[1:j+1], a[::-1])
        gamma   = -g / epsi[j-1]
        a_ext   = np.append(a, 0.0)
        a       = a_ext + gamma * np.conj(a_ext[::-1])
        epsi[j] = epsi[j-1] * (1.0 - abs(gamma)**2)
        delta   = np.dot(r[1:j+1], x[::-1])
        q       = (b[j] - delta) / epsi[j]
        x_ext   = np.append(x, 0.0)
        x       = x_ext + q * np.conj(a[::-1])
    return x


# =============================================================================
# Signal generation
# =============================================================================
k    = 200
nois = np.random.randn(k)
d    = AutoRnoise([1, -0.6, -0.3], nois, k)   # desired signal

noi  = np.random.randn(k)
v1   = AutoRnoise([1, -0.8], noi, k)           # primary-sensor noise
v2   = AutoRnoise([1,  0.6], noi, k)           # secondary-sensor noise
x    = d + v1

rv1  = TimeAC(v1)
rv2  = TimeAC(v2)
rxv2 = TimeCC(x, v2)

# =============================================================================
# 1. Conventional (full-rank) Wiener solution
# =============================================================================
p    = 12
Rv2  = toeplitz(rv2[:p])
w_opt = np.linalg.inv(Rv2) @ rxv2[:p]

conv_rv2         = convolmtx(v2, p)
estimated_noise  = convol(v2, w_opt)
recovered_signal = x - estimated_noise[:k]
MSE              = rv1[0] - w_opt @ Rv2 @ w_opt

plt.figure()
plt.plot(range(1, k+1), d,                'r', label='Original')
plt.plot(range(1, k+1), recovered_signal, 'b', label='Estimated')
plt.legend(); plt.title('Full-rank Wiener Filter')
plt.xlabel('Sample Index'); plt.ylabel('Amplitude')

# =============================================================================
# 2. Eigen Decomposition rank reduction
# =============================================================================
temp             = rxv2[:p].copy()
eig_vals, V_eig  = np.linalg.eigh(Rv2)   # ascending eigenvalues (same as MATLAB)

desired_signal_eigen = np.zeros((p, k))
MSE_eig              = np.zeros(p)

for idx in range(p):
    T     = V_eig[:, idx:]                 # columns idx..p-1
    I_sub = np.diag(eig_vals[idx:])
    w_e   = np.linalg.inv(T.T @ Rv2 @ T) @ (T.T @ temp)
    est   = (conv_rv2 @ T) @ w_e
    desired_signal_eigen[idx, :] = x - est[:k]
    MSE_eig[idx] = rv1[0] - temp @ T @ np.linalg.inv(I_sub) @ T.T @ temp

rank_idx = 5   # rank-6 (0-based)

# =============================================================================
# 3. Cross-Spectral (CS) metric rank reduction
# =============================================================================
o         = 6
combi     = list(combinations(range(p), o))
CS_metric = np.array([
    temp @ V_eig[:, list(c)]
         @ np.linalg.inv(np.diag(eig_vals[list(c)]))
         @ V_eig[:, list(c)].T @ temp
    for c in combi
])

best_cols = list(combi[np.argmax(CS_metric)])
T_cs      = V_eig[:, best_cols]
lamda_CS  = np.diag(eig_vals[best_cols])
w_cs      = np.linalg.inv(T_cs.T @ Rv2 @ T_cs) @ (T_cs.T @ temp)
est_cs    = (conv_rv2 @ T_cs) @ w_cs
desired_signal_CS = x - est_cs[:k]
MSE_CS    = rv1[0] - temp @ T_cs @ np.linalg.inv(lamda_CS) @ T_cs.T @ temp

plt.figure()
plt.plot(range(1,k+1), recovered_signal,              label='Full rank Wiener')
plt.plot(range(1,k+1), desired_signal_eigen[rank_idx],label='Eigen Decomp (6)')
plt.plot(range(1,k+1), desired_signal_CS,             label='CS (6)')
plt.legend(); plt.title('Wiener vs Eigen Decomp vs CS')
plt.xlabel('Sample Index'); plt.ylabel('Amplitude')

plt.figure()
plt.plot(range(1,k+1), d,                              label='Original')
plt.plot(range(1,k+1), desired_signal_eigen[rank_idx], label='Eigen Decomp (6)')
plt.plot(range(1,k+1), desired_signal_CS,              label='CS (6)')
plt.legend(); plt.title('Original vs ED vs CS')
plt.xlabel('Sample Index'); plt.ylabel('Amplitude')

plt.figure()
plt.plot(range(1,k+1), d,                 label='Original')
plt.plot(range(1,k+1), recovered_signal,  label='Full rank Wiener')
plt.plot(range(1,k+1), desired_signal_CS, label='CS (6)')
plt.legend(); plt.title('Original vs Full-rank Wiener vs CS')
plt.xlabel('Sample Index'); plt.ylabel('Amplitude')

print(f"Full-rank Wiener MSE : {MSE:.6f}")
print(f"Eigen Decomp MSE (6) : {MSE_eig[rank_idx]:.6f}")
print(f"CS MSE               : {MSE_CS:.6f}")

# =============================================================================
# 4. Lanczos-based MSNWF
# =============================================================================
# Index philosophy: keep MATLAB's 1-based indices intact.
# Every MATLAB access Cfirst(r, c) becomes Cfirst[r, c] (no offset arithmetic).
# Arrays are allocated with size D+2 so indices 2..D+1 are all valid.
#
#   MATLAB              Python
#   T(:,2)              T_L[:, 2]
#   Cfirst(2,2)         Cfirst[2, 2]
#   Cfirst(2:i, i)      Cfirst[2:i+1, i]   ← Python slice upper bound is exclusive
#   Clast(2, i-1)       Clast[2, i-1]
# =============================================================================

D      = 6
sz     = D + 2                          # indices 0..D+1; active range 2..D+1

Cfirst = np.zeros((sz, sz))
Clast  = np.zeros((sz, sz))
T_L    = np.zeros((p, sz))
r_L    = np.zeros((sz, sz))
b_L    = np.zeros(sz)
MSE_RR = np.zeros(D)

no          = np.linalg.norm(rxv2[:p])
T_L[:, 2]   = rxv2[:p] / no             # T(:,1)=zeros by default; T(:,2)=rxv2/no

r_L[2, 2]   = T_L[:, 2] @ Rv2 @ T_L[:, 2]
b_L[2]      = r_L[2, 2]
Cfirst[2,2] = 1.0 / r_L[2, 2]
Clast [2,2] = 1.0 / r_L[2, 2]
MSE_RR[0]   = rv1[0] - no**2 * Cfirst[2, 2]

for i in range(3, D + 2):               # i = 3, 4, ..., D+1  (exact MATLAB values)

    u            = (Rv2 @ T_L[:, i-1]
                    - r_L[i-1, i-1] * T_L[:, i-1]
                    - r_L[i-2, i-1] * T_L[:, i-2])

    r_L[i-1, i]  = np.linalg.norm(u)
    T_L[:, i]    = u / r_L[i-1, i]
    r_L[i,   i]  = T_L[:, i] @ Rv2 @ T_L[:, i]
    b_L[i]       = r_L[i, i] - r_L[i-1, i]**2 / b_L[i-1]

    # Cfirst(2:i, i) — MATLAB rows 2..i  →  Python slice [2 : i+1]
    # Cfirst(2:i-1, i-1) — MATLAB rows 2..i-1  →  Python slice [2 : i]
    Cfirst[2:i+1, i] = (
        np.append(Cfirst[2:i, i-1], 0.0)
        + (1.0 / b_L[i]) * Clast[2, i-1]
        * np.append(r_L[i-1, i]**2 * Clast[2:i, i-1], -r_L[i-1, i])
    )

    Clast[2:i+1, i] = (1.0 / b_L[i]) * np.append(
        -r_L[i-1, i] * Clast[2:i, i-1], 1.0
    )

    MSE_RR[i-2] = rv1[0] - no**2 * Cfirst[2, i]   # MSE_RR(i-1) in MATLAB

# V = T(:, 2:D+1) in MATLAB  →  columns 2..D+1 in Python
V_lanc = T_L[:, 2:D+2]                 # shape (p, D)

print(f"MSNWF MSE (stage 6)  : {MSE_RR[D-1]:.6f}")

# C_first = Cfirst(2:end, 2:end)  →  Cfirst[2:D+2, 2:D+2]
C_first_sub = Cfirst[2:D+2, 2:D+2]    # (D × D)
RR_wopt     = C_first_sub * no          # MATLAB: C_first * no  is scalar multiply

estimated_RR = np.zeros((D, k + p - 1))
recovered_RR = np.zeros((D, k))

for m in range(1, D + 1):
    conv_d = conv_rv2 @ V_lanc[:, :m]
    est    = conv_d @ RR_wopt[:m, m-1]       # RR_wopt(1:m, m) in MATLAB
    estimated_RR[m-1, :len(est)] = est
    recovered_RR[m-1, :]         = x - estimated_RR[m-1, :k]

plt.figure()
plt.plot(range(1,k+1), recovered_signal,   'r',
         label=f'Full rank Wiener MSE: {MSE:.4f}')
plt.plot(range(1,k+1), recovered_RR[5,:], 'k--',
         label=f'MSNWF MSE: {MSE_RR[5]:.4f}', markersize=3)
plt.plot(range(1,k+1), desired_signal_CS,  'b--',
         label=f'CS MSE: {MSE_CS:.4f}', markersize=3)
plt.title('Full rank Wiener and MSNWF, CF (r=6) Comparison', fontsize=13)
plt.xlabel('Sample Index', fontsize=13); plt.ylabel('Amplitude', fontsize=13)
plt.legend()

# =============================================================================
# 5. Filter order selection (Levinson over increasing order)
# =============================================================================
desired_point = 100
MSE_new = np.zeros(desired_point)

for newvar in range(1, desired_point + 1):
    Rv2_sub          = toeplitz(rv2[:newvar])
    w_lev            = Levinson(rv2[:newvar], rxv2[:newvar])
    MSE_new[newvar-1] = rv1[0] - w_lev @ Rv2_sub @ w_lev

plt.figure()
plt.plot(range(1, desired_point+1), MSE_new)
plt.title('MSE vs Filter Order (Levinson)')
plt.xlabel('Filter Order'); plt.ylabel('MSE')

plt.tight_layout()
plt.show()
