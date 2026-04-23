"""
MSNWF for Adaptive Noise Cancellation - Python Implementation
Based on Lanczos algorithm for reduced-rank Wiener filtering
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz, solve

def auto_correlation(x, max_lag=None):
    """
    Compute autocorrelation function
    """
    n = len(x)
    if max_lag is None:
        max_lag = n - 1
    
    r = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        r[lag] = np.mean(x[:n-lag] * x[lag:])
    return r

def cross_correlation(x, y, max_lag=None):
    """
    Compute cross-correlation function
    """
    n = min(len(x), len(y))
    if max_lag is None:
        max_lag = n - 1
    
    r = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        r[lag] = np.mean(x[:n-lag] * y[lag:])
    return r

def generate_ar_noise(coeffs, noise, length):
    """
    Generate AR noise: similar to AutoRnoise in MATLAB
    coeffs: AR coefficients [1, a1, a2, ...]
    noise: white noise input
    """
    # Using filter: y[n] = x[n] + a1*y[n-1] + a2*y[n-2] + ...
    b = np.array([1.0])
    a = np.array(coeffs)
    
    # Initialize output
    y = np.zeros(length)
    ar_order = len(a) - 1
    
    for n in range(length):
        y[n] = noise[n]
        for k in range(1, min(n+1, ar_order+1)):
            y[n] -= a[k] * y[n-k]
    
    return y

def levinson_durbin(r):
    """
    Levinson-Durbin recursion for Toeplitz systems
    Solves R*w = r where R is Toeplitz
    """
    n = len(r)
    w = np.zeros(n)
    
    # Use numpy solve for Toeplitz matrix
    R = toeplitz(r)
    w = solve(R, r, assume_a='pos')
    
    return w

def msnwf_lanczos(Rx, rxd, sigma_d2, D):
    """
    Multi-Stage Nested Wiener Filter using Lanczos algorithm
    
    Parameters:
    -----------
    Rx : ndarray (p x p)
        Autocorrelation matrix of observation
    rxd : ndarray (p,)
        Cross-correlation vector between observation and desired signal
    sigma_d2 : float
        Variance of desired signal
    D : int
        Number of stages (reduced rank)
    
    Returns:
    --------
    V : ndarray (p x D)
        Krylov subspace basis vectors
    RR_wopt : ndarray (D x D)
        Reduced rank optimal weights for each stage
    MSE_RR : ndarray (D,)
        MSE at each stage
    """
    p = Rx.shape[0]
    
    # Initialize matrices
    Cfirst = np.zeros((D+1, D+1))
    Clast = np.zeros((D+1, D+1))
    T = np.zeros((p, D+2))
    r = np.zeros((D+2, D+2))
    b = np.zeros(D+2)
    MSE_RR = np.zeros(D)
    
    # Initialization
    no = np.linalg.norm(rxd)
    T[:, 1] = rxd / no
    
    r[0, 1] = 0
    r[1, 1] = T[:, 1].T @ Rx @ T[:, 1]
    b[1] = r[1, 1]
    Cfirst[1, 1] = 1.0 / r[1, 1]
    Clast[1, 1] = 1.0 / r[1, 1]
    
    MSE_RR[0] = sigma_d2 - (no**2) * Cfirst[1, 1]
    
    # Lanczos iteration
    for i in range(2, D+1):
        u = Rx @ T[:, i-1] - r[i-1, i-1] * T[:, i-1] - r[i-2, i-1] * T[:, i-2]
        r[i-1, i] = np.linalg.norm(u)
        T[:, i] = u / r[i-1, i]
        r[i, i] = T[:, i].T @ Rx @ T[:, i]
        b[i] = r[i, i] - r[i-1, i]**2 / b[i-1]
        
        # Update Cfirst
        temp1 = np.concatenate([Cfirst[1:i, i-1], [0]])
        temp2 = np.concatenate([r[i-1, i]**2 * Clast[1:i, i-1], [-r[i-1, i]]])
        Cfirst[1:i+1, i] = temp1 + (1.0/b[i]) * Clast[1, i-1] * temp2
        
        # Update Clast
        temp3 = np.concatenate([-r[i-1, i] * Clast[1:i, i-1], [1]])
        Clast[1:i+1, i] = (1.0/b[i]) * temp3
        
        MSE_RR[i-1] = sigma_d2 - (no**2) * Cfirst[1, i]
    
    # Extract results
    V = T[:, 1:D+1]
    C_first = Cfirst[1:D+1, 1:D+1]
    RR_wopt = C_first * no
    
    return V, RR_wopt, MSE_RR

def main():
    """
    Main function: Adaptive noise cancellation using MSNWF
    """
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    #=========================================================================
    # 1. GENERATE SIGNALS
    #=========================================================================
    
    # Desired signal (sinusoid)
    w0 = 0.05 * np.pi
    n = np.arange(200)
    phi = np.random.rand() * 2 * np.pi
    d = np.sin(w0 * n + phi)
    k = len(d)
    
    # Noise source
    noi = np.random.randn(k)
    
    # Primary sensor noise (correlated with source)
    v1 = generate_ar_noise([1, -0.8], noi, k)
    
    # Secondary sensor noise (correlated with source, different AR model)
    v2 = generate_ar_noise([1, 0.6], noi, k)
    
    # Observation (signal + noise)
    x = d + v1
    
    #=========================================================================
    # 2. COMPUTE STATISTICS
    #=========================================================================
    
    # Autocorrelation of v2
    rv2 = auto_correlation(v2)
    
    # Cross-correlation between x and v2
    rxv2 = cross_correlation(x, v2)
    
    # Autocorrelation of v1 (for MSE computation)
    rv1 = auto_correlation(v1)
    
    #=========================================================================
    # 3. CONVENTIONAL WIENER SOLUTION (Full Rank)
    #=========================================================================
    
    p = 12  # Order of the noise estimator
    
    # Form autocorrelation matrix
    Rv2 = toeplitz(rv2[:p])
    
    # Optimal Wiener filter
    w_opt = solve(Rv2, rxv2[:p], assume_a='pos')
    
    # Estimate noise
    estimated_noise = np.convolve(v2, w_opt, mode='full')
    
    # Recover signal
    recovered_signal = x - estimated_noise[:k]
    
    # MSE for full rank Wiener
    MSE_full = rv1[0] - w_opt.T @ Rv2 @ w_opt
    
    #=========================================================================
    # 4. MSNWF (Lanczos-based)
    #=========================================================================
    
    D = 6  # Number of stages
    
    # Variance of v1
    sigma_v1_2 = rv1[0]
    
    # Apply MSNWF
    V, RR_wopt, MSE_RR = msnwf_lanczos(Rv2, rxv2[:p], sigma_v1_2, D)
    
    # Estimate noise using MSNWF at different stages
    recovered_signal_RR = np.zeros((D, k))
    for m in range(D):
        w_rr = V[:, :m+1] @ RR_wopt[:m+1, m]
        estimated_RR = np.convolve(v2, w_rr, mode='full')
        recovered_signal_RR[m, :] = x - estimated_RR[:k]
    
    #=========================================================================
    # 5. CROSS-SPECTRAL (CS) METRIC (for comparison)
    #=========================================================================
    
    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(Rv2)
    
    # Sort in descending order
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # For rank 6 CS metric
    o = 6
    V_CS = eigenvectors[:, :o]
    D_CS = np.diag(eigenvalues[:o])
    
    # CS metric optimal solution
    w_CS = V_CS @ np.linalg.inv(D_CS) @ V_CS.T @ rxv2[:p]
    estimated_noise_CS = np.convolve(v2, w_CS, mode='full')
    desired_signal_CS = x - estimated_noise_CS[:k]
    
    MSE_CS = rv1[0] - rxv2[:p].T @ V_CS @ np.linalg.inv(D_CS) @ V_CS.T @ rxv2[:p]
    
    #=========================================================================
    # 6. VISUALIZATION
    #=========================================================================
    
    # Plot 1: Original vs Recovered Signal
    plt.figure(figsize=(15, 10))
    
    plt.subplot(3, 2, 1)
    plt.plot(n, d, 'r', label='Original Signal', linewidth=2)
    plt.plot(n, recovered_signal, 'b--', label='Full Rank Wiener', linewidth=1.5)
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')
    plt.title('Original vs Full Rank Wiener')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 2, 2)
    plt.plot(n, d, 'r', label='Original Signal', linewidth=2)
    plt.plot(n, recovered_signal_RR[D-1, :], 'g--', label=f'MSNWF (D={D})', linewidth=1.5)
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')
    plt.title(f'Original vs MSNWF (D={D})')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 2, 3)
    plt.plot(n, recovered_signal, 'b', label='Full Rank Wiener', linewidth=1.5)
    plt.plot(n, recovered_signal_RR[D-1, :], 'g--', label=f'MSNWF (D={D})', linewidth=1.5)
    plt.plot(n, desired_signal_CS, 'k:', label='CS (6)', linewidth=1.5)
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')
    plt.title('Comparison: Full Rank vs MSNWF vs CS')
    plt.legend()
    plt.grid(True)
    
    # Plot 2: MSE vs Rank
    plt.subplot(3, 2, 4)
    plt.plot(range(1, D+1), MSE_RR, 'go-', linewidth=2, markersize=8, label='MSNWF')
    plt.axhline(y=MSE_full, color='b', linestyle='--', linewidth=2, label='Full Rank Wiener')
    plt.axhline(y=MSE_CS, color='k', linestyle=':', linewidth=2, label='CS (6)')
    plt.xlabel('Number of Stages (Rank)')
    plt.ylabel('Mean Squared Error')
    plt.title('MSE vs Rank')
    plt.legend()
    plt.grid(True)
    
    # Plot 3: Filter weights comparison
    plt.subplot(3, 2, 5)
    plt.stem(w_opt, linefmt='b-', markerfmt='bo', basefmt='b-', label='Full Rank')
    plt.stem(V[:, :D] @ RR_wopt[:D, D-1], linefmt='g-', markerfmt='gs', 
             basefmt='g-', label=f'MSNWF (D={D})')
    plt.xlabel('Tap Index')
    plt.ylabel('Filter Weight')
    plt.title('Filter Weights Comparison')
    plt.legend()
    plt.grid(True)
    
    # Plot 4: Noise estimates
    plt.subplot(3, 2, 6)
    plt.plot(n, v1, 'r', label='True Noise (v1)', linewidth=1.5)
    plt.plot(n, estimated_noise[:k], 'b--', label='Estimated (Full Rank)', linewidth=1.5)
    w_rr_final = V[:, :D] @ RR_wopt[:D, D-1]
    estimated_RR_final = np.convolve(v2, w_rr_final, mode='full')
    plt.plot(n, estimated_RR_final[:k], 'g:', label=f'Estimated (MSNWF D={D})', linewidth=1.5)
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')
    plt.title('True vs Estimated Noise')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    
    #=========================================================================
    # 7. PRINT RESULTS
    #=========================================================================
    
    print("\n" + "="*60)
    print("ADAPTIVE NOISE CANCELLATION - RESULTS")
    print("="*60)
    print(f"Filter Order (p): {p}")
    print(f"Number of Samples: {k}")
    print(f"\nMEAN SQUARED ERROR:")
    print(f"  Full Rank Wiener:  {MSE_full:.6f}")
    print(f"  MSNWF (D={D}):      {MSE_RR[D-1]:.6f}")
    print(f"  CS Metric (6):     {MSE_CS:.6f}")
    print(f"\nMSE at each MSNWF stage:")
    for i, mse in enumerate(MSE_RR):
        print(f"  Stage {i+1}: {mse:.6f}")
    print("="*60 + "\n")
    
    plt.show()

if __name__ == "__main__":
    main()
