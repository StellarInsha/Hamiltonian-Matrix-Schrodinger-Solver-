
import numpy as np
import matplotlib.pyplot as plt
 
 
# ── Solver (identical to the Harmonic Oscillator one) ─────────────────────────────────────────────────
 
x_min, x_max, N = -8.0, 8.0, 1000
x = np.linspace(x_min, x_max, N)
h = x[1] - x[0]
 
def solve(V_func, n_states=6):
    """Build H matrix for given V(x), return eigenvalues and wavefunctions."""
    diag = 1.0/h**2 + V_func(x)
    offd = -0.5/h**2 * np.ones(N - 1)
    H    = np.diag(diag) + np.diag(offd, 1) + np.diag(offd, -1)
 
    vals, vecs = np.linalg.eigh(H)
 
    psis = []
    for n in range(n_states):
        psi  = vecs[:, n]
        norm = np.sqrt(np.trapezoid(psi**2, x))
        psis.append(psi / norm)
 
    return vals, psis
 
 
# ── Three potentials ──────────────────────────────────────────────────────────
 
# 1. Anharmonic oscillator
#    Harmonic oscillator + small x^4 correction.
#    Real molecules are never perfect parabolas so this is the first correction
#    The x^4 term pushes the walls outward, so energy levels shift up and slopes are steep
#    The equal spacing of the HO is broken.
def V_anharmonic(x):
    return 0.5*x**2 + 0.1*x**4
 
# 2. Double well
#    Two bowls separated by a central hill.
#    The two lowest states are nearly degenerate (almost same energy).
#    That tiny gap between E_0 and E_1 IS quantum tunneling —
#    the particle can tunnel through the barrier even with E < hill height.
def V_double_well(x):
    return (x**2 - 2.0)**2
 
# 3. Finite square well
#    You know this from Griffiths. Griffiths solves it graphically
#    using a transcendental equation. Your code gets the same answer numerically.
#    This is a direct check: does your solver reproduce Griffiths?
def V_finite_well(x):
    return np.where(np.abs(x) < 2.0, -10.0, 0.0)
 
 
# ── Solve all three ───────────────────────────────────────────────────────────
 
E_anharmonic,  psi_anharmonic  = solve(V_anharmonic)
E_double_well, psi_double_well = solve(V_double_well)
E_finite_well, psi_finite_well = solve(V_finite_well)
 
 
# ── Print results ─────────────────────────────────────────────────────────────
 
print("── Anharmonic oscillator ─────────────────────────────")
print("   Compare with HO exact (n+0.5): levels shift UP, spacing grows")
print(f"  {'n':>3}  {'Anharmonic E':>14}  {'HO exact':>10}  {'shift':>8}")
for n in range(5):
    ho = n + 0.5
    ah = E_anharmonic[n]
    print(f"  {n:>3}  {ah:>14.6f}  {ho:>10.4f}  {ah-ho:>+8.4f}")
 
print()
print("── Double well ───────────────────────────────────────")
print("   E_0 and E_1 are nearly equal — the gap = tunneling splitting")
for n in range(5):
    print(f"  n={n}: E = {E_double_well[n]:.6f}")
print(f"  Tunneling splitting ΔE = {E_double_well[1] - E_double_well[0]:.6f}")
 
print()
print("── Finite square well ────────────────────────────────")
print("   Negative energies = bound states inside the well (depth = 10)")
for n in range(5):
    print(f"  n={n}: E = {E_finite_well[n]:.6f}")
 
 
# ── Plotting them all ──────────────────────────────────────────────────────────────────────
 
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle("Three Potentials With No Analytical Solution", fontsize=14, fontweight='bold')
 
configs = [
    (V_anharmonic,  E_anharmonic,  psi_anharmonic,
     "Anharmonic  V = ½x² + 0.1x⁴", (-5,5), (-0.5, 8),  4),
    (V_double_well, E_double_well, psi_double_well,
     "Double well  V = (x²−2)²",     (-4,4), (-0.5, 6),  4),
    (V_finite_well, E_finite_well, psi_finite_well,
     "Finite square well  V₀ = 10",  (-5,5), (-11,  2),  3),
]
 
for col, (Vf, evals, wavefns, title, xlim, ylim, n_show) in enumerate(configs):
    Vx = Vf(x)
 
    # Top row: potential + energy levels
    ax = axes[0, col]
    ax.plot(x, Vx, 'k-', lw=2.5)
    ax.fill_between(x, Vx, alpha=0.07, color='gray')
    for n in range(n_show):
        E = evals[n]
        if ylim[0] < E < ylim[1]:
            ax.axhline(y=E, color=f'C{n}', ls='--', lw=1.6,
                       label=f'E_{n} = {E:.3f}')
    ax.set_xlim(*xlim);  ax.set_ylim(*ylim)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('x');  ax.set_ylabel('Energy')
    ax.legend(fontsize=8);  ax.grid(alpha=0.2)
 
    # Bottom row: wavefunctions raised to their energy levels
    ax = axes[1, col]
    ax.plot(x, Vx, 'k-', lw=1.5, alpha=0.3)
    ax.fill_between(x, Vx, alpha=0.07, color='gray')
    for n in range(n_show):
        E = evals[n]
        if ylim[0] < E < ylim[1]:
            ax.plot(x, wavefns[n]*0.5 + E, color=f'C{n}', lw=1.8,
                    label=f'ψ_{n}')
            ax.axhline(y=E, color=f'C{n}', ls=':', lw=0.8, alpha=0.4)
    ax.set_xlim(*xlim);  ax.set_ylim(*ylim)
    ax.set_title("Wavefunctions", fontsize=10)
    ax.set_xlabel('x')
    ax.legend(fontsize=8);  ax.grid(alpha=0.2)
 
plt.tight_layout()
plt.savefig('02_unsolvable_potentials.png', dpi=150, bbox_inches='tight')
plt.show()
