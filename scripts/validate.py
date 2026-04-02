"""
validate.py
===========
Validate the DVM simulation against Theodorsen's analytical solution
for a flat plate in pure sinusoidal plunge.

THEODORSEN THEORY (pure heave, h positive downward):
  L = π ρ b² (−ḧ)  +  2π ρ U b C(k) ḣ
  where b = c/2, k = ωc/(2U), C(k) = F(k) + iG(k) (Theodorsen function)

Expanding for h(t) = h₀ sin(ωt):
  Cl(t) = 2π (h₀/c) · [ k² sin(ωt)  +  2k (F cos(ωt) − G sin(ωt)) ]

Quasi-steady limit (k → 0, C → 1):
  Cl(t) ≈ 4π kh cos(ωt)   where kh = k · h₀/c

Theodorsen function:
  C(k) = H₁⁽²⁾(k) / [H₁⁽²⁾(k) + i H₀⁽²⁾(k)]
  H_n⁽²⁾ = J_n − i Y_n   (Hankel function of the second kind)
  F(k) = Re[C(k)],   G(k) = −Im[C(k)]

SIGN CONVENTION NOTE:
  The Fortran code uses h positive DOWNWARD, matching Theodorsen.
  Downward motion → positive Cl → positive circulation.
  The validation comparison is therefore direct with no sign flip needed.

Usage:
  python validate.py [--config config.nml] [--output_dir output]
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, yv  # Bessel functions J_n and Y_n

# ─────────────────────────────────────────────────────────────
# Theodorsen function
# ─────────────────────────────────────────────────────────────


def theodorsen_FG(k):
    """
    Compute F(k) and G(k) from the Theodorsen function C(k) = F + iG.

    C(k) = H₁⁽²⁾(k) / [H₁⁽²⁾(k) + i H₀⁽²⁾(k)]

    Using the explicit formula:
      F =  (J₁(J₁+Y₀) + Y₁(J₀−Y₁)) / D
        =  (J₁² + Y₁² + J₁Y₀ − Y₁J₀) / D
      G =  (J₁J₀ + Y₁Y₀) / D
      D =  (J₁ + Y₀)² + (J₀ − Y₁)²
    """
    J0 = jv(0, k)
    J1 = jv(1, k)
    Y0 = yv(0, k)
    Y1 = yv(1, k)

    D = (J1 + Y0) ** 2 + (J0 - Y1) ** 2

    F = (J1**2 + Y1**2 + J1 * Y0 - Y1 * J0) / D
    G = (J1 * J0 + Y1 * Y0) / D
    return F, G


def theodorsen_Cl(t_arr, k, h0_over_c, omega):
    """
    Analytical Cl for pure sinusoidal plunge (h positive downward).
    Cl(t) = 2π A_h [k² sin(ωt) + 2k (F cos(ωt) − G sin(ωt))]
    where A_h = h₀/c.
    """
    F, G = theodorsen_FG(k)
    A_h = h0_over_c
    phase = omega * t_arr
    Cl = (
        2
        * np.pi
        * A_h
        * (k**2 * np.sin(phase) + 2 * k * (F * np.cos(phase) - G * np.sin(phase)))
    )
    return Cl


def theodorsen_Cl_amplitude_phase(k, h0_over_c):
    """
    Return the complex amplitude of Cl as a function of reduced frequency.
    Cl(t) = Re[Cl_amp * exp(iωt)]
    """
    F, G = theodorsen_FG(k)
    A_h = h0_over_c
    # Cl(t) = 2π A_h [ (2kF) cos(ωt) + (k² − 2kG) sin(ωt) ]
    # = Re[ Cl_amp * e^{iωt} ] with Cl_amp = 2π A_h [ 2kF − i(k²−2kG) ]
    #   ... but easier to just return |Cl| and phase_lag.
    # Using phasor: Cl(t) = Cl_max cos(ωt + φ)
    cos_coeff = 2 * np.pi * A_h * 2 * k * F
    sin_coeff = 2 * np.pi * A_h * (k**2 - 2 * k * G)
    amplitude = np.sqrt(cos_coeff**2 + sin_coeff**2)
    phase_lag = np.arctan2(-sin_coeff, cos_coeff)  # lag behind plunge velocity
    return amplitude, np.degrees(phase_lag)


# ─────────────────────────────────────────────────────────────
# Namelist parser (minimal, for reading config.nml)
# ─────────────────────────────────────────────────────────────


def parse_nml(filename):
    """
    Parse a Fortran namelist file into a flat dict.
    Handles:  key = value  and  key = 'string'
    Returns all keys in lower-case.
    """
    params = {}
    try:
        with open(filename) as f:
            for raw in f:
                line = raw.strip()
                if (
                    not line
                    or line.startswith("!")
                    or line.startswith("&")
                    or line.startswith("/")
                ):
                    continue
                if "=" not in line:
                    continue
                key, _, rest = line.partition("=")
                key = key.strip().lower()
                val = rest.split("!")[0].strip().rstrip(",").strip()
                # Remove Fortran quotes
                val = val.strip("'\"")
                # Try numeric conversion
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val
    except FileNotFoundError:
        print(f"  WARNING: {filename} not found; using hard-coded defaults.")
    return params


# ─────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────


def load_forces(filepath):
    """Load forces.dat → numpy array.  Skips lines starting with '#'."""
    if not os.path.exists(filepath):
        sys.exit(
            f"  ERROR: Force file not found: {filepath}\n"
            f"  Run the Fortran simulation first: make run"
        )
    data = np.loadtxt(filepath, comments="#")
    # Columns: t  Cl  Ct  Cm_qc  Gamma_total
    return data


def last_cycle(data, n_cycles, steps_per_cycle, n_steady=2):
    """
    Return the last n_steady cycles of data for steady-state comparison.
    The first few cycles contain the impulsive-start transient.
    """
    rows_per_cycle = steps_per_cycle
    start = max(0, len(data) - n_steady * rows_per_cycle)
    return data[start:]


# ─────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────


def plot_time_series(t, Cl_dvm, Cl_theo, omega, k, kh, ax):
    cycles = t * omega / (2 * np.pi)
    ax.plot(cycles, Cl_dvm, lw=2.0, color="tab:blue", label="DVM")
    ax.plot(
        cycles, Cl_theo, lw=1.5, color="tab:orange", linestyle="--", label="Theodorsen"
    )
    ax.set_xlabel("Number of cycles")
    ax.set_ylabel(r"$C_L$")
    ax.set_title(rf"Lift history  ($k={k:.3f}$,  $kh={kh:.3f}$)")
    ax.legend()
    ax.grid(True, alpha=0.4)


def plot_steady_comparison(t_ss, Cl_dvm_ss, Cl_theo_ss, omega, k, kh, ax):
    cycles = t_ss * omega / (2 * np.pi)
    ax.plot(
        cycles - cycles[0],
        Cl_dvm_ss,
        lw=2.0,
        color="tab:blue",
        label=f"DVM   (amp={np.max(np.abs(Cl_dvm_ss)):.4f})",
    )
    ax.plot(
        cycles - cycles[0],
        Cl_theo_ss,
        lw=1.5,
        color="tab:orange",
        linestyle="--",
        label=f"Theo  (amp={np.max(np.abs(Cl_theo_ss)):.4f})",
    )
    ax.set_xlabel("Cycles (steady-state window)")
    ax.set_ylabel(r"$C_L$")
    ax.set_title(rf"Steady-state comparison  ($k={k:.3f}$,  $kh={kh:.3f}$)")
    ax.legend()
    ax.grid(True, alpha=0.4)


def plot_FG_vs_k(ax):
    """Plot Theodorsen function F(k) and G(k) as a reference."""
    k_arr = np.linspace(0.01, 2.0, 300)
    F_arr = np.array([theodorsen_FG(k)[0] for k in k_arr])
    G_arr = np.array([theodorsen_FG(k)[1] for k in k_arr])
    ax.plot(k_arr, F_arr, label=r"$F(k)$", color="tab:green", lw=2)
    ax.plot(k_arr, G_arr, label=r"$G(k)$", color="tab:red", lw=2)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"Reduced frequency $k$")
    ax.set_ylabel(r"$F$,  $G$")
    ax.set_title("Theodorsen function components")
    ax.legend()
    ax.grid(True, alpha=0.4)


def plot_Cl_amplitude_vs_k(ax, h0_over_c=0.1):
    """
    Normalised Cl amplitude vs k for pure plunge.
    Cl_max / (4π kh)  →  1 in the quasi-steady limit.
    """
    k_arr = np.linspace(0.01, 2.0, 300)
    amps = np.array([theodorsen_Cl_amplitude_phase(k, h0_over_c * k)[0] for k in k_arr])
    # Normalise by quasi-steady amplitude 4π kh
    norm = 4 * np.pi * k_arr * h0_over_c
    ax.plot(
        k_arr,
        amps / norm,
        color="tab:purple",
        lw=2,
        label=r"$C_{L,\max}\,/\,(4\pi kh)$",
    )
    ax.axhline(1.0, color="k", lw=0.8, linestyle="--", label="Quasi-steady")
    ax.set_xlabel(r"Reduced frequency $k$")
    ax.set_ylabel(r"Normalised $C_{L}$ amplitude")
    ax.set_title("Amplitude attenuation by wake effects")
    ax.legend()
    ax.grid(True, alpha=0.4)


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Validate DVM against Theodorsen")
    parser.add_argument("--config", default="config.nml", help="Namelist config file")
    parser.add_argument("--output_dir", default=None, help="Override output directory")
    args = parser.parse_args()

    # ── Read config ──────────────────────────────────────────
    cfg = parse_nml(args.config)
    k = float(cfg.get("reduced_freq", 0.5))
    h0_over_c = float(cfg.get("h0", 0.1))  # h0/chord since chord=1 default
    chord = float(cfg.get("chord", 1.0))
    U_inf = float(cfg.get("u_inf", 1.0))
    N_cycles = int(cfg.get("n_cycles", 6))
    spc = int(cfg.get("steps_per_cycle", 100))
    out_dir = args.output_dir or str(cfg.get("output_dir", "output"))

    h0_over_c /= chord  # normalise if chord ≠ 1
    omega = 2.0 * k * U_inf / chord
    kh = k * h0_over_c

    print(f"\n  Validation parameters")
    print(f"  ── k       = {k}")
    print(f"  ── h0/c    = {h0_over_c}")
    print(f"  ── kh      = {kh:.4f}")
    print(f"  ── omega   = {omega:.4f}")

    # ── Load simulation output ────────────────────────────────
    force_path = cfg.get("force_file", "forces.dat")
    force_path = os.path.join(out_dir, force_path)
    data = load_forces(force_path)

    t = data[:, 0]
    Cl_dvm = data[:, 1]

    # ── Theodorsen Cl at same time points ────────────────────
    Cl_theo = theodorsen_Cl(t, k, h0_over_c, omega)

    # ── Steady-state window (last 2 cycles) ──────────────────
    ss_data = last_cycle(data, N_cycles, spc, n_steady=2)
    t_ss = ss_data[:, 0]
    Cl_dvm_ss = ss_data[:, 1]
    Cl_theo_ss = theodorsen_Cl(t_ss, k, h0_over_c, omega)

    # ── Quantitative comparison ───────────────────────────────
    amp_dvm = np.max(np.abs(Cl_dvm_ss))
    amp_theo = np.max(np.abs(Cl_theo_ss))
    amp_err = abs(amp_dvm - amp_theo) / amp_theo * 100
    qs_amp = 4 * np.pi * kh
    print(f"\n  ── Steady-state Cl amplitude ──────────────────")
    print(f"     DVM:              {amp_dvm:.5f}")
    print(f"     Theodorsen:       {amp_theo:.5f}")
    print(f"     Quasi-steady:     {qs_amp:.5f}")
    print(f"     Error DVM/Theo:   {amp_err:.2f}%")
    F, G = theodorsen_FG(k)
    print(f"     F(k)={F:.4f}   G(k)={G:.4f}")

    # ── Plots ─────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"DVM Validation vs Theodorsen  (k={k}, kh={kh:.3f})",
        fontsize=14,
        fontweight="bold",
    )

    # (0,0) – Full time history
    plot_time_series(t, Cl_dvm, Cl_theo, omega, k, kh, axes[0, 0])

    # (0,1) – Steady-state zoom
    plot_steady_comparison(t_ss, Cl_dvm_ss, Cl_theo_ss, omega, k, kh, axes[0, 1])

    # (1,0) – F and G vs k
    plot_FG_vs_k(axes[1, 0])
    axes[1, 0].axvline(k, color="k", lw=0.8, linestyle=":", label=f"current k={k}")
    axes[1, 0].legend()

    # (1,1) – Cl amplitude vs k
    plot_Cl_amplitude_vs_k(axes[1, 1], h0_over_c)
    axes[1, 1].axvline(k, color="k", lw=0.8, linestyle=":", label=f"current k={k}")
    axes[1, 1].legend()

    plt.tight_layout()
    out_png = os.path.join(out_dir, "validation.png")
    fig.savefig(out_png, dpi=150)
    print(f"\n  Saved: {out_png}")
    plt.show()


if __name__ == "__main__":
    main()
