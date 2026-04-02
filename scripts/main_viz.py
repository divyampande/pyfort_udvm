"""
main_viz.py
===========
Comprehensive post-processing and visualisation for the Unsteady DVM.

Capabilities
────────────
1. Force history plot  (Cl, Ct, Cm vs non-dimensional time)
2. Wake pattern animation  (matplotlib)
3. Wake vorticity colourmap at a chosen snapshot
4. Parametric kh sweep  (requires re-running the Fortran code)
5. Cl / Ct vs kh  (requires sweep)
6. Interactive Plotly wake animation  (optional, writes HTML)

Usage
─────
  # Plot a single run:
  python main_viz.py --mode single

  # Run a kh sweep (re-runs ./udvm for each kh value):
  python main_viz.py --mode sweep --kh_values 0.05 0.1 0.2 0.3 0.5

  # Animate wake from snapshot files:
  python main_viz.py --mode animate

  # All of the above:
  python main_viz.py --mode all

Dependencies
────────────
  numpy, matplotlib, scipy (optional: plotly)
"""

import argparse
import os
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

try:
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

# ─────────────────────────────────────────────────────────────
# Config / namelist helpers
# ─────────────────────────────────────────────────────────────


def parse_nml(filename):
    """Minimal Fortran namelist parser → flat dict (lowercase keys)."""
    params = {}
    try:
        with open(filename) as f:
            for raw in f:
                line = raw.strip()
                if not line or line[0] in ("!", "&", "/"):
                    continue
                if "=" not in line:
                    continue
                key, _, rest = line.partition("=")
                key = key.strip().lower()
                val = rest.split("!")[0].strip().rstrip(",").strip().strip("'\"")
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val
    except FileNotFoundError:
        print(f"  WARNING: config file {filename} not found")
    return params


def write_nml_override(base_nml, out_nml, overrides):
    """
    Write a new namelist file based on base_nml with key=value overrides.
    Used for the kh sweep to vary reduced_freq and h0.
    """
    with open(base_nml) as f:
        lines = f.readlines()
    with open(out_nml, "w") as f:
        for line in lines:
            written = False
            for key, val in overrides.items():
                stripped = line.strip().lower()
                if stripped.startswith(key + " ") or stripped.startswith(key + "="):
                    if isinstance(val, str):
                        f.write(f"  {key} = '{val}'\n")
                    else:
                        f.write(f"  {key} = {val}\n")
                    written = True
                    break
            if not written:
                f.write(line)


# ─────────────────────────────────────────────────────────────
# Data loaders
# ─────────────────────────────────────────────────────────────


def load_forces(output_dir, force_file="forces.dat"):
    path = os.path.join(output_dir, force_file)
    if not os.path.exists(path):
        sys.exit(f"  ERROR: {path} not found.  Run simulation first.")
    data = np.loadtxt(path, comments="#")
    # Columns: t  Cl  Ct  Cm_qc  Gamma_total
    cols = {"t": 0, "Cl": 1, "Ct": 2, "Cm_qc": 3, "Gamma_total": 4}
    return {name: data[:, idx] for name, idx in cols.items()}


def load_wake_snapshot(filepath):
    """
    Load one wake snapshot file.
    Returns dict with keys: type (0=bound,1=wake), x, y, Gamma.
    """
    if not os.path.exists(filepath):
        return None
    raw = np.loadtxt(filepath, comments="#")
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    return {
        "type": raw[:, 0].astype(int),
        "x": raw[:, 1],
        "y": raw[:, 2],
        "Gamma": raw[:, 3],
    }


def load_all_snapshots(output_dir, base="wake"):
    """Load all wake snapshots sorted by step number."""
    pattern = os.path.join(output_dir, f"{base}_*.dat")
    files = sorted(glob.glob(pattern))
    snaps = []
    for f in files:
        s = load_wake_snapshot(f)
        if s is not None:
            step = int(os.path.basename(f).split("_")[-1].split(".")[0])
            snaps.append((step, s))
    return snaps


# ─────────────────────────────────────────────────────────────
# 1. Force history plot
# ─────────────────────────────────────────────────────────────


def plot_forces(cfg, output_dir="output", save=True):
    """Plot Cl, Ct, Cm vs time and vs number of cycles."""
    forces = load_forces(output_dir)
    k = float(cfg.get("reduced_freq", 0.5))
    h0_c = float(cfg.get("h0", 0.1))
    chord = float(cfg.get("chord", 1.0))
    U_inf = float(cfg.get("u_inf", 1.0))
    omega = 2.0 * k * U_inf / chord
    kh = k * h0_c / chord

    t = forces["t"]
    cycles = t * omega / (2 * np.pi)
    Cl = forces["Cl"]
    Ct = forces["Ct"]
    Cm = forces["Cm_qc"]

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle(rf"Force History  $k={k:.3f}$,  $kh={kh:.4f}$", fontsize=14)

    axes[0].plot(cycles, Cl, lw=1.5, color="tab:blue")
    axes[0].set_ylabel(r"$C_L$")
    axes[0].grid(True, alpha=0.3)
    axes[0].axhline(0, color="k", lw=0.5)

    axes[1].plot(cycles, Ct, lw=1.5, color="tab:red")
    axes[1].set_ylabel(r"$C_T$")
    axes[1].grid(True, alpha=0.3)
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_title("(Ct = 0 until leading-edge suction is implemented)", fontsize=9)

    axes[2].plot(cycles, Cm, lw=1.5, color="tab:green")
    axes[2].set_ylabel(r"$C_{m,c/4}$")
    axes[2].set_xlabel("Number of cycles")
    axes[2].grid(True, alpha=0.3)
    axes[2].axhline(0, color="k", lw=0.5)

    plt.tight_layout()
    if save:
        path = os.path.join(output_dir, "forces.png")
        fig.savefig(path, dpi=150)
        print(f"  Saved: {path}")
    plt.show()


# ─────────────────────────────────────────────────────────────
# 2. Wake pattern at a single snapshot
# ─────────────────────────────────────────────────────────────


def plot_wake_snapshot(snap, cfg, step=None, output_dir="output", save=True):
    """
    Plot bound + wake vortex positions, coloured by Gamma.
    """
    chord = float(cfg.get("chord", 1.0))
    mask_bound = snap["type"] == 0
    mask_wake = snap["type"] == 1

    xb, yb, Gb = snap["x"][mask_bound], snap["y"][mask_bound], snap["Gamma"][mask_bound]
    xw, yw, Gw = snap["x"][mask_wake], snap["y"][mask_wake], snap["Gamma"][mask_wake]

    fig, ax = plt.subplots(figsize=(14, 4))

    # Colour scale: symmetric around 0
    all_G = np.concatenate([Gb, Gw]) if len(Gw) > 0 else Gb
    vmax = np.max(np.abs(all_G)) if len(all_G) > 0 else 1.0

    cmap = cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Bound vortices
    sc_b = ax.scatter(
        xb, yb, c=Gb, cmap=cmap, norm=norm, s=60, marker="o", zorder=5, label="Bound Γ"
    )
    # Wake vortices
    if len(xw) > 0:
        sc_w = ax.scatter(
            xw,
            yw,
            c=Gw,
            cmap=cmap,
            norm=norm,
            s=18,
            marker="^",
            alpha=0.85,
            zorder=4,
            label="Wake Γ",
        )

    # Flat plate outline
    ax.plot([0, chord], [0, 0], "k-", lw=3, label="Flat plate")

    plt.colorbar(sc_b, ax=ax, label=r"$\Gamma$", shrink=0.8)
    ax.set_xlabel("x / chord")
    ax.set_ylabel("y / chord")
    title_str = f"Wake pattern (step {step})" if step is not None else "Wake pattern"
    ax.set_title(title_str)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()

    if save and step is not None:
        path = os.path.join(output_dir, f"wake_plot_{step:06d}.png")
        fig.savefig(path, dpi=120)
    plt.show()


# ─────────────────────────────────────────────────────────────
# 3. Animated wake (matplotlib)
# ─────────────────────────────────────────────────────────────


def animate_wake(
    cfg, output_dir="output", snap_base="wake", save_gif=True, interval=80
):
    """
    Animate wake evolution from snapshot files.
    Saves an animated GIF (or uses FFMpeg if available).
    """
    chord = float(cfg.get("chord", 1.0))
    snaps = load_all_snapshots(output_dir, snap_base)
    if not snaps:
        print("  No wake snapshots found.  Check snap_interval in config.nml.")
        return

    # Find global colour limits
    all_G = []
    for _, s in snaps:
        all_G.extend(s["Gamma"].tolist())
    vmax = np.max(np.abs(all_G)) if all_G else 1.0

    cmap = cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Find x-limits
    all_x = []
    for _, s in snaps:
        all_x.extend(s["x"].tolist())
    x_max = max(max(all_x) if all_x else chord * 5, chord * 5)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.set_xlim(-0.5 * chord, x_max * 1.05)
    ax.set_ylim(-3 * chord, 3 * chord)
    ax.set_xlabel("x / chord")
    ax.set_ylabel("y / chord")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25)
    ax.plot([0, chord], [0, 0], "k-", lw=3)

    sc_bound = ax.scatter(
        [], [], c=[], cmap=cmap, norm=norm, s=60, marker="o", zorder=5
    )
    sc_wake = ax.scatter(
        [], [], c=[], cmap=cmap, norm=norm, s=15, marker="^", alpha=0.8, zorder=4
    )
    title_obj = ax.set_title("")
    plt.colorbar(sc_bound, ax=ax, label=r"$\Gamma$", shrink=0.8)
    plt.tight_layout()

    def _update(frame_idx):
        step, snap = snaps[frame_idx]
        mb = snap["type"] == 0
        mw = snap["type"] == 1
        # Update bound
        xy_b = (
            np.column_stack([snap["x"][mb], snap["y"][mb]])
            if mb.any()
            else np.empty((0, 2))
        )
        sc_bound.set_offsets(xy_b)
        sc_bound.set_array(snap["Gamma"][mb])
        # Update wake
        xy_w = (
            np.column_stack([snap["x"][mw], snap["y"][mw]])
            if mw.any()
            else np.empty((0, 2))
        )
        sc_wake.set_offsets(xy_w)
        sc_wake.set_array(snap["Gamma"][mw])
        title_obj.set_text(f"Wake pattern  step={step}  n_wake={mw.sum()}")
        return sc_bound, sc_wake, title_obj

    ani = FuncAnimation(fig, _update, frames=len(snaps), interval=interval, blit=False)

    if save_gif:
        gif_path = os.path.join(output_dir, "wake_animation.gif")
        writer = PillowWriter(fps=max(1, 1000 // interval))
        ani.save(gif_path, writer=writer)
        print(f"  Saved: {gif_path}")

    plt.show()
    return ani


# ─────────────────────────────────────────────────────────────
# 4. Optional: Plotly interactive wake animation
# ─────────────────────────────────────────────────────────────


def plotly_wake_animation(cfg, output_dir="output", snap_base="wake"):
    """
    Build an interactive HTML animation using Plotly.
    Each frame shows the bound + wake vortices coloured by Gamma.
    """
    if not HAS_PLOTLY:
        print("  plotly not installed.  Install with: pip install plotly")
        return

    chord = float(cfg.get("chord", 1.0))
    snaps = load_all_snapshots(output_dir, snap_base)
    if not snaps:
        print("  No snapshots for Plotly animation.")
        return

    all_G = np.concatenate([s["Gamma"] for _, s in snaps])
    vmax = np.max(np.abs(all_G))

    frames = []
    for step, snap in snaps:
        mb = snap["type"] == 0
        mw = snap["type"] == 1
        data = [
            go.Scatter(
                x=snap["x"][mb],
                y=snap["y"][mb],
                mode="markers",
                marker=dict(
                    color=snap["Gamma"][mb],
                    colorscale="RdBu",
                    cmin=-vmax,
                    cmax=vmax,
                    size=8,
                ),
                name="Bound Γ",
                showlegend=False,
            ),
            go.Scatter(
                x=snap["x"][mw],
                y=snap["y"][mw],
                mode="markers",
                marker=dict(
                    color=snap["Gamma"][mw],
                    colorscale="RdBu",
                    cmin=-vmax,
                    cmax=vmax,
                    size=5,
                    symbol="triangle-up",
                ),
                name="Wake Γ",
                showlegend=False,
            ),
        ]
        frames.append(go.Frame(data=data, name=str(step)))

    # Initial data (first frame)
    s0 = snaps[0][1]
    mb0, mw0 = (s0["type"] == 0), (s0["type"] == 1)
    fig = go.Figure(
        data=[
            go.Scatter(
                x=s0["x"][mb0],
                y=s0["y"][mb0],
                mode="markers",
                marker=dict(
                    color=s0["Gamma"][mb0],
                    colorscale="RdBu",
                    cmin=-vmax,
                    cmax=vmax,
                    size=8,
                    colorbar=dict(title="Γ"),
                ),
            ),
            go.Scatter(
                x=s0["x"][mw0],
                y=s0["y"][mw0],
                mode="markers",
                marker=dict(
                    color=s0["Gamma"][mw0],
                    colorscale="RdBu",
                    cmin=-vmax,
                    cmax=vmax,
                    size=5,
                    symbol="triangle-up",
                ),
            ),
            go.Scatter(
                x=[0, chord],
                y=[0, 0],
                mode="lines",
                line=dict(color="black", width=4),
                name="Plate",
            ),
        ],
        frames=frames,
        layout=go.Layout(
            title="Unsteady DVM – Wake Evolution",
            xaxis=dict(title="x / chord"),
            yaxis=dict(title="y / chord", scaleanchor="x", scaleratio=1),
            updatemenus=[
                dict(
                    type="buttons",
                    showactive=False,
                    buttons=[
                        dict(
                            label="▶ Play",
                            method="animate",
                            args=[
                                None,
                                dict(frame=dict(duration=80), fromcurrent=True),
                            ],
                        )
                    ],
                )
            ],
            sliders=[
                dict(
                    steps=[
                        dict(
                            method="animate",
                            args=[[str(step)], dict(mode="immediate")],
                            label=str(step),
                        )
                        for step, _ in snaps
                    ]
                )
            ],
        ),
    )
    html_path = os.path.join(output_dir, "wake_interactive.html")
    fig.write_html(html_path)
    print(f"  Saved: {html_path}")


# ─────────────────────────────────────────────────────────────
# 5. kh parametric sweep
# ─────────────────────────────────────────────────────────────


def run_sweep(
    kh_values, base_config="config.nml", executable="./udvm", output_base="output_sweep"
):
    """
    Run the Fortran solver for a range of kh values.
    For each value, writes a temporary config and collects steady-state Cl amplitude.

    kh = k * h0/c.  We fix h0/c = 0.1 and vary k to achieve the desired kh.
    Returns: list of (kh, Cl_amp, Ct_mean) tuples.
    """
    if not os.path.exists(executable):
        print(f"  ERROR: executable {executable} not found.  Run 'make' first.")
        return []

    cfg_base = parse_nml(base_config)
    h0_over_c = 0.1  # Fix amplitude, vary k
    results = []
    spc = int(cfg_base.get("steps_per_cycle", 100))
    n_cycles = int(cfg_base.get("n_cycles", 6))

    os.makedirs(output_base, exist_ok=True)

    for kh in kh_values:
        k = kh / h0_over_c  # k = kh / (h0/c)
        print(f"  Sweep: kh={kh:.4f}  k={k:.4f}")

        run_dir = os.path.join(output_base, f"kh_{kh:.4f}")
        os.makedirs(run_dir, exist_ok=True)

        # Write modified config
        tmp_cfg = os.path.join(run_dir, "config_tmp.nml")
        overrides = {
            "reduced_freq": k,
            "h0": h0_over_c,  # chord = 1 so h0/c = h0
            "output_dir": f"'{run_dir}'",
        }
        write_nml_override(base_config, tmp_cfg, overrides)

        # Run simulation
        result = subprocess.run([executable, tmp_cfg], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  WARNING: Simulation failed for kh={kh}\n{result.stderr}")
            continue

        # Load results
        force_file = str(cfg_base.get("force_file", "forces.dat"))
        try:
            forces = load_forces(run_dir, force_file)
        except SystemExit:
            print(f"  WARNING: No output for kh={kh}")
            continue

        # Steady-state window
        n_ss = 2 * spc
        Cl_ss = forces["Cl"][-n_ss:]
        Ct_ss = forces["Ct"][-n_ss:]
        Cl_amp = np.max(np.abs(Cl_ss))
        Ct_mean = np.mean(Ct_ss)  # Mean thrust (should be ~0 until Ct implemented)

        results.append((kh, Cl_amp, Ct_mean))
        print(f"    → Cl_amp = {Cl_amp:.5f},  Ct_mean = {Ct_mean:.5f}")

    return results


def plot_sweep_results(results, output_dir="output_sweep"):
    """
    Plot Cl amplitude and Ct vs kh from the sweep results.
    Also overlays the Theodorsen prediction for Cl.
    """
    from validate import theodorsen_Cl_amplitude_phase  # reuse from validate.py

    if not results:
        print("  No sweep results to plot.")
        return

    kh_arr = np.array([r[0] for r in results])
    Cl_arr = np.array([r[1] for r in results])
    Ct_arr = np.array([r[2] for r in results])
    h0_over_c = 0.1

    # Theodorsen Cl amplitude at same kh values
    k_arr = kh_arr / h0_over_c
    Cl_theo = np.array(
        [theodorsen_Cl_amplitude_phase(k, kh)[0] for k, kh in zip(k_arr, kh_arr)]
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(r"Parametric sweep vs $kh$", fontsize=14)

    axes[0].plot(kh_arr, Cl_arr, "o-", lw=2, color="tab:blue", label="DVM")
    axes[0].plot(kh_arr, Cl_theo, "s--", lw=1.5, color="tab:orange", label="Theodorsen")
    axes[0].set_xlabel(r"$kh = k \cdot h_0/c$")
    axes[0].set_ylabel(r"$C_{L,\max}$ (steady-state amplitude)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.35)
    axes[0].set_title(r"Lift coefficient amplitude vs $kh$")

    axes[1].plot(kh_arr, Ct_arr, "o-", lw=2, color="tab:red")
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_xlabel(r"$kh = k \cdot h_0/c$")
    axes[1].set_ylabel(r"$\overline{C_T}$ (mean thrust)")
    axes[1].grid(True, alpha=0.35)
    axes[1].set_title(r"Mean thrust coefficient vs $kh$  (TODO: LES)")

    plt.tight_layout()
    path = os.path.join(output_dir, "sweep_results.png")
    fig.savefig(path, dpi=150)
    print(f"  Saved: {path}")
    plt.show()


# ─────────────────────────────────────────────────────────────
# Main CLI
# ─────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Unsteady DVM – visualisation and post-processing"
    )
    parser.add_argument(
        "--mode",
        default="single",
        choices=["single", "animate", "sweep", "all"],
        help="Visualisation mode",
    )
    parser.add_argument("--config", default="config.nml")
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--snap_base", default="wake")
    parser.add_argument(
        "--kh_values",
        nargs="+",
        type=float,
        default=[0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50],
        help="kh values for sweep mode",
    )
    parser.add_argument(
        "--plotly", action="store_true", help="Also generate Plotly HTML animation"
    )
    parser.add_argument(
        "--no_save", action="store_true", help="Do not save any files, just display"
    )
    args = parser.parse_args()

    cfg = parse_nml(args.config)
    out = args.output_dir or str(cfg.get("output_dir", "output"))

    save = not args.no_save

    if args.mode in ("single", "all"):
        print("\n  ─── Force history ────────────────────────────")
        plot_forces(cfg, output_dir=out, save=save)

    if args.mode in ("animate", "all"):
        print("\n  ─── Wake animation ────────────────────────────")
        snap_base = str(cfg.get("wake_snap_base", args.snap_base))
        animate_wake(cfg, output_dir=out, snap_base=snap_base, save_gif=save)
        if args.plotly:
            plotly_wake_animation(cfg, output_dir=out, snap_base=snap_base)

    if args.mode in ("single", "all"):
        # Plot the last available wake snapshot
        snap_base = str(cfg.get("wake_snap_base", args.snap_base))
        snaps = load_all_snapshots(out, snap_base)
        if snaps:
            step, snap = snaps[-1]
            print(f"\n  ─── Wake pattern (last snapshot, step {step}) ────")
            plot_wake_snapshot(snap, cfg, step=step, output_dir=out, save=save)

    if args.mode in ("sweep", "all"):
        print("\n  ─── kh parametric sweep ───────────────────────")
        sweep_dir = out + "_sweep"
        results = run_sweep(
            args.kh_values, base_config=args.config, output_base=sweep_dir
        )
        if results:
            plot_sweep_results(results, output_dir=sweep_dir)


if __name__ == "__main__":
    main()
