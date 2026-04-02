# Unsteady Vortex Solver (UDVM) 🌪️

A computational aerodynamics project implementing the **Unsteady Discrete Vortex Method (UDVM)** to simulate an airfoil undergoing plunging motion. This solver computes the unsteady wake evolution, calculates aerodynamic forces over time, and validates the results against Theodorsen's classical analytical theory.

![Wake Evolution](assets/wake_placeholder.gif)
> Placeholder: Animation of the shed vortex wake rolling up behind the plunging plate.

## 🎯 Project Overview
This project was developed as a term project for an Unsteady Aerodynamics course. The primary objective is to investigate the relationship between the non-dimensional plunge velocity ($kh$) and the resulting wake patterns and lift/thrust coefficients $(C_l$ / $C_t)$. 

In this simulation, $kh$ is defined as the non-dimensional amplitude scaled into the reduced frequency $(kh = k \cdot \frac{h_0}{c})$.

The architecture is split into two primary components:
1. **The Core Solver (Fortran 90):** Handles the heavy lifting—influence matrices via Biot-Savart, time-stepping, Kutta condition enforcement through Kelvin's circulation theorem, Gaussian elimination, and wake convection.
2. **The Post-Processor (Python):** Handles validation, data parsing, and high-quality visualizations using Matplotlib, Pandas, and Plotly.

## 🚀 Features
* **Modular Fortran Architecture:** Cleanly separated modules for Kinematics, Geometry, Wake, and Influence matrices to avoid namespace clutter.
* **Classical Flat Plate Geometry:** Automatic meshing using the standard 1/4 - 3/4 panel rule.
* **Wake Management:** Shed wake modeling utilizing a "Frozen Wake" assumption (wake convects horizontally at freestream velocity).
* **Force Computation:** Calculates Lift ($C_l$) using the unsteady Bernoulli equation (combining quasi-steady Kutta-Joukowski with added-mass terms) and quarter-chord moments ($C_m$).
* **Analytical Validation:** Dynamically computes Theodorsen's function $C(k) = F(k) + iG(k)$ utilizing Hankel/Bessel functions from `scipy.special` to benchmark the Fortran $C_l$ output.
* **Parametric Sweep Capability:** Automates the testing of varied $kh$ values to analyze steady-state lift amplitude behavior.
* **Interactive Visualization:** Capable of generating static force history graphs, Matplotlib `.gif` wake animations, and Plotly interactive HTML files.

## 📂 Repository Structure
``` plaintext
├── Makefile                # Target-mapped build instructions (parallel-make safe)
├── config.nml              # Main namelist configuration file (generated/required)
├── src/
│   ├── parameters.f90      # Simulation parameters and precision types
│   ├── geometry.f90        # Airfoil panel discretization
│   ├── kinematics.f90      # Motion equations (sinusoidal plunge/pitch)
│   ├── influence.f90       # Biot-Savart velocity calculations & AIC Matrix
│   ├── wake.f90            # Wake storage and advection handling
│   ├── solver.f90          # Gaussian linear system assembly and solver
│   ├── forces.f90          # Unsteady Bernoulli lift calculations
│   ├── io.f90              # Input config parsing and output dumping
│   └── main.f90            # Main execution loop
├── scripts/
│   ├── validate.py         # Python script for Theodorsen analytical comparison
│   └── main_viz.py         # Matplotlib/Plotly visualization hub
└── README.md
```

## 🛠️ Building and Running

### Prerequisites
* **gfortran** (GNU Fortran compiler)
* **GNU Make**
* **Python 3.x** (with `numpy`, `matplotlib`, `scipy`, `pandas`)
* *Optional:* `plotly` for interactive HTML animations.

### Execution Steps
1. **Compile the Solver:**
   ```bash
   make
   ```

2. **Run the Simulation:**
   ```bash
   make run
   # Or run directly: ./udvm config.nml
   ```
   *This populates the `output/` directory with `forces.dat` and `wake_NNNNNN.dat` snapshot files.*

3. **Visualize the Results:**
   Generate force histories, wake snapshots, and animations in one command:
   ```bash
   python scripts/main_viz.py --mode all
   ```

4. **Run Analytical Validation:**
   Compare the computed simulation directly against Theodorsen theory:
   ```bash
   python scripts/validate.py
   ```

## 📊 Results & Validation

### Validation against Theodorsen's Theory
![Validation Plot](assets/validation_placeholder.png)
> Placeholder: Dual plot showing time-history $C_l$ convergence overlaying the analytical Theodorsen solution, alongside amplitude attenuation $F(k)$ and $G(k)$ vs reduced frequency.

### Aerodynamic Forces ($C_l$ & $C_t$)
![Forces Plot](assets/forces_placeholder.png)
> Placeholder: Dual plot showing Lift and Thrust coefficients stabilizing over multiple plunging cycles.

## 📝 Future Improvements (TODOs)
- [ ] Implement free wake roll-up model (relaxing the frozen wake assumption via Biot-Savart self-induction).
- [ ] Implement full leading-edge suction calculations for proper thrust ($C_t$) extraction.
- [ ] Add active pitching kinematics and user-defined combined pitch-plunge phasing.
- [ ] Implement cambered NACA airfoil geometries via a combined source-vortex panel method.
