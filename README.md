# Unsteady Vortex Solver (UDVM) 🌪️

A computational aerodynamics project implementing the **Unsteady Discrete Vortex Method (UDVM)** to simulate an airfoil undergoing plunging motion. This solver computes the unsteady wake evolution, calculates aerodynamic forces over time, and validates the results against Theodorsen's classical analytical theory.

![Wake Evolution](assets/wake_placeholder.gif)
> *Placeholder: Animation of the shed vortex wake rolling up behind the plunging plate.*

## 🎯 Project Overview
This project was developed as a term project for an Unsteady Aerodynamics course. The primary objective is to investigate the relationship between the non-dimensional plunge velocity ($kh$) and the resulting wake patterns and lift/thrust coefficients ($C_l$ / $C_t$).

The architecture is split into two parts:
1. **The Core Solver (Fortran 90):** Handles the heavy lifting—influence matrices, time-stepping, Kutta condition enforcement, and wake convection.
2. **The Post-Processor (Python):** Handles validation, data parsing, and high-quality visualizations using Matplotlib and Pandas.

## 🚀 Features
* **Modular Fortran Architecture:** Cleanly separated modules for Kinematics, Geometry, Wake, and Influence matrices.
* **Frozen Wake Assumption:** (Base implementation) Shed vortices convect at freestream velocity.
* **Force Computation:** Quasi-steady approximations with infrastructure ready for full Unsteady Bernoulli integration.
* **Analytical Validation:** Compares simulated $C_l$ against the Theodorsen function $C(k) = F(k) + iG(k)$.

## 📂 Repository Structure

├── Makefile                # Build instructions
├── src/
│   ├── parameters.f90      # Simulation parameters and constants
│   ├── geometry.f90        # Airfoil panel discretization (1/4 - 3/4 rule)
│   ├── kinematics.f90      # Plunge motion equations
│   ├── influence.f90       # Biot-Savart velocity calculations
│   ├── wake.f90            # Wake management and convection
│   ├── solver.f90          # Linear system assembly and solver
│   ├── forces.f90          # Cl and Ct calculations
│   ├── io.f90              # Data export to .dat files
│   └── main.f90            # Main execution loop
├── scripts/
│   ├── validate.py         # Python script for Theodorsen comparison
│   └── main_viz.py         # Visualization for wake and forces
└── README.md

## 🛠️ Building and Running

### Prerequisites
* **gfortran** (GNU Fortran compiler)
* **make**
* **Python 3.x** (with `numpy`, `matplotlib`, `scipy`, `pandas`)

### Execution Steps
1. **Compile the Solver:**
   make

2. **Run the Simulation:**
   ./udvm_solver
   
   *This will generate `output_forces.dat` and `output_wake.dat`.*

3. **Visualize the Results:**
   python scripts/main_viz.py

4. **Run Analytical Validation:**
   python scripts/validate.py

## 📊 Results & Validation

### Validation against Theodorsen's Theory
![Validation Plot](assets/validation_placeholder.png)
> *Placeholder: Plot showing Fortran output $C_l$ overlaying the analytical Theodorsen solution.*

### Aerodynamic Forces ($C_l$ & $C_t$)
![Forces Plot](assets/forces_placeholder.png)
> *Placeholder: Dual plot showing Lift and Thrust coefficients over multiple cycles.*

## 📝 Future Improvements (TODOs)
- [ ] Implement free wake roll-up model (relaxing the frozen wake assumption).
- [ ] Integrate full unsteady Bernoulli equation for pressure distribution.
- [ ] Add pitching kinematics and combined pitch-plunge motion.
- [ ] Implement cambered NACA airfoil geometry via source-vortex panel method.