# BaoabMD

A lightweight Molecular Dynamics (MD) engine written in JavaScript with GPU acceleration via OpenCL/C++. 

The engine processes PDB (Protein Data Bank) files and outputs full atomic trajectories into XYZ format for visualization.

---

## Math & Physics

* **Integrator:** BAOAB Langevin thermostat
* **Long-Range / Short-Range Forces:** Particle Mesh Ewald (PME) & Lennard-Jones 12-6 potential
* **Constraints:** SHAKE / RATTLE algorithms for rigid hydrogen bonds
* **Minimization:** L-BFGS geometry optimization *(larger files may take over 1 minute)*
* **GPU Acceleration:** C++ native add-ons (`.gyp`) and OpenCL offloading put baseline runtime from ~24.7s down to ~3.3s.

---

## Prerequisites

* [Node.js](https://nodejs.org/) (v16+ recommended)
* *OpenCL dependencies are pre-packaged in `BaoabMD/CL`.*
* *Intructions are in the `BaoabMD/INSTRUCTIONS.txt` file.
---

## Quickstart

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ZanderChampion/BaoabMD
   cd BaoabMD