# Misfolded Protein Propagation on 82‑Node Brain Network

Project illustrating the spread of misfolded proteins on a standardized 82‑region human brain network using two computational models.

## Folder Overview

- **Data/**
  - Contains raw and processed connectivity/region data for running simulations and input to the BrainNet Toolbox.

- **Project 82 nodes/**
  - **fisher_kolmogorov/**: Scripts implementing the Fisher–Kolmogorov reaction‑diffusion model for aggregate growth and spread.
  - **heterodimer/**: Code for heterodimer interaction dynamics between misfolded and native proteins.

- **Images/**
  - **results/**: Output plots and figures from simulations.
  - **toolbox/**: Example network visualizations generated with BrainNet Toolbox.

---

## What This Project Does

1. **Load & preprocess** structural connectivity and region definitions.
2. **Simulate** protein propagation with two models:
   - *Fisher–Kolmogorov*: Diffusion + local growth over the connectome.
   - *Heterodimer*: Interaction and seeding dynamics.
3. **Visualize** results on 3D brain meshes via BrainNet Toolbox screenshots and generated figures.

Explore, modify parameters, and extend analyses to study patterns of proteinopathy spread in neurodegenerative disease contexts.

