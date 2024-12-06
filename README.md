# Brain Network Project

This repository is part of a collaborative project to study and simulate dynamics on brain networks using the **Brain Connectivity Toolbox (BCT)**. Our goal is to analyze network properties, simulate dynamics (e.g., Susceptible-Infected models), and visualize 3D brain networks based on coactivation data.

---

## Overview
Brain networks are modeled as graphs where:
- **Nodes** represent brain regions.
- **Edges** represent connections (e.g., coactivation or structural links).

We leverage tools from the **BCT** to process, analyze, and simulate dynamics on these networks. This project primarily uses MATLAB for implementation.

---

## Key Features
1. **3D Visualization**: Visualize brain networks in 3D using spatial coordinates of nodes.
2. **Dynamic Simulations**: Simulate processes like signal propagation and infection spread using graph-based models.
3. **Generative Models**: Generate synthetic brain networks to compare with real data.

---

## Important Files
### 1. **Data Files** (in `data_and_demos/`)
These files store the core data used for analysis and simulations:
- **`Coactivation_matrix.mat`**:
  - Contains the adjacency matrix representing coactivation between brain regions.
  - Includes `Coord`, a matrix of 3D spatial coordinates for nodes.

- **`GroupAverage_rsfMRI_matrix.mat`**:
  - Stores an averaged resting-state functional connectivity matrix.

- **`macaque47.mat` / `macaque71.mat`**:
  - Adjacency matrices for structural brain networks of macaque monkeys.

### 2. **Key Scripts**
These scripts are central to our project:
- **`Dynamic_Brain.m`**:
  - Simulates dynamic processes (e.g., infection spread) on the brain network.
  - Visualizes the network evolution in 3D over multiple time steps.

- **`evaluate_generative_model.m`**:
  - Generates synthetic networks using generative models (e.g., spatial or clustering-based).
  - Evaluates the similarity of synthetic networks to real networks.

- **`resource_efficiency_bin.m`**:
  - Calculates resource efficiency and shortest-path probabilities between nodes.

- **`demo_efficiency_measures.m`**:
  - Demonstrates how to compute efficiency measures on example data.

### 3. **Visualization Tools**
- **`Brain_Map.m`**:
  - Visualizes the brain network in 3D, highlighting specific nodes (e.g., initial states).

- **`adjacency_plot_und.m`**:
  - Plots adjacency matrices for undirected graphs.

---

## How to Run
### Prerequisites
- MATLAB installed on your system.
- Ensure the `data_and_demos` folder and BCT toolbox are included in MATLAB’s path:
  ```matlab
  addpath(genpath('/path/to/BCT/2019_03_03_BCT'));
