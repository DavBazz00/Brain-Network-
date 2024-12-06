# Brain Network Project

This repository is part of a collaborative project to study and simulate dynamics on brain networks using the **Brain Connectivity Toolbox (BCT)** and custom implementations. Our goal is to analyze network properties, simulate dynamics (e.g., Susceptible-Infected models), and visualize 3D brain networks based on coactivation data.

---

## Overview
Brain networks are modeled as graphs where:
- **Nodes** represent brain regions.
- **Edges** represent connections (e.g., coactivation or structural links).

We use tools from the **BCT** alongside custom MATLAB scripts to process, analyze, and simulate dynamics on these networks.

---

## Key Features
1. **3D Visualization**: Visualize brain networks in 3D using spatial coordinates of nodes.
2. **Dynamic Simulations**: Simulate processes like signal propagation and infection spread using graph-based models.
3. **Custom Implementations**: Expand the BCT with tailored scripts for visualization and dynamics.

---

## Important Files
### **Data Files** (in `data_and_demos/`)
These files store the core data used for analysis and simulations:
- **`Coactivation_matrix.mat`**:
  - Contains the adjacency matrix representing coactivation between brain regions.
  - Includes `Coord`, a matrix of 3D spatial coordinates for nodes.

---

### **Brain Connectivity Toolbox Scripts** (Original Files)
These files are part of the BCT and provide the foundational functions for our analysis:
- **`adjacency_plot_und.m`**:
  - Plots adjacency matrices for undirected graphs.
  - 
- **`evaluate_generative_model.m`**:
  - Generates synthetic networks using generative models (e.g., spatial or clustering-based).
  - Evaluates the similarity of synthetic networks to real networks.

- **`resource_efficiency_bin.m`**:
  - Calculates resource efficiency and shortest-path probabilities between nodes.

---

### **Custom Scripts** (Created by Us)
These scripts were developed by the team to extend and apply BCT functionality. **These are the files to modify and implement as needed:**
- **`Dynamic_Brain.m`**:
  - Simulates dynamic processes (e.g., infection spread) on the brain network.
  - Visualizes the network evolution in 3D over multiple time steps.
  - Implements Susceptible-Infected (SI) dynamics with probabilistic transitions.

- **`Brain_Map.m`**:
  - Visualizes the brain network in 3D, highlighting specific nodes (e.g., initial states).
  - Useful for understanding the structure and initial setup of the network.

- **`Visualization_Brain.m`**:
  - Visualizes the connectivity matrix using the `imagesc` command. This command shows the connectivity matrix as a  grid.
  - The strongest connections are assigned warm colors while non-existent or weak connections are assigned cooler colors.
 
- **`Community.m`**:
   - Groups nodes in communities.
   - Uses an implemented funtction in the *Brain Connectivity Toolbox* as the `community_louvain` function.
   - The final community assignments for each node is made after averaging over multiple runs since the `community_louvain` function is stochastic.

---

## How to Run
### Prerequisites
- MATLAB installed on your system.
- Ensure the `data_and_demos` folder and BCT toolbox are included in MATLAB’s path:
  ```matlab
  addpath(genpath('/path/to/BCT/2019_03_03_BCT'));
