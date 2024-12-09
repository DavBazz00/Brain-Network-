# Brain Network Project

This repository include an ongoing project for the exam in Biomedicine at the Politecnico di Torino. 

This is a collaborative project to study and simulate dynamics on brain networks using the **Brain Connectivity Toolbox (BCT)** and custom implementations. Our goal is to analyze network properties, simulate dynamics (e.g., Susceptible-Infected models), and visualize 3D brain networks based on coactivation data.

---

## Table of Contents
1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Important Files](#important-files)
    - [Data Files](#data-files-in-data_and_demos)
    - [Predefined Functions from the Brain Connectivity Toolbox](#predefined-functions-from-the-brain-connectivity-toolbox)
    - [Custom Scripts (Created by Us)](#custom-scripts-created-by-us)
4. [How to Run](#how-to-run)

---

## Overview
Brain networks are modeled as graphs where:
- **Nodes** represent brain regions.
- **Edges** represent connections (e.g., coactivation or structural links).

We use tools from the **BCT** alongside custom MATLAB scripts to process, analyze, and simulate dynamics on these networks. 

We also suggest giving a quick read at this [repository](https://github.com/brain-networks/PSY-P457) on how to use the **BCT toolbox**.

---

## Key Features
1. **3D Visualization**: Visualize brain networks in 3D using spatial coordinates of nodes.
2. **Dynamic Simulations**: Simulate processes like signal propagation and infection spread using graph-based models.
3. **Custom Implementations**: Expand the BCT with tailored scripts for visualization and dynamics.

---

## Important Files

### Data Files (in `data_and_demos/`)
These files store the core data used for analysis and simulations:
- **`Coactivation_matrix.mat`**:
  - Contains the adjacency matrix representing coactivation between brain regions.
  - Includes `Coord`, a matrix of 3D spatial coordinates for nodes.

---

### Predefined Functions from the Brain Connectivity Toolbox

The Brain Connectivity Toolbox (BCT) offers a wide array of predefined functions to analyze brain networks. Below is an overview of some important functions used in this project. Each description starts with a formal explanation, followed by a simplified, informal explanation.

#### 1. Visualization in 3-D and Preprocess
- **`adjacency_plot_und`**:
  - It takes:
    -  *adjacency matrix AIJ*
    -  *node spatial coordinates COOR*
    generates three vectors that can be used for quickly plotting the edges

- **`threshold_absolute`**:
  -  This function thresholds the connectivity matrix by absolute weight magnitude;
  -  All weights below the given threshold are set to 0.

- **`threshold_proportional`**:
  - This function "thresholds" the connectivity matrix by preserving a proportion p (0<p<1) of the strongest;
  - All other weights are set to 0.


#### 2. Community and Clustering Analysis
- **`modularity_und`**:
  -  It makes a subdivision of the network into nonoverlapping groups of nodes


---

### Custom Scripts (Created by Us)
These scripts were developed by the team to extend and apply BCT functionality. **These are the files to modify and implement as needed:**
- **`Dynamic_Brain.m`**:
  - Simulates dynamic processes (e.g., infection spread) on the brain network.
  - Visualizes the network evolution in 3D over multiple time steps.
  - Implements all the other scripts of the files below.

- **`Brain_Map.m`**:
  - Visualizes the brain network in 3D, highlighting specific nodes (e.g., initial states).
  - Useful for understanding the structure and initial setup of the network.

- **`Cluster_modularities.m`**:
   - Groups nodes in communities.
   - Uses an implemented function in the *Brain Connectivity Toolbox* such as the `modularity_und` function.
     
- **`Visualization_Brain.m`**:
  - Visualizes the connectivity matrix using the `imagesc` command. This command shows the connectivity matrix as a grid.
  - The strongest connections are assigned warm colors while non-existent or weak connections are assigned cooler colors.
 

---

## How to Run

### Prerequisites
- MATLAB installed on your system.
- Ensure the `data_and_demos` folder and BCT toolbox are included in MATLAB’s path:
  ```matlab
  addpath(genpath('/path/to/BCT/2019_03_03_BCT'));
