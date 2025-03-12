# Brain Network Project

This repository include an ongoing project for the exam in Biomedicine at the Politecnico di Torino. 

This is a collaborative project to study and simulate dynamics on brain networks using the **Brain Connectivity Toolbox (BCT)** and custom implementations. Our goal is to analyze network properties, simulate dynamics (e.g., Susceptible-Infected models), and visualize 3D brain networks based on coactivation data.

---

## Table of Contents
1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Predefined Functions from the Brain Connectivity Toolbox (BCT)](#predefined-functions-from-the-brain-connectivity-toolbox)
4. [Custom Scripts (Created by Us)](#custom-scripts-created-by-us)
5. [How to Run](#how-to-run)

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

---
### Predefined Functions from the Brain Connectivity Toolbox (BCT)

The Brain Connectivity Toolbox (BCT) offers a wide array of predefined functions to analyze brain networks. Below is an overview of some important functions used in this project. Each description starts with a formal explanation, followed by a simplified, informal explanation.

#### 1. Visualization in 3-D and Preprocess
- **`adjacency_plot_und`**:
  - It takes:
    -  *adjacency matrix AIJ*
    -  *node spatial coordinates COOR*
  - It generates three vectors that can be used for quickly plotting the edges

- **`threshold_absolute`**:
  -  This function thresholds the connectivity matrix by absolute weight magnitude;
  -  All weights below the given threshold are set to 0.

---

### Custom Scripts (Created by Us)
These scripts were developed by the team to extend and apply BCT functionality. **These are the files to modify and implement as needed:**
- **`Project.m`**:
  - Simulates dynamic processes (e.g., infection spread) on the brain network.
  - Visualizes the network evolution in 3D over multiple time steps.
  - Implements all the other scripts of the files below.


