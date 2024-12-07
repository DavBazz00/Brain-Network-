# Brain Network Project

This repository include a project for the exam in Biomedicine at the Politecnico di Torino. 

This is a collaborative project to study and simulate dynamics on brain networks using the **Brain Connectivity Toolbox (BCT)** and custom implementations. Our goal is to analyze network properties, simulate dynamics (e.g., Susceptible-Infected models), and visualize 3D brain networks based on coactivation data.

---

## Table of Contents
1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Important Files](#important-files)
    - [Data Files](#data-files-in-data_and_demos)
    - [Predefined Functions from the Brain Connectivity Toolbox](#predefined-functions-from-the-brain-connectivity-toolbox)
    - [Quick Reference](#quick-reference)
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

#### 1. Community and Clustering Analysis
- **`agreement.m`**: 
    - This function calculates the agreement matrix from multiple community partitions. The agreement matrix indicates how frequently pairs of nodes are assigned to the same community across different partitions.
  - It counts how often two nodes belong to the same group when we run community detection multiple times.

- **`clique_communities.m`**:
  - It detects overlapping communities in binary undirected networks using the clique percolation method. It finds groups of nodes that form tightly connected subgraphs (cliques).
  - This function finds small, super-tight groups of nodes that overlap, like cliques at a party.

- **`link_communities.m`**:
  - It identifies overlapping communities by clustering links (edges) instead of nodes. It is generalized for directed and weighted networks.
  - Instead of grouping people (nodes), this function groups relationships (edges) between them.


#### 2. Network Backbone and Synthetic Networks
- **`backbone_wu.m`**:
  - It extracts the backbone of a weighted undirected network using a minimum-spanning-tree algorithm. This reduces the network to its most important connections while preserving its structure.
  - It trims the network down to just the essential connections, making it easier to see the big picture.

- **`evaluate_generative_model.m`**:
  - It generates synthetic networks and evaluates their energy function using specified generative rules, such as spatial models or clustering-based models.
  - This function builds fake networks based on rules we choose and tells us how "good" they are compared to the real network.

- **`generative_model.m`**: 
  - It runs simulations to create synthetic networks using different generative models, such as spatial proximity or node similarity.
  - It builds fake networks based on specific recipes, like connecting nearby nodes or similar ones.


#### 3. Centrality Measures
- **`pagerank_centrality.m`**: 
  - It calculates the PageRank centrality, which measures a node's influence based on how often it is visited in a random walk with restarts.
  - Let's think of this as how popular a node is, like Google's PageRank algorithm for websites.

- **`subgraph_centrality.m`**: 
  - It computes the subgraph centrality for each node, which is a weighted sum of closed walks of various lengths starting and ending at that node.
  - It tells us how well-connected a node is to itself through loops.

#### 4. Functional and Structural Network Analysis
- **`generate_fc.m`**: 
  - It creates synthetic functional connectivity matrices based on structural connectivity data and network measures like shortest paths or search information.
  - It predicts how nodes communicate based on their physical connections.

- **`get_components.m`**: 
  - It identifies connected components in an undirected graph. Each component consists of nodes that are directly or indirectly connected.
  - It groups isolated islands of nodes in the network.

---

### Quick Reference
| **Function**             | **Purpose**                             | **Key Input**         | **Key Output**            |
|---------------------------|-----------------------------------------|-----------------------|---------------------------|
| `agreement.m`            | Compute agreement matrix               | Partitions (`Ci`)     | Agreement matrix (`D`)   |
| `clique_communities.m`   | Detect overlapping communities          | Binary adjacency (`A`)| Community affiliation (`M`) |
| `backbone_wu.m`          | Extract network backbone                | Weighted matrix (`CIJ`)| Backbone matrices         |
| `pagerank_centrality.m`  | Compute PageRank centrality             | Adjacency (`A`)       | PageRank vector (`r`)    |
| `subgraph_centrality.m`  | Compute subgraph centrality             | Adjacency (`CIJ`)     | Centrality vector (`Cs`) |
| `generate_fc.m`          | Generate synthetic FC matrices          | SC matrix, predictors | Predicted FC matrix       |
| `get_components.m`       | Identify connected components           | Binary matrix (`adj`) | Components and sizes      |

---

### Custom Scripts (Created by Us)
These scripts were developed by the team to extend and apply BCT functionality. **These are the files to modify and implement as needed:**
- **`Dynamic_Brain.m`**:
  - Simulates dynamic processes (e.g., infection spread) on the brain network.
  - Visualizes the network evolution in 3D over multiple time steps.
  - Implements Susceptible-Infected (SI) dynamics with probabilistic transitions.

- **`Brain_Map.m`**:
  - Visualizes the brain network in 3D, highlighting specific nodes (e.g., initial states).
  - Useful for understanding the structure and initial setup of the network.

- **`Visualization_Brain.m`**:
  - Visualizes the connectivity matrix using the `imagesc` command. This command shows the connectivity matrix as a grid.
  - The strongest connections are assigned warm colors while non-existent or weak connections are assigned cooler colors.
 
- **`Community.m`**:
   - Groups nodes in communities.
   - Uses an implemented function in the *Brain Connectivity Toolbox* such as the `community_louvain` function.
   - The final community assignments for each node are made after averaging over multiple runs since the `community_louvain` function is stochastic.

---

## How to Run

### Prerequisites
- MATLAB installed on your system.
- Ensure the `data_and_demos` folder and BCT toolbox are included in MATLAB’s path:
  ```matlab
  addpath(genpath('/path/to/BCT/2019_03_03_BCT'));
