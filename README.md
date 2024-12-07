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

### **Important Predefined Functions**

The Brain Connectivity Toolbox provides a variety of powerful functions for analyzing brain networks. Below are some important predefined functions and their roles:

#### **1. Community and Clustering Analysis**
- **`agreement.m`**:
  - **Purpose**: Computes an agreement matrix from multiple partitions, showing how often pairs of nodes are assigned to the same community.
  - **Usage**:
    ```matlab
    D = agreement(Ci);
    ```
    - `Ci`: A matrix where each column represents a community partition.
    - `D`: The resulting agreement matrix.

- **`clique_communities.m`**:
  - **Purpose**: Detects overlapping community structures in a binary undirected network using the clique percolation method.
  - **Usage**:
    ```matlab
    M = clique_communities(A, cq_thr);
    ```
    - `A`: Binary adjacency matrix.
    - `cq_thr`: Clique size threshold.

- **`link_communities.m`**:
  - **Purpose**: Identifies overlapping communities using hierarchical clustering of network links.
  - **Usage**:
    ```matlab
    M = link_communities(W, 'complete');
    ```
    - `W`: Directed (binary or weighted) connection matrix.
    - `'complete'`: Optional clustering type (default is `'single'`).

---

#### **2. Network Backbone and Synthetic Networks**
- **`backbone_wu.m`**:
  - **Purpose**: Extracts the backbone of a weighted undirected network using a minimum-spanning-tree algorithm.
  - **Usage**:
    ```matlab
    [CIJtree, CIJclus] = backbone_wu(CIJ, avgdeg);
    ```
    - `CIJ`: Weighted undirected connection matrix.
    - `avgdeg`: Desired average degree of the backbone.

- **`evaluate_generative_model.m`**:
  - **Purpose**: Generates synthetic networks and evaluates their energy functions using specified generative rules.
  - **Usage**:
    ```matlab
    [B, E, K] = evaluate_generative_model(A, Atgt, D, 'sptl', 'powerlaw', params);
    ```
    - `A`: Binary seed connections.
    - `D`: Euclidean distance matrix.
    - `'sptl'`: Generative model type (e.g., spatial model).
    - `'powerlaw'`: Generative model variant.
    - `params`: Model parameters.

- **`generative_model.m`**:
  - **Purpose**: Runs generative model simulations to create synthetic networks.
  - **Usage**:
    ```matlab
    B = generative_model(A, D, m, 'neighbors', {'powerlaw'}, params);
    ```

---

#### **3. Centrality Measures**
- **`pagerank_centrality.m`**:
  - **Purpose**: Computes the PageRank centrality, a variant of eigenvector centrality, to rank nodes in the network.
  - **Usage**:
    ```matlab
    r = pagerank_centrality(A, 0.85, falff);
    ```
    - `A`: Adjacency matrix.
    - `0.85`: Damping factor (typical value).
    - `falff`: Initial PageRank probabilities.

- **`subgraph_centrality.m`**:
  - **Purpose**: Calculates the subgraph centrality, a measure of the weighted sum of closed walks for each node.
  - **Usage**:
    ```matlab
    Cs = subgraph_centrality(CIJ);
    ```

---

#### **4. Functional and Structural Network Analysis**
- **`generate_fc.m`**:
  - **Purpose**: Generates synthetic functional connectivity (FC) matrices using structural connectivity (SC) as predictors.
  - **Usage**:
    ```matlab
    [FCpre, pred_data, Fcorr] = generate_fc(SC, beta, ED, {'SPLwei_log', 'SIwei_log'}, FC);
    ```
    - `SC`: Structural connectivity matrix.
    - `beta`: Regression coefficients.
    - `ED`: Euclidean distance matrix (optional).

- **`get_components.m`**:
  - **Purpose**: Identifies connected components of a binary undirected network.
  - **Usage**:
    ```matlab
    [comps, comp_sizes] = get_components(adj);
    ```
    - `adj`: Binary undirected adjacency matrix.

---

### **Quick Reference**
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
