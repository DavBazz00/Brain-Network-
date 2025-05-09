# Misfolded Protein Propagation on an 82-Node Human Brain Network

This repository implements two complementary computational models of misfolded-protein spread on a standardized 82-region human connectome. Explore healthy, aging, and treatment scenarios to study Alzheimer-type proteinopathy propagation dynamics.

---

## Project Structure

### **Data/**  
- **Edge.csv** and **Node.csv**  
  - Define the 82-node network: `Edge.csv` contains connection weights (adjacency matrix), and `Node.csv` lists region labels.  
  - Used by all simulation scripts to drive propagation along the connectome.

- **Data for Brain Network/**  
  - Contains 3D mesh files (`.nv`, `.in`, masks, etc.) and region definitions for the BrainNet Toolbox.  
  - These files let you render brain surfaces and overlay simulation data.  
  - Learn more and download the toolbox [**BrainNet Viewer**](https://www.nitrc.org/projects/bnv/)

---

### **Project 82 nodes/**  
Contains MATLAB functions for simulating and plotting each model.  

#### **fisher_kolmogorov/**  
- **FK_propagation.m**  
  Runs the baseline Fisher–Kolmogorov reaction–diffusion model (diffusion + local growth) on the 82-node connectome.  

- **FK_propagation_dynamic.m**  
  Simulates an aging brain by stochastically weakening connection strengths over time.  

- **FK_propagation_treatment.m**  
  Models a treatment scenario: random connection weakening **and** reduced local misfolded-protein growth rate.  

- **FK_propagation_treatment_aging.m**  
  Combines aging and treatment.

#### **heterodimer/**  
- **HeterodimerInfection.m**  
  Implements the heterodimer interaction model between misfolded and native proteins on the network.  

- **HeterodimerInfection_dynamic.m**  
  Applies aging-related stochastic connection weakening to the heterodimer model.  

- **HeterodimerInfection_treatment.m**  
  Applies treatment: connection weakening **and** reduced local growth in the heterodimer framework.  

- **HeterodimerInfection_treatment_aging.m**  
  Combines aging and treatment.

---

### **Images/**  
There are the most important plots and result in a `.png` format.

---

## Getting Started

If you’re new to GitHub, follow these steps to download and run the project:

1. **Install prerequisites**  
   - **MATLAB** (R2018b or later recommended)  
   - **BrainNet Toolbox**  
     1. Download from [https://www.nitrc.org/projects/bnv/](https://www.nitrc.org/projects/bnv/). 
     2. Add the downloaded toolbox folder to your MATLAB path (e.g., using `addpath`).

2. **Clone the repository**  
   - Open a terminal or Git Bash.  
   - Run:
     ```bash
     git clone https://github.com/DavBazz00/Brain-Network-.git
     ```
   - Or, navigate to the repository page (https://github.com/DavBazz00/Brain-Network-/tree/main) in your web browser and click the green **Code** button, then **Download ZIP**.

---



Happy exploring! Feel free to open an issue on GitHub if you encounter problems or have suggestions.  
