# 🎯 CC-Docking Analysis   

This repository is dedicated to the **analysis of molecular docking results**, focusing on **ligand-protein distances** and **scoring** to facilitate **predictive modeling**. 📊🔬  

## 🛠️ Prerequisites  
Before running molecular docking, **ligand** and **receptor** preparation is necessary:  

- **Ligand Preparation**:  
  - Requires a **list of compounds** in **SMILES format**.  
  - Additional data is necessary for output analysis.  
-  **Receptor Preparation**:  
  - Done using **AutoDockTools** following **traditional methodologies**.  

## Docking Execution  
- Docking was performed using **AutoDock Vina**   
- Output **analysis scripts** are tailored for **Vina-generated results**.  
- **⚠️ Docking execution scripts are NOT included**, as they depend on your specific **computing setup**.  

## 📂 Output Processing  
- The system generates an **`.out` file** containing docking results for all compounds.  
- Initial preprocessing splits this file into **individual `compound_name.pdbqt` files**.  

## 🔬 Analysis Workflow  
After preprocessing, the provided **sequence of scripts** allows:  
✅ **Extraction and treatment of docking results**.  
✅ **Evaluation of ligand-protein interactions**.  
✅ **Development of a predictive model** using **binary classification** (Active 🔵 / Inactive ⚪).  

 
