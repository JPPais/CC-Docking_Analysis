from pubchempy import get_compounds
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
from rdkit.Chem.rdMolDescriptors import CalcTPSA
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# List of compounds to convert to SMILES
compounds = [
    "Rifampin", "Isoniazid", "Pyrazinamide", "Ethambutol",
    "Kanamycin", "Streptomycin", "Capreomycin", "Amikacin", "Levofloxacin",
    "Moxifloxacin", "Gatifloxacin", "Bedaquiline", "Delamanid", "Linezolid",
    "Pretomanid", "Delpazolid", "Sutezolid", "SQ109", "PBTZ169", "Q203", "Clofazimine"
]

def get_smiles(name):
    """Fetches the SMILES representation of a compound by name from PubChem."""
    try:
        compounds = get_compounds(name, 'name')
        return compounds[0].canonical_smiles if compounds else None
    except Exception as e:
        print(f"Error fetching {name}: {e}")
        return None

def get_molecule_properties(name, smiles):
    """Computes molecular properties using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    properties = {
        "Name": name,
        "SMILES": smiles,
        "Molecular Weight": Descriptors.MolWt(mol),
        "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "TPSA": CalcTPSA(mol),
        "HB Donors": rdMolDescriptors.CalcNumHBD(mol),
        "HB Acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "LogP": Crippen.MolLogP(mol),
        "Heavy Atom Count": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "Fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "Heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "Formal Charge": Chem.GetFormalCharge(mol),
        "Refractivity": Crippen.MolMR(mol),
        "Lipinski Violations": sum([
            Descriptors.MolWt(mol) > 500,
            rdMolDescriptors.CalcNumHBD(mol) > 5,
            rdMolDescriptors.CalcNumHBA(mol) > 10,
            Crippen.MolLogP(mol) > 5
        ]), 
    }
    return properties

# Fetch SMILES and compute properties
data = []
for name in compounds:
    smiles = get_smiles(name)
    if smiles:
        props = get_molecule_properties(name, smiles)
        if props:
            data.append(props)

# Convert to DataFrame and save
properties_df = pd.DataFrame(data)
output_file = "molecular_properties.csv"
properties_df.to_csv(output_file, index=False)
print(f"Molecular properties saved to {output_file}")

# Visualization
sns.set(style="whitegrid")

# Save directory
output_dir = "output_graphics/"
os.makedirs(output_dir, exist_ok=True)

# Histograms for property distributions
properties_df.hist(figsize=(12, 8), bins=20)
plt.tight_layout()
plt.savefig(f"{output_dir}histograms.png")
plt.close()

# Correlation Heatmap (excluding non-numeric columns)
numeric_df = properties_df.select_dtypes(include=["number"])
plt.figure(figsize=(10, 8))
sns.heatmap(numeric_df.corr(), annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Heatmap of Molecular Properties")
plt.savefig(f"{output_dir}correlation_heatmap.png")
plt.close()

# Scatter Plots
plt.figure(figsize=(8, 6))
sns.scatterplot(x=properties_df["Molecular Weight"], y=properties_df["LogP"], 
                hue=properties_df["Lipinski Violations"], palette="coolwarm", s=100)
plt.xlabel("Molecular Weight")
plt.ylabel("LogP (Hydrophobicity)")
plt.title("Molecular Weight vs. LogP")
plt.savefig(f"{output_dir}molecular_weight_vs_logp.png")
plt.close()

plt.figure(figsize=(8, 6))
sns.scatterplot(x=properties_df["TPSA"], y=properties_df["HB Donors"], 
                hue=properties_df["Aromatic Rings"], palette="viridis", s=100)
plt.xlabel("Topological Polar Surface Area (TPSA)")
plt.ylabel("H-Bond Donors")
plt.title("TPSA vs. H-Bond Donors")
plt.savefig(f"{output_dir}tpsa_vs_hb_donors.png")
plt.close()

# Pairplot for multivariate relationships
pairplot = sns.pairplot(numeric_df, diag_kind="kde", corner=True)
pairplot.savefig(f"{output_dir}pairplot.png")
plt.close()

print(f"All plots saved in {output_dir}")
