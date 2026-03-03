from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Create water molecule from SMILES (need to learn the notation)
mol = Chem.MolFromSmiles("C")

# Add explicit hydrogens
mol = Chem.AddHs(mol)
    
# Generate 3D coordinates
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

# Convert molecule to MOL block format
mol_block = Chem.MolToMolBlock(mol)

# Send to Py3Dmol viewer
view = py3Dmol.view(width=400, height=400)
view.addModel(mol_block, "mol")
view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
view.zoomTo()
view.show()
with open('templates/html_visualisation.html', 'w') as f:
    view.write_html(f, fullpage=True)
