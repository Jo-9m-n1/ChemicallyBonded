from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def userMolecule():
    chemName = list(input("Input your molecule (In this format:" \
    " NaOH or NH3): ").strip())
    i = 0
    chemDict = {}
    while len(chemName) > i:
        atom = chemName[i]
        num = ""
        i += 1

        if len(chemName) > i and chemName[i].islower():
            atom += chemName[i]
            i += 1

        while i < len(chemName) and chemName[i].isdigit():
            num += chemName[i]
            i += 1 
        
        if num:
            count = int(num)
        else:
            count = 1
        chemDict[atom] = count
    
    return chemDict

# Create water molecule from SMILES (need to learn the notation)
mol = Chem.MolFromSmiles("[N+](=O)([O-])[O-]")

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

print(userMolecule())