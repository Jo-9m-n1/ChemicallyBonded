from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json

def userMolecule():
    chemName = list(input("Input your molecule (In this format:" \
    " BeCl2 or H2O): ").strip())
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

def findCentral(chemDict):
    """
    Find the central atom in the lewis structure by finding which
    one is the least electronegative.

    Args:
        chemDict (dict): The dictionary that includes the amount
        of atoms in the molecule
    
    Returns:
        centralAtom
    """
    file_path = "static/periodic_table.json"
    try:
        with open(file_path, "r") as file:
            data = json.load(file)
    except:
        print(f"Error: The file '{file_path}' was not found.")
    
    electroDict = {}
    for key in chemDict:
        atomElectro = data[key]["electronegativity"]
        electroDict[key] = atomElectro

    centralAtom = min(electroDict, key=electroDict.get)
    
    return centralAtom
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
#view.show()
with open('templates/html_visualisation.html', 'w') as f:
    view.write_html(f, fullpage=True)

molecule = userMolecule()
print(findCentral(molecule))