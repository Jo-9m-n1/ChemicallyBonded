from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json

# (bonding_pairs, lone_pairs): geometry
VSEPR_TABLE = {
    (2, 0): "linear",
    (3, 0): "trigonal planar",
    (2, 1): "bent",
    (4, 0): "tetrahedral",
    (3, 1): "trigonal pyramidal",
    (2, 2): "bent",
    (5, 0): "trigonal bipyramidal",
    (4, 1): "seesaw",
    (3, 2): "T-shaped",
    (2, 3): "linear",
    (6, 0): "octahedral",
    (5, 1): "square pyramidal",
    (4, 2): "square planar",
}

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

def findLonePairs(chemDict, centralAtom):
    file_path = "static/periodic_table.json"
    try:
        with open(file_path, "r") as file:
            data = json.load(file)
    except:
        print(f"Error: The file '{file_path}' was not found.")

    totalValence = 0

    for key, value in chemDict.items():
        atomElectro = data[key]["valence_electrons"]
        totalValence += (int(atomElectro) * value)
    """
    i just realized im doing alex's code
    """ 
    return totalValence

def findSMILES(structure3D, centralAtom, chemDict):
    """
    given: VSEPR strcture, centralAtom
    should return: SMILES format
    """
    containsHydrogen = False

    for key in chemDict: # going to change to accomodate for 2 or more surrounding atoms
        if centralAtom != key:
            surroundingAtom = key

    for key in chemDict:
        if key == "H":
            containsHydrogen = True

    if structure3D == "linear":
        smilesForm = f"{surroundingAtom}" + f"[{centralAtom}]" + f"{surroundingAtom}"
    
    if not surroundingAtom:
        smilesForm = centralAtom

    if containsHydrogen:
        smilesForm = f"[{centralAtom}" + f"{surroundingAtom}" + "2]"

    return smilesForm, containsHydrogen

molecule = userMolecule()
central = findCentral(molecule)
smiles, containsHydrogen = findSMILES("linear", central, molecule)

print(smiles, containsHydrogen)
# Create water molecule from SMILES (need to learn the notation)
mol = Chem.MolFromSmiles(smiles)

# Add explicit hydrogens
if containsHydrogen:
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