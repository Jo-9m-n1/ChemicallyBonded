from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json

from prototype_lewis_structure import assign_lone_pairs

# (bonding_pairs, lone_pairs): geometry
VSEPR_TABLE = {
    (2, 0): "linear2",
    (3, 0): "trigonal planar",
    (2, 1): "bent",
    (4, 0): "tetrahedral",
    (3, 1): "trigonal pyramidal",
    (2, 2): "bent",
    (5, 0): "trigonal bipyramidal",
    (4, 1): "seesaw",
    (3, 2): "T-shaped",
    (2, 3): "linear5",
    (6, 0): "octahedral",
    (5, 1): "square pyramidal",
    (4, 2): "square planar",
    (3, 3): "T-shaped",
    (2, 4): "linear6"
}

def userMolecule():
    """
    Parse a chemical formula string into a dictionary containing
    each atom and its quantity in the molecule. (placeholder from alex's-feature)

    Args:
        formula (str): The chemical formula (e.g., "H2O", "CO2")

    Returns:
        dict: A dictionary mapping each atom to its count
    """
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

# TODO: need to be able to compute this from a function
#       this is a temp dictionary for now
bonds = [
    {"to": "O", "bond_order": 2},
    {"to": "O", "bond_order": 2}
]

def findSMILES(structure3D, centralAtom, chemDict, bond_order):
    """
    Convert the user molecule to SMILES format with the VSEPR structure, central atom,
    and bond order.
    
    Args:
        VSEPR strcture (string)
        centralAtom (string)
        chemDict (dict)
        bond_order (list of dict)

    Returns:
        smilesForm (string)
        containsHydrogen (bool)

    TODO: going to change to accommodate for 2 or more surrounding atoms
    """

    containsHydrogen = False
    bond_symbol = {
        1: "",
        2: "=",
        3: "#"
    }

    for key in chemDict:
        if key == "H":
            containsHydrogen = True
    print(structure3D, centralAtom, chemDict, bond_order)

    if structure3D == "linear2":
        for bond in bonds:
            if bond["to"] == "H":
                bond["to"] = "[H]"

        leftForm = bond_order[0]["to"] + bond_symbol[bond_order[0]["bond_order"]]
        rightForm = bond_symbol[bond_order[1]["bond_order"]] + bond_order[1]["to"]
        
        smilesForm = f"{leftForm}[{centralAtom}]{rightForm}"
    elif structure3D == "trigonal planar":
        smilesForm = f"{centralAtom}({bond_order[0]['to']})({bond_order[1]['to']}){bond_order[2]['to']}"


    return smilesForm, containsHydrogen

molecule = userMolecule()
central = findCentral(molecule)
#bonds, centralLonePairs, structure3D =  buildBondListAndVSEPR(molecule, central)
smiles, containsHydrogen = findSMILES("linear2", central, molecule, bonds)

# print(smiles, containsHydrogen)
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