from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json
import prototype_lewis_structure

# (bonding_pairs, lone_pairs): geometry
VSEPR_TABLE = {
    (2, 0): "linear2",  # done
    (3, 0): "trigonal planar",  # done
    (2, 1): "bent1",  # done
    (4, 0): "tetrahedral",  # done
    (3, 1): "trigonal pyramidal",  # done
    (2, 2): "bent2",  # done
    (5, 0): "trigonal bipyramidal",  # the rest has to be done manually
    (4, 1): "seesaw",
    (3, 2): "T-shaped",
    (2, 3): "linear5",
    (6, 0): "octahedral",
    (5, 1): "square pyramidal",
    (4, 2): "square planar",
    (3, 3): "T-shaped",
    (2, 4): "linear6"
}

# added: simple colors for legend
ELEMENT_COLORS = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#90E050",
    "Cl": "#1FF01F",
    "Br": "#A62929",
    "I": "#940094",
    "S": "#FFFF30",
    "P": "#FF8000",
    "B": "#FFB5B5",
    "Si": "#F0C8A0",
    "Se": "#FFA100",
    "Xe": "#429EB0",
    "Be": "#C2FF00"
}

HIGH_DOMAIN_STRUCTURES = {
    "trigonal bipyramidal",  # the rest has to be done manually
    "seesaw",
    "T-shaped",
    "linear5",
    "octahedral",
    "square pyramidal",
    "square planar",
    "T-shaped",
    "linear6"
}

def add_legend_to_html(output_file, atom_symbols):
    items = ""
    for symbol in sorted(atom_symbols):
        color = ELEMENT_COLORS.get(symbol, "#BBBBBB")
        border = "border:2px solid #000000;"

        items += f'''
        <div>
            <span style="
                display:inline-block;
                width:12px;
                height:12px;
                border-radius:50%;
                background:{color};
                {border}
                margin-right:6px;
            "></span>{symbol}
        </div>
        '''

    legend_html = f'''
    <div style="
        position:absolute;
        top:5px;
        right:5px;
        background:#E5F2FF;
        padding:10px;
        border:2px solid #000000;
        border-radius:15px;
        z-index:9999;
        font-family:Trebuchet MS, sans-serif;
        font-size:15px;
        border: 0px;
        overflow: hidden;
    ">
        <b>Legend:</b>
        {items}
    </div>
    '''

    with open(output_file, "r", encoding="utf-8") as f:
        html = f.read()

    body_start = html.find("<body")
    if body_start == -1:
        return

    body_end = html.find(">", body_start)
    if body_end == -1:
        return

    html = html[:body_end + 1] + legend_html + html[body_end + 1:]

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html)


def findBonds(G):
    bonds = []

    for atom1, atom2, bond_data in G.edges(data=True):
        bond = {}
        stringAtom2 = ''.join(char for char in atom2 if not char.isdigit())
        bond["to"] = stringAtom2  # this includes the general atom e.g "O"
        bond["toString"] = atom2  # this includes the specific atom e.g "O1"
        bond["from"] = atom1
        bond["bond_order"] = bond_data["order"]
        bonds.append(bond)
    return bonds


def findVSEPR(G, central):
    bonding_pairs = len(G.edges())
    lone_pairs = G.nodes[central]["lone_pairs"]
    structureVSEPR = ""
    for key, values in VSEPR_TABLE.items():
        if (bonding_pairs, lone_pairs) == key:
            structureVSEPR = values

    if not structureVSEPR:
        structureVSEPR = "linear1"

    return structureVSEPR


def findFormalCharges(atom, G, bond_order, central):
    with open("static/periodic_table.json") as file:
        elements = json.load(file)

    print(atom, G, bond_order, central)

    if atom != central:
        atomDict = None
        for i in range(len(bond_order)):
            if atom == bond_order[i]['to']:
                atomDict = bond_order[i]
                break

        if atomDict is None:
            return ""

        atom = atomDict['to']
        specificAtom = atomDict["toString"]

        valenceElectrons = elements[atom]["valence_electrons"]
        bondingOrder = int(atomDict["bond_order"])
        lonePairs = int(G.nodes[specificAtom]["lone_pairs"]) * 2

        formalCharges = valenceElectrons - lonePairs - bondingOrder
    else:
        valenceElectrons = elements[central]["valence_electrons"]
        lonePairs = int(G.nodes[central]["lone_pairs"]) * 2

        bondingOrder = 0
        for i in range(len(bond_order)):
            if bond_order[i]["from"] == central:
                bondingOrder += int(bond_order[i]["bond_order"])

        formalCharges = valenceElectrons - lonePairs - bondingOrder
    print(valenceElectrons, lonePairs, bondingOrder)
    if formalCharges > 0:
        stringFormalCharge = "+" + str(formalCharges)
    elif formalCharges < 0:
        stringFormalCharge = str(formalCharges)
    else:
        stringFormalCharge = ""
    return stringFormalCharge.replace("1", "")


def findSMILES(structure3D, centralAtom, bond_order, containsHydrogen, formalChargeDict):
    """
    Convert the user molecule to SMILES format with the VSEPR structure, central atom,
    and bond order.
    """
    bond_symbol = {
        1: "",
        2: "=",
        3: "#",
    }

    for bond in bond_order:
        if bond["to"] == "H":
            bond["to"] = "[H]"

    print(structure3D, centralAtom, bond_order, containsHydrogen, formalChargeDict)

    if structure3D == "linear1":
        if formalChargeDict:
            sAtom = bond_order[0]["to"]
            smilesForm = (
                f"[{centralAtom}{formalChargeDict.get(centralAtom, '')}]"
                + bond_symbol[bond_order[0]["bond_order"]]
                + f"[{bond_order[0]['to']}{formalChargeDict.get(sAtom, '')}]"
            )
        else:
            stringCentralAtom = ''.join(char for char in centralAtom if not char.isdigit())
            smilesForm = (
                f"{stringCentralAtom}"
                + bond_symbol[bond_order[0]["bond_order"]]
                + f"{bond_order[0]['to']}"
            )

    elif structure3D in ("linear2", "bent1", "bent2"):
        leftForm = bond_order[0]["to"] + bond_symbol[bond_order[0]["bond_order"]]
        rightForm = bond_symbol[bond_order[1]["bond_order"]] + bond_order[1]["to"]

        smilesForm = f"{leftForm}[{centralAtom}]{rightForm}"

    elif structure3D in ("trigonal planar", "trigonal pyramidal", "trigonal bipyramidal",
                         "seesaw", "T-shaped", "linear5", "octahedral", "square pyramidal",
                         "square planar", "linear6"):
        smilesForm = f"{centralAtom}"
        for i in range(len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"

    elif structure3D == "tetrahedral":
        smilesForm = bond_order[0]['to'] + bond_symbol[bond_order[0]["bond_order"]] + f"[{centralAtom}]"
        for i in range(1, len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"

        # special case with HClO4
        if centralAtom == "Cl" and len(bond_order) == 5:
            smilesForm = "[OH][Cl](=O)(=O)=O"

    else:
        smilesForm = centralAtom

    return smilesForm, containsHydrogen


bond_order = [
    {'to': 'N', 'toString': 'N1', 'from': 'N', 'bond_order': 3},
]
"""
list of molecules that doesn't work with alex's code but works with mines:
SO3 (computes all O as bond order 1, instead of 2)
SO2 (computes all O as bond order 1, instead of 2, and wrong VSEPR format)
SeO2Cl2, SeO2F2 (computes O as bond order 1, instead of 2, and wrong VSEPR format)
SeCl2 (Se is very silly)

SeO2, XeO4 doesn't work properly for mines.
"""

def expanded_vsepr_structures(structure):
    return

def main(formula):
    molecule, central, networkxGraph = prototype_lewis_structure.main(formula)

    bonds = findBonds(networkxGraph)
    structureVSEPR = findVSEPR(networkxGraph, central)

    formalChargeDict = {}
    containsHydrogen = False
    for key in molecule:
        if key == "H":
            containsHydrogen = True
        formalCharge = findFormalCharges(key, networkxGraph, bonds, central)
        if formalCharge:
            formalChargeDict[key] = formalCharge

    smiles, containsHydrogen = findSMILES(structureVSEPR, central, bonds, containsHydrogen, formalChargeDict)
    print(smiles)
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
    view = py3Dmol.view(width=200, height=200)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    view.zoomTo()

    output_file = "static/html_visualisation.html"
    with open(output_file, "w", encoding="utf-8") as f:
        view.write_html(f, fullpage=True)

    # added: get atoms from molecule and add legend
    atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    add_legend_to_html(output_file, atom_symbols)