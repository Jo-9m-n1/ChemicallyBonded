from rdkit import Chem
from rdkit.Chem import AllChem, RWMol, Conformer
import py3Dmol
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import json
import math

import prototype_lewis_structure
# (bonding_pairs, lone_pairs): geometry
VSEPR_TABLE = {
    (2, 0): "Linear2",  # done
    (3, 0): "Trigonal planar",  # done
    (2, 1): "Bent1",  # done
    (4, 0): "Tetrahedral",  # done
    (3, 1): "Trigonal pyramidal",  # done
    (2, 2): "Bent2",  # done
    (5, 0): "Trigonal bipyramidal",  # done
    (4, 1): "Seesaw", # done
    (3, 2): "T-shaped2", # done
    (2, 3): "Linear5", # done
    (6, 0): "Octahedral", # done
    (5, 1): "Square pyramidal", # done
    (4, 2): "Square planar", # done
    (3, 3): "T-shaped3", # done
    (2, 4): "Linear6" # done
}

ELECTRON_PAIR_TABLE = {
    2: {"geometry": "Linear", "hybridization": "sp"},
    3: {"geometry": "Trigonal planar", "hybridization": "sp2"},
    4: {"geometry": "Tetrahedral", "hybridization": "sp3"},
    5: {"geometry": "Trigonal bipyramidal", "hybridization": "sp3d"},
    6: {"geometry": "Octahedral", "hybridization": "sp3d2"}
}

# added: simple colors for legend
ELEMENT_COLORS = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#D9A420",
    "Cl": "#1FF01F",
    "Br": "#A62929",
    "I": "#940094",
    "S": "#fac431",
    "P": "#FF8000",
    "B": "#1FF01F",
    "Si": "#F0C8A0",
    "Xe": "#429EB0",
    "As": "#EB1388",
}

HIGH_DOMAIN_STRUCTURES = {
    "Trigonal bipyramidal",
    "Seesaw",
    "T-shaped2",
    "Linear5",
    "Octahedral",
    "Square pyramidal",
    "Square planar",
    "T-shaped3",
    "Linear6"
}

def add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structure3D, 
                       electronGeometry, molecule, hybridization, lone_pairs):
    items = ""
    stringVSEPR = ''.join(char for char in structure3D if not char.isdigit())
    for symbol in sorted(atom_symbols):
        color = ELEMENT_COLORS.get(symbol, "#EB1388")
        if expandedVSEPR:
            if symbol == "F":
                color = "#D9A420"
        border = "border:2px solid #000000;"

        items += f'''
        <div>
            <span style="
                display:inline-block;
                width:12px;
                height:12px;
                border-radius:100%;
                background:{color};
                {border}
                margin-right:6px;
            "></span>{symbol}
        </div>
        '''
    if lone_pairs:
        items += f'''
        <div>
            <span style="
                display:inline-block;
                width:12px;
                height:12px;
                border-radius:50%;
                background:blue;
                border:0.5px solid #000000;
                margin-right:6px;
            "></span>Lone pairs
        </div>
        '''

    legend_html = f'''
    <div style="
        position:absolute;
        width: 96%;
        background: rgba(112, 165, 225, 0.25);
        padding:6px;
        border:none;
        border-radius:25px;
        z-index:9999;
        display: flex;
        font-size:15px;
        gap: 2px;
    ">
        <img src="/static/molecule_chart.png" alt="Distribution Chart" style="height: 80px; border-radius: 10px;">

        <ul style="font-size: 8px; color: black;">
            <li><b>Molecule Name:</b> {molecule}</li>
            <li><b>Electron-Pair Geometry:</b> {electronGeometry}</li>
            <li><b>Molecular Geometry:</b> {stringVSEPR}</li>
            <li><b>Hybridization:</b> {hybridization}</li>
            <li style="display:flex; gap: 3px;"> {items}</li>
        </ul>
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
    total_bonds = bonding_pairs + lone_pairs

    structureVSEPR = ""
    electronGeometry = ""
    hybridization = ""

    for key, values in VSEPR_TABLE.items():
        if (bonding_pairs, lone_pairs) == key:
            structureVSEPR = values

    for key, values in ELECTRON_PAIR_TABLE.items():
        if total_bonds == key:
            electronGeometry = values["geometry"]
            hybridization = values["hybridization"]

    if not structureVSEPR:
        structureVSEPR = "linear1"

    if not electronGeometry:
        structureVSEPR = "No electron geometry found"

    showPieChart(bonding_pairs, lone_pairs)

    return structureVSEPR, electronGeometry, hybridization, lone_pairs


def showPieChart(bonding_pairs, lone_pairs):
    total = bonding_pairs + lone_pairs
    
    sizes = np.array([
        (bonding_pairs / total) * 100, 
        (lone_pairs / total) * 100
    ])

    labels = np.array(["Bonding Pairs", "Lone Pairs"])
    mask = sizes > 0

    plt.pie(sizes[mask], labels=labels[mask], autopct='%1.1f%%', shadow=True, 
            startangle=140, textprops={'fontsize': 14, 'weight': 'bold'})
    plt.axis("equal")
    plt.title("Bonding Distribution", fontsize=16, weight='bold', pad=20)
    plt.savefig("static/molecule_chart.png", bbox_inches='tight', transparent=True)
    return


def generatingLonePairs(view, structure, mol, central):
    if mol.GetNumConformers() == 0:
        try:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except:
            print("Could not generate 3D coordinates for this molecule.")
            return

    conf = mol.GetConformer()

    positions_central = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        if symbol == central:
            positions_central.append(pos.x)
            positions_central.append(pos.y)
            positions_central.append(pos.z)
    
    trigonal_bipyramidal_positions = [
        {'x': -0.5 * (math.sqrt(3) / 2.0), 'y': -0.5/2, 'z': 0.0},
        {'x': 0, 'y': 0.5, 'z': 0.0},
        {'x': 0.5 * (math.sqrt(3) / 2.0), 'y': -0.5/2, 'z': 0.0},
    ]

    if structure == "Bent1":
        view.addSphere({
            'center': {'x': positions_central[0], 'y': positions_central[1] + 0.5, 'z': positions_central[2]},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Trigonal pyramidal":
        view.addSphere({
            'center': {'x': positions_central[0], 'y': positions_central[1], 'z': positions_central[2] + 0.55},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Bent2":
        view.addSphere({
            'center': {'x': positions_central[0], 'y': positions_central[1] + 0.35, 'z': positions_central[2] + 0.35},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
        view.addSphere({
            'center': {'x': positions_central[0], 'y': positions_central[1] + 0.35, 'z': positions_central[2] - 0.35},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Seesaw":
        view.addSphere({
            'center': trigonal_bipyramidal_positions[0],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "T-shaped2":
        view.addSphere({
            'center': trigonal_bipyramidal_positions[0],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
        view.addSphere({
            'center': trigonal_bipyramidal_positions[1],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Linear5":
        view.addSphere({
            'center': trigonal_bipyramidal_positions[0],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
        view.addSphere({
            'center': trigonal_bipyramidal_positions[1],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
        view.addSphere({
            'center': trigonal_bipyramidal_positions[2],
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Square pyramidal":
        view.addSphere({
            'center': {'x': 0, 'y': 0, 'z': 0.5},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    elif structure == "Square Planar":
        view.addSphere({
            'center': {'x': 0, 'y': 0, 'z': 0.5},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
        view.addSphere({
            'center': {'x': 0, 'y': 0, 'z': -0.5},
            'radius': 0.3,
            'color': 'blue',
            'opacity': '0.90'
        })
    return


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

    expandedVSEPR = False

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

    elif structure3D in ("Linear2", "Bent1", "Bent2"):
        leftForm = bond_order[0]["to"] + bond_symbol[bond_order[0]["bond_order"]]
        rightForm = bond_symbol[bond_order[1]["bond_order"]] + bond_order[1]["to"]

        smilesForm = f"{leftForm}[{centralAtom}]{rightForm}"

    elif structure3D in ("Trigonal planar", "Trigonal pyramidal"):
        smilesForm = f"{centralAtom}"
        for i in range(len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"

    elif structure3D == "Tetrahedral":
        smilesForm = bond_order[0]['to'] + bond_symbol[bond_order[0]["bond_order"]] + f"[{centralAtom}]"
        for i in range(1, len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"

        # special case with HClO4
        if centralAtom == "Cl" and len(bond_order) == 5:
            smilesForm = "[OH][Cl](=O)(=O)=O"
    
    elif structure3D in (HIGH_DOMAIN_STRUCTURES):
        expandedVSEPR = True

        pt = Chem.GetPeriodicTable()

        central_num = pt.GetAtomicNumber(centralAtom)
        surroundings_num = pt.GetAtomicNumber(bond_order[0]['to'])

        mol = RWMol()
        s_idx = mol.AddAtom(Chem.Atom(central_num))

        f_indices = []
        for i in range(len(bond_order)):
            f_indices.append(mol.AddAtom(Chem.Atom(surroundings_num)))
            
        for f_idx in f_indices:
            mol.AddBond(s_idx, f_idx, Chem.BondType.SINGLE)

        smilesForm = mol.GetMol()

        coords = expanded_vsepr_structures(structure3D)

        conf = Conformer(smilesForm.GetNumAtoms())
        for i in range(smilesForm.GetNumAtoms()):
            conf.SetAtomPosition(i, coords[i])
        smilesForm.AddConformer(conf)

    else:
        smilesForm = centralAtom

    return smilesForm, containsHydrogen, expandedVSEPR


"""
list of molecules that doesn't work with alex's code but works with mines:
SeOCl2, SeOF2 (computes O as bond order 1, instead of 2, and wrong VSEPR format)
SeO2, XeO4 doesn't work properly for mines.
"""

def expanded_vsepr_structures(structure):
    d = 1.5
    sqrt3_over_2 = math.sqrt(3) / 2.0
    # Standard octahedral vertices: (+-d, 0, 0), (0, +-d, 0), (0, 0, +-d)
    octahedral = np.array([
        [0, 0, 0],
        [d, 0, 0], # right
        [-d, 0, 0], # left
        [0, d, 0], # forward
        [0, -d, 0], # back
        [0, 0, d], # top
        [0, 0, -d] # bottom
    ])

    trigonal_bipyramidal = np.array([
        [0, 0, 0],
        [d * sqrt3_over_2, -d/2, 0], # right
        [-d * sqrt3_over_2, -d/2, 0], # left
        [0, d, 0], # forward
        [0, 0, d], # top
        [0, 0, -d] # bottom
    ])

    if structure == "Trigonal bipyramidal":
        return trigonal_bipyramidal
    
    if structure == "Seesaw":
        return [
            trigonal_bipyramidal[0],
            trigonal_bipyramidal[1],
            trigonal_bipyramidal[3],
            trigonal_bipyramidal[4],
            trigonal_bipyramidal[5],
        ]
    
    if structure == "T-shaped2":
        return [
            trigonal_bipyramidal[0],
            trigonal_bipyramidal[1],
            trigonal_bipyramidal[4],
            trigonal_bipyramidal[5],
        ]
    
    if structure == "Linear5":
        return [
            trigonal_bipyramidal[0],
            trigonal_bipyramidal[4],
            trigonal_bipyramidal[5],
        ]
    
    if structure == "Octahedral":
        return octahedral
    
    if structure == "Square pyramidal":
        return [
            octahedral[0],
            octahedral[1],
            octahedral[2],
            octahedral[3],
            octahedral[4],
            octahedral[6],            
        ]
    
    if structure == "Square planar":
        return [
            octahedral[0],
            octahedral[1],
            octahedral[2],
            octahedral[3],
            octahedral[4],
        ]
    
    if structure == "T-shaped3":
        return [
            octahedral[0],
            octahedral[4],
            octahedral[5],
            octahedral[6],
        ]
    
    if structure == "Linear6":
        return [
            octahedral[0],
            octahedral[5],
            octahedral[6],
        ]
    return []

def main(formula):

    molecule, central, networkxGraph = prototype_lewis_structure.main(formula)

    if central == None:
        return molecule

    bonds = findBonds(networkxGraph)
    structureVSEPR, electronVSEPR, hybridization, lone_pairs = findVSEPR(networkxGraph, central)

    formalChargeDict = {}
    containsHydrogen = False
    for key in molecule:
        if key == "H":
            containsHydrogen = True
        formalCharge = findFormalCharges(key, networkxGraph, bonds, central)
        if formalCharge:
            formalChargeDict[key] = formalCharge

    smiles, containsHydrogen, expandedVSEPR = findSMILES(structureVSEPR, central, bonds, containsHydrogen, formalChargeDict)
    print(smiles, structureVSEPR)

    if not expandedVSEPR:
        mol = Chem.MolFromSmiles(smiles)

        # Add explicit hydrogens
        if containsHydrogen:
            mol = Chem.AddHs(mol)
        
        try:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except:
            print("Failed to parse molecule to SMILES.")
            return
    else:
        mol = smiles
    
    # Convert molecule to MOL block format
    mol_block = Chem.MolToMolBlock(mol)

    # Send to Py3Dmol viewer
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})

    generatingLonePairs(view, structureVSEPR, mol, central)

    view.zoomTo({"sphere": {"scale": 0.5}}, 2000, True)

    output_file = "static/html_visualisation.html"
    with open(output_file, "w", encoding="utf-8") as f:
        view.write_html(f, fullpage=True)

    # added: get atoms from molecule and add legend
    atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structureVSEPR, 
                       electronVSEPR, formula, hybridization, lone_pairs)