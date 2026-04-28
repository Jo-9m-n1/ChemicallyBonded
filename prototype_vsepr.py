from rdkit import Chem
from rdkit.Chem import AllChem, RWMol, Conformer, rdMolTransforms
import py3Dmol
import numpy as np
import matplotlib.pyplot as plt

import json
import math

import prototype_lewis_structure

with open("static/vsepr.json", "r") as file:
    data = json.load(file)

VSEPR_TABLE = {
    tuple(map(int, key.split(","))): value
    for key, value in data["VSEPR_TABLE"].items()
}

ELECTRON_PAIR_TABLE = {
    key: value
    for key, value in data["ELECTRON_PAIR_TABLE"].items()
}

ELEMENT_COLORS = data["ELEMENT_COLORS"]

HIGH_DOMAIN_STRUCTURES = set(data["HIGH_DOMAIN_STRUCTURES"])

def add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structure3D, 
                       electronGeometry, molecule, hybridization, lone_pairs):
    items = ""
    stringVSEPR = ''.join(char for char in structure3D if not char.isdigit())
    
    for symbol in sorted(atom_symbols):
        color = ELEMENT_COLORS.get(symbol, "#EB1388")
        if expandedVSEPR and symbol == "F":
            color = "#D9A420"
        
        items += f'''
        <div style="display:flex; align-items:center; margin-bottom:2px;">
            <span style="display:inline-block; width:9px; height:9px; border-radius:50%; 
                         background:{color}; border:0.5px solid #000; margin-right:6px;"></span>
            <span style="font-size:10px; font-family:sans-serif;">{symbol}</span>
        </div>
        '''
        
    if lone_pairs:
        items += f'''
        <div style="display:flex; align-items:center; margin-bottom:2px;">
            <span style="display:inline-block; width:9px; height:9px; border-radius:50%;
            background:blue; border:0.5px solid #000; margin-right:6px;"></span>
            <span style="font-size:10px; font-family:sans-serif;">Lone pairs</span>
        </div>
        '''

    legend_html = f'''
    <div style="
        position: absolute; 
        bottom: 10px; 
        left: 8px;
        width: 190px; 
        background: rgba(255, 255, 255, 0.9);
        padding: 10px;
        border-radius: 10px;
        box-shadow: 0 4px 10px rgba(0,0,0,0.1);
        z-index: 9999;
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif;
        box-sizing: border-box;
        line-height: 1;
    ">
        <div style="font-size: 13px; font-weight: bold; margin-bottom: 6px; color: #000;">{molecule}</div>
        
        <div style="font-size: 11px; color: #444; margin-bottom: 8px;">
            <b>E-Geom:</b> {electronGeometry}<br>
            <b>M-Geom:</b> {stringVSEPR}<br>
            <b>Hybrid:</b> {hybridization}
        </div>

        <div style="border-top: 0.5px solid #eee; padding-top: 8px;">
            {items}
        </div>
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
        if total_bonds == int(key):
            electronGeometry = values["geometry"]
            hybridization = values["hybridization"]

    if not structureVSEPR:
        structureVSEPR = "linear1"

    if not electronGeometry:
        structureVSEPR = "No electron geometry found"

    #showPieChart(bonding_pairs, lone_pairs)

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


def lonePairPositions(view, positions):
    view.addSphere({
        'center': positions,
        'radius': 0.3,
        'color': 'blue',
        'opacity': '0.90'
    })
    return 


def generatingLonePairs(view, structure, mol, central, lone_pairs):
    if mol.GetNumConformers() == 0:
        try:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except:
            print("Could not generate 3D coordinates for this molecule.")
            return

    conf = mol.GetConformer()

    c_pos = [] # positions of central atom to put the lone pairs

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        if symbol == central:
            c_pos.append(pos.x)
            c_pos.append(pos.y)
            c_pos.append(pos.z)
    
    angle = rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
    print(f"Angle: {angle:.2f}°")

    trigonal_bipyramidal_pos = [
        {'x': -0.5 * (math.sqrt(3) / 2.0), 'y': -0.5/2, 'z': 0.0},
        {'x': 0, 'y': 0.5, 'z': 0.0},
        {'x': 0.5 * (math.sqrt(3) / 2.0), 'y': -0.5/2, 'z': 0.0},
    ]

    octahedral_pos = [
        {'x': 0.0, 'y': 0.0, 'z': 0.5},
        {'x': 0.0, 'y': 0.0, 'z': -0.5},
        {'x': 0.0, 'y': -0.5, 'z': 0.0},
        {'x': 0.0, 'y': 0.5, 'z': 0.0},
    ]

    if structure == "Bent1":
        lonePairPositions(view, 
                          {'x': c_pos[0], 'y': c_pos[1] + 0.5, 'z': c_pos[2]})
    elif structure == "Trigonal pyramidal":
        lonePairPositions(view, 
                          {'x': c_pos[0], 'y': c_pos[1], 'z': c_pos[2] + 0.55})
    elif structure == "Bent2":
        for i in range(2):
            lonePairPositions(view, 
                              {'x': c_pos[0], 'y': c_pos[1] + 0.35, 'z': c_pos[2] + 0.35 * (-1)**i})
    elif structure in ("Seesaw", "T-shaped2", "Linear5"):
        for i in range(lone_pairs):
            lonePairPositions(view,
                              trigonal_bipyramidal_pos[i])    
    elif structure in ("Square pyramidal", "Square planar", "T-shaped3", "Linear6"):
        for i in range(lone_pairs):
            lonePairPositions(view,
                              octahedral_pos[i])
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
    # standard octahedral vertices: (+-d, 0, 0), (0, +-d, 0), (0, 0, +-d)
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
            octahedral[1],
            octahedral[2],
            octahedral[3],
        ]
    
    if structure == "Linear6":
        return [
            octahedral[0],
            octahedral[1],
            octahedral[2],
        ]
    return []


def main(formula):
    try:
        with open(f"static/generated/vsepr_{formula}.html", "r") as f:
            content = f.read()
            return
    except FileNotFoundError:
        pass
    molecule, central, networkxGraph = prototype_lewis_structure.main(formula)

    if central == None:
        return molecule
    
    try:

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
                error_message = "Failed to parse molecule to SMILES."
                return None, error_message
        else:
            mol = smiles
        
        # Convert molecule to MOL block format
        mol_block = Chem.MolToMolBlock(mol)

        # Send to Py3Dmol viewer
        view = py3Dmol.view(width=400, height=400)
        view.addModel(mol_block, "mol")
        view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})

        generatingLonePairs(view, structureVSEPR, mol, central, lone_pairs)

        view.zoomTo({"sphere": {"scale": 0.5}}, 2000, True)

        output_file = f"static/generated/vsepr_{formula}.html"
        with open(output_file, "w", encoding="utf-8") as f:
            view.write_html(f, fullpage=True)

        # added: get atoms from molecule and add legend
        atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
        add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structureVSEPR, 
                        electronVSEPR, formula, hybridization, lone_pairs)
        
    except:
        # James!!!! Change this later!!
        # When there is an error in your function, make sure to return None and then the error message!
        return None, "This is an error :)"


if __name__ == "__main__":
    main("CO2")