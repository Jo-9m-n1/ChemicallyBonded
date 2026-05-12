from rdkit import Chem
from rdkit.Chem import AllChem, RWMol, Conformer, rdMolTransforms
import py3Dmol
import numpy as np
import matplotlib.pyplot as plt

import json
import math

from itertools import combinations

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

BOND_ANGLES_TABLE = data["BOND_ANGLE_TABLE"]


def findBonds(G):
    """
    Creates a list of bonding between atoms including important properties. The
    bonding between atoms is represented as a dictionary in this form: 
    {'to': 'H', 'toString': 'H1', 'from': 'O', 'bond_order': 1} 
    'toString' represents the unique id of atom itself. 
    
    Args: 
        G (networkx.Graph): Graph representing the molecule 
    
    Returns: 
        bonds (list of dictionaries): Each dictionary is the bond between 
        central atom and the surrounding atom. 
    """
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


def showPieChart(bonding_pairs, lone_pairs):
    """
    This if for metric purposes. It creates a pie chart that distribution on the
    number of bonding pairs and lone pairs in a comprehensive pie chart. This
    shows which is type of domain doing the most repulsion.

    Args:
        bonding_pairs (int): Any surrounding that is connected to the central
        atom.
        lone_pairs (int): Lone pairs of valence electrons to the central.

    Returns:
        None
    """
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


def findVSEPR(G, central):
    """
    Finds the VSEPR / molecular geometry of the molecule from the number of
    bonding pairs and lone pairs. Each molecular structure is associated with
    a specific pattern. For example, "tetrahedral" is associated with 4 bonding
    pairs and 0 lone pairs, "trigonal pyramidal" is associated with 3 bonding
    pairs and 1 lone pairs, etc.

    Args:
        G (networkx.Graph): Graph representing the molecule 
        central (str): Name of central atom.

    Returns:
        structureVSEPR (str): Molecular Geometry
        electronGeometry (str): Electron Geometry is the molecular geometry
        without the lone pairs.
        hybridization (str): The orbital mixing type (e.g., 'sp3')
        based on the number of domain.
        lone_pairs (int): Number of lone pairs.
    """
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
        electronGeometry = "No electron geometry found."

    #showPieChart(bonding_pairs, lone_pairs) # uncomment this for metric purposes

    return structureVSEPR, electronGeometry, hybridization, lone_pairs


def findFormalCharges(atom2, G, bond, central):
    """
    Finds formal charges for a specific surrounding atom based on its bonding.

    Args:
        atom2 (str): The unique id of atom itself.
        G (networkx.Graph): G (networkx.Graph): Graph representing the molecule
        bond_order (dict): Includes the properties of the molecule.
        central (str): Name of central atom.

    Returns:
        stringFormalCharge (str): The charge (e.g. "2+", "-", "2-")
    """
    with open("static/periodic_table.json") as file:
        elements = json.load(file)

    atom = ''.join(char for char in atom2 if not char.isdigit())

    if atom != central:
        if bond is None:
            return ""

        atom = bond['to']
        specificAtom = bond["toString"]

        valenceElectrons = elements[atom]["valence_electrons"]
        bondingOrder = int(bond["bond_order"])
        lonePairs = int(G.nodes[specificAtom]["lone_pairs"]) * 2

        formalCharges = valenceElectrons - lonePairs - bondingOrder
    else:
        valenceElectrons = elements[central]["valence_electrons"]
        lonePairs = int(G.nodes[central]["lone_pairs"]) * 2

        bondingOrder = int(bond["bond_order"])

        formalCharges = valenceElectrons - lonePairs - bondingOrder

    if formalCharges > 0:
        stringFormalCharge = "+" + str(formalCharges)
    elif formalCharges < 0:
        stringFormalCharge = str(formalCharges)
    else:
        stringFormalCharge = ""
    return stringFormalCharge.replace("1", "")


def findSMILES(structure3D, centralAtom, bond_order, formalChargeDict):
    """
    Parses the molecule to proper SMILES formatting that RDKit can read.

    Args:
        structureVSEPR (str): Molecular Geometry
        centralAtom (str): Name of central atom.
        bonds (list of dictionaries): Each dictionary is the bond between 
        central atom and the surrounding atom. 
        formalChargeDict (dict): It contains each atom's specific formal charge.

    Returns:
        smilesForm (str): SMILES formatted string.
        expandedVSEPR (bool): Bool that checks if the VSEPR structure is more
        than 4 domains and has to be created manually.
    """
    bond_symbol = {
        1: "",
        2: "=",
        3: "#",
    }

    for bond in bond_order:
        if bond["to"] == "H":
            bond["to"] = "[H]"
    if centralAtom == "Al":
        centralAtom = "[Al]"

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
            if stringCentralAtom == "H":
                stringCentralAtom = "[H]"
            smilesForm = (
                f"{stringCentralAtom}"
                + bond_symbol[bond_order[0]["bond_order"]]
                + f"{bond_order[0]['to']}"
            )

    elif structure3D in ("Linear2", "Bent1", "Bent2"):
        leftForm = bond_order[0]["to"] + bond_symbol[bond_order[0]["bond_order"]]
        rightForm = bond_symbol[bond_order[1]["bond_order"]] + bond_order[1]["to"]

        smilesForm = f"{leftForm}[{centralAtom}]{rightForm}"

        if formalChargeDict:
            smilesForm = f"{leftForm}[{centralAtom}]{rightForm}"

    elif structure3D in ("Trigonal planar", "Trigonal pyramidal"):
        if not centralAtom == "[Al]":
            smilesForm = f"[{centralAtom}]"
        else:
            smilesForm = f"{centralAtom}"
        for i in range(len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"

    elif structure3D == "Tetrahedral":
        smilesForm = bond_order[0]['to'] + bond_symbol[bond_order[0]["bond_order"]] + f"[{centralAtom}]"
        for i in range(1, len(bond_order)):
            smilesForm += "(" + bond_symbol[bond_order[i]["bond_order"]] + bond_order[i]['to'] + ")"
    
    elif structure3D in (HIGH_DOMAIN_STRUCTURES):
        expandedVSEPR = True

        pt = Chem.GetPeriodicTable()

        central_num = pt.GetAtomicNumber(centralAtom)
        surroundings_list = []

        for i in range(len(bond_order)):
            surroundings_list.append(pt.GetAtomicNumber(bond_order[i]['to']))

        mol = RWMol()
        s_idx = mol.AddAtom(Chem.Atom(central_num))

        f_indices = []

        for surrounding_atom in surroundings_list:
            f_indices.append(mol.AddAtom(Chem.Atom(surrounding_atom)))
            
        for index, f_idx in enumerate(f_indices):
            if bond_order[index]['bond_order'] == 2:
                mol.AddBond(s_idx, f_idx, Chem.BondType.DOUBLE)
            elif bond_order[index]['bond_order'] == 3:
                mol.AddBond(s_idx, f_idx, Chem.BondType.TRIPLE)
            else:    
                mol.AddBond(s_idx, f_idx, Chem.BondType.SINGLE)

        smilesForm = mol.GetMol()

        coords = expanded_vsepr_structures(structure3D)

        conf = Conformer(smilesForm.GetNumAtoms())
        for i in range(smilesForm.GetNumAtoms()):
            conf.SetAtomPosition(i, coords[i])
        smilesForm.AddConformer(conf)

    else:
        smilesForm = centralAtom

    return smilesForm, expandedVSEPR


def expanded_vsepr_structures(structure):
    """
    Creates the atom in specific position in the py3Dmol viewer to mimic the
    correct VSEPR structure.

    Args:
        structureVSEPR (str): Molecular Geometry

    Returns:
        coords (np.array): Coordinates of the position of surrounding atoms.
    """
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


def lonePairPositions(view, positions):
    """
    Helper function that adds a sphere to the position in the py3Dmol viewer.

    Args:
        view (py3Dmol.view): The 3D py3Dmol viewer where the lone pair sphere 
        is added.
        positions (dict): The x, y, and z coordinates where the sphere will be 
        placed.

    Returns:
        None
    """
    view.addSphere({
        'center': positions,
        'radius': 0.3,
        'color': 'blue',
        'opacity': '0.90'
    })
    return 


def generatingLonePairs(view, structure, mol, central, lone_pairs):
    """
    Creates the lone pairs associated with the central atom of the molecule.

    Args:
        view (py3Dmol.view): The 3D py3Dmol viewer where the lone pair sphere 
        is added.
        structure (str): Molecular Geometry
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object containing the 
        atoms, bonds, and 3D coordinates of the molecule.
        central (str): Name of central atom.
        lone_pairs (int): Number of lone pairs.

    Returns:
        None
    """
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
                              {'x': c_pos[0], 'y': c_pos[1] + 0.35, 
                               'z': c_pos[2] + 0.35 * (-1)**i})
    elif structure in ("Seesaw", "T-shaped2", "Linear5"):
        for i in range(lone_pairs):
            lonePairPositions(view,
                              trigonal_bipyramidal_pos[i])    
    elif structure in ("Square pyramidal", "Square planar", "T-shaped3", 
                       "Linear6"):
        for i in range(lone_pairs):
            lonePairPositions(view,
                              octahedral_pos[i])
    return


def findBondAngles(mol, central, structure3D):
    """
    Calculates the bond angles of the constructed molecule using RDKit's
    conformer. Given the 3D position of each atom, it calculates the angle
    between the two in 2D space. 

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object containing the 
        atoms, bonds, and 3D coordinates of the molecule.
        central (str): Name of central atom.
        structure3D (str): Molecular Geometry

    Returns:
        angle_data (list of dictionaries): Stores the data associated with
        bond angle, such as the names of the two atoms calculated.
        angles (set): Avoid duplicates of bond angles and this variable will
        be the one to show to the user.
    """
    conf = mol.GetConformer()
    center_atom_obj = None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == central:
            center_atom_obj = atom

    if center_atom_obj is None:
        return [], None

    center_atom_index = center_atom_obj.GetIdx()

    neighbors = []
    for neighbor in center_atom_obj.GetNeighbors():
        neighbors.append(neighbor.GetIdx())

    if len(neighbors) < 2:
        return [], None
    
    angle_data = []
    angles = set()

    for atom1_index, atom3_index in combinations(neighbors, 2):
        angle = rdMolTransforms.GetAngleDeg(
            conf,
            atom1_index,
            center_atom_index,
            atom3_index
        )
        angle = round(angle, 2)
    
        if structure3D not in ("Linear2", "Linear5", "Linear6"):
            if angle == 180.0:
                continue

        if angle in angles:
            continue
        
        angles.add(angle)

        angle_data.append({
            "atom1": atom1_index,
            "center": center_atom_index,
            "atom3": atom3_index,
            "angle": angle
        })

    return angle_data, angles


def bondAngleGraph(angles, structure3D):
    """
    For metrics purposes, creates a graph that compares the bond angles and
    the awaited bond angles associated with its VSEPR structure.

    Args:
        angles (set): Description here.
        structure3D (str): Molecular Geometry

    Returns:
        None
    """
    if not angles:
        return

    angles = sorted(list(angles))
    ideal_angles = BOND_ANGLES_TABLE[structure3D]

    # making sure that both lists are the same minimum length
    min_len = min(len(angles), len(ideal_angles))
    angles = angles[:min_len]
    ideal_angles = ideal_angles[:min_len]

    angles_label = []
    for i in range(min_len):
        angles_label.append(f"Angles {i + 1}")

    x = np.arange(min_len)
    w = 0.4
    
    plt.figure(figsize=(7,5))

    bar1 = plt.bar(x-w/2, ideal_angles, w, label="Ideal Angles")
    bar2 = plt.bar(x+w/2, angles, w, label="Angles Calculated")
    
    plt.bar_label(bar1, fmt="%.3f")
    plt.bar_label(bar2, fmt="%.3f")

    all_values = ideal_angles + angles
    y_min = min(all_values) - 5
    y_max = max(all_values) + 5
    plt.ylim(y_min, y_max)
    plt.ylim(y_min, y_max)

    plt.xticks(x, angles_label)
    plt.ylabel("Degree")
    plt.title(f"Ideal Bond Angle vs. Calculated Bond Angles: {structure3D}")
    plt.legend()

    plt.savefig("double_bar_graph.png", dpi=300, bbox_inches="tight")
    plt.close()
    return


def xyz(pos):
    """
    Helper function that converts a 3D position into a dictionary format that
    py3Dmol can use.

    Args:
        pos (np.array): The x, y, and z coordinates of a position.

    Returns:
        position_dict (dict): Contains the x, y, and z coordinates.
    """
    return {
        'x': float(pos[0]),
        'y': float(pos[1]),
        'z': float(pos[2])
    }


def get_pos(conf, atom_idx):
    """
    Gets the 3D coordinates from the RDKit conformer.

    Args:
        conf (rdkit.Chem.rdchem.Conformer): The RDKit conformer containing the
        3D coordinates of the molecule.
        atom_idx (int): The index of the atom.

    Returns:
        position_array (np.array): Array containing the x, y, and z coordinates
        of the atom.
    """
    p = conf.GetAtomPosition(atom_idx)
    return np.array([p.x, p.y, p.z])


def addBondAngleArc(view, mol, atom1_idx, center_idx, atom3_idx,
                          angle):
    """
    Adds a visual bond angle arc between two bonded atoms around the central
    atom in the py3Dmol viewer.

    Args:
        view (py3Dmol.view): The 3D py3Dmol viewer where the bond angle arc is
        added.
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object containing the 
        atoms, bonds, and 3D coordinates of the molecule.
        atom1_idx (int): Index of the first surrounding atom.
        center_idx (int): Index of the central atom.
        atom3_idx (int): Index of the second surrounding atom.
        angle (set): The calculated bond angles in degrees.

    Returns:
        None
    """
    size = 0.75
    color = "yellow"
    
    conf = mol.GetConformer()
    
    p1 = get_pos(conf, atom1_idx)
    pc = get_pos(conf, center_idx)
    p3 = get_pos(conf, atom3_idx)

    # find positions from where to start the arc and where to end
    start = pc + (p1 - pc) * size
    end = pc + (p3 - pc) * size

    middle = (start + end) / 2

    push_direction = middle - pc

    if np.linalg.norm(push_direction) != 0:
        push_direction = push_direction / np.linalg.norm(push_direction)
        middle = middle + push_direction * 0.15
    else:
        middle = pc + np.array([0, 0, 0.45])

    view.addCylinder({
        'start': xyz(start),
        'end': xyz(middle),
        'radius': 0.025,
        'color': color,
        'opacity': 0.9
    })

    view.addCylinder({
        'start': xyz(middle),
        'end': xyz(end),
        'radius': 0.025,
        'color': color,
        'opacity': 0.9
    })

    view.addLabel(f"{angle}°", {
        'position': xyz(middle),
        'fontColor': color,
        'backgroundColor': 'black',
        'backgroundOpacity': 0.5,
        'fontSize': 14
    })


def add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structure3D, 
                       electronGeometry, molecule, hybridization, lone_pairs,
                        angle):
    """
    Overwrites the outputted HTML file with a HTML tab that shows extra 
    information about the VSEPR structure and molecule.
    """
    
    
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
            <b>Hybrid:</b> {hybridization}<br>
            <b>Bond Angles:</b> {angle}
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


def main(formula):
    """
    Main code logic and error detection. It creates the HTML file with the
    interactive 3D molecule.

    Args:
        formula (str): The user inputted string taken from the app.py and the
        chatbot recognition.

    Returns:
        None
    """
    try:
        with open(f"static/generated/vsepr_{formula}.html", "r") as f:
            content = f.read()
            return None, None
    except FileNotFoundError:
        pass

    molecule, central, networkxGraph = prototype_lewis_structure.main(formula)

    if central == None:
        return molecule, molecule

    bonds = findBonds(networkxGraph)
    structureVSEPR, electronVSEPR, hybridization, lone_pairs = findVSEPR(networkxGraph, central)
    
    formalChargeDict = {}
    containsHydrogen = False

    for key in molecule:
        if key == "H":
            containsHydrogen = True
    
    for bond in bonds:
        key = bond["toString"]
        formalCharge = findFormalCharges(key, networkxGraph, bond, central)
        if formalCharge:
            formalChargeDict[key] = formalCharge

    try:
        smiles, expandedVSEPR = findSMILES(structureVSEPR, central, bonds, formalChargeDict)
    except:
        error_message = "Even though I have a PhD in chemistry, I cannot " \
        "generate the correct SMILES for this molecule."
        return None, error_message

    if not expandedVSEPR:
        try:
            mol = Chem.MolFromSmiles(smiles)

            # Add explicit hydrogens
            if containsHydrogen:
                mol = Chem.AddHs(mol)

            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except:
            error_message = "Even though I have a PhD in chemistry, I cannot " \
            "generate the correct SMILES for this molecule."
            return None, error_message
    else:
        mol = smiles
    
    # Convert molecule to MOL block format
    mol_block = Chem.MolToMolBlock(mol)

    # Send to Py3Dmol viewer
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})

    try:
        generatingLonePairs(view, structureVSEPR, mol, central, lone_pairs)
    except:
        error_message = "Even though I have a PhD in chemistry, I cannot " \
        "generate this molecule's lone pairs."
        return None, error_message

    try:
        angle_data, angle = findBondAngles(mol, central, structureVSEPR)

        for item in angle_data:
            addBondAngleArc(
                view,
                mol,
                item["atom1"],
                item["center"],
                item["atom3"],
                item["angle"]
            )
        #bondAngleGraph(angle, structureVSEPR)
    except:
        angle = "No VSEPR angles associated with this molecule."

    view.zoomTo({"sphere": {"scale": 0.5}}, 2000, True)

    output_file = f"static/generated/VSEPR_{formula}.html"
    with open(output_file, "w", encoding="utf-8") as f:
        view.write_html(f, fullpage=True)

    # added: get atoms from molecule and add legend
    atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    add_legend_to_html(output_file, atom_symbols, expandedVSEPR, structureVSEPR, 
                       electronVSEPR, formula, hybridization, lone_pairs,
                       angle)
    
    return None, None