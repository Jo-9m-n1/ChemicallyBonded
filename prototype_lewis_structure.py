import networkx as nx
import numpy as np
import svgwrite
import json

with open("static/periodic_table.json") as file:
    elements = json.load(file)

def parse_formula(formula):
    """
    Parse a chemical formula string into a dictionary containing
    each atom and its quantity in the molecule.

    Args:
        formula (str): The chemical formula (e.g., "H2O", "CO2")

    Returns:
        dict: A dictionary mapping each atom to its count
    """
    result = {}
    i = 0

    while i < len(formula):
        atom = formula[i]
        i += 1

        if i < len(formula) and formula[i].islower():
            atom += formula[i]
            i += 1
        
        num = ""
        while i < len(formula) and formula[i].isdigit():
            num += formula[i]
            i += 1
        
        count = 0
        if num:
            count = int(num)
        else:
            count = 1
        
        result[atom] = result.get(atom, 0) + count

    return result

def validate_molecule(chemDict):
    """
    Validate if the molecule is supported:
    - Reject ionic (metal + nonmetal)
    - Reject impossible valence combinations

    Returns:
        (bool, str): (is_valid, message)
    """
    has_metal = any(elements[a]["type"] == "metal" for a in chemDict)
    has_nonmetal = any(elements[a]["type"] == "nonmetal" for a in chemDict)

    if has_metal and has_nonmetal:
        return False, "Ionic compounds are not supported."

    for atom in chemDict:
        if elements[atom]["max_bonds"] == 0:
            return False, f"{atom} does not form bonds."

    total_bonds_possible = sum(
        elements[a]["max_bonds"] * count
        for a, count in chemDict.items()
    )

    total_atoms = sum(chemDict.values())

    if total_bonds_possible < (total_atoms - 1) * 2:
        return False, "Molecule is not chemically valid (insufficient bonding)."

    return True, "Valid molecule"

def is_diatomic(chemDict):
    """
    Check if the molecule is diatomic (consists of exactly two identical atoms).

    Args:
        chemDict (dict): Dictionary containing atoms and their counts

    Returns:
        bool: True if the molecule is diatomic, False otherwise
    """
    return len(chemDict) == 1 and list(chemDict.values())[0] == 2

def get_valence_electrons(atom):
    """
    Retrieve the number of valence electrons for a given atom
    from the periodic table data.

    Args:
        atom (str): The chemical symbol of the atom

    Returns:
        int: Number of valence electrons for the atom
    """
    return elements[atom]["valence_electrons"]

def find_central_atom(chemDict):
    """
    Determine the central atom of the molecule using a priority list
    of common central atoms.

    Args:
        chemDict (dict): Dictionary containing atoms and their counts

    Returns:
        str: The selected central atom of the molecule
    """
    priority = ["C", "N", "P", "S", "O", "Br", "Cl", "I", "F"]

    candidates = [atom for atom in chemDict if atom != "H"]

    for p in priority:
        if p in candidates:
            return p

    return candidates[0]

def assign_lone_pairs(G):
    """
    Calculate and assign the number of lone pairs for each atom
    based on remaining non-bonding electrons.

    Args:
        G (networkx.Graph): Graph representing the molecule

    Returns:
        None
    """
    for node in G.nodes():
        element = G.nodes[node]["element"]

        if element == "H":
            G.nodes[node]["lone_pairs"] = 0
            continue

        electrons = G.nodes[node]["lone_electrons"]
        G.nodes[node]["lone_pairs"] = electrons // 2

def build_diatomic_graph(chemDict):
    """
    Construct a graph for a diatomic molecule consisting of two identical atoms
    connected by a bond.

    Args:
        chemDict (dict): Dictionary of atoms and their counts

    Returns:
        tuple: (Graph of the molecule, one atom label used as reference)
    """
    G = nx.Graph()

    atom = list(chemDict.keys())[0]

    atom1 = f"{atom}1"
    atom2 = f"{atom}2"

    valence = get_valence_electrons(atom)

    G.add_node(atom1, element = atom, valence = valence, lone_electrons = valence, bonds = 0)
    G.add_node(atom2, element = atom, valence = valence, lone_electrons = valence, bonds = 0)

    G.add_edge(atom1, atom2, order=1)

    G.nodes[atom1]["bonds"] += 1
    G.nodes[atom2]["bonds"] += 1

    G.nodes[atom1]["lone_electrons"] -= 1
    G.nodes[atom2]["lone_electrons"] -= 1

    return G, atom1

def build_graph(chemDict, central):
    """
    Construct a graph representation of the molecule where atoms
    are nodes and bonds are edges, initially with single bonds.

    Args:
        chemDict (dict): Dictionary of atoms and their counts
        central (str): The central atom of the molecule

    Returns:
        networkx.Graph: Graph representing the molecule structure
    """
    G = nx.Graph()

    G.add_node(central, element = central, valence = get_valence_electrons(central), lone_electrons = get_valence_electrons(central), bonds = 0)

    for atom, count in chemDict.items():
        if atom == central:
            continue

        for i in range(count):
            atom_name = f"{atom}{i+1}"

            G.add_node(atom_name, element = atom, valence = get_valence_electrons(atom), lone_electrons = get_valence_electrons(atom), bonds = 0)
            G.add_edge(central, atom_name, order = 1)

            G.nodes[central]["lone_electrons"] -= 1
            G.nodes[atom_name]["lone_electrons"] -= 1

            G.nodes[central]["bonds"] += 1
            G.nodes[atom_name]["bonds"] += 1
    return G

def get_max_bonds(element):
    """
    Determine the maximum number of bonds an atom can form
    based on simplified chemical rules.

    Args:
        element (str): The chemical symbol of the atom

    Returns:
        int: Maximum number of bonds allowed for the atom
    """
    if element == "H":
        return 1
    elif element in ["F", "Cl", "Br", "I"]:
        return 1
    elif element == "O":
        return 2
    elif element == "N":
        return 3
    elif element == "C":
        return 4
    else:
        return 1

def upgrade_bonds(G, central):
    """
    Upgrade bonds by increasing bond order based on which atom
    needs more bonds. Works for both normal and diatomic molecules.

    Args:
        G (networkx.Graph): Graph representing the molecule
        central (str): Reference atom (central or first atom in diatomic)

    Returns:
        None
    """
    neighbors = list(G.neighbors(central))

    # diatomic case
    if len(neighbors) == 1:
        neighbor = neighbors[0]
        element = G.nodes[neighbor]["element"]

        max_bonds = get_max_bonds(element)

        while G[central][neighbor]["order"] < max_bonds:
            G[central][neighbor]["order"] += 1

            G.nodes[central]["lone_electrons"] -= 1
            G.nodes[neighbor]["lone_electrons"] -= 1

            G.nodes[central]["bonds"] += 1
            G.nodes[neighbor]["bonds"] += 1

        return

    # normal cases
    while True:
        best_neighbor = None
        best_need = 0

        for neighbor in neighbors:
            element = G.nodes[neighbor]["element"]
            current_bonds = G.nodes[neighbor]["bonds"]
            max_bonds = get_max_bonds(element)

            need = max_bonds - current_bonds

            # prevent hydrogen from upgrading
            if element == "H":
                continue

            if need > best_need and G[central][neighbor]["order"] < 3:
                best_need = need
                best_neighbor = neighbor

        if best_neighbor is None:
            break

        G[central][best_neighbor]["order"] += 1

        G.nodes[central]["lone_electrons"] -= 1
        G.nodes[best_neighbor]["lone_electrons"] -= 1

        G.nodes[central]["bonds"] += 1
        G.nodes[best_neighbor]["bonds"] += 1

def generate_positions(G, central, bond_length = 1.0):
    """
    Generate 2D positions for each atom by placing surrounding atoms
    evenly around the central atom in a circular layout.

    Args:
        G (networkx.Graph): Graph representing the molecule
        central (str): The central atom
        bond_length (float): Distance from central atom to neighbors

    Returns:
        dict: Dictionary mapping each node to its 2D position
    """
    positions = {}
    positions[central] = np.array([0.0, 0.0])
    
    neighbors = list(G.neighbors(central))
    n = len(neighbors)
    
    if n == 2 and central == "O":
        angle = np.radians(104.5)
        positions[neighbors[0]] = np.array([-0.5, -0.5])
        
        positions[neighbors[1]] = np.array([0.5, -0.5])

    else:
        if n == 0:
            return positions
        angle_step = 2 * np.pi / n

        for i, atom in enumerate(neighbors):
            angle = i * angle_step
            positions[atom] = np.array([
                bond_length * np.cos(angle),
                bond_length * np.sin(angle)
            ])

    return positions

def angle_diff(a, b):
    """
    Compute the smallest difference between two angles, taking into
    account circular wrapping.

    Args:
        a (float): First angle in radians
        b (float): Second angle in radians

    Returns:
        float: Absolute angular difference
    """
    return abs((a - b + np.pi) % (2 * np.pi) - np.pi)

def draw_molecule(G, positions):
    """
    Render the molecule as an SVG image, including atoms, bonds,
    and lone pairs positioned symmetrically around each atom.
    """
    position = {node: tuple(positions[node]) for node in positions}
    dwg = svgwrite.Drawing("molecule.svg", size=("400px", "400px"))

    scale = 150
    canvas_offset_x = 200
    canvas_offset_y = 200

    for atom1, atom2, bond_data in G.edges(data=True):
        order = bond_data.get("order", 1)
        x1, y1 = position[atom1]
        x2, y2 = position[atom2]

        direction_x = x2 - x1
        direction_y = y2 - y1
        length = np.sqrt(direction_x**2 + direction_y**2)

        perp_x = -direction_y / length
        perp_y = direction_x / length

        for i in range(order):
            offset = (i - (order - 1) / 2) * 10
            dwg.add(dwg.line(
                start=(
                    x1 * scale + canvas_offset_x + perp_x * offset,
                    -y1 * scale + canvas_offset_y + perp_y * offset
                ),
                end=(
                    x2 * scale + canvas_offset_x + perp_x * offset,
                    -y2 * scale + canvas_offset_y + perp_y * offset
                ),
                stroke="white"
            ))

    for node in G.nodes():
        x, y = position[node]
        neighbors = list(G.neighbors(node))
        n = len(neighbors)

        bond_angles = []
        for neighbor in neighbors:
            nb_x, nb_y = position[neighbor]
            angle = np.arctan2(nb_y - y, nb_x - x)
            bond_angles.append(angle)
        
        dwg.add(dwg.circle(
            center=(x * scale + canvas_offset_x, -y * scale + canvas_offset_y),
            r=20,
            fill="white",
            stroke="black"
        ))

        lone_pairs = G.nodes[node].get("lone_pairs", 0)
        if lone_pairs > 0:
            if n == 2:
                avg_bond_angle = np.arctan2(
                    np.sum([np.sin(a) for a in bond_angles]),
                    np.sum([np.cos(a) for a in bond_angles])
                )
                opposite_angle = avg_bond_angle + np.pi
                
                spread = np.radians(50) 
                if lone_pairs == 1:
                    pair_angles = [opposite_angle]
                else:
                    pair_angles = [opposite_angle - spread/2, opposite_angle + spread/2]
            else:
                all_angles = [i * (np.pi / 2) for i in range(4)]
                pair_angles = [a for a in all_angles if all(angle_diff(a, b) > 0.5 for b in bond_angles)]

            for i in range(min(lone_pairs, len(pair_angles))):
                angle = pair_angles[i]
                dist = 32
                dot_sep = 4
                
                px, py = -np.sin(angle) * dot_sep, np.cos(angle) * dot_sep

                for offset in [-1, 1]:
                    dwg.add(dwg.circle(
                        center=(
                            x * scale + canvas_offset_x + dist * np.cos(angle) + px * offset,
                            -y * scale + canvas_offset_y - dist * np.sin(angle) - py * offset
                        ),
                        r=2,
                        fill="white"
                    ))

        dwg.add(dwg.text(
            G.nodes[node]["element"],
            insert=(x * scale + canvas_offset_x - 8, -y * scale + canvas_offset_y + 5),
            fill="black",
            font_size="16px",
            font_family="Arial"
        ))

    dwg.save()

def main(formula):
    parsed = parse_formula(formula)

    valid, message = validate_molecule(parsed)

    if not valid:
        print(message)
        exit()

    print("Valence electrons per atom:")

    for atom, count in parsed.items():
        valence = get_valence_electrons(atom)
        print(f"{atom}: {valence} × {count} = {valence * count}")

    total_valence = sum(
        get_valence_electrons(atom) * count
        for atom, count in parsed.items()
    )

    print("Total valence electrons:", total_valence)

    if is_diatomic(parsed):
        G, central = build_diatomic_graph(parsed)
        upgrade_bonds(G, central)

    else:
        central = find_central_atom(parsed)
        G = build_graph(parsed, central)
        upgrade_bonds(G, central)
        
    assign_lone_pairs(G)
    positions = generate_positions(G, central)
    draw_molecule(G, positions)
    return parsed, central, G