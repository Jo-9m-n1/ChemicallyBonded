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
    Find the central atom in the Lewis structure by selecting
    the atom with the lowest electronegativity (excluding hydrogen).

    Args:
        chemDict (dict): Dictionary containing atoms and their counts

    Returns:
        str: The central atom of the molecule
    """
    candidates = {}
    for atom in chemDict:
        if atom == "H":
            continue
        candidates[atom] = elements[atom]["electronegativity"]
    
    central_atom = None
    min_value = None

    for atom, value in candidates.items():
        if min_value is None or value < min_value:
            min_value = value
            central_atom = atom

    return central_atom

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
    based on simple valence rules.

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
    Upgrade single bonds to double or triple bonds in order to
    satisfy bonding requirements of the central atom.

    Args:
        G (networkx.Graph): Graph representing the molecule
        central (str): The central atom of the molecule

    Returns:
        None
    """
    neighbors = list(G.neighbors(central))

    limit = get_max_bonds(G.nodes[central]["element"])

    while G.nodes[central]["bonds"] < limit:
        updated = False

        for neighbor in neighbors:
            if G.nodes[neighbor]["bonds"] < 4 and G.nodes[central]["lone_electrons"] > 0 and G[central][neighbor]["order"] < 3:
                G[central][neighbor]["order"] += 1

                G.nodes[central]["lone_electrons"] -= 1
                G.nodes[neighbor]["lone_electrons"] -= 1

                G.nodes[central]["bonds"] += 1
                G.nodes[neighbor]["bonds"] += 1

                updated = True

        if not updated:
            break

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
    and lone pairs positioned around each atom.

    Args:
        G (networkx.Graph): Graph representing the molecule
        positions (dict): Dictionary of node positions

    Returns:
        None
    """
    position = {node: tuple(positions[node])for node in positions}
    dwg = svgwrite.Drawing("static/molecule.svg", size=("400px", "400px"))

    scale = 150
    canvas_offset_x = 200
    canvas_offset_y = 200

    # Draw bonds
    for atom1, atom2, bond_data in G.edges(data = True):
        order = bond_data.get("order", 1)

        x1, y1 = position[atom1]
        x2, y2 = position[atom2]

        # direction vector
        direction_x = x2 - x1
        direction_y = y2 - y1
        length = np.sqrt(direction_x**2 + direction_y**2)

        # perpendicular vector
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
                stroke="black"
            ))

        print(atom1, "-", atom2, "order:", bond_data["order"])

    # Draw atoms
    for node in G.nodes():
        x, y = position[node]

        neighbors = list(G.neighbors(node))

        bond_angles = []
        for neighbor in neighbors:
            neighbor_x, neighbor_y = position[neighbor]
            angle = np.arctan2(neighbor_y - y, neighbor_x - x)
            bond_angles.append(angle)
        
        all_angles = [i * (np.pi / 2) for i in range(4)]  # 0°, 90°, 180°, 270°

        free_angles = []

        for a in all_angles:
            if all(angle_diff(a, b) > 0.5 for b in bond_angles):
                free_angles.append(a)

        dwg.add(dwg.circle(
            center=(x * scale + canvas_offset_x, -y * scale + canvas_offset_y),
            r=20,
            fill="white",
            stroke="black"
        ))

        lone_pairs = G.nodes[node].get("lone_pairs", 0)
        for i in range(min(lone_pairs, len(free_angles))):
            angle = free_angles[i]

            direction_x = 30 * np.cos(angle)
            direction_y = 30 * np.sin(angle)

            # draw 2 dots per pair
            dwg.add(dwg.circle(
                center=(x * scale + canvas_offset_x + direction_x - 3,
                        -y * scale + canvas_offset_y + direction_y),
                r=2,
                fill="black"
            ))

            dwg.add(dwg.circle(
                center=(x * scale + canvas_offset_x + direction_x + 3,
                        -y * scale + canvas_offset_y + direction_y),
                r=2,
                fill="black"
            ))

        dwg.add(dwg.text(
            G.nodes[node]["element"],
            insert=(x * scale + canvas_offset_x - 8, -y * scale + canvas_offset_y + 5)
        ))

    dwg.save()

def main(formula):
    parsed = parse_formula(formula)

    print("Valence electrons per atom:")

    for atom, count in parsed.items():
        valence = get_valence_electrons(atom)
        print(f"{atom}: {valence} × {count} = {valence * count}")

    total_valence = sum(
        get_valence_electrons(atom) * count
        for atom, count in parsed.items()
    )

    print("Total valence electrons:", total_valence)

    central = find_central_atom(parsed)
    G = build_graph(parsed, central)
    upgrade_bonds(G, central)
    assign_lone_pairs(G)
    positions = generate_positions(G, central)
    draw_molecule(G, positions)

    return parsed, central, G