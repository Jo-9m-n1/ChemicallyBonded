import networkx as nx
import numpy as np
import svgwrite
import json

from pathlib import Path

with open("static/periodic_table.json") as file:
    elements = json.load(file)

covalent_metals = ["Be", "Al", "Sn"]
expanded_octet_atoms = ["P", "S", "Se", "Te", "Cl", "Br", "I", "Xe"]
central_atom_priority = ["C", "N", "P", "S", "Be", "Al", "Sn", "Xe", "Cl", "Br", "O"]
halogen_atoms = ["F", "Cl", "Br", "I"]
terminal_atoms = ["H", "F", "Cl", "Br", "I"]

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
    - Reject atoms that don't exist
    - Reject ionic (metal + nonmetal)
    - Reject impossible valence combinations

    Returns:
        (bool, str): (is_valid, message)
    """
    has_metal = False
    has_nonmetal = False
    has_expanded_octet_atom = False

    for atom in chemDict:
        if atom not in elements:
            return False, f"{atom} is not a valid atom symbol in my knowledge! Could you please try with another molecule?"

        atom_type = elements[atom]["type"]

        if atom_type == "metal" and atom not in covalent_metals:
            has_metal = True
        elif atom_type == "nonmetal":
            has_nonmetal = True

    if has_metal and has_nonmetal:
        return False, "I may have a PhD in chemistry but I do not know how to draw ionic compounds."

    for atom in chemDict:
        if elements[atom]["max_bonds"] == 0:
            return False, f"Hmmmmm, {atom} does not form bonds."

    if not is_supported_structure(chemDict):
        return False, "I may have a PhD in chemistry but I have no clue how to generate it."

    total_bonds_possible = 0
    for atom, count in chemDict.items():
        max_bonds = elements[atom]["max_bonds"]
        total_bonds_possible += max_bonds * count

    total_atoms = sum(chemDict.values())
    for atom in chemDict:
        if atom in expanded_octet_atoms:
            has_expanded_octet_atom = True
            break

    if not has_expanded_octet_atom:
        minimum_bonds_needed = (total_atoms - 1) * 2
        if total_bonds_possible < minimum_bonds_needed:
            return False, "This molecule is chemically invalid (insufficient bonding)... Try again, unless you want discover a new molecule."

    return True, "Valid molecule"

def is_diatomic(chemDict):
    """
    Check if the molecule is diatomic (consists of exactly two identical atoms).

    Args:
        chemDict (dict): Dictionary containing atoms and their counts

    Returns:
        bool: True if the molecule is diatomic, False otherwise
    """
    if len(chemDict) != 1:
        return False
    
    return list(chemDict.values())[0] == 2

def is_supported_structure(chemDict):
    """
    Check whether the molecule fits the structural limits of this program.

    Supported structures are:
    - diatomic molecules
    - single-element molecules
    - molecules with up to 7 total atoms
    - molecules with 1 central atom and up to 6 surrounding atoms

    Args:
        chemDict (dict): Dictionary containing atoms and their counts

    Returns:
        bool: True if the structure is supported, False otherwise
    """
    if is_diatomic(chemDict):
        return True

    total_atoms = sum(chemDict.values())
    if total_atoms > 7:
        return False

    if len(chemDict) == 1:
        only_atom_count = list(chemDict.values())[0]
        return only_atom_count <= 3

    central = find_central_atom(chemDict)
    outer_atoms = total_atoms - chemDict[central]

    if outer_atoms > 6:
        return False

    return True

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
    if len(chemDict) == 1:
        return list(chemDict.keys())[0]

    non_hydrogen_atoms = [atom for atom in chemDict if atom != "H"]

    non_halogen_atoms = [atom for atom in non_hydrogen_atoms if atom not in halogen_atoms]

    if non_halogen_atoms:
        candidate_atoms = non_halogen_atoms
    else:
        candidate_atoms = non_hydrogen_atoms

    for atom in central_atom_priority:
        if atom in candidate_atoms:
            return atom
    
    return candidate_atoms[0]

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

def add_atom_node(G, atom_name, atom):
    """
    Add an atom node to the graph with its default properties.

    Args:
        G (networkx.Graph): Graph representing the molecule
        atom_name (str): Unique node name
        atom (str): Atom symbol

    Returns:
        None
    """
    valence = get_valence_electrons(atom)
    G.add_node(atom_name, element=atom, valence=valence, lone_electrons=valence, bonds=0)

def connect_atoms(G, atom1, atom2):
    """
    Connect two atoms with a single bond and update their
    bond and lone electron counts.

    Args:
        G (networkx.Graph): Graph representing the molecule
        atom1 (str): First atom label
        atom2 (str): Second atom label

    Returns:
        None
    """
    G.add_edge(atom1, atom2, order=1)
    G.nodes[atom1]["bonds"] += 1
    G.nodes[atom2]["bonds"] += 1
    G.nodes[atom1]["lone_electrons"] -= 1
    G.nodes[atom2]["lone_electrons"] -= 1

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

    add_atom_node(G, atom1, atom)
    add_atom_node(G, atom2, atom)

    connect_atoms(G, atom1, atom2)

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
    add_atom_node(G, central, central)

    for atom, count in chemDict.items():
        if atom == central:
            atoms_to_add = count - 1
        else:
            atoms_to_add = count

        for i in range(atoms_to_add):
            atom_name = f"{atom}{i+1}"
            add_atom_node(G, atom_name, atom)
            connect_atoms(G, central, atom_name)

    return G

def increase_bond_order(G, atom1, atom2):
    """
    Increase the bond order between two atoms by 1 and update
    their bond and lone electron counts.

    Args:
        G (networkx.Graph): Graph representing the molecule
        atom1 (str): First atom label
        atom2 (str): Second atom label

    Returns:
        None
    """
    G[atom1][atom2]["order"] += 1
    G.nodes[atom1]["lone_electrons"] -= 1
    G.nodes[atom2]["lone_electrons"] -= 1
    G.nodes[atom1]["bonds"] += 1
    G.nodes[atom2]["bonds"] += 1

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

    if len(G.nodes()) == 3 and len(neighbors) == 2:
        central_element = G.nodes[central]["element"]

        all_neighbors_match_central = True
        for neighbor in neighbors:
            if G.nodes[neighbor]["element"] != central_element:
                all_neighbors_match_central = False
                break

        if all_neighbors_match_central:
            chosen_neighbor = neighbors[0]
            increase_bond_order(G, central, chosen_neighbor)
            return

    # diatomic case
    if len(neighbors) == 1:
        neighbor = neighbors[0]
        element = G.nodes[neighbor]["element"]

        max_bonds = elements[element]["max_bonds"]

        while G[central][neighbor]["order"] < max_bonds:
            increase_bond_order(G, central, neighbor)

        return

    # normal cases
    while True:
        best_neighbor = None
        best_need = 0

        for neighbor in neighbors:
            element = G.nodes[neighbor]["element"]
            current_bonds = G.nodes[neighbor]["bonds"]

            if element in terminal_atoms:
                continue

            max_bonds = elements[element]["max_bonds"]
            need = max_bonds - current_bonds

            if need > best_need and G[central][neighbor]["order"] < 3:
                best_need = need
                best_neighbor = neighbor

        if best_neighbor is None:
            break

        increase_bond_order(G, central, best_neighbor)

def generate_positions(G, central):
    """
    Generate 2D positions for each atom by placing surrounding atoms
    around the central atom.

    Args:
        G (networkx.Graph): Graph representing the molecule
        central (str): The central atom

    Returns:
        dict: Dictionary mapping each node to its 2D position
    """
    positions = {central: np.array([0.0, 0.0])}

    neighbors = list(G.neighbors(central))
    n = len(neighbors)
    lone_pairs = G.nodes[central].get("lone_pairs", 0)

    if n == 2:
        if lone_pairs == 3:
            positions[neighbors[0]] = np.array([-1.0, 0.0])
            positions[neighbors[1]] = np.array([1.0, 0.0])
            return positions

        if lone_pairs >= 1:
            bond_angle = np.radians(104.5)
            half = bond_angle / 2

            positions[neighbors[0]] = np.array([-np.cos(half), -np.sin(half)])
            positions[neighbors[1]] = np.array([np.cos(half), -np.sin(half)])
            return positions

        positions[neighbors[0]] = np.array([-1.0, 0.0])
        positions[neighbors[1]] = np.array([1.0, 0.0])
        return positions

    if n == 0:
        return positions

    angle_step = 2 * np.pi / n

    for i, atom in enumerate(neighbors):
        angle = i * angle_step
        positions[atom] = np.array([np.cos(angle), np.sin(angle)])

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

def to_canvas_coords(x, y, scale, canvas_offset_x, canvas_offset_y):
    """
    Convert molecule coordinates to SVG canvas coordinates.

    Args:
        x (float): X position in molecule space
        y (float): Y position in molecule space
        scale (float): Scale factor for drawing
        canvas_offset_x (float): Horizontal canvas offset
        canvas_offset_y (float): Vertical canvas offset

    Returns:
        tuple: (canvas_x, canvas_y)
    """
    canvas_x = x * scale + canvas_offset_x
    canvas_y = -y * scale + canvas_offset_y
    return canvas_x, canvas_y

def get_lone_pair_angles(bond_angles, lone_pairs, neighbor_count):
    """
    Determine the angles where lone pairs should be drawn.

    Args:
        bond_angles (list): List of bond angles in radians
        lone_pairs (int): Number of lone pairs
        neighbor_count (int): Number of bonded neighbors

    Returns:
        list: Angles for lone pair placement
    """
    if lone_pairs <= 0:
        return []

    if neighbor_count == 2:
        sin_sum = 0
        cos_sum = 0

        for angle in bond_angles:
            sin_sum += np.sin(angle)
            cos_sum += np.cos(angle)
        avg_bond_angle = np.arctan2(sin_sum, cos_sum)

        opposite_angle = avg_bond_angle + np.pi

        if lone_pairs == 1:
            return [opposite_angle]

        if lone_pairs == 2:
            spread = np.radians(50)
            left_angle = opposite_angle - spread / 2
            right_angle = opposite_angle + spread / 2
            return [left_angle, right_angle]

        if lone_pairs == 3:
            upper_angle = opposite_angle
            left_angle = opposite_angle + (2 * np.pi / 3)
            right_angle = opposite_angle - (2 * np.pi / 3)
            return [upper_angle, left_angle, right_angle]

        return []

    primary_angles = [i * (np.pi / 2) for i in range(4)]
    secondary_angles = [np.pi / 4 + i * (np.pi / 2) for i in range(4)]
    all_angles = primary_angles + secondary_angles

    pair_angles = []
    for angle in all_angles:
        far_from_all_bonds = True
        for bond_angle in bond_angles:
            if angle_diff(angle, bond_angle) <= 0.5:
                far_from_all_bonds = False
                break

        if far_from_all_bonds:
            pair_angles.append(angle)

    return pair_angles

def draw_molecule(G, positions, molecule):
    """
    Render the molecule as an SVG image, including atoms, bonds,
    and lone pairs positioned symmetrically around each atom.

    Args:
        G (networkx.Graph): Graph representing the molecule
        positions (dict): Dictionary mapping each node to its 2D position

    Returns:
        None
    """
    file_path = Path(f"static/generated/molecule_{molecule}.svg")
    if file_path.is_file():
        return
    
    dwg = svgwrite.Drawing(f"static/generated/molecule_{molecule}.svg", size=("400px", "400px"))

    scale = 150
    canvas_offset_x = 200
    canvas_offset_y = 200

    for atom1, atom2, bond_data in G.edges(data=True):
        order = bond_data.get("order", 1)
        x1, y1 = positions[atom1]
        x2, y2 = positions[atom2]

        direction_x = x2 - x1
        direction_y = y2 - y1
        length = np.sqrt(direction_x**2 + direction_y**2)

        perp_x = -direction_y / length
        perp_y = direction_x / length

        for i in range(order):
            offset = (i - (order - 1) / 2) * 10

            start_x, start_y = to_canvas_coords(
                x1 + (perp_x * offset) / scale,
                y1 - (perp_y * offset) / scale,
                scale,
                canvas_offset_x,
                canvas_offset_y
            )
            end_x, end_y = to_canvas_coords(
                x2 + (perp_x * offset) / scale,
                y2 - (perp_y * offset) / scale,
                scale,
                canvas_offset_x,
                canvas_offset_y
            )

            dwg.add(dwg.line(
                start=(start_x, start_y),
                end=(end_x, end_y),
                stroke="black"
            ))

    for node in G.nodes():
        x, y = positions[node]
        neighbors = list(G.neighbors(node))

        bond_angles = []
        for neighbor in neighbors:
            neighbor_x, neighbor_y = positions[neighbor]
            angle = np.arctan2(neighbor_y - y, neighbor_x - x)
            bond_angles.append(angle)

        center_x, center_y = to_canvas_coords(
            x, y, scale, canvas_offset_x, canvas_offset_y
        )

        dwg.add(dwg.circle(
            center=(center_x, center_y),
            r=20,
            fill="white",
            stroke="black"
        ))

        lone_pairs = G.nodes[node].get("lone_pairs", 0)
        pair_angles = get_lone_pair_angles(bond_angles, lone_pairs, len(neighbors))

        for i in range(min(lone_pairs, len(pair_angles))):
            angle = pair_angles[i]
            distance = 32
            dot_separation = 4

            offset_x = -np.sin(angle) * dot_separation
            offset_y = np.cos(angle) * dot_separation

            for dot_direction in [-1, 1]:
                dot_x = x * scale + canvas_offset_x + distance * np.cos(angle) + offset_x * dot_direction
                dot_y = -y * scale + canvas_offset_y - distance * np.sin(angle) - offset_y * dot_direction

                dwg.add(dwg.circle(
                    center=(dot_x, dot_y),
                    r=2,
                    fill="black"
                ))

        dwg.add(dwg.text(
            G.nodes[node]["element"],
            insert=(center_x - 8, center_y + 5),
            fill="black",
            font_size="16px",
            font_family="Arial"
        ))

    dwg.save()

def main(formula):
    parsed = parse_formula(formula)

    valid, message = validate_molecule(parsed)

    if not valid:
        return message, None, None

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
    draw_molecule(G, positions, formula)
    return parsed, central, G