import networkx as nx
import numpy as np
import svgwrite
import json

with open("periodic_table.json") as file:
    elements = json.load(file)

def parse_formula(formula):
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
    return elements[atom]["valence_electrons"]

def find_central_atom(chemDict):
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
    for node in G.nodes():
        element = G.nodes[node]["element"]

        if element == "H":
            G.nodes[node]["lone_pairs"] = 0
            continue

        electrons = G.nodes[node]["lone_electrons"]
        G.nodes[node]["lone_pairs"] = electrons // 2

def build_graph(chemDict, central):
    # https://networkx.org/documentation/stable/tutorial.html
    G = nx.Graph()

    G.add_node(central, element = central, valence = get_valence_electrons(central), lone_electrons = get_valence_electrons(central), bonds = 0)

    for atom, count in chemDict.items():
        if atom == central:
            continue

        for i in range(count):
            name = f"{atom}{i+1}"

            G.add_node(name, element = atom, valence = get_valence_electrons(atom), lone_electrons = get_valence_electrons(atom), bonds = 0)
            G.add_edge(central, name, order = 1)

            G.nodes[central]["lone_electrons"] -= 1
            G.nodes[name]["lone_electrons"] -= 1

            G.nodes[central]["bonds"] += 1
            G.nodes[name]["bonds"] += 1
    return G

def get_max_bonds(element):
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
    return abs((a - b + np.pi) % (2 * np.pi) - np.pi)

def draw_molecule(G, positions):
    position = {node: tuple(positions[node])for node in positions}
    dwg = svgwrite.Drawing("molecule.svg", size=("400px", "400px"))

    scale = 150
    offset_x = 200
    offset_y = 200

    # Draw bonds
    for u, v, data in G.edges(data = True):
        order = data.get("order", 1)

        x1, y1 = position[u]
        x2, y2 = position[v]

        # direction vector
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)

        # perpendicular vector
        px = -dy / length
        py = dx / length

        for i in range(order):
            offset = (i - (order - 1) / 2) * 10

            dwg.add(dwg.line(
                start=(
                    x1 * scale + offset_x + px * offset,
                    -y1 * scale + offset_y + py * offset
                ),
                end=(
                    x2 * scale + offset_x + px * offset,
                    -y2 * scale + offset_y + py * offset
                ),
                stroke="black"
            ))

        print(u, "-", v, "order:", data["order"])

    # Draw atoms
    for node in G.nodes():
        x, y = position[node]

        neighbors = list(G.neighbors(node))

        bond_angles = []
        for neighbor in neighbors:
            nx, ny = position[neighbor]
            angle = np.arctan2(ny - y, nx - x)
            bond_angles.append(angle)
        
        all_angles = [i * (np.pi / 2) for i in range(4)]  # 0°, 90°, 180°, 270°

        free_angles = []

        for a in all_angles:
            if all(angle_diff(a, b) > 0.5 for b in bond_angles):
                free_angles.append(a)

        dwg.add(dwg.circle(
            center=(x * scale + offset_x, -y * scale + offset_y),
            r=20,
            fill="white",
            stroke="black"
        ))

        lone_pairs = G.nodes[node].get("lone_pairs", 0)
        for i in range(min(lone_pairs, len(free_angles))):
            angle = free_angles[i]

            dx = 30 * np.cos(angle)
            dy = 30 * np.sin(angle)

            # draw 2 dots per pair
            dwg.add(dwg.circle(
                center=(x * scale + offset_x + dx - 3,
                        -y * scale + offset_y + dy),
                r=2,
                fill="black"
            ))

            dwg.add(dwg.circle(
                center=(x * scale + offset_x + dx + 3,
                        -y * scale + offset_y + dy),
                r=2,
                fill="black"
            ))

        dwg.add(dwg.text(
            G.nodes[node]["element"],
            insert=(x * scale + offset_x - 8, -y * scale + offset_y + 5)
        ))

    dwg.save()

formula = input("Enter a chemical formula: ")

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