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

def build_graph(chemDict, central):
    # https://networkx.org/documentation/stable/tutorial.html
    G = nx.Graph()

    G.add_node(central, element = central, valence = elements[central]["valence_electrons"])

    for atom, count in chemDict.items():
        if atom == central:
            continue

        for i in range(count):
            name = f"{atom}{i+1}"

            G.add_node(name, element = atom, valence = elements[atom]["valence_electrons"])
            G.add_edge(central, name, order = 1)
    return G

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

def draw_molecule(G, positions):
    position = {node: tuple(positions[node])for node in positions}
    dwg = svgwrite.Drawing("molecule.svg", size=("400px", "400px"))

    scale = 150
    offset_x = 200
    offset_y = 200

    # Draw bonds
    for u, v in G.edges():
        x1, y1 = position[u]
        x2, y2 = position[v]

        dwg.add(dwg.line(
            start=(x1 * scale + offset_x, -y1 * scale + offset_y),
            end=(x2 * scale + offset_x, -y2 * scale + offset_y),
            stroke="black"
        ))

    # Draw atoms
    for node in G.nodes():
        x, y = position[node]

        dwg.add(dwg.circle(
            center=(x * scale + offset_x, -y * scale + offset_y),
            r=20,
            fill="white",
            stroke="black"
        ))

        dwg.add(dwg.text(
            node,
            insert=(x * scale + offset_x - 8, -y * scale + offset_y + 5)
        ))

    dwg.save()

formula = "CH4"

parsed = parse_formula(formula)
central = find_central_atom(parsed)
G = build_graph(parsed, central)
positions = generate_positions(G, central)
draw_molecule(G, positions)