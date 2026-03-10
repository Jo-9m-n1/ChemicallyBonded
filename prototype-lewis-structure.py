import networkx as nx
import numpy as np
import svgwrite

# https://networkx.org/documentation/stable/tutorial.html
G = nx.Graph()

# this makes atoms into nodes (later will be connected)
G.add_node("O", element="O", valence=6, lone_pairs=2)
G.add_node("H1", element="H", valence=1)
G.add_node("H2", element="H", valence=1)

# make edges to connect nodes
G.add_edge("O", "H1", order=1)
G.add_edge("O", "H2", order=1)

# this is a basic bond angle with bond length for matplotlib
bond_angle = 109.5
bond_length = 1.0

positions = {}

# this makes oxygen the origin. it lets me place it in matplotlib
positions["O"] = np.array([0.0, 0.0])

# this makes bond angles into radians.
# https://numpy.org/doc/2.2/reference/generated/numpy.radians.html
# we are splitting the total angle in 2 so that both hydrogen molecules
# are symmetric.
theta = np.radians(bond_angle / 2)

# these puts the hydrogen atoms in place
# it is a basic circle that we are then doing the equation r * cos(angle)
# and r * sin(angle)
positions["H1"] = np.array([
    bond_length * np.cos(theta),
    bond_length * np.sin(theta)
])

positions["H2"] = np.array([
    bond_length * np.cos(theta),
    # flip the y coordinate, so H2 can be on the bottom instead
    -bond_length * np.sin(theta)
])

# converts Numpy arrays to tuples so NetworkX can draw them
position = {}
for node in positions:
    position[node] = tuple(positions[node])

print(positions)

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