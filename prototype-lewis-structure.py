import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

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
bond_angle = 180
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

# Render
nx.draw(G, position, with_labels=True, node_size=2000)
plt.title("H2O Lewis Structure Prototype")
plt.axis('equal')
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.show()