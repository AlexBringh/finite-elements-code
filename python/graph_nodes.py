import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib

nodes = [
    (0,  0.0, 0.0),
    (1,  0.0, 0.5),
    (2,  0.0, 1.0),
    (3,  0.5, 0.0),
    (4,  0.5, 0.5),
    (5,  0.5, 1.0),
    (6,  1.0, 0.0),
    (7,  1.0, 0.5),
    (8,  1.0, 1.0),
    (9,  1.5, 0.0),
    (10, 1.5, 0.5),
    (11, 1.5, 1.0),
    (12, 2.0, 0.0),
    (13, 2.0, 0.5),
    (14, 2.0, 1.0)
]

# Define connections (pairs of node indices)
connections = [
    (0, 1), (1, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (3, 6), (4, 7), (5, 8), (6, 7), (7, 8), (6, 9), (7, 10), (8, 11), (9, 10), (10, 11), (9, 12), (10, 13), (11, 14), (12, 13), (13, 14)  
]

# Extract coordinates for plotting
x = [n[1] for n in nodes]
y = [n[2] for n in nodes]

plt.figure(figsize=(6, 3))
plt.scatter(x, y, color='red')

# Draw connections
for start, end in connections:
    x_coords = [nodes[start][1], nodes[end][1]]
    y_coords = [nodes[start][2], nodes[end][2]]
    plt.plot(x_coords, y_coords, 'b-')


# Annotate node numbers
for nodeid, xi, yi in nodes:
    plt.text(xi, yi, str(nodeid), fontsize=12, ha='right', va='bottom')


plt.gca().set_aspect('equal')
plt.xlim(-0.5, 2.5)
plt.ylim(-0.5, 1.5)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Quad Elements Nodes and Connections')
plt.grid(False)
plt.show()


ux = np.array([0.00, 0.00, 0.00, 0.000693924, 0.000085971, -0.000506973, 0.001179199, 0.000185465, -0.000805690, 0.001434859, 0.000222447, -0.000790254, 0.001637106, 0.000234011, -0.000406566]) # x -  Displacement results
uy = np.array([0.00, 0.00, 0.00, 0.000692527, 0.000510651, 0.000604994, 0.001716460, 0.001627596, 0.001647243, 0.002969674, 0.002998450, 0.002976607, 0.004721342, 0.004284395, 0.004036121]) # y - Displacement results

# Plot deformed mesh in a new figure
x_def = [n[1] + ux[n[0]] for n in nodes]
y_def = [n[2] + uy[n[0]] for n in nodes]

plt.figure(figsize=(6, 3))
plt.scatter(x_def, y_def, color='green')

# Draw connections for deformed mesh
for start, end in connections:
    x_coords = [x_def[start], x_def[end]]
    y_coords = [y_def[start], y_def[end]]
    plt.plot(x_coords, y_coords, 'g-')

# Annotate node numbers
for nodeid, xi, yi in zip(range(len(nodes)), x_def, y_def):
    plt.text(xi, yi, str(nodeid), fontsize=12, ha='right', va='bottom')

plt.gca().set_aspect('equal')
plt.xlim(-0.5, 2.5)
plt.ylim(-0.5, 1.5)
plt.xlabel('X (deformed)')
plt.ylabel('Y (deformed)')
plt.title('2D Quad Elements Nodes and Connections (Deformed)')
plt.grid(False)
plt.show()

# Plot deformed mesh in a new figure (scaled and unscaled, side by side)
fig, axes = plt.subplots(2, 1, figsize=(6, 6))

# Unscaled deformed mesh
axes[0].scatter(x_def, y_def, color='green')
for start, end in connections:
    axes[0].plot([x_def[start], x_def[end]], [y_def[start], y_def[end]], 'g-')
for nodeid, xi, yi in zip(range(len(nodes)), x_def, y_def):
    axes[0].text(xi, yi, str(nodeid), fontsize=12, ha='right', va='bottom')
axes[0].set_aspect('equal')
axes[0].set_xlim(-0.5, 2.5)
axes[0].set_ylim(-0.5, 1.5)
axes[0].set_xlabel('X (deformed)')
axes[0].set_ylabel('Y (deformed)')
axes[0].set_title('Deformed (scale=1)')
axes[0].grid(True)

# Scaled deformed mesh (deformation * 10)
x_def_scaled = [n[1] + 20*ux[n[0]] for n in nodes]
y_def_scaled = [n[2] + 20*uy[n[0]] for n in nodes]
axes[1].scatter(x_def_scaled, y_def_scaled, color='orange')
for start, end in connections:
    axes[1].plot([x_def_scaled[start], x_def_scaled[end]], [y_def_scaled[start], y_def_scaled[end]], 'orange')
for nodeid, xi, yi in zip(range(len(nodes)), x_def_scaled, y_def_scaled):
    axes[1].text(xi, yi, str(nodeid), fontsize=12, ha='right', va='bottom')
axes[1].set_aspect('equal')
axes[1].set_xlim(-0.5, 2.5)
axes[1].set_ylim(-0.5, 1.5)
axes[1].set_xlabel('X (deformed, x20)')
axes[1].set_ylabel('Y (deformed, x20)')
axes[1].set_title('Deformed (scale=20)')
axes[1].grid(True)

plt.tight_layout()
plt.show()