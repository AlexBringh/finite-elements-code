import numpy as np

# Material properties
E = 210e9         # Elastic modulus [Pa]
nu = 0.3          # Poisson's ratio
sigma_y = 250e6   # Initial yield stress [Pa]
H = 1e9           # Isotropic hardening modulus [Pa]

# Geometry
L = 1.0
H_elem = 1.0

# Nodes: square element
nodes = np.array([
    [0.0, 0.0],   # Node 0
    [L, 0.0],     # Node 1
    [L, H_elem],  # Node 2
    [0.0, H_elem] # Node 3
])

elements = [[0, 1, 2, 3]]

# Gauss points and weights for 2x2 integration
gauss_pts = [(-1/np.sqrt(3), -1/np.sqrt(3)),
             ( 1/np.sqrt(3), -1/np.sqrt(3)),
             ( 1/np.sqrt(3),  1/np.sqrt(3)),
             (-1/np.sqrt(3),  1/np.sqrt(3))]

weights = [1.0, 1.0, 1.0, 1.0]

# Strain-displacement matrix B for Q4
def B_matrix(xi, eta, coords):
    dN_dxi = 0.25 * np.array([
        [-(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta)],
        [-(1 - xi),  -(1 + xi),   (1 + xi),   (1 - xi)]
    ])
    
    J = dN_dxi @ coords
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)
    
    dN_dx = invJ @ dN_dxi

    B = np.zeros((3, 8))
    for i in range(4):
        B[0, 2*i]     = dN_dx[0, i]
        B[1, 2*i+1]   = dN_dx[1, i]
        B[2, 2*i]     = dN_dx[1, i]
        B[2, 2*i+1]   = dN_dx[0, i]
    return B, detJ

# Constitutive law: Von Mises with isotropic hardening
def constitutive_update(strain, ep_prev, alpha_prev):
    # Elastic stiffness matrix (plane strain)
    C = E / ((1 + nu)*(1 - 2*nu)) * np.array([
        [1 - nu,     nu,       0],
        [nu,       1 - nu,     0],
        [0,         0,   0.5 - nu]
    ])

    trial_stress = C @ (strain - ep_prev)
    s_trial = trial_stress - np.mean(trial_stress[:2]) * np.array([1, 1, 0])
    seq = np.sqrt(1.5 * np.dot(s_trial, s_trial))
    f_trial = seq - (sigma_y + H * alpha_prev)

    if f_trial <= 0:
        # Elastic step
        return trial_stress, ep_prev, alpha_prev, C
    else:
        # Plastic correction
        dgamma = f_trial / (1.5 * np.dot(s_trial, C @ s_trial) + H)
        n = s_trial / seq
        delta_ep = dgamma * 1.5 * n
        stress = trial_stress - C @ delta_ep
        ep_new = ep_prev + delta_ep
        alpha_new = alpha_prev + dgamma
        # Consistent tangent modulus (optional: full derivation)
        return stress, ep_new, alpha_new, C  # Simplified: use elastic C for now

# DOF mapping
n_nodes = nodes.shape[0]
n_dofs = n_nodes * 2
u = np.zeros(n_dofs)

# History variables at 4 Gauss points
plastic_strain = np.zeros((4, 3))   # 3 components of strain tensor
alpha = np.zeros(4)                 # isotropic hardening var

# Prescribed vertical displacement
u_top = -0.01  # 1 cm downward

# Boundary conditions: bottom nodes fixed, top nodes moved
bc = {0: 0.0, 1: 0.0, 6: 0.0, 7: u_top}

# Newton-Raphson loop
for it in range(20):
    K = np.zeros((n_dofs, n_dofs))
    R = np.zeros(n_dofs)

    for e in elements:
        coords = nodes[e]
        Ke = np.zeros((8, 8))
        Re = np.zeros(8)

        for gp, w in zip(gauss_pts, weights):
            xi, eta = gp
            print(f"Gauss point: {gp}")
            B, detJ = B_matrix(xi, eta, coords)
            gp_id = gauss_pts.index(gp)

            strain = B @ u[np.ix_([2*i for i in e] + [2*i+1 for i in e])]
            stress, ep, a, D = constitutive_update(strain, plastic_strain[gp_id], alpha[gp_id])
            plastic_strain[gp_id] = ep
            alpha[gp_id] = a

            Ke += B.T @ D @ B * detJ * w
            Re += B.T @ stress * detJ * w

        edofs = []
        for i in e:
            edofs += [2*i, 2*i+1]
        K[np.ix_(edofs, edofs)] += Ke
        R[edofs] += Re

    # External force is zero; residual is -R
    residual = -R

    # Apply BCs
    K_mod = K.copy()
    res_mod = residual.copy()
    for dof, val in bc.items():
        K_mod[dof, :] = 0
        K_mod[:, dof] = 0
        K_mod[dof, dof] = 1
        res_mod[dof] = val - u[dof]

    du = np.linalg.solve(K_mod, res_mod)
    u += du

    norm_res = np.linalg.norm(res_mod)
    print(f"Iteration {it+1}: Residual = {norm_res:.3e}")
    if norm_res < 1e-6:
        break

# Print results
print("\nDisplacements (u):")
print(u.reshape(-1, 2))
