import numpy as np
import sys
sys.path.append(r"C:\Documents\Grad_School\!AE_Masters_Project\AE_Masters_Project\Lib")
from KeCalc import KeCalc
from StressStrainCalc import SSEleCalc


def run_q4_patch_test_solve():
    print("\n============================")
    print(" Q4 PATCH TEST (KeCalc + SSEleCalc)")
    print("============================\n")

    # -------------------------------------------------------
    # 1. Distorted Q4 geometry (edit to test other shapes)
    # -------------------------------------------------------
    # Node order: [N1, N2, N3, N4]
    Points = np.array([
        [0.00, 0.00, 0.0],   # N1 - left bottom
        [1.00, 0.00, 0.0],   # N2 - right bottom
        [1.00, 0.25, 0.0],   # N3 - right top
        [0.80, 0.25, 0.0],   # N4 - left top (distorted)
    ])

    # -------------------------------------------------------
    # 2. Material properties & plane type
    # -------------------------------------------------------
    E = 30e6     # psi
    v = 0.3
    t = 0.25     # thickness
    planeType = "PlaneStress"

    # -------------------------------------------------------
    # 3. Element stiffness via your KeCalc (with selective Q4)
    # -------------------------------------------------------
    Ke = KeCalc(Points, planeType, E, v, t, ID=1)  # 8x8 matrix
    ndof = Ke.shape[0]

    # -------------------------------------------------------
    # 4. Build load vector for uniform axial traction on right edge
    # -------------------------------------------------------
    # Desired uniform axial stress sigma_x:
    sigma_x = 30_000.0  # psi

    # Right edge is between node 2 and 3
    P2 = Points[1, :2]  # [x2, y2]
    P3 = Points[2, :2]  # [x3, y3]
    edge_vec = P3 - P2
    edge_length = np.linalg.norm(edge_vec)

    # Uniform traction tx = sigma_x in +x direction, constant along this edge.
    # Equivalent nodal forces for a 2-noded linear edge:
    #   Fx2 = Fx3 = tx * t * L / 2
    Fx_edge = sigma_x * t * edge_length / 2.0

    F = np.zeros((ndof, 1))
    # DOFs: node i -> [2*i, 2*i+1] => [Ux_i, Uy_i]
    # Node 2 = index 1
    F[2 * 1, 0] += Fx_edge  # Ux at node 2
    # Node 3 = index 2
    F[2 * 2, 0] += Fx_edge  # Ux at node 3

    # -------------------------------------------------------
    # 5. Apply boundary conditions
    # -------------------------------------------------------
    # Fix Ux, Uy at nodes 1 and 4 (left edge)
    fixed_dofs = [
        0, 1,   # node 1: Ux, Uy
        6, 7    # node 4: Ux, Uy
    ]

    all_dofs = np.arange(ndof)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    # Partitioned system: solve K_ff * u_f = F_f
    K_ff = Ke[np.ix_(free_dofs, free_dofs)]
    F_f = F[free_dofs, :]

    # Solve for free displacements
    u_f = np.linalg.solve(K_ff, F_f)

    # Build full displacement vector
    u = np.zeros((ndof, 1))
    u[free_dofs, :] = u_f
    # fixed DOFs remain 0

    # -------------------------------------------------------
    # 6. Use SSEleCalc to compute strains/stresses
    # -------------------------------------------------------
    strain, stress = SSEleCalc(Points, u, planeType, E, v, t, ID=1)

    # -------------------------------------------------------
    # 7. Report
    # -------------------------------------------------------
    print("Node-by-node STRAIN (E1, E2, E12):")
    for i in range(4):
        print(f"N{i+1}:", strain[i])

    print("\nNode-by-node STRESS (S1, S2, S12):")
    for i in range(4):
        print(f"N{i+1}:", stress[i])

    s1_vals = stress[:, 0]
    s2_vals = stress[:, 1]
    s12_vals = stress[:, 2]

    print("\nExpected S1 ≈", sigma_x)
    print("\nPatch Test Summary:")
    print("-------------------")
    print(f"Max deviation in S1:  {np.max(np.abs(s1_vals - sigma_x)):.6f}")
    print(f"Max |S2|:             {np.max(np.abs(s2_vals)):.6f}")
    print(f"Max |S12|:            {np.max(np.abs(s12_vals)):.6f}")

    tol = 0.01 * sigma_x   # 1% of S1

    if (np.max(np.abs(s1_vals - sigma_x)) < tol and
        np.max(np.abs(s2_vals)) < tol and
        np.max(np.abs(s12_vals)) < tol):
        print("➡️  PASS: Element reproduces near-constant axial stress.")
    else:
        print("❌ FAIL: Significant parasitic S2/S12 or S1 deviation detected.")


if __name__ == "__main__":
    run_q4_patch_test_solve()
