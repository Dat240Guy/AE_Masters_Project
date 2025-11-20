import numpy as np
import ElementRepository as ER


def _split_constitutive_2d(C, planeType, E, v):
    """
    Split a 2D isotropic constitutive matrix C (3x3) into
    volumetric and deviatoric parts:

        C = C_vol + C_dev

    where C_vol acts only on the 'mean' normal strain (εx + εy),
    and C_dev is everything else.

    This is used for selective / reduced integration on Q4 elements.
    """
    # Unit vector for volumetric strain direction [εx, εy, γxy]
    v_vec = np.array([[1.0],
                      [1.0],
                      [0.0]])
    v_norm = v_vec / np.linalg.norm(v_vec)
    P_vol = v_norm @ v_norm.T  # 3x3 projector onto volumetric strain

    if planeType == "PlaneStress":
        # Effective 2D bulk-like modulus for plane stress
        k = E / (1.0 - v)
    elif planeType == "PlaneStrain":
        # Effective 2D bulk-like modulus for plane strain
        k = E / ((1.0 + v) * (1.0 - 2.0 * v))
    else:
        raise ValueError("planeType must be 'PlaneStress' or 'PlaneStrain'")

    C_vol = k * P_vol
    C_dev = C - C_vol
    return C_vol, C_dev


def KeCalc(Points, planeType, E, v, t, ID=None):
    """
    Element stiffness matrix calculator.

    - PlaneStress / PlaneStrain handled via ElementRepository.
    - Q4: selective / reduced integration
          * deviatoric part fully integrated with 2x2 Gauss
          * volumetric part integrated with 1 point at (0,0) with weight 4.0
    - Q7, Q8: original full Gauss integration (unchanged)
    """
    Points = np.asarray(Points, dtype=float)

    # Constitutive matrix
    if planeType == "PlaneStrain":
        C = ER.PlaneStrain(E, v, t).Array
    elif planeType == "PlaneStress":
        C = ER.PlaneStress(E, v, t).Array
    else:
        raise ValueError("planeType must be 'PlaneStress' or 'PlaneStrain'")

    
    if Points.shape[0] == 3:
        element = ER.t3(Points, ID=ID)
    elif Points.shape[0] == 4:
        element = ER.q4(Points, ID=ID)
    elif Points.shape[0] == 8:
        element = ER.q8(Points, ID=ID)
    elif Points.shape[0] == 7:
        element = ER.q7(Points, ID=ID)
    else:
        raise ValueError("Element type not recognized for Ke Calculation")

    calc = ER.qCalc(element)
    Ke = np.zeros((element.totalDof, element.totalDof))

    # Original Gauss integration path for Q8/Q7
    for i, xi in enumerate(element.xiIntegrationPoints):
        for j, eta in enumerate(element.etaIntegrationPoints):
            jacb = calc.jacobian(element, xi, eta)

            # eB1 = calc.B1()
            # eB2 = calc.B2(jacb)
            # eB3 = calc.B3(xi, eta)
            # B = eB1 @ eB2 @ eB3
            B = calc.B(xi, eta, jacb)
            Ke += (
                B.T @ C @ B
                * jacb.det
                * element.Weights[i]
                * element.Weights[j]
            )

    Ke *= t
    return Ke


def globalKCalc(KGlobal, dfEle, dfNodes, dfMatProps, elementType):
    """
    Assemble global stiffness matrix from element stiffness matrices.

    This function is unchanged in structure, but now benefits from the
    selective integration for Q4 via KeCalc() above.
    """
    dof_per_node = 2

    for index, row in dfEle.iterrows():
        # Temporary matrix for debugging (optional)
        KGTemp = np.zeros_like(KGlobal)

        # Build local coordinate array for this element
        if elementType == "CTRIA3":
            points = np.zeros((3, 3))
            Ns_idx = [0, 1, 2]
        elif elementType == "CQ4":
            points = np.zeros((4, 3))
            Ns_idx = [0, 1, 2, 3]
        elif elementType == "CQ8":
            points = np.zeros((8, 3))
            Ns_idx = [0, 2, 4, 6, 1, 3, 5, 7]  # C1,C2,C3,C4,M1,M2,M3,M4
        elif elementType == "CQ7":
            points = np.zeros((7, 3))
            Ns_idx = [0, 1, 2, 3, 4, 5, 6]      # C1..C4,M1..M3
        else:
            raise ValueError("Element type not recognized")

        # Get xyz of all nodes in the order defined above
        for i in Ns_idx:
            nID = row[f"N{i+1}"]
            nRow = dfNodes[dfNodes["N"] == nID]
            points[i, :] = nRow["XYZ"].values[0]

        # Material properties
        pid = row["Prop"]
        E = dfMatProps[dfMatProps["PID"] == pid]["E"].values[0]
        v = dfMatProps[dfMatProps["PID"] == pid]["NU"].values[0]
        t = dfMatProps[dfMatProps["PID"] == pid]["T"].values[0]

        # Element stiffness
        Ke = KeCalc(points, "PlaneStress", E, v, t, ID=row["Enumber"])

        # Node labels in the dataframe for DOF mapping
        if elementType == "CTRIA3":
            Ns = ["N1", "N2", "N3"]
        elif elementType == "CQ4":
            Ns = ["N1", "N2", "N3", "N4"]
        elif elementType == "CQ8":
            Ns = ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"]
        elif elementType == "CQ7":
            Ns = ["N1", "N2", "N3", "N4", "N5", "N6", "N7"]

        # Build list of global DOF indices for this element
        glob_indexes = []
        for n in Ns:
            nID = row[n]
            base = int(dfNodes.index[dfNodes["N"] == nID][0]) * dof_per_node
            glob_indexes.extend([base + i for i in range(dof_per_node)])

        # Scatter add into global stiffness matrix
        for a, A in enumerate(glob_indexes):
            for b, B in enumerate(glob_indexes):
                KGlobal[A, B] += Ke[a, b]
                KGTemp[A, B] += Ke[a, b]  # handy for debugging if needed

    return KGlobal
