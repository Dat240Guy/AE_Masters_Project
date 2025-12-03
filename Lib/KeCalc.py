import numpy as np
import ElementRepository as ER



def KeCalc(Points, planeType, E, v, t, ID=None):

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
    elif Points.shape[0] == 6:
        element = ER.q6(Points, ID=ID)
    else:
        raise ValueError("Element type not recognized for Ke Calculation")

    calc = ER.qCalc(element)
    Ke = np.zeros((element.totalDof, element.totalDof))

    for i, xi in enumerate(element.xiIntegrationPoints):
        for j, eta in enumerate(element.etaIntegrationPoints):
            jacb = calc.jacobian(element, xi, eta)

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
        elif elementType == "CQ6":
            points = np.zeros((6, 3))
            Ns_idx = [0, 1, 2, 3, 5, 4]      # C1..C4,M1,M2
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
        elif elementType == "CQ6":
            Ns = ["N1", "N2", "N3", "N4", "N5", "N6"]

        # Build list of global DOF indices for this element
        glob_indexes = []
        for n in Ns:
            nID = row[n]
            base = int(dfNodes.index[dfNodes["N"] == nID][0]) * dof_per_node
            glob_indexes.extend([base + i for i in range(dof_per_node)])

        # Scatter add into global stiffness matrix
        for a, A in enumerate(glob_indexes):
            for b, B in enumerate(glob_indexes):
                if b >= a:
                    KGlobal[A, B] += Ke[a, b]
                    # KGTemp[A, B] += Ke[a, b]  # handy for debugging if needed

    return KGlobal
