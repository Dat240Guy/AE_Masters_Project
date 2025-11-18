import pandas as pd
import numpy as np
import re
from KeCalc import globalKCalc  

def DispCalc(dfNodes, dfEle4, dfEle8, dfEle7, dfForces, dfConstraints, dfMatProps):
    dof = len(dfNodes)*2 #Total degrees of freedom assuming x, y displacement at each node
    KGlobal = np.zeros((dof, dof)) # Blank global stiffness matrix
    uGlobal = np.zeros((dof, 1))+1 # Blank displacement vector
    fGlobal = np.zeros((dof, 1)) # Blank force vector

    # Generating force vector
    for index, row in dfForces.iterrows():
        nID = row["NID"]
        dofId = int(dfNodes.index[dfNodes["N"] == nID][0])*2 # multiply by 2 for 2 degrees of freedom
        fGlobal[dofId, 0] = row["Factor"] * row["F1"] # x force
        fGlobal[dofId+1, 0] = row["Factor"] * row["F2"] # y force

    # Generating Constrain vector
    for index, row in dfConstraints.iterrows():
        nID = row["NID"]
        dofId = int(dfNodes.index[dfNodes["N"] == nID][0])*2
        cDOF = str(row["Comp"])
        # dofId = dfNodes.index[dfNodes["N"] == nID]*2
        for char in range(len(cDOF)): # finding the x and y constraints only for the 2d plate situation 
            if cDOF[char] == "1":
                uGlobal[dofId, 0] = 0
            elif cDOF[char] == "2":
                uGlobal[dofId+1, 0] = 0
            else:
                pass

    #Calculate the Global Stiffness Matrix
    if isinstance(dfEle4, pd.DataFrame):
        KGlobal = globalKCalc(KGlobal, dfEle4, dfNodes, dfMatProps, "CQ4")
    if isinstance(dfEle8, pd.DataFrame):
        KGlobal = globalKCalc(KGlobal, dfEle8, dfNodes, dfMatProps, "CQ8")
    if isinstance(dfEle7, pd.DataFrame):
        KGlobal = globalKCalc(KGlobal, dfEle7, dfNodes, dfMatProps, "CQ7")
    # KGlobal = KGlobal + KGlobal.T - np.diag(KGlobal.diagonal())
    
    #Applying boundary Conditions

    #Solving for displacements
    print("Debug Point")
    # Zeroing rows and cols of K and f where the dof is fixed (constrained)
    for u in range(len(uGlobal)):
        if uGlobal[u] == 0:
            KGlobal[u, :] = 0
            KGlobal[:, u] = 0
            fGlobal[u] = 0

    # Removing fixed rows and cols 
    KGlobal_clean = KGlobal[~np.all(KGlobal == 0, axis=1)][:, ~np.all(KGlobal == 0, axis=0)]
    fGlobalClean = fGlobal[~np.all(KGlobal == 0, axis=1)]
    
    # Solving for global displacements
    uGlob = np.linalg.solve(KGlobal_clean, fGlobalClean)
    
    # Creating results dataframes    
    j = 0
    for i in range(len(uGlobal)):
        if uGlobal[i] != 0:
            uGlobal[i] = uGlob[j]
            j += 1
        else:
            pass
        
    results = []
    for index, row in dfNodes.iterrows():
        nID = row["N"]
        results.append([nID, uGlobal[index*2, 0], uGlobal[index*2+1, 0]])
        
    dfDisp = pd.DataFrame(results, columns = ["NID", "U1", "U2"])
    dfDisp = dfDisp.astype({"NID":int, "U1":float, "U2":float}).sort_values(by="NID")
    return dfDisp

import numpy as np
import pandas as pd

def build_deformed_nodes(dfNodes, dfDisp, scale=1.0):
    """
    dfNodes : DataFrame with columns:
              N  - node ID (int)
              XYZ - np.array([x,y,z])
    dfDisp  : DataFrame with columns:
              NID - node ID (int)
              U1, U2[, U3] - displacements in global x,y[,z]
    scale   : visualization scale factor for displacements

    Returns:
        dfDef : DataFrame like dfNodes but with an extra column 'XYZ_def'
                containing the displaced coordinates.
    """
    # Merge nodes with displacements
    df = dfNodes.merge(dfDisp, how="left", left_on="N", right_on="NID")

    # Assume 2D plate: z displacement = 0 if not present
    def make_xyz_def(row):
        x, y, z = row["XYZ"]
        u1 = row.get("U1", 0.0)
        u2 = row.get("U2", 0.0)
        u3 = row.get("U3", 0.0) if "U3" in row else 0.0
        return np.array([x + scale*u1, y + scale*u2, z + scale*u3])

    df["XYZ_def"] = df.apply(make_xyz_def, axis=1)

    # Keep same columns as dfNodes plus XYZ_def
    dfDef = df[["N", "XYZ", "XYZ_def"]].copy()
    return dfDef


import matplotlib.pyplot as plt
import matplotlib.tri as tri

def plot_deformed_mesh(dfDefNodes, dfEles, show_undeformed=True):
    # Use XYZ_def for deformed coordinates
    coords_def = np.array(dfDefNodes["XYZ_def"].tolist())[:, :2]
    nids = dfDefNodes["N"].astype(int).to_numpy()
    nid_to_idx = {nid: i for i, nid in enumerate(nids)}

    # Build triangles from connectivity (same logic as in ContourMesh)
    tris = []
    for dfE in dfEles:
        node_cols = [c for c in dfE.columns if c.startswith("N")]
        for _, row in dfE[node_cols].iterrows():
            nodes_all = [int(row[c]) for c in node_cols if int(row[c]) != 0]
            nnode = len(nodes_all)

            if nnode in (3, 6):       # tri / ctria6
                corners = nodes_all[:3]
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3 = (nid_to_idx[n] for n in corners)
                tris.append([i1, i2, i3])

            elif nnode in (4, 7, 8):  # quad family
                corners = nodes_all[:4]
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3, i4 = (nid_to_idx[n] for n in corners)
                tris.append([i1, i2, i3])
                tris.append([i1, i3, i4])

    tris = np.array(tris, dtype=int)

    tri_def = tri.Triangulation(coords_def[:,0], coords_def[:,1], triangles=tris)

    plt.figure()
    if show_undeformed:
        coords_und = np.array(dfDefNodes["XYZ"].tolist())[:, :2]
        tri_und = tri.Triangulation(coords_und[:,0], coords_und[:,1], triangles=tris)
        plt.triplot(tri_und, linewidth=0.5, alpha=0.4, label="undeformed")

    plt.triplot(tri_def, linewidth=0.8, label="deformed")
    plt.axis("equal")
    plt.legend()
    plt.title("Undeformed vs Deformed Mesh")
    plt.show()
