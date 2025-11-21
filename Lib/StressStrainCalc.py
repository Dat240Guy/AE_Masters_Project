import numpy as np
import pandas as pd
import ElementRepository as ER

def EleLocalCoordCalc(Points):
    zDir = np.array([0, 0, 1])
    xDir = (Points[1] - Points[0]) / np.linalg.norm((Points[1] - Points[0]))
    yDir = np.cross(xDir, zDir)
    yDir = yDir/np.linalg.norm(yDir)
    return xDir, yDir, zDir

def RotMatrix(theta, ds):
    c = np.cos(theta)
    s = np.sin(theta)
    if ds == "3d":
        rot =  np.array([[c, -s, 0],
                        [s, c, 0],
                        [0, 0, 1]])
    else: 
        rot =  np.array([[c, -s],
                        [s, c]])
    return rot

def stressRotMatrix(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    c2, s2, cs = c*c, s*s, c*s
    return np.array([
        [ c2,   s2,   2*cs],
        [ s2,   c2,  -2*cs],
        [-cs,   cs,  c2 - s2],
    ])

def maxPrinCalc(df, Type):
    if Type == "Stress":
        Prefix = "S"
    else:
        Prefix = "E"
        
    S1 = df[Prefix + "1"].values
    S2 = df[Prefix + "2"].values
    S12 = df[Prefix + "12"].values   # engineering shear

    # Convert engineering shear → tensor shear τxy
    tau = S12 / 2.0

    # Principal value formula
    avg = 0.5 * (S1 + S2)
    rad = np.sqrt(((S1 - S2) / 2.0)**2 + tau**2)

    df[Prefix + "_max"] = avg + rad
    df[Prefix + "_min"] = avg - rad  # optional, but helpful

    return df

def SSEleCalc(Points, disp, planeType, E, v, t, ID = None):
    if planeType == "PlaneStrain":
        C = ER.PlaneStrain(E, v, t).Array
    elif planeType == "PlaneStress":
        C = ER.PlaneStress(E, v, t).Array
    
    eleX, eleY, eleZ = EleLocalCoordCalc(Points)
    theta = np.arctan2(eleX[1], eleX[0])
    # print("Theta Rotation in Degrees: ", np.degrees(theta))
    PointsLocal = np.empty_like(Points)
    for i, point in enumerate(Points):
        PointsLocal[i, :] = RotMatrix(-theta, "3d") @ point
    theta = 0
    dispLocal = np.empty_like(disp)
    for i in range(int(len(dispLocal)/2)):
        dispLocal[i*2:i*2+2] = RotMatrix(-theta, "2d") @ disp[i*2:i*2+2]
        
    if len(Points) == 3:
        element = ER.t3(PointsLocal, ID = ID)
    elif len(Points) == 4:
        element = ER.q4(PointsLocal, ID = ID)
    elif len(Points) == 8:
        element = ER.q8(PointsLocal, ID = ID)
    elif len(Points) == 7:
        element = ER.q7(PointsLocal, ID = ID)
    elif len(Points) == 6:
        element = ER.q6(PointsLocal, ID = ID)
    else:
        raise ValueError("Element type not recognized for Ke Calculation")
    
    calc = ER.qCalc(element)
    n_nodes = element.nodeCount
    nodalStrain = np.zeros((n_nodes, 3)) #Empty arrays to fill
    nodalStress = np.zeros((n_nodes, 3)) #Empty arrays to fill
    nodalStrainGlobal = np.zeros((n_nodes, 3)) #Empty arrays to fill
    nodalStressGlobal = np.zeros((n_nodes, 3)) #Empty arrays to fill
    
    
    #looping through the xi, eta points in the natrual coord sys to calculate strain/stress 
    # each xi, eta point corresponds to the respective node in the global system
    for  p, (xiCalc, etaCalc) in enumerate(element.localCoord):
        # print("xiCalc", xiCalc, "etaCalc", etaCalc)
        jacb = calc.jacobian(element, xiCalc, etaCalc)
        B = calc.B(xiCalc, etaCalc, jacb)

        eps = (B @ dispLocal).reshape(-1)        # [εx, εy, γxy]
        sig = (C @ (B @ dispLocal)).reshape(-1)  # [σx, σy, τxy]

        nodalStrain[p, :] += eps
        nodalStress[p, :] += sig

        # Rotation to global if you decide to re-enable theta:
        T = stressRotMatrix(-theta)
        nodalStrainGlobal[p, :] += T @ nodalStrain[p, :]
        nodalStressGlobal[p, :] += T @ nodalStress[p, :]
    return nodalStrainGlobal, nodalStressGlobal

def StressStrainCalc(dfEles, eTypes, dfDisp, dfNodes, planeType, dfMatProps):
    # dfStress = pd.DataFrame(columns=["Element", "Node", "S1", "S2", "S12"])
    # dfStrain = pd.DataFrame(columns=["Element", "Node", "E1", "E2", "E12"])
    
    for df, elementType in zip(dfEles, eTypes):
        for index, row in df.iterrows():
            if elementType == "CTRIA3":
                points = np.zeros((3,3))
                Ns = [0, 1, 2]
            elif elementType == "CQ4":
                points = np.zeros((4, 3))
                Ns = [0, 1, 2, 3]
            elif elementType == "CQ8":
                points = np.zeros((8, 3))
                Ns = [0, 2, 4, 6, 1, 3, 5, 7]
            elif elementType == "CQ7":
                points = np.zeros((7, 3))
                Ns = [0, 1, 2, 3, 4, 5, 6]
            elif elementType == "CQ6":
                points = np.zeros((6, 3))
                Ns = [0, 1, 2, 3, 4, 5]
            else:
                raise ValueError("Element type not recognized")
            uElement = np.zeros((len(Ns)*2, 1)) # Blank array for elemental displacements
            for i in Ns:
                nID = row[f"N{i+1}"] #finding the node Id for each node column of the element
                nRow = dfNodes[dfNodes["N"] == nID] # Returning the whole row from dfNodes of the node of interest from the element
                points[i, :] = nRow["XYZ"].values[0] #Getting the YXZ values from dfNodes for the node of interest
                uElement[i*2, 0] = dfDisp[dfDisp["NID"] == nID]["U1"].values[0]
                uElement[i*2+1, 0] = dfDisp[dfDisp["NID"] == nID]["U2"].values[0]
            E = dfMatProps[dfMatProps["PID"] == row["Prop"]]["E"].values[0]
            v = dfMatProps[dfMatProps["PID"] == row["Prop"]]["NU"].values[0]
            t = dfMatProps[dfMatProps["PID"] == row["Prop"]]["T"].values[0]
            if row["Enumber"] == 15484:
                pass
            nStrain, nStress = SSEleCalc(points, uElement, planeType, E, v, t, row["Enumber"])
            for i in range(len(Ns)):
                if i == 0 and eTypes.index(elementType) == 0 and index == 0:
                    dfStress = pd.DataFrame([[row["Enumber"], row[f"N{i+1}"], nStress[i,0], nStress[i,1], nStress[i,2]]], columns=["Element", "NID", "S1", "S2", "S12"])
                    dfStrain = pd.DataFrame([[row["Enumber"], row[f"N{i+1}"], nStrain[i,0], nStrain[i,1], nStrain[i,2]]], columns=["Element", "NID", "E1", "E2", "E12"])
                else:
                    dfnStress = pd.DataFrame([[row["Enumber"], row[f"N{i+1}"], nStress[i,0], nStress[i,1], nStress[i,2]]], columns=["Element", "NID", "S1", "S2", "S12"])
                    dfStress = pd.concat([dfStress, dfnStress], ignore_index=True)
                    dfnStrain = pd.DataFrame([[row["Enumber"], row[f"N{i+1}"], nStrain[i,0], nStrain[i,1], nStrain[i,2]]], columns=["Element", "NID", "E1", "E2", "E12"])
                    dfStrain = pd.concat([dfStrain, dfnStrain], ignore_index=True) 
          
    dfStress = dfStress.astype({"Element":int, "NID":int, "S1":float, "S2":float, "S12":float}).sort_values(by=["Element"])
    dfStrain = dfStrain.astype({"Element":int, "NID":float, "E1":float, "E2":float, "E12":float}).sort_values(by=["Element"])            
    dfStress = maxPrinCalc(dfStress, "Stress") 
    dfStrain = maxPrinCalc(dfStrain, "Strain")
    return dfStress, dfStrain

import numpy as np
import pandas as pd

def rThetaStress(dfNodes, dfStress, hole_center, hole_node_ids):
    """
    Compute σ_rr, σ_θθ, τ_rθ for a set of nodes around a hole.

    Parameters
    ----------
    dfNodes : DataFrame
        Columns: ["N", "XYZ"]
    dfStress : DataFrame
        Columns: ["NID", "S1", "S2", "S12"]  (global σx, σy, τxy)
    hole_center : tuple (xc, yc)
        Coordinates of the hole center
    hole_node_ids : list
        List of node IDs along the hole boundary

    Returns
    -------
    dfOut : DataFrame
        Columns: ["NID", "theta_deg", "sigma_rr", "sigma_tt", "tau_rt"]
        Sorted in increasing theta
    """

    xc, yc = hole_center
    out_rows = []

    for nid in hole_node_ids:
        # Get node coordinates
        node_row = dfNodes[dfNodes["N"] == nid]
        x, y, _ = node_row["XYZ"].values[0]

        # Vector from hole center → node
        dx = x - xc
        dy = y - yc
        r = np.hypot(dx, dy)

        if r == 0:
            continue  # ignore exactly-at-center points

        # Angle of radial vector (global)
        theta = np.arctan2(dy, dx)  # radians
        c = np.cos(theta)
        s = np.sin(theta)

        # Get FE stress at node (global components)
        s_row = dfStress[dfStress["NID"] == nid].iloc[0]
        sx = s_row["S1"]
        sy = s_row["S2"]
        txy = s_row["S12"]

        # Transform stresses to (r,θ)
        sigma_rr = sx*c*c + sy*s*s + 2*txy*s*c
        sigma_tt = sx*s*s + sy*c*c - 2*txy*s*c
        tau_rt   = (sy - sx)*s*c + txy*(c*c - s*s)

        out_rows.append([nid, np.degrees(theta), sigma_rr, sigma_tt, tau_rt])

    dfOut = pd.DataFrame(out_rows, columns=["NID", "theta_deg", "sigma_rr", "sigma_tt", "tau_rt"])
    dfOut = dfOut.sort_values("theta_deg").reset_index(drop=True)
    return dfOut
