import numpy as np
import ElementRepository as ER


def KeCalc(Points, planeType, E, v, t, ID = None):
    if planeType == "PlaneStrain":
        C = ER.PlaneStrain(E, v, t).Array
    elif planeType == "PlaneStress":
        C = ER.PlaneStress(E, v, t).Array
    
    if len(Points) == 4:
        element = ER.q4(Points, ID = ID)
    elif len(Points) == 8:
        element = ER.q8(Points, ID = ID)
    elif len(Points) == 7:
        element = ER.q7(Points, ID = ID)
    else:
        raise ValueError("Element type not recognized for Ke Calculation")
    
    xiIntegrationPoints = element.xiIntegrationPoints
    etaIntegrationPoints = element.etaIntegrationPoints
    
    calc = ER.qCalc(element)    
    Ke = np.zeros([element.totalDof, element.totalDof]) # Blank Ke matrix
    # Gauss Integration
    for i, xi in enumerate(xiIntegrationPoints):
        for j, eta in enumerate(etaIntegrationPoints):
            jacb = calc.jacobian(element, xi, eta)
            print("Jacobian of element ", ID, "is: \n", jacb.J)
            print("det of the Jacobnian for element ", ID, " is= ", jacb.det)
            eB1 = calc.B1() # strain displacement relation, relating the derivatives of displacement with resepct to the dofs to the components of strain
            eB2 = calc.B2(jacb) # Scalling matrix containing the inverse jacobian transposed
            eB3 = calc.B3(xi, eta) # Matrix of derivates of the shape functions with respect to xi and eta
            B = eB1 @ eB2 @ eB3 # matrix mulptiplication into the total B matrix
            Ke += B.transpose() @ C @ B * jacb.det * element.Weights[i] * element.Weights[j] # Gauss Integration Statement
    
    Ke = Ke * t
      
    return Ke 


'''
START HERE TO CONFRIM THAT THE BELOW IS CORRECT AND AND COMMENTS!!!
'''
def globalKCalc(KGlobal, dfEle, dfNodes, dfMatProps, elementType):
    for index, row in dfEle.iterrows():
        KGTemp = np.zeros_like(KGlobal) # making a blank array of matching size to the gloabl stiffness matrix
        if elementType == "CQ4":
            points = np.zeros((4, 3))
            Ns = [0, 1, 2, 3]
        elif elementType == "CQ8":
            points = np.zeros((8, 3))
            Ns = [0, 2, 4, 6, 1, 3, 5, 7] # Ordered C1, C2, C3, C4, M1, M2, M3, M4
        elif elementType == "CQ7":
            points = np.zeros((7, 3))
            Ns = [0, 1, 2, 3, 4, 5, 6] # Ordered C1, C2, C3, C4, M1, M2, M3 NO M4!
        else:
            raise ValueError("Element type not recognized")
        #Getting an array of xyz values for all the nodes in the order defined above
        for i in Ns:
            nID = row[f"N{i+1}"]
            nRow = dfNodes[dfNodes["N"] == nID]
            points[i, :] = nRow["XYZ"].values[0]
        E = dfMatProps[dfMatProps["PID"] == row["Prop"]]["E"].values[0]
        v = dfMatProps[dfMatProps["PID"] == row["Prop"]]["NU"].values[0]
        t = dfMatProps[dfMatProps["PID"] == row["Prop"]]["T"].values[0]
        points = points
        Ke = KeCalc(points, "PlaneStress", E, v, t, ID = row["Enumber"]) # Executing the elemental stiffness calc function 
        # print(Ke)
        # Column headers used to match with dfEle dataframe
        if elementType == "CQ4":
            Ns = ["N1", "N2", "N3", "N4"]
        elif elementType == "CQ8":   
            Ns = ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"]
        elif elementType == "CQ7":
            Ns = ["N1", "N2", "N3", "N4", "N5", "N6", "N7"]
        
        # NsPrime = [Ns[x-1] for x in nOrder]s
        dof_per_node = 2
        glob_indexes = []
        for n in Ns: 
            nID = row[n] # Getting the node ID
            base = int(dfNodes.index[dfNodes["N"] == nID][0])*dof_per_node # Converting node ID into a KGlobal index location
            glob_indexes.extend([base + i for i in range(dof_per_node)]) # expanding for the DOFs
        for a, A in enumerate(glob_indexes):
            for b, B in enumerate(glob_indexes):
                # if B < A: # Filling only the upper right trianlge of the matrix
                #     pass
                # else:
                #     KGlobal[A, B] += Ke[a, b]
                #     # KGTemp[A, B] += Ke[a, b] # Used for debugging 
                KGlobal[A, B] += Ke[a, b]
                KGTemp[A, B] += Ke[a, b]
        print("Done Ke Calc")

    # KGlobal = KGlobal + KGlobal.T - np.diag(KGlobal.diagonal()) # Used for debugging populates the entire matrix
    return KGlobal