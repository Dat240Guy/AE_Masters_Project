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
        fGlobal[dofId, 0] = row["F1"] # x force
        fGlobal[dofId+1, 0] = row["F2"] # y force

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