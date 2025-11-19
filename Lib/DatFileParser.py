
import pandas as pd
import numpy as np
import re

def nastranFloat(s):
    if re.search(r'(?<=\d)([+-]\d+)', s) != None:
        s = re.sub(r'(?<=\d)([+-]\d+)', r'E\1', s)
        return float(s)
    else:
        return s

def wrap(s, width):
    chunks = [s[i:i+width] for i in range(0, len(s), width)]
    trimmed_chunks = [chunk.strip() for chunk in chunks]
    trimmed_chunks = [nastranFloat(chunk) for chunk in trimmed_chunks if chunk != ""]  # Remove empty strings
    return trimmed_chunks

def transitionEleParsing(dfNodes, dfEle4, dfEle8, e8):
    dfEle7 = None
    dfEle6 = None
    
    # Creating a dicitonary of nodes and the quantity of occurances that nodes has in CQ4 and CQ8 elements
    ndict = {}
    for node in dfNodes.iterrows():
        n4Count = (dfEle4[['N1', 'N2', 'N3', 'N4']] == node[1]['N']).sum().sum()
        n8Count = (dfEle8[['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']] == node[1]['N']).sum().sum()
        ndict[node[1]['N']] = {"n4Count":n4Count, "n8Count":n8Count, "totalCount":n4Count + n8Count}
    
    nCommon = {x: ndict[x] for x in ndict.keys() if ndict[x]["n8Count"] > 0 and ndict[x]["n4Count"] > 0} # Finding the common nodes between cq4 and cq8 elements
    e7Common = [ele for ele in e8 if any(str(n) in ele[3:11] for n in nCommon.keys())] # Finding the cq8 elements with the nCommon nodes
    
    #routine to find the e7 nodes
    e7, nFree = [], []
    for i, e in enumerate(e7Common):
        efull = e
        e = e[3:11]
        freeEdgeNodes, freeEdgeMidNodes = [], []
        #Determining if the element has any free edge where free midside nodes are expected
        for j, n in enumerate(e):
            if j <= 3: # iterating through the first 4 nodes of the cq* ele. These are the 4 corner nodes
                if ndict[int(n)]["totalCount"] == 2: #Coner node that only attaches two 2 elements is at a free edge
                    freeEdgeNodes.append(n)
        #Determinig the two corner nodes that lie on the free edge to find the corresponding free midside node
        if len(freeEdgeNodes) > 0:
            if e.index(freeEdgeNodes[0]) == 0 and e.index(freeEdgeNodes[1]) == 1:
                freeEdgeMidNodes = e[4]
            elif e.index(freeEdgeNodes[0]) == 1 and e.index(freeEdgeNodes[1]) == 2:
                freeEdgeMidNodes = e[5]
            elif e.index(freeEdgeNodes[0]) == 2 and e.index(freeEdgeNodes[1]) == 3:
                freeEdgeMidNodes = e[6]
            elif e.index(freeEdgeNodes[0]) == 3 and e.index(freeEdgeNodes[1]) == 0:
                freeEdgeMidNodes = e[7]
            elif e.index(freeEdgeNodes[0]) == 0 and e.index(freeEdgeNodes[1]) == 3:
                freeEdgeMidNodes = e[7]
        '''
        Determining the extra node of the CQ8 element that is not connected to any other elements
        and not on a free edge. This node is deleted and the 7 remaining nodes are the 7 nodes
        of the CQ7 element 
        '''
        for j, n in enumerate(e):
            if ndict[int(n)]["totalCount"] == 1 and n not in freeEdgeMidNodes: #If the node is only found on a single element aka free and not a free edge node
                e[e.index(n)] = "Blank"
                # re--ordering as necessary to ensure the free "8th" node lies between the 1 and 4 nodes as is assumed by the shape functions
                # re-odering does not have any impact on the down stream caclulations as the Jacobian maps global to natural coordinate systesm
                if e.index("Blank") == 7:
                    e.remove("Blank")
                elif e.index("Blank") == 4:
                    e = [e[1], e[2], e[3], e[0], e[5], e[6], e[7]]
                elif e.index("Blank") == 5:
                    e = [e[2], e[3], e[0], e[1], e[6], e[7], e[4]]
                elif e.index("Blank") == 6:
                    e = [e[3], e[0], e[1], e[2], e[7], e[4], e[5]]
                e7.append(e)
                nFree.append(n)
    e7Final = []
    if len(e7Common) != len(e7):
        raise ValueError("6 noded transition elements not yet implemented. Check for duplicate cq7 entries indicating a possible 6 noded transition element")
    for eA, eB in zip(e7Common, e7):
        e7Final.append(["CQUAD7", eA[1], eA[2]] + eB) #Re-assembling the full CQ7 Element Card
    if e7 != []:
        dfEle7 = pd.DataFrame(e7Final, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3", "N4", "N5", "N6", "N7"])
        dfEle7 = dfEle7.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int, "N4": int,
                            "N5":int, "N6":int, "N7":int})
    lenDFNodes = len(dfNodes)
    dfNodes = dfNodes[~dfNodes["N"].isin([nFree])]
    if lenDFNodes == len(dfNodes) and isinstance(dfEle7, pd.DataFrame):
        RuntimeWarning("I do not believe the free nodes have been removed properly ")
    
    return dfEle7, dfEle6, dfNodes

def DatFileParsing(dat):
    with open(dat) as datFile:
        contents = datFile.readlines()
    datFile.close()
    for c in contents:
        if c.startswith("+"): # appending multiline cards into a single line 
            contents[contents.index(c)-1] = contents[contents.index(c)-1].strip() + c.strip()[1:]
            contents.remove(c)
            
    nodes = [wrap(x, 8) for x in contents if x.startswith("GRID")]
    nArray = [[x[0], int(x[1]), int(x[2]),
            np.array([float(x[3]), float(x[4]), float(x[5])]), int(x[6])] for x in nodes]    

    dfNodes = pd.DataFrame(nArray, columns = ["Type", "N", "CP", "XYZ", "CD"])
    dfNodes = dfNodes[["N", "XYZ"]]
    
    e4 = [wrap(x, 8) for x in contents if x.startswith("CQUAD4")]
    e8 = [wrap(x, 8) for x in contents if x.startswith("CQUAD8")]
    e8 = [[y for y in x if y != "+"] for x in e8]
    
    dfEle4 = pd.DataFrame(e4, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3", "N4"])
    dfEle4 = dfEle4.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int, "N4": int})

    dfEle8 = pd.DataFrame(e8, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"])
    dfEle8 = dfEle8.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int, "N4": int,
                            "N5":int, "N6":int, "N7":int, "N8":int})
    if len(dfEle8) == 0:
        dfEle8, dfEle7, dfEle6 = None, None, None
    else:
        dfEle7, dfEle6, dfNodes = transitionEleParsing(dfNodes, dfEle4, dfEle8, e8)
    
    #processing materials
    materials = [wrap(x, 8)[0:7] for x in contents if x.startswith("MAT1")]
    print(materials)
    mTypes = {"Type": str, "MID":int, "E":float, 
                "NU":float, "RHO":float, "A":float,
                "TREF":float}#, "Geo":float}

    dfMaterials = pd.DataFrame(materials, columns = mTypes).astype(mTypes)
    dfMaterials = dfMaterials[["MID", "E", "NU", "RHO"]]
    
    #processing property cards
    properties = [wrap(x, 8)[0:4] for x in contents if x.startswith("PSHELL")]
    pTypes = {"Type": str, "PID":int, "MID1":int, "T":float}#,
            # "MID2":int, "MID3":int, "NSM":float}
    dfProperties = pd.DataFrame(properties, columns = pTypes).astype(pTypes)
    dfProperties = dfProperties[["PID", "MID1", "T"]]
    
    #combining into a single dataframe
    dfMatProps = pd.merge(dfProperties, dfMaterials, left_on = "MID1", right_on = "MID", suffixes = ("_Prop", "_Mat"))

    #Processing forces and constraints
    forces = [wrap(x, 8) for x in contents if x.startswith("FORCE")]
    constraints = [wrap(x, 8) for x in contents if x.startswith("SPC1")]
    fTypes = {"Type": str, "SID":int, "NID":int, "CID":int,
            "Factor":float, "F1":float, "F2":float, "F3":float}
    dfForces = pd.DataFrame(forces, columns = fTypes).astype(fTypes)

    cTypes = {"Type": str, "SID":int, "Comp":int, "NID":int}
    dfConstraints = pd.DataFrame(constraints, columns = cTypes).astype(cTypes)
    
    if len(dfEle4) == 0:
        dfEle4 = None
    
    dOut = [dfNodes, dfEle4, dfEle6, dfEle7, dfEle8, 
            dfMatProps, dfForces, dfConstraints]
    return dOut