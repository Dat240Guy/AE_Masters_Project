
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
    import pandas as pd

    dfEle7 = None
    dfEle6 = None

    # ----------------------------------------------------
    # 1. BUILD NODE USAGE DICTIONARY (counts in CQ4 and CQ8)
    # ----------------------------------------------------
    ndict = {}
    for _, node in dfNodes.iterrows():
        nid = node["N"]
        n4Count = (dfEle4[['N1', 'N2', 'N3', 'N4']] == nid).sum().sum()
        n8Count = (dfEle8[['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']] == nid).sum().sum()
        ndict[nid] = {
            "n4Count": n4Count,
            "n8Count": n8Count,
            "totalCount": n4Count + n8Count
        }

    # ----------------------------------------------------
    # 2. IDENTIFY CQ8 ELEMENTS THAT TOUCH CQ4 (TRANSITION REGION)
    # ----------------------------------------------------
    nCommon = {nid: info for nid, info in ndict.items()
               if info["n8Count"] > 0 and info["n4Count"] > 0}

    # Filter CQ8 list: Those touching common nodes
    e8_candidates = [
        ele for ele in e8
        if any(str(nid) in ele[3:11] for nid in nCommon.keys())
    ]

    # Further filter: must have < 3 "strong CQ4-adjacent" corner nodes
    e8_transition = []
    for ele in e8_candidates:
        cnt = 0
        for nid in ele[3:7]:  # only corner nodes
            if ndict[int(nid)]["n4Count"] >= 2:
                cnt += 1
        if cnt <= 3:
            e8_transition.append(ele)

    # ----------------------------------------------------
    # 3. PROCESS EACH TRANSITION CQ8 ELEMENT → CQ7 or CQ6
    # ----------------------------------------------------
    e7_list = []    # list of tuples: (CQ8_full_card, [7 nodes])
    e6_list = []    # list of tuples: (CQ8_full_card, [6 nodes])
    nFree = []      # nodes to delete later
    blankCountList = []
    for ele in e8_transition:
        efull = ele                        # original CQ8 structure
        nodes = ele[3:11].copy()           # list of 8 nodes [C1,C2,C3,C4,MB,MR,MT,ML]

        # ---------------------------------------------
        # 3a. Identify free-edge corner nodes
        # ---------------------------------------------
        freeEdgeCorners = []
        for j, nid in enumerate(nodes[:4]):
            if ndict[int(nid)]["totalCount"] == 2:
                freeEdgeCorners.append(nid)

        # Determine the free-edge midside node (if present)
        freeEdgeMid = None
        if len(freeEdgeCorners) >= 2:
            i1 = nodes.index(freeEdgeCorners[0])
            i2 = nodes.index(freeEdgeCorners[1])

            # Corner index pairs → matching midside
            if (i1, i2) in [(0, 1), (1, 0)]:
                freeEdgeMid = nodes[4]   # bottom mid
            elif (i1, i2) in [(1, 2), (2, 1)]:
                freeEdgeMid = nodes[5]   # right mid
            elif (i1, i2) in [(2, 3), (3, 2)]:
                freeEdgeMid = nodes[6]   # top mid
            elif (i1, i2) in [(3, 0), (0, 3)]:
                freeEdgeMid = nodes[7]   # left mid

        # ---------------------------------------------
        # 3b. Blank nodes: those only appearing once & not free-edge-mid
        # ---------------------------------------------
        blankCount = 0
        
        for idx, nid in enumerate(nodes):
            if ndict[int(nid)]["totalCount"] == 1 and nid != freeEdgeMid:
                nodes[idx] = "Blank"
                blankCount += 1
                blankCountList.append(nid)
                nFree.append(nid)

        # ---------------------------------------------
        # 3c. CASE 1: CQ7 (1 blank)
        # ---------------------------------------------
        if blankCount == 1:
            bidx = nodes.index("Blank")

            # remove Blank or rotate appropriately
            if bidx == 7:
                # already between N1 and N4
                nodes.remove("Blank")
            elif bidx == 4:
                # Left shift corners
                nodes = [nodes[1], nodes[2], nodes[3], nodes[0],
                         nodes[5], nodes[6], nodes[7]]
            elif bidx == 5:
                nodes = [nodes[2], nodes[3], nodes[0], nodes[1],
                         nodes[6], nodes[7], nodes[4]]
            elif bidx == 6:
                nodes = [nodes[3], nodes[0], nodes[1], nodes[2],
                         nodes[7], nodes[4], nodes[5]]

            e7_list.append((efull, nodes))

        # ---------------------------------------------
        # 3d. CASE 2: CQ6 (2 blanks)
        # ---------------------------------------------
        elif blankCount == 2:

            # find blank positions
            blank_indices = sorted(i for i, nid in enumerate(nodes) if nid == "Blank")

            # rotate corners until blanks sit at MT(6) and ML(7)
            def rotate_element(e):
                """
                Rotate CQ8 node list CCW while preserving
                the relative mapping of midside nodes.
                """
                C1, C2, C3, C4, MB, MR, MT, ML = e

                return [
                    C2, C3, C4, C1,   # rotated corners
                    MR, MT, ML, MB    # rotated midsides
                ]

            for _ in range(4):
                blank_indices = [i for i, n in enumerate(nodes) if n == "Blank"]
                if blank_indices == [6, 7]:
                    break
                nodes = rotate_element(nodes)

            if blank_indices != [6, 7]:
                raise ValueError("Failed to align blanks to MT/ML for CQ6 element")

            # Remove MT and ML → keep first 6 nodes
            reduced = nodes[:6]
            e6_list.append((efull, reduced))

        else:
            raise ValueError("Invalid number of free nodes: expected 1 or 2.")

    # ----------------------------------------------------
    # 4. Assemble dfEle7 dataframe
    # ----------------------------------------------------
    if e7_list:
        e7Final = []
        for efull, nodes7 in e7_list:
            e7Final.append(["CQUAD7", efull[1], efull[2]] + nodes7)

        dfEle7 = pd.DataFrame(
            e7Final,
            columns=["Type", "Enumber", "Prop",
                     "N1", "N2", "N3", "N4", "N5", "N6", "N7"]
        ).astype({
            "Type": str, "Enumber": int, "Prop": int,
            "N1": int, "N2": int, "N3": int, "N4": int,
            "N5": int, "N6": int, "N7": int
        })

    # ----------------------------------------------------
    # 5. Assemble dfEle6 dataframe
    # ----------------------------------------------------
    if e6_list:
        e6Final = []
        for efull, nodes6 in e6_list:
            e6Final.append(["CQUAD6", efull[1], efull[2]] + nodes6)

        dfEle6 = pd.DataFrame(
            e6Final,
            columns=["Type", "Enumber", "Prop",
                     "N1", "N2", "N3", "N4", "N5", "N6"]
        ).astype({
            "Type": str, "Enumber": int, "Prop": int,
            "N1": int, "N2": int, "N3": int,
            "N4": int, "N5": int, "N6": int
        })

    # ----------------------------------------------------
    # 6. Remove free nodes from DFNodes
    # ----------------------------------------------------
    lenBefore = len(dfNodes)
    dfNodes = dfNodes[~dfNodes["N"].isin(nFree)]

    if lenBefore == len(dfNodes) and (dfEle7 is not None or dfEle6 is not None):
        print("WARNING: free nodes may not have been removed properly")

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
    
    t3 = [wrap(x, 8) for x in contents if x.startswith("CTRIA3")]
    e4 = [wrap(x, 8) for x in contents if x.startswith("CQUAD4")]
    e8 = [wrap(x, 8) for x in contents if x.startswith("CQUAD8")]
    e8 = [[y for y in x if y != "+"] for x in e8]
    
    dfEle3 = pd.DataFrame(t3, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3"])
    dfEle3 = dfEle3.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int})
    
    dfEle4 = pd.DataFrame(e4, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3", "N4"])
    dfEle4 = dfEle4.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int, "N4": int})

    dfEle8 = pd.DataFrame(e8, columns= ["Type", "Enumber", "Prop", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8"])
    dfEle8 = dfEle8.astype({"Type": str, "Enumber":int, "Prop":int, "N1":int, "N2":int, "N3":int, "N4": int,
                            "N5":int, "N6":int, "N7":int, "N8":int})
    if len(dfEle8) == 0:
        dfEle8, dfEle7, dfEle6 = None, None, None
    else:
        dfEle7, dfEle6, dfNodes = transitionEleParsing(dfNodes, dfEle4, dfEle8, e8)
    
    if len(dfEle3) == 0: dfEle3 = None
    
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
    
    dOut = [dfNodes, dfEle4, dfEle3, dfEle6, dfEle7, dfEle8, 
            dfMatProps, dfForces, dfConstraints]
    return dOut