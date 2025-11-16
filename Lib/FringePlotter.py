import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pandas as pd
import numpy as np

def Contour(dfNodes, dfValues, Components):
    
    try: 
        df = pd.merge(dfNodes, dfValues, how = "inner",
                      left_on = "N", right_on = "NID")
    except:
        raise KeyError("Merging columns not found")
    for comp in Components:
        df[comp + "_Avg"] = df.groupby("N")[comp].transform("mean")
    df = df.drop_duplicates("N")
    coords = np.array(df["XYZ"].tolist()) # Getting xyz locations of each node into a list
    triang = tri.Triangulation(coords[:,0], coords[:,1])
    for comp in Components:
        plt.tricontourf(triang, df[comp+"_Avg"].to_numpy(), cmap = "jet", levels = 12)
        plt.title(comp + "_Avg")
        plt.gca().ticklabel_format(style='plain')

        plt.colorbar()
        plt.show(block = True)
    
    