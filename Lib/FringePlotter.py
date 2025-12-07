
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri


def build_triangles_from_elements(df_nodes_merged, dfEles):

    nids = df_nodes_merged["N"].astype(int).to_numpy()
    nid_to_idx = {nid: i for i, nid in enumerate(nids)}

    coords = np.array(df_nodes_merged["XYZ"].tolist())[:, :2]  # x, y

    triangles = []

    for dfE in dfEles:
        node_cols = [c for c in dfE.columns if c.startswith("N")]
        if not node_cols:
            continue

        for _, row in dfE[node_cols].iterrows():
            nodes_all = [int(row[c]) for c in node_cols if int(row[c]) != 0]

            nnode = len(nodes_all)

            # CTRIAS
            if nnode in (2, 3):
                
                corners = nodes_all[:3]
                
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3 = (nid_to_idx[n] for n in corners)
                triangles.append([i1, i2, i3])

            # CQUADS
            elif nnode in (4, 6, 7, 8):
                corners = nodes_all[:4]
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3, i4 = (nid_to_idx[n] for n in corners)

                triangles.append([i1, i2, i3])
                triangles.append([i1, i3, i4])

            else:
                continue

    triangles = np.array(triangles, dtype=int)
    return coords, triangles


def Contour(dfNodes, dfValues, Components, dfEles, Averaging="Nodal"):

    coords, triangles = build_triangles_from_elements(dfNodes, dfEles)

    XY = np.vstack(dfNodes["XYZ"].values)
    X = XY[:, 0]
    Y = XY[:, 1]

    # Averaging
    if Averaging == "Nodal":
        df = dfNodes.merge(dfValues, how="inner", left_on="N", right_on="NID")

        for comp in Components:
            df[comp + "_Plot"] = df.groupby("N")[comp].transform("mean")

        df = df.drop_duplicates("N")  # one row per node, with averaged values
        Zfields = {comp: df[comp + "_Plot"].to_numpy() for comp in Components}

        triang = tri.Triangulation(X, Y, triangles=triangles)

    elif Averaging == "Elemental":

        ele_centroids = []
        ele_fields = {comp: [] for comp in Components}

        for dfEle in dfEles:
            for _, row in dfEle.iterrows():
                # Node IDs for this element
                nid_list = []
                for c in row.index:
                    if c.startswith("N") and c[1:].isdigit():
                        nid_list.append(row[c])

                # Coordinates
                pts = np.vstack(dfNodes[dfNodes["N"].isin(nid_list)]["XYZ"])

                # Element centroid
                centroid = pts[:, :2].mean(axis=0)
                ele_centroids.append(centroid)

                # Element-averaged values
                df_sub = dfValues[dfValues["NID"].isin(nid_list)]
                for comp in Components:
                    ele_fields[comp].append(df_sub[comp].mean())

        ele_centroids = np.array(ele_centroids)
        Xc = ele_centroids[:, 0]
        Yc = ele_centroids[:, 1]
        
        triang = tri.Triangulation(Xc, Yc)

        Zfields = {comp: np.array(vals) for comp, vals in ele_fields.items()}

    else:
        raise ValueError('Averaging must be "Nodal" or "Elemental"')

    for comp in Components:
        Z = Zfields[comp]

        plt.figure()
        plt.tricontourf(triang, Z, levels=14, cmap="jet")
        # plt.tricontourf(triang, Z, levels = 10, cmap= "jet", 
        #                 vmin = 0, vmax = 20000)
        plt.title(f"{comp} ({Averaging} Avg)")
        plt.gca().ticklabel_format(style="plain")
        plt.colorbar()
        plt.axis("equal")
        plt.show(block=True)

        # Z = Zfields[comp]

        # plt.figure()

        # contour = plt.tricontourf(
        #     triang, Z,
        #     levels=10,
        #     cmap="jet",
        #     vmin=0, vmax=20000
        # )

        # plt.title(f"{comp} ({Averaging} Avg)")
        # plt.gca().ticklabel_format(style="plain")

        # cbar = plt.colorbar(contour)

        # # Disable scientific notation offsets
        # cbar.formatter.set_useOffset(False)
        # cbar.update_ticks()

        # plt.axis("equal")
        # plt.show(block=True)




