# import matplotlib.pyplot as plt
# import matplotlib.tri as tri
# import pandas as pd
# import numpy as np

# def Contour(dfNodes, dfValues, Components):
    
#     try: 
#         df = pd.merge(dfNodes, dfValues, how = "inner",
#                       left_on = "N", right_on = "NID")
#     except:
#         raise KeyError("Merging columns not found")
#     for comp in Components:
#         df[comp + "_Avg"] = df.groupby("N")[comp].transform("mean")
#     df = df.drop_duplicates("N")
#     coords = np.array(df["XYZ"].tolist()) # Getting xyz locations of each node into a list
#     triang = tri.Triangulation(coords[:,0], coords[:,1])
#     for comp in Components:
#         plt.tricontourf(triang, df[comp+"_Avg"].to_numpy(), cmap = "jet", levels = 12)
#         plt.title(comp + "_Avg")
#         plt.gca().ticklabel_format(style='plain')

#         plt.colorbar()
#         plt.show(block = True)



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# ------------------------------------------------------------
# Build triangles directly from FE connectivity
# ------------------------------------------------------------

def build_triangles_from_elements(df_nodes_merged, dfEles):
    """
    df_nodes_merged : DataFrame after merging dfNodes + dfValues and
                      dropping duplicates. Must contain columns "N" and "XYZ".
    dfEles          : list of element DataFrames [dfEle4, dfEle8, dfEle7, dfEle6, ...]

    Returns:
        coords    : (Nnode, 2) array of x,y
        triangles : (Ntri, 3) array of integer indices into coords
    """

    # Node ID -> local index in coords
    nids = df_nodes_merged["N"].astype(int).to_numpy()
    nid_to_idx = {nid: i for i, nid in enumerate(nids)}

    coords = np.array(df_nodes_merged["XYZ"].tolist())[:, :2]  # just x,y

    triangles = []

    for dfE in dfEles:
        # Node columns: N1, N2, ...
        node_cols = [c for c in dfE.columns if c.startswith("N")]
        if not node_cols:
            continue

        for _, row in dfE[node_cols].iterrows():
            # full node list for this element (strip zeros/padding)
            nodes_all = [int(row[c]) for c in node_cols if int(row[c]) != 0]

            nnode = len(nodes_all)

            # ---- Triangular elements (CTRIA3 / CTRIA6 style) ----
            if nnode in (3, 6):
                # Use first 3 nodes as corners (N1,N2,N3)
                corners = nodes_all[:3]
                # skip if this element lies outside the region in df (e.g. different part)
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3 = (nid_to_idx[n] for n in corners)
                triangles.append([i1, i2, i3])

            # ---- Quadrilateral elements (CQUAD4 / 7 / 8 style) ----
            elif nnode in (4, 7, 8):
                # For Nastran/Femap quads, N1..N4 are the corners.
                corners = nodes_all[:4]
                if not all(n in nid_to_idx for n in corners):
                    continue
                i1, i2, i3, i4 = (nid_to_idx[n] for n in corners)

                # Split quad into two triangles: (1,2,3) and (1,3,4)
                triangles.append([i1, i2, i3])
                triangles.append([i1, i3, i4])

            else:
                # Unknown topology – ignore or raise if you prefer
                # raise ValueError(f"Unsupported element with {nnode} nodes: {nodes_all}")
                continue

    triangles = np.array(triangles, dtype=int)
    return coords, triangles


# ------------------------------------------------------------
# Contour function: triangulation only over actual elements
# ------------------------------------------------------------

def Contour(dfNodes, dfValues, Components, dfEles):
    """
    dfNodes   : DataFrame with columns N, XYZ
    dfValues  : DataFrame with columns NID, component fields (e.g. S1, S2, S12, ...)
    Components: list of strings, names of fields to contour (e.g. ["S1", "S_max"])
    dfEles    : [dfEle4, dfEle8, dfEle7, dfEle6, ...] element dataframes

    This builds a triangulation from the real FE connectivity, so
    contours are drawn only over meshed regions (holes stay empty).
    """

    # Merge values onto nodes
    df = dfNodes.merge(dfValues, how="inner", left_on="N", right_on="NID")

    # Nodal averaging (in case multiple element contributions per node)
    for comp in Components:
        df[comp + "_Avg"] = df.groupby("N")[comp].transform("mean")

    # One row per node
    df = df.drop_duplicates("N")

    # Build triangles directly from element connectivity
    coords, triangles = build_triangles_from_elements(df, dfEles)

    # Custom triangulation (no Delaunay fill)
    triang = tri.Triangulation(coords[:, 0], coords[:, 1], triangles=triangles)

    # Plot for each requested component
    for comp in Components:
        z = df[comp + "_Avg"].to_numpy()

        plt.figure()
        plt.tricontourf(triang, z, levels=14, cmap="jet")
        plt.title(comp + "_Avg")
        plt.gca().ticklabel_format(style="plain")
        plt.colorbar()
        plt.axis("equal")
        plt.show(block=True)



































'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from collections import defaultdict
import shapely.geometry as geom


# ============================================================
# 0. Build connectivity from your element dataframes
#    dfEles = [dfEle4, dfEle8, dfEle7, dfEle6, ...]
# ============================================================

def build_connectivity(dfEles):
    """
    dfEles : list of DataFrames for each element type.
             Each df must have columns N1, N2, ..., Nn (corner+mid nodes).

    Returns:
        connectivity : list of lists of node IDs (per element)
    """
    connectivity = []

    for df in dfEles:
        node_cols = [c for c in df.columns if c.startswith("N")]
        elems = df[node_cols].astype(int).values.tolist()

        for e in elems:
            # Remove zeros / padding if present
            e_clean = [int(n) for n in e if int(n) != 0]
            connectivity.append(e_clean)

    return connectivity


# ============================================================
# 1. TRUE geometric edges from an element (corner-to-corner)
#    Handles 4, 6, 7, 8-noded elems (CQUAD4, CTRIA6, CQUAD7, CQUAD8)
# ============================================================

def element_edges(elem):
    n = elem

    if len(n) == 4:      # CQUAD4
        c = [n[0], n[1], n[2], n[3]]

    elif len(n) == 8:    # CQUAD8 (Nastran / Femap)
        # corners at N1, N3, N5, N7
        c = [n[0], n[2], n[4], n[6]]

    elif len(n) == 7:    # CQUAD7 transition in your mesh
        # in your model the CQUAD7 is a CQUAD8 with one midside collapsed,
        # so it shares the SAME corner positions as CQUAD8
        c = [n[0], n[2], n[4], n[6]]

    elif len(n) == 6:    # CTRIA6
        # corners at N1, N3, N5
        c = [n[0], n[2], n[4]]

    else:
        raise ValueError(f"Unsupported element topology: {elem}")

    return [(c[i], c[(i+1) % len(c)]) for i in range(len(c))]





# ============================================================
# 2. Extract edges + boundary edges
# ============================================================

def extract_edges_from_elements(connectivity):
    edge_count = defaultdict(int)
    for elem in connectivity:
        for e in element_edges(elem):
            edge_count[tuple(sorted(e))] += 1
    return edge_count


def get_boundary_edges(edge_count):
    # Edges used by exactly one element are on the boundary
    return [e for e, c in edge_count.items() if c == 1]


# ============================================================
# 3. Build ordered loops from boundary edges (bullet-proof)
# ============================================================

def build_loops(boundary_edges):
    adj = defaultdict(list)
    for a, b in boundary_edges:
        adj[a].append(b)
        adj[b].append(a)

    # All boundary nodes must have degree 2 for a proper loop
    for node, neigh in adj.items():
        if len(neigh) != 2:
            raise ValueError(
                f"Boundary node {node} has {len(neigh)} neighbors → invalid boundary (check connectivity/midside mapping)."
            )

    loops = []
    visited_nodes = set()

    for start in adj:
        if start in visited_nodes:
            continue

        loop = []
        current = start
        prev = None

        while True:
            loop.append(current)
            visited_nodes.add(current)

            neighbors = adj[current]
            # next = neighbor that is not the previous node
            nxt = neighbors[0] if neighbors[1] == prev else neighbors[1]

            prev, current = current, nxt
            if current == start:
                break

        loops.append(loop)

    return loops


# ============================================================
# 4. Polygon area + separate outer vs holes
# ============================================================

def polygon_area(coords):
    x = coords[:, 0]
    y = coords[:, 1]
    return 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])


def separate_outer_and_holes(loops, node_coords):
    """
    loops       : list of loops (node ID lists)
    node_coords : dict {nodeID: (x,y)}

    Returns:
        outer_pts     : Nx2 array (CCW)
        hole_boundaries : list of arrays (each CW)
    """
    loop_info = []

    for loop in loops:
        pts = np.array([node_coords[n] for n in loop])

        # close loop
        if not np.allclose(pts[0], pts[-1]):
            pts = np.vstack([pts, pts[0]])

        A = polygon_area(pts)
        loop_info.append((pts, A))

    # Outer boundary = loop with largest |area|
    outer_pts, outer_A = max(loop_info, key=lambda x: abs(x[1]))

    # Ensure CCW for outer boundary
    if outer_A < 0:
        outer_pts = outer_pts[::-1]

    hole_boundaries = []
    for pts, A in loop_info:
        if np.array_equal(pts, outer_pts):
            continue

        # Ensure CW for holes
        if A > 0:
            pts = pts[::-1]

        hole_boundaries.append(pts)

    return outer_pts, hole_boundaries


# ============================================================
# 5. Master: extract mesh boundaries from dfNodes + connectivity
# ============================================================

def extract_mesh_boundaries(dfNodes, connectivity):
    """
    dfNodes: DataFrame with columns "N" (node ID) and "XYZ" (np.array([x,y,z]))
    connectivity: list of elements (list of node IDs)

    Returns:
        outer_boundary : Nx2 array (CCW)
        hole_boundaries : list of arrays (CW)
    """
    node_coords = {int(row["N"]): tuple(row["XYZ"][:2])
                   for _, row in dfNodes.iterrows()}

    edge_count = extract_edges_from_elements(connectivity)
    boundary_edges = get_boundary_edges(edge_count)
    loops = build_loops(boundary_edges)
    outer_boundary, hole_boundaries = separate_outer_and_holes(loops, node_coords)

    return outer_boundary, hole_boundaries


# ============================================================
# 6. Smart, domain-aware contour plotting
# ============================================================

def ContourSmart(dfNodes, dfValues, Components, dfEles):
    """
    dfNodes    : DataFrame with columns N, XYZ
    dfValues   : DataFrame with columns NID, component fields
    Components : list of column names in dfValues to contour (e.g. ["S1", "S2"])
    dfEles     : list of element DataFrames [dfEle4, dfEle8, dfEle7, dfEle6, ...]

    This function:
      - builds connectivity from dfEles
      - auto-detects mesh boundary + holes
      - masks triangulation outside the domain
      - plots nodal-averaged contours for each component
    """

    # Build mixed-element connectivity
    connectivity = build_connectivity(dfEles)

    # Merge nodal coords with values
    df = dfNodes.merge(dfValues, how="inner", left_on="N", right_on="NID")

    # Nodal averaging of duplicate contributions
    for comp in Components:
        df[comp + "_Avg"] = df.groupby("N")[comp].transform("mean")

    df = df.drop_duplicates("N")

    coords = np.array(df["XYZ"].tolist())
    x = coords[:, 0]
    y = coords[:, 1]

    # Extract mesh boundaries (outer + holes)
    outer_boundary, hole_boundaries = extract_mesh_boundaries(dfNodes, connectivity)

    # Build polygon domain
    domain_polygon = geom.Polygon(outer_boundary, holes=hole_boundaries)

    # Triangulate all nodes
    triang = tri.Triangulation(x, y)

    # Mask triangles whose centroids lie outside the domain
    tri_nodes = triang.triangles
    xc = x[tri_nodes].mean(axis=1)
    yc = y[tri_nodes].mean(axis=1)

    mask = []
    for xi, yi in zip(xc, yc):
        p = geom.Point(xi, yi)
        mask.append(not domain_polygon.contains(p))

    triang.set_mask(mask)

    # Plot each component
    for comp in Components:
        z = df[comp + "_Avg"].to_numpy()

        plt.figure()
        plt.tricontourf(triang, z, levels=14, cmap="jet")
        plt.title(comp + "_Avg")
        plt.gca().ticklabel_format(style="plain")
        plt.colorbar()
        plt.axis("equal")
        plt.show(block=True)
'''

