import numpy as np

class PlaneStress:
    def __init__(self, E, v, t):
        self.E = E
        self.v = v
        self.t = t
        self.Array = (E / (1 - v**2)) * np.array([[1, v, 0],
                                              [v, 1, 0],
                                              [0, 0, (1-v) / 2]])

class PlaneStrain:
    def __init__(self, E, v, t):
        self.E = E
        self.v = v
        self.t = t
        self.Array = (E / ((1 + v)*(1-2*v))) * np.array([[1-v, v, 0],
                                                         [v, 1-v, 0],
                                                         [0, 0, (1 - 2*v) / 2]])
        
# Universal class for calculating the element Jacobian and B Matrix        
class qCalc: 
    def __init__(self, element):
        self.element = element
    
    class jacobian:
        def __init__(self, element, xi, eta):
            self.J = np.array([
                [sum(element.dN_dxi[i](xi, eta) * element.globalCoord[i][0] for i in range(element.nodeCount)),
                 sum(element.dN_dxi[i](xi, eta) * element.globalCoord[i][1] for i in range(element.nodeCount))],
                [sum(element.dN_deta[i](xi, eta) * element.globalCoord[i][0] for i in range(element.nodeCount)),
                 sum(element.dN_deta[i](xi, eta) * element.globalCoord[i][1] for i in range(element.nodeCount))]
            ])
            self.det  = np.linalg.det(self.J)
            if self.det <= 0:
                print("WARNING: detJ <= 0 at element", element.ID,"xi,eta", xi,eta, "detJ=", self.det)
            elif self.det < 1e-6:
                print("Small detJ at element", element.ID, "detJ=", self.det)

            self.invT = np.linalg.inv(self.J).T
            
            self.inv = np.linalg.inv(self.J)
            
    def Nvals(self, element, xi, eta):
        return np.array([N(xi, eta) for N in element.N])

    def dN_dxi_vals(self, element, xi, eta):
        return np.array([dN(xi, eta) for dN in element.dN_dxi])

    def dN_deta_vals(self, element, xi, eta):
        return np.array([dN(xi, eta) for dN in element.dN_deta]) 
    

    def B(self, xi, eta, jacb):

            e = self.element
            n = e.nodeCount
            totalDof = e.totalDof

            # Derivatives 
            dN_dxi  = self.dN_dxi_vals(e, xi, eta)   # (n,)
            dN_deta = self.dN_deta_vals(e, xi, eta)  # (n,)

            invJ = jacb.invT.T                       
            
            grads_nat = np.vstack((dN_dxi, dN_deta))  # derivatives in natural coords
            grads_xy  = invJ @ grads_nat # multivariable chain rule for derivatives in global coords

            dN_dx = grads_xy[0, :]
            dN_dy = grads_xy[1, :]

            B = np.zeros((3, totalDof))
            for i in range(n):
                col_u = 2 * i
                col_v = 2 * i + 1

                B[0, col_u] = dN_dx[i]  
                B[1, col_v] = dN_dy[i]  
                B[2, col_u] = dN_dy[i]  
                B[2, col_v] = dN_dx[i]

            return B

class q4:
    def __init__(self, globalCoord, ID = None):
        self.nodeCount = 4
        self.dimensions = 2
        self.dofPerNode = 2
        self.totalDof = self.nodeCount * self.dofPerNode
        self.ID = ID
        self.globalCoord = globalCoord
        # self.globalCoord, self.nOrder = self.globalCoordReOrder(nodeOrder=True) #THIS HAS BEEN COMMENTED OUT DO NOT THINK IT IS NEEDED
        
        # Node points of the natural element (xi, eta)
        self.localCoord = np.array([[-1, -1],
                                    [1, -1],
                                    [1, 1],
                                    [-1, 1]])

        #Shape functions 1 at their corner location and zero everywhere else See Sympy
        self.N = [lambda xi, eta: 0.25 * (1 - xi) * (1 - eta),
                  lambda xi, eta: 0.25 * (1 + xi) * (1 - eta),
                  lambda xi, eta: 0.25 * (1 + xi) * (1 + eta),
                  lambda xi, eta: 0.25 * (1 - xi) * (1 + eta)]

        self.dN_dxi = [lambda xi, eta: -0.25 * (1 - eta),
                       lambda xi, eta: 0.25 * (1 - eta),
                       lambda xi, eta: 0.25 * (1 + eta),
                       lambda xi, eta: -0.25 * (1 + eta)]

        self.dN_deta = [lambda xi, eta: -0.25 * (1 - xi),
                        lambda xi, eta: -0.25 * (1 + xi),
                        lambda xi, eta: 0.25 * (1 + xi),
                        lambda xi, eta: 0.25 * (1 - xi)]
        # Integration points
        self.xiIntegrationPoints = [-1/np.sqrt(3), 1/np.sqrt(3)]
        self.etaIntegrationPoints = [-1/np.sqrt(3), 1/np.sqrt(3)]
        self.Weights = [1, 1]
        
class q8:
    def __init__(self, globalCoord, ID = None):
        self.nodeCount = 8
        self.dimensions = 2
        self.dofPerNode = 2
        self.ID = ID
        self.totalDof = self.nodeCount * self.dofPerNode
        self.globalCoord = [x[0:2] for x in globalCoord]

        self.localCoord = np.array([[-1, -1], #LL
                                    [1, -1],  #LR
                                    [1, 1],   #UR
                                    [-1, 1],  #UL
                                    [0, -1],  #MB
                                    [1, 0],   #MR
                                    [0, 1],   #MU
                                    [-1, 0]]) #ML

        self.N = [lambda xi, eta: -1/4 * (1-eta) * (1-xi) * (1+xi+eta), # N1 = -1/4 * (1-eta) * (1-xi) * (1+xi+eta)
                  lambda xi, eta: -1/4 * (1-eta) * (1+xi) * (1-xi+eta), # N2 = -1/4 * (1-eta) * (1+xi) * (1-xi+eta)
                  lambda xi, eta: -1/4 * (1+eta) * (1+xi) * (1-xi-eta), # N3 = -1/4 * (1+eta) * (1+xi) * (1-xi-eta)
                  lambda xi, eta: -1/4 * (1+eta) * (1-xi) * (1+xi-eta), # N4 = -1/4 * (1+eta) * (1-xi) * (1+xi-eta)
                  lambda xi, eta:  1/2 * (1-eta) * (1-xi) * (1+xi),     # N5 =  1/2 * (1-eta) * (1-xi) * (1+xi)
                  lambda xi, eta:  1/2 * (1-eta) * (1+xi) * (1+eta),    # N6 =  1/2 * (1-eta) * (1+xi) * (1+eta)
                  lambda xi, eta:  1/2 * (1+eta) * (1-xi) * (1+xi),     # N7 =  1/2 * (1+eta) * (1-xi) * (1+xi)
                  lambda xi, eta:  1/2 * (1-eta) * (1-xi) * (1+eta)]    # N8 =  1/2 * (1-eta) * (1-xi) * (1+eta)

        self.dN_dxi = [lambda xi, eta: 0.25*(-eta-2*xi)*(eta-1),
                       lambda xi, eta: 0.25*(eta-1)*(eta-2*xi),
                       lambda xi, eta: 0.25*(eta+1)*(eta+2*xi),
                       lambda xi, eta: 0.25*(-eta+2*xi)*(eta+1),
                       lambda xi, eta: xi*(eta-1),
                       lambda xi, eta: 0.5-0.5*eta**2,
                       lambda xi, eta: xi*(-eta-1),
                       lambda xi, eta: 0.5*eta**2 - 0.5]


        self.dN_deta = [lambda xi, eta: 0.25*(-2*eta-xi)*(xi-1),
                        lambda xi, eta: 0.25*(2*eta-xi)*(xi+1),
                        lambda xi, eta: 0.25*(2*eta+xi)*(xi+1),
                        lambda xi, eta: 0.25*(-2*eta+xi)*(xi-1),
                        lambda xi, eta: 0.5*xi**2-0.5,
                        lambda xi, eta: eta*(-xi-1),
                        lambda xi, eta: 0.5-0.5*xi**2,
                        lambda xi, eta: eta*(xi-1)]
        
        #Gauss points and weights
        self.xiIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.etaIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.Weights = [5/9, 8/9, 5/9]
        
class q7:
    def __init__(self, globalCoord, ID = None):
        self.nodeCount = 7
        self.dimensions = 2
        self.dofPerNode = 2
        self.ID = ID
        self.totalDof = self.nodeCount * self.dofPerNode
        self.globalCoord = [x[0:2] for x in globalCoord]

        self.localCoord = np.array([[-1, -1],  #LL
                                    [1, -1],   #LR
                                    [1, 1],    #UR
                                    [-1, 1],   #UL
                                    [0, -1],   #MB
                                    [1, 0],    #MR
                                    [0, 1]])   #MU

        self.N = [lambda xi, eta: 0.25*xi*(1-xi)*(eta-1),
                  lambda xi, eta: 0.25*(eta**2*xi + eta**2 - eta*xi**2 - eta*xi + xi**2 -1),
                  lambda xi, eta: 0.25*(eta**2*xi + eta**2 + eta*xi**2 + eta*xi + xi**2 -1),
                  lambda xi, eta: 0.25*xi*(xi-1)*(eta+1),
                  lambda xi, eta: 0.5*(xi**2-1)*(eta-1),
                  lambda xi, eta: -0.5*(xi+1)*(eta**2-1),
                  lambda xi, eta: -0.5*(xi**2-1)*(eta+1)]

        self.dN_dxi = [lambda xi, eta: 0.25*eta-0.5*xi*(eta-1)-0.25,
                       lambda xi, eta: 0.25*eta**2-0.25*eta-0.5*xi*(eta-1),
                       lambda xi, eta: 0.25*eta**2+0.25*eta+0.5*xi*(eta+1),
                       lambda xi, eta: -0.25*eta+0.5*xi*(eta+1)-0.25,
                       lambda xi, eta: xi*(eta-1),
                       lambda xi, eta: 0.5-0.5*eta**2,
                       lambda xi, eta: -xi*(eta+1)]


        self.dN_deta = [lambda xi, eta: 0.25*xi*(1-xi),
                        lambda xi, eta: 0.5*eta*(xi+1)-0.25*xi**2-0.25*xi,
                        lambda xi, eta: 0.5*eta*(xi+1)+0.25*xi**2+0.25*xi,
                        lambda xi, eta: 0.25*xi*(xi-1),
                        lambda xi, eta: 0.5*xi**2-0.5,
                        lambda xi, eta: -eta*(xi+1),
                        lambda xi, eta: 0.5-0.5*xi**2]
        #Gauss points and weights 
        self.xiIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.etaIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.Weights = [5/9, 8/9, 5/9]
        
class q6:

    def __init__(self, globalCoord, ID=None):
        self.nodeCount   = 6
        self.dimensions  = 2
        self.dofPerNode  = 2
        self.ID          = ID
        self.totalDof    = self.nodeCount * self.dofPerNode
        self.globalCoord = [x[0:2] for x in globalCoord]

        self.localCoord = np.array([
            [-1.0, -1.0],  
            [ 1.0, -1.0],  
            [ 1.0,  1.0],  
            [-1.0,  1.0],  
            [ 0.0, -1.0],  
            [ 1.0,  0.0],  
        ])

        self.N = [
            lambda xi, eta: 0.25*xi*(1-xi)*(eta-1),
            lambda xi, eta: 0.25*eta**2*xi + 0.25*eta**2 - 0.25*eta*(xi**2) - 0.25*eta*xi + 0.25*xi**2 - 0.25,
            lambda xi, eta: 0.25*eta*(eta+1)*(xi+1),
            lambda xi, eta: -0.25*(eta+1)*(xi-1),
            lambda xi, eta: 0.5*(eta-1)*(xi**2-1),
            lambda xi, eta: -0.5*(eta**2-1)*(xi+1),
        ]

        self.dN_dxi = [
            lambda xi, eta: 0.25*eta-0.5*xi*(eta-1)-0.25,
            lambda xi, eta: 0.25*eta**2-0.25*eta-0.5*xi*(eta-1),
            lambda xi, eta: 0.25*eta*(eta+1),
            lambda xi, eta: -0.25*eta-0.25,
            lambda xi, eta: xi*(eta-1),
            lambda xi, eta: 0.5-0.5*eta**2,
        ]

        self.dN_deta = [
            lambda xi, eta: 0.25*xi*(1-xi),
            lambda xi, eta: 0.5*eta*(xi+1)-0.25*xi**2-0.25*xi,
            lambda xi, eta: 0.5*eta*(xi+1)+0.25*xi+0.25,
            lambda xi, eta: 0.25-0.25*xi,
            lambda xi, eta: 0.5*xi**2-0.5,
            lambda xi, eta: -eta*(xi+1),
        ]

        # self.xiIntegrationPoints  = [-np.sqrt(1/3), np.sqrt(1/3)]
        # self.etaIntegrationPoints = [-np.sqrt(1/3), np.sqrt(1/3)]
        # self.Weights              = [1, 1]

        self.xiIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.etaIntegrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        self.Weights = [5/9, 8/9, 5/9]
        
class t3:
    def __init__(self, globalCoord, ID = None):
        self.ID = ID
        self.nodeCount = 3
        self.totalDof = 2 * self.nodeCount
 
        self.globalCoord = [x[0:2] for x in globalCoord]


        self.xiIntegrationPoints  = [1.0 / 3.0]
        self.etaIntegrationPoints = [1.0 / 3.0]
        
        self.Weights = [0.5]

        self.localCoord = np.array([
            [0.0, 0.0],  # N1
            [1.0, 0.0],  # N2
            [0.0, 1.0],  # N3
        ])

        self.N = [lambda xi, eta: 1 - xi - eta,
                  lambda xi, eta: xi,
                  lambda xi, eta: eta]

        self.dN_dxi = [lambda xi, eta: -1,
                       lambda xi, eta: 1,
                       lambda xi, eta: 0]
        
        self.dN_deta = [lambda xi, eta: -1,
                       lambda xi, eta: 0,
                       lambda xi, eta: 1]


# =====================================================================
# Shape function validation tests (can be used with pytest or __main__)
# =====================================================================

def _check_kronecker_delta(element, tol=1e-12):
    """
    Check: N_a(xi_b, eta_b) = delta_ab at each nodal point.
    """
    for a, (xi_a, eta_a) in enumerate(element.localCoord):
        vals = [N(xi_a, eta_a) for N in element.N]
        # Own node should be 1
        assert abs(vals[a] - 1.0) < tol, (
            f"{element.__class__.__name__}: N[{a}] at its own node is {vals[a]}, not 1"
        )
        # All other nodes should be 0
        for b, val in enumerate(vals):
            if b == a:
                continue
            assert abs(val) < tol, (
                f"{element.__class__.__name__}: N[{b}] at node {a} is {val}, not 0"
            )


def _check_partition_of_unity(element, sample_points, tol=1e-12):
    """
    Check: sum_a N_a(xi, eta) = 1 at given sample points.
    """
    for (xi, eta) in sample_points:
        vals = [N(xi, eta) for N in element.N]
        s = sum(vals)
        assert abs(s - 1.0) < tol, (
            f"{element.__class__.__name__}: sum(N) at ({xi},{eta}) = {s}, not 1"
        )


def _sample_points_quad():
    """
    Some points inside [-1,1]x[-1,1] for quads.
    """
    return [
        (0.0, 0.0),
        (-0.5, -0.3),
        (0.7, -0.2),
        (-0.2, 0.8),
        (0.3, -0.9),
    ]


def _sample_points_t3():
    """
    Some points inside the reference triangle for t3: xi>=0, eta>=0, xi+eta<=1
    """
    return [
        (1.0 / 3.0, 1.0 / 3.0),
        (0.2, 0.2),
        (0.6, 0.2),
        (0.1, 0.7),
    ]


def _check_q7_left_edge_linearity(element, tol=1e-12):
    """
    For q7 (and similarly for q6) the left edge (xi = -1) is linear, shared with Q4.
    On that edge:
      - Only corner nodes 1 (LL) and 4 (UL) should be active.
      - N1 = 0.5 * (1 - eta), N4 = 0.5 * (1 + eta)
    This ensures compatibility with a Q4 neighbour.
    """
    for eta in [-1.0, -0.5, 0.0, 0.5, 1.0]:
        xi = -1.0
        vals = [N(xi, eta) for N in element.N]

        # Other nodes should be zero
        for i, v in enumerate(vals):
            if i not in (0, 3):  # nodes 1 and 4 in Python indexing
                assert abs(v) < tol, (
                    f"{element.__class__.__name__}: N[{i}] nonzero on left edge, "
                    f"value={v} at eta={eta}"
                )

        N1, N4 = vals[0], vals[3]
        N1_expected = 0.5 * (1.0 - eta)
        N4_expected = 0.5 * (1.0 + eta)

        assert abs(N1 - N1_expected) < tol, (
            f"{element.__class__.__name__}: N1({xi},{eta})={N1}, "
            f"expected {N1_expected}"
        )
        assert abs(N4 - N4_expected) < tol, (
            f"{element.__class__.__name__}: N4({xi},{eta})={N4}, "
            f"expected {N4_expected}"
        )
        # Also check they sum to 1 on the edge
        assert abs(N1 + N4 - 1.0) < tol, (
            f"{element.__class__.__name__}: N1+N4={N1+N4} on left edge, not 1"
        )


# ---------------------- Pytest-style test functions ---------------------- #

def test_q4_shape_functions():
    coords = [[0.0, 0.0]] * 4
    e = q4(coords)
    _check_kronecker_delta(e)
    _check_partition_of_unity(e, _sample_points_quad())


def test_q8_shape_functions():
    coords = [[0.0, 0.0]] * 8
    e = q8(coords)
    _check_kronecker_delta(e)
    _check_partition_of_unity(e, _sample_points_quad())


def test_q7_shape_functions():
    coords = [[0.0, 0.0]] * 7
    e = q7(coords)
    _check_kronecker_delta(e)
    _check_partition_of_unity(e, _sample_points_quad())
    _check_q7_left_edge_linearity(e)


def test_q6_shape_functions():
    coords = [[0.0, 0.0]] * 6
    e = q6(coords)
    _check_kronecker_delta(e)
    _check_partition_of_unity(e, _sample_points_quad())
    _check_q7_left_edge_linearity(e)  # same left-edge behaviour as q7


def test_t3_shape_functions():
    coords = [[0.0, 0.0]] * 3
    e = t3(coords)
    _check_kronecker_delta(e)
    _check_partition_of_unity(e, _sample_points_t3())


# =====================================================================
# Derivative validation tests (numerical check against N)
# =====================================================================

def _numeric_derivative(f, xi, eta, var="xi", h=1e-6):
    """
    Central finite difference derivative of f(xi, eta)
    with respect to var ('xi' or 'eta').
    """
    if var == "xi":
        return (f(xi + h, eta) - f(xi - h, eta)) / (2.0 * h)
    elif var == "eta":
        return (f(xi, eta + h) - f(xi, eta - h)) / (2.0 * h)
    else:
        raise ValueError("var must be 'xi' or 'eta'")


def _sample_points_quad_for_deriv():
    # A few interior-ish points in [-1,1]x[-1,1]
    return [
        (0.0, 0.0),
        (-0.3, -0.2),
        (0.4, -0.5),
        (-0.6, 0.3),
        (0.7, 0.6),
    ]


def _sample_points_t3_for_deriv():
    # Points inside reference triangle xi>=0, eta>=0, xi+eta<=1
    return [
        (0.2, 0.2),
        (0.3, 0.1),
        (0.1, 0.4),
        (0.4, 0.2),
    ]


def _check_numeric_derivatives(element, sample_points, tol=1e-7):
    """
    For each shape function N_i of an element, check that:

        dN_dxi_i(xi,eta)  ≈ ∂N_i/∂xi (central finite diff)
        dN_deta_i(xi,eta) ≈ ∂N_i/∂eta (central finite diff)

    at a bunch of sample points.
    """
    for xi, eta in sample_points:
        for i, N in enumerate(element.N):
            dNdxi_impl = element.dN_dxi[i](xi, eta)
            dNdeta_impl = element.dN_deta[i](xi, eta)

            dNdxi_num = _numeric_derivative(N, xi, eta, var="xi")
            dNdeta_num = _numeric_derivative(N, xi, eta, var="eta")

            assert abs(dNdxi_impl - dNdxi_num) < tol, (
                f"{element.__class__.__name__}: dN_dxi[{i}] mismatch at "
                f"(xi,eta)=({xi},{eta}). impl={dNdxi_impl}, num={dNdxi_num}"
            )
            assert abs(dNdeta_impl - dNdeta_num) < tol, (
                f"{element.__class__.__name__}: dN_deta[{i}] mismatch at "
                f"(xi,eta)=({xi},{eta}). impl={dNdeta_impl}, num={dNdeta_num}"
            )


# ---------------------- Pytest-style derivative tests ---------------------- #

def test_q4_derivatives():
    coords = [[0.0, 0.0]] * 4
    e = q4(coords)
    _check_numeric_derivatives(e, _sample_points_quad_for_deriv())


def test_q8_derivatives():
    coords = [[0.0, 0.0]] * 8
    e = q8(coords)
    _check_numeric_derivatives(e, _sample_points_quad_for_deriv())


def test_q7_derivatives():
    coords = [[0.0, 0.0]] * 7
    e = q7(coords)
    _check_numeric_derivatives(e, _sample_points_quad_for_deriv())


def test_q6_derivatives():
    coords = [[0.0, 0.0]] * 6
    e = q6(coords)
    _check_numeric_derivatives(e, _sample_points_quad_for_deriv())


def test_t3_derivatives():
    coords = [[0.0, 0.0]] * 3
    e = t3(coords)
    _check_numeric_derivatives(e, _sample_points_t3_for_deriv())


# Optionally update your __main__ runner to include derivative checks too:
def _run_all_tests():
    # existing shape-function tests you already had:
    test_q4_shape_functions()
    test_q8_shape_functions()
    test_q7_shape_functions()
    test_q6_shape_functions()
    test_t3_shape_functions()

    # new derivative tests:
    test_q4_derivatives()
    test_q8_derivatives()
    test_q7_derivatives()
    test_q6_derivatives()
    test_t3_derivatives()

    print("All shape-function and derivative tests passed.")


if __name__ == "__main__":
    _run_all_tests()
