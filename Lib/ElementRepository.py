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
            
    def Nvals(self, element, xi, eta):
        return np.array([N(xi, eta) for N in element.N])

    def dN_dxi_vals(self, element, xi, eta):
        return np.array([dN(xi, eta) for dN in element.dN_dxi])

    def dN_deta_vals(self, element, xi, eta):
        return np.array([dN(xi, eta) for dN in element.dN_deta]) 
    
    # strain displacement relation, relating the derivatives of displacement with resepct to the dofs to the components of strain
    def B1(self): 
        B1 = np.array([[1, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 1, 1, 0]])
        return B1
    
    # Scalling matrix containing the inverse jacobian transposed
    def B2(self, jacobian):
        # computting the inverse transposed of the jacobian matrix
        jinvT = jacobian.invT # responsible for all scalling and rotations
        B2 = np.zeros([4,4]) # this is the same regardles of element order
        # Filling in the scalling matrix
        B2[0:2, 0:2] = jinvT
        B2[2:4, 2:4] = jinvT
        return B2
    
    def B3(self, xi, eta):
        e = self.element
        B3 = np.zeros([4, e.totalDof])
        B3[0, 0::2] = self.dN_dxi_vals(e, xi, eta)
        B3[1, 0::2] = self.dN_deta_vals(e, xi, eta)
        B3[2, 1::2] = self.dN_dxi_vals(e, xi, eta)
        B3[3, 1::2] = self.dN_deta_vals(e, xi, eta)
        return B3

    def B(self, xi, eta, jacb):
            """
            Build the standard 3 x (2*n) B-matrix for 2D plane stress/strain:

                [ εx ]   [ dN1/dx  0   ... dNn/dx   0   ] [u1]
                [ εy ] = [ 0      dN1/dy ... 0     dNn/dy] [v1]
                [ γxy]   [ dN1/dy dN1/dx ... dNn/dy dNn/dx] [vn]
            """
            e = self.element
            n = e.nodeCount
            totalDof = e.totalDof

            # Derivatives w.r.t natural coords
            dN_dxi  = self.dN_dxi_vals(e, xi, eta)   # (n,)
            dN_deta = self.dN_deta_vals(e, xi, eta)  # (n,)

            # J^{-1} from stored invT (invT = (J^{-1})^T)
            invJ = jacb.invT.T                       # 2x2

            # [dN/dx; dN/dy] = invJ @ [dN/dξ; dN/dη]
            grads_nat = np.vstack((dN_dxi, dN_deta))  # 2 x n
            grads_xy  = invJ @ grads_nat             # 2 x n

            dN_dx = grads_xy[0, :]
            dN_dy = grads_xy[1, :]

            B = np.zeros((3, totalDof))
            for i in range(n):
                col_u = 2 * i
                col_v = 2 * i + 1

                B[0, col_u] = dN_dx[i]  # εx = Σ dNi/dx * ui
                B[1, col_v] = dN_dy[i]  # εy = Σ dNi/dy * vi
                B[2, col_u] = dN_dy[i]  # γxy = Σ dNi/dy * ui + dNi/dx * vi
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
                  lambda xi, eta: 0.25*(xi**2*eta + xi**2 + xi*eta**2 - xi*eta + eta**2 - 1),
                  lambda xi, eta: 0.25*(xi**2*eta + xi**2 + xi*eta**2 + xi*eta + eta**2 - 1),
                  lambda xi, eta: 0.25*xi*(xi-1)*(eta+1),
                  lambda xi, eta: 0.5*(xi**2-1)*(eta-1),
                  lambda xi, eta: -0.5*(xi+1)*(eta**2-1),
                  lambda xi, eta: -0.5*(xi**2-1)*(eta+1)]

        self.dN_dxi = [lambda xi, eta: 0.25*eta-0.5*xi*(eta-1)-0.25,
                       lambda xi, eta: 0.25*eta**2 - 0.25*eta-0.5*xi*(eta-1),
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