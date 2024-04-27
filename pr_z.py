import numpy as np

# Analytic solution for maximum real root of cubic polynomial
# a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
def max_root(a):
    if a[0] != 1:
        a = np.array(a) / a[0] # Normalize to unity exponent for Z^3
    p = (3 * a[2]- a[1]**2)/3
    q = (2 * a[1]**3 - 9 * a[1] * a[2] + 27 * a[3])/27
    root_diagnostic = q**2/4 + p**3/27

    if root_diagnostic < 0:
        m = 2*np.sqrt(-p/3)
        qpm = 3*q/p/m
        theta1 = np.arccos(qpm)/3
        roots = np.array([m*np.cos(theta1), m*np.cos(theta1+4*np.pi/3), m*np.cos(theta1+2*np.pi/3)])
        Zs = roots - a[1] / 3
    else:
        P = (-q/2 + np.sqrt(root_diagnostic))
        if P >= 0:
            P = P **(1/3)
        else:
            P = -(-P)**(1/3)
            
        Q = (-q/2 - np.sqrt(root_diagnostic))
        if Q >=0:
            Q = Q **(1/3)
        else:
            Q = -(-Q)**(1/3)
        Zs = np.array([P + Q]) - a[1] / 3
    return max(Zs)
    
R, mwAir, degF2R = 10.731577089016, 28.97, 459.67

# Modified Sutton Tc and Pc
def tc_pc(sg_hc):
    x = [-4.01841583e-03, -2.75078362e+00, 7.07644217e+02, -8.86950457e-02, 1.29170902e+01, 1.57873612e+02]
    coefic_pc = x[:3]
    coefic_tc = x[3:]
    mw_gas = mwAir * sg_hc
    ppc_hc = coefic_pc[0] * mw_gas ** 2 + coefic_pc[1] * mw_gas + coefic_pc[2]
    tpc_hc = coefic_tc[0] * mw_gas ** 2 + coefic_tc[1] * mw_gas + coefic_tc[2] 
    return (tpc_hc, ppc_hc)    


names = ['CO2', 'H2S', 'N2', 'Hydrocarbon']
mws = np.array([44.01, 34.082, 28.014, 0])
tcs = np.array([547.416, 672.120, 227.160, 1])
pcs = np.array([1069.51, 1299.97, 492.84, 1])
ACF = np.array([0.22500, 0.09000, 0.03700, -0.04896])
VSHIFT = np.array([-0.27025, -0.11575, -0.20976, -0.39490])
OmegaA = np.array([0.441272, 0.457236, 0.457236, 0.429188])
OmegaB = np.array([0.0701277, 0.0777961, 0.0777961, 0.0692551]) 

kij = np.array([[0.0,       0.14345,   0.80000,   0.04506],
                [0.14345,   0.0,      -0.50085,   0.14667],
                [0.80000,  -0.50085,   0.0,      -0.24376],
                [0.04506,   0.14667,  -0.24376,   0.0]])   
                
def peng_robinson_z(degf, psia, sg, co2 = 0, h2s = 0, n2 = 0):
    if co2 + h2s + n2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 or sg < 0:
        return None
    degR = degf + degF2R
    z = np.array([co2, h2s, n2, 1 - co2 - h2s - n2])
    if n2 + co2 + h2s < 1:
        sg_hc = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2]) / mwAir) / (1 - co2 - h2s - n2)
    else:
        sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
    tcs[-1], pcs[-1] = tc_pc(sg_hc) # Hydrocarbon Tc and Pc from SG using Modified Sutton correlation
    trs = degR / tcs
    prs = psia / pcs
    
    m = 0.37464 + 1.54226 * ACF - 0.26992 * ACF**2
    alpha = (1 + m * (1 - np.sqrt(trs)))**2    
     
    Ai, Bi = OmegaA * alpha * prs / trs**2, OmegaB * prs / trs
    A, B = np.sum(z[:, None] * z * np.sqrt(np.outer(Ai, Ai)) * (1 - kij)), np.sum(z * Bi)

    # Coefficients of Cubic: a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
    coeffics = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]

    return max_root(coeffics) - np.sum(z * VSHIFT * Bi)  # Volume translated Z
    
    