import numpy as np

# Analytic solution for real root(s) of cubic polynomial
# a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
# Flag = 1 return Max root, = -1 returns minimum root, = 0 returns all real roots
def cubic_root(a, flag = 0):
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
        
    if flag == -1:      # Return minimum root
        return min(Zs)
    if flag == 1:       # Return maximum root
        return max(Zs)
    else:               # Return all roots
        return Zs    


# Modified Sutton Tc and Pc
R, mwAir, degF2R = 10.731577089016, 28.97, 459.67
def tc_pc(sg_hc):
    x = [-3.66048577e-03, -2.73719309e+00, 7.01192731e+02, -9.91619530e-02, 1.18144133e+01, 1.78887597e+02]
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
VCVIS = np.array([1.35480, 1.37424, 1.25466, 0]) # cuft/lbmol

def calc_bips(hc_mw, degf):
    degR = degf + degF2R
    intcpts = np.array([-4.427770E-01, -2.305710E-01, 2.657560E+00, 9.215660E-01, 8.314670E-03, -5.405930E+00])
    degR_slopes = np.array([3.741680E+02, 1.330530E+02, -1.601900E+03, 1.892310E+01, 2.093930E+01, 2.136390E+03])
    mw_slopes = np.array([0.000000E+00, 0.000000E+00, 0.000000E+00, -5.247030E-02, 5.168520E-03, 7.717860E-02])
    bips = intcpts + degR_slopes/degR + mw_slopes * hc_mw
    #                 CO2        H2S        N2         Gas
    return np.array([[0.0,      bips[0],   bips[1],   bips[3]],
                    [bips[0],   0.0,       bips[2],   bips[4]],
                    [bips[1],   bips[2],   0.0,       bips[5]],
                    [bips[3],   bips[4],   bips[5],   0.0]]) 
                
def peng_robinson_z(degf, psia, sg, co2 = 0, h2s = 0, n2 = 0):
    if co2 + h2s + n2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 :
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
    
    kij = calc_bips(sg_hc*mwAir, degf)
    
    m = 0.37464 + 1.54226 * ACF - 0.26992 * ACF**2
    alpha = (1 + m * (1 - np.sqrt(trs)))**2    
     
    Ai, Bi = OmegaA * alpha * prs / trs**2, OmegaB * prs / trs
    A, B = np.sum(z[:, None] * z * np.sqrt(np.outer(Ai, Ai)) * (1 - kij)), np.sum(z * Bi)

    # Coefficients of Cubic: a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
    coeffics = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]

    #      Return maximum real root
    return cubic_root(coeffics, flag = 1) - np.sum(z * VSHIFT * Bi)  # Volume translated Z
    
# From https://wiki.whitson.com/bopvt/visc_correlations/
def lbc(Z, degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0):
    
    if co2 + h2s + n2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 :
        return None
    degR = degf + degF2R
    zi = np.array([co2, h2s, n2, 1 - co2 - h2s - n2])
    if n2 + co2 + h2s < 1:
        sg_hc = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2]) / mwAir) / (1 - co2 - h2s - n2)
    else:
        sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
    
    hc_gas_mw = sg_hc * mwAir
        
    def vcvis_hc(mw): # Returns hydrocarbon gas VcVis for LBC viscosity calculations
        return  0.0532116692724726 * mw + 0.487478676470591 # ft3/lbmol

    
    mws[-1]  = hc_gas_mw
    tcs[-1], pcs[-1] = tc_pc(hc_gas_mw/mwAir)
    VCVIS[-1] = vcvis_hc(hc_gas_mw)
    degR = degf + degF2R

    def stiel_thodos(degR, mws):
        #Calculate the viscosity of a pure component using the Stiel-Thodos correlation.
        Tr = degR / tcs
        ui = []
        Tc = tcs * 5/9 # (deg K)
        Pc = pcs / 14.696
        eta = Tc**(1/6) / (mws**(1/2) * Pc**(2/3)) # Tc and Pc must be in degK and Atm respectively
        
        for i in range(len(Tr)):
            if Tr[i] <= 1.5:
                ui.append(34e-5 * Tr[i]**0.94 / eta[i])
            else:
                ui.append(17.78e-5 * (4.58 * Tr[i] - 1.67)**(5/8) / eta[i])
        return np.array(ui)
    
    def u0(zi, ui, mws, Z):  # dilute gas mixture viscosity from Herning and Zippener
        sqrt_mws = np.sqrt(mws)
        return np.sum(zi * ui * sqrt_mws)/np.sum(zi * sqrt_mws)
    
    a = [0.1023, 0.023364, 0.058533, -3.47874e-02, 7.74541e-03] # P3 and P4 have been modified
    # Calculate the viscosity of the mixture using the Lorenz-Bray-Clark method.
    rhoc = 1/np.sum(VCVIS*zi)
    Tc = tcs * 5/9    # (deg K)
    Pc = pcs / 14.696 # (Atm)
    
    eta = np.abs(np.sum(zi*Tc)**(1/6)) / (np.abs(np.sum(zi * mws))**0.5 * np.abs(np.sum(zi * Pc))**(2/3)) # Note 0.5 exponent from original paper
    mw = np.sum(zi * mws)
    rhor = psia / (Z * R * degR) / rhoc
    lhs = a[0] + a[1]*rhor + a[2]*rhor**2 + a[3]*rhor**3 + a[4]*rhor**4
    ui = stiel_thodos(degR, mws)
    vis = (lhs**4 - 1e-4)/eta + u0(zi, ui, mws, Z)
    return vis    