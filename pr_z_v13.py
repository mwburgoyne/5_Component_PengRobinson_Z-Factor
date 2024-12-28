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

R, mwAir, degF2R = 10.731577089016, 28.97, 459.67

# Custom Tc and Pc
def tc_pc(sg_hc):
    v1, v2, v3 = 0.000229975, 0.186415901, 2.470903632
    offset, pl, vl = -0.032901049, 42.6061669, 3007.108548
 
    mw_gas = mwAir * sg_hc
    vc_on_zc = v1 * mw_gas**2 + v2*mw_gas + v3
    tpc_hc = (offset + vc_on_zc) * vl / (offset + vc_on_zc + pl) 
    ppc_hc = tpc_hc * R / vc_on_zc
    
    return (tpc_hc, ppc_hc)    
    
names = ['CO2', 'H2S', 'N2', 'H2', 'Hydrocarbon']
mws = np.array([44.01, 34.082, 28.014, 2.016, 0])
tcs = np.array([547.416, 672.120, 227.160, 47.430, 1]) # H2 Tc has been modified
pcs = np.array([1069.51, 1299.97, 492.84, 187.5300, 1])
ACF = np.array([0.12256, 0.04916, 0.037, -0.21700, -0.03899])
VSHIFT = np.array([-0.27593, -0.22896, -0.21066, -0.32400, -0.19076])
OmegaA = np.array([0.427705, 0.436743, 0.457236, 0.457236, 0.457236])
OmegaB = np.array([0.0696460, 0.0724373, 0.0777961, 0.0777961, 0.0777961]) 
VCVIS = np.array([1.46020, 1.46460, 1.35422, 0.67967, 0]) # cuft/lbmol    
    
def calc_bips(hc_mw, degf):                                                                     
    degR = degf + degF2R  

    # Hydrocarbon-Inert BIPS (Regressed to Wichert & Synthetic GERG Data)
    # BIP = intcpt + degR_slope/degR + mw_slope * hc_mw
    #                      CO2      H2S        N2        H2 
    intcpts = np.array([0.386557, 0.267007, 0.486589, 0.776917])
    mw_slopes = np.array([-0.00219806, -0.00396541, -0.00316789, 0.0106061])
    degR_slopes = np.array([-158.333, -58.611, -226.239, -474.283])

    hc_bips = list(intcpts + degR_slopes/degR + mw_slopes * hc_mw)
    
    # Inert:Inert BIP Pairs
    #            CO2:H2S       CO2:N2     H2S:N2    CO2:H2  H2S:H2  N2:H2
    inert_bips = [0.0600319, -0.229807, -0.18346, 0.646796, 0.65, 0.369087]
    bips = np.array(hc_bips + inert_bips)
    bip_pairs = [(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]  
    bip_matrix = np.zeros((5, 5))
    
    for i in range(5):
        for j in range(5):
            for p, pair in enumerate(bip_pairs):
                if (i, j) == pair:
                    bip_matrix[i, j] = bips[p]
                    bip_matrix[j, i] = bips[p]
                    continue
    return bip_matrix
                
def peng_robinson_z(degf, psia, sg, co2 = 0, h2s = 0, n2 = 0, h2 = 0):
    if co2 + h2s + n2 + h2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 or h2 < 0:
        return None
    degR = degf + degF2R
    z = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
    if n2 + co2 + h2s + h2 < 1:
        sg_hc = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2] + h2 * mws[3]) / mwAir) / (1 - co2 - h2s - n2 - h2)
    else:
        sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
    sg_hc = max(sg_hc, 0.553779772) # Methane is lower limit
    hc_mw = sg_hc*mwAir
    
    tcs[-1], pcs[-1] = tc_pc(sg_hc) # Hydrocarbon Tc and Pc from SG using custom correlation
    trs = degR / tcs
    prs = psia / pcs
    
    kij = calc_bips(hc_mw, degf)
    
    m = 0.37464 + 1.54226 * ACF - 0.26992 * ACF**2
    alpha = (1 + m * (1 - np.sqrt(trs)))**2    
     
    Ai, Bi = OmegaA * alpha * prs / trs**2, OmegaB * prs / trs
    A, B = np.sum(z[:, None] * z * np.sqrt(np.outer(Ai, Ai)) * (1 - kij)), np.sum(z * Bi)

    # Coefficients of Cubic: a[0] * Z**3 + a[1]*Z**2 + a[2]*Z + a[3] = 0
    coeffics = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]
    
    #      Return maximum real root
    return cubic_root(coeffics, flag = 1) - np.sum(z * VSHIFT * Bi)  # Volume translated Z


# From https://wiki.whitson.com/bopvt/visc_correlations/
def lbc(Z, degf, psia, sg, co2=0.0, h2s=0.0, n2=0.0, h2 = 0.0):
    
    if co2 + h2s + n2 + h2 > 1 or co2 < 0 or h2s < 0 or n2 < 0 or h2 < 0:
        return None
    degR = degf + degF2R
    zi = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
    if n2 + co2 + h2s + h2 < 1:
        sg_hc = (sg - (co2 * mws[0] + h2s * mws[1] + n2 * mws[2] + h2 * mws[3]) / mwAir) / (1 - co2 - h2s - n2 - h2)
    else:
        sg_hc = 0.75 # Irrelevant, since hydrocarbon fraction = 0
    
    sg_hc = max(sg_hc, 0.553779772) # Methane is lower limit
    
    hc_gas_mw = sg_hc * mwAir
        
    def vcvis_hc(mw): # Returns hydrocarbon gas VcVis for LBC viscosity calculations   
        return  0.057511062 *  mw + 0.478400158 # ft3/lbmol      
                                            
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

    a = [0.1023, 0.023364, 0.058533, -3.92835e-02,  9.28591e-03] # P3 and P4 have been modified
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
