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

    coefic_pc = [-1.18331,1104.13,0.00702049875,1.50107866]
    coefic_tc = [4.76601,45.4709,0.00244823,0.3219414605]
    mw_gas = mwAir * sg_hc
    ppc_hc = (coefic_pc[0] * mw_gas + coefic_pc[1])/(coefic_pc[2] * mw_gas + coefic_pc[3])
    tpc_hc = (coefic_tc[0] * mw_gas + coefic_tc[1])/(coefic_tc[2] * mw_gas + coefic_tc[3])
    return (tpc_hc, ppc_hc)    

names = ['CO2', 'H2S', 'N2', 'H2', 'Hydrocarbon']
mws = np.array([44.01, 34.082, 28.014, 2.016, 0])
tcs = np.array([547.416, 672.120, 227.160, 59.3640, 1])
pcs = np.array([1069.51, 1299.97, 492.84, 187.5300, 1])
ACF = np.array([0.10560, 0.01839, 0.037, -0.21700, -0.04051])
VSHIFT = np.array([-0.19240, -0.18906, -0.21067, -0.12899, -0.19185])
OmegaA = np.array([0.440853, 0.441796, 0.457236, 0.457236, 0.457236])
OmegaB = np.array([0.0730166, 0.0739200, 0.0777961, 0.0777961, 0.0777961]) 
VCVIS = np.array([1.55342, 1.57142, 1.41785, 1.18168, 0]) # cuft/lbmol    
                                                                                           
def calc_bips(hc_mw, degf):                                                                     
    degR = degf + degF2R  

    # Hydrocarbon-Inert BIPS (Regressed to Wichert & Synthetic GERG Data)
    # BIP = intcpt + degR_slope/degR + mw_slope * hc_mw
    #                      CO2           H2S          N2          H2 
    intcpts = np.array([0.356872, 0.13333, -0.0429868, 0.712742])
    mw_slopes = np.array([-0.00128858, -0.00141921, 0.000728009, 0.00110919])
    degR_slopes = np.array([-179.465, -24.7713, -3.39882, -349.25])
    hc_bips = list(intcpts + degR_slopes/degR + mw_slopes * hc_mw)
   
    # Inert:Inert BIP Pairs                                                                        
    slopes = [[-0.000811630909090909, -0.00112389, 0.00113879071428571, -3.949E-05, 0.00377893, 0.000227832857142857],                                
              [0.000199935714285714, 0.000648015753811406, 0.00193186545454545, 0.00116141530208922, 0.00119598402597403, 0.000141143636363636]]                                    
    intcpts = [[0.244400581818182, -0.1389305, 0.251688807142857, 0.661337, 0.5633018, 0.199822661904762],                                   
               [0.0843990523809524, -0.294767026538679, 0.0794874545454544, 0.530504334274421, 0.831735004329004, 0.216575909090909]]                               
    evals = [np.max, np.max, np.max, np.max, np.min, np.min]                                                                  
    pairs = ['H2S with CO2', 'CO2 with N2', 'H2S WITH N2', 'CO2 with H2', 'H2S WITH H2', 'N2 WITH H2']  
                        
    inert_bips = [evals[i]([slopes[0][i]*degf + intcpts[0][i], slopes[1][i]*degf + intcpts[1][i]]) for i in range(len(evals))]
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
        return  0.070294775 * mw + 0.333582429 # ft3/lbmol                                

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

    a = [0.1023, 0.023364, 0.058533, -4.01411e-02, 9.00617e-03] # P3 and P4 have been modified
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
