' bns_VBA_Module1.bas
'
' Module1 exported from '5-Component Universal EOS.xlsm', with the Peneloux volume
' shift carried into the caloric path (see the c_mix, H_vshift_Btu and V_shifted
' additions below). Exported so the VBA is reviewable and diffable in git rather
' than living only inside the workbook binary.
'
' This file and the workbook are in step as of 15 August 2026. Excel owns the copy
' inside the .xlsm, so keep them in step by hand: edit in the VBE, then export
' Module1 over this file. Patching the workbook from outside Excel is not safe,
' because Excel can keep running the cached p-code and ignore patched source.
'
Attribute VB_Name = "Module1"
Option Explicit
'===================================================================
' Purpose: Computes Peng–Robinson Z, thermo properties (H, Cp, Cv, mu_JT),
'          and LBC viscosity.
' Author: Mark Burgoyne, June 2025
'
' Content described in ADIPEC paper SPE-229932-MS, 2025 titled:
' "A Universal, EOS-Based Correlation for Z-Factor, Viscosity and Enthalpy For Hydrocarbon and H2, N2, CO2, H2S Gas Mixtures"
' M. W. Burgoyne, Santos; M. H. Nielsen, Whitson AS; M. Stanko, Whitson AS
'===================================================================

'--- Module-level constants (field units) ---
Private Const R_field       As Double = 10.731577089016    ' ft³·psia/(lb-mol·°R)
Private Const R_THERMO      As Double = 1.98588            ' Btu/(lb-mol·°R)
Private Const MW_AIR        As Double = 28.97              ' lb/lb-mol
Private Const DEG_F_TO_R    As Double = 459.67             ' °F to °R shift
Private Const FT3_PSIA_TO_BTU As Double = 5.403            ' conversion ft³·psia - Btu

'--- Pure-fluid constants (1=CO2, 2=H2S, 3=N2, 4=H2, 5=Gas) ---
Private mws(1 To 5)      As Double
Private tcs(1 To 5)      As Double
Private pcs(1 To 5)      As Double
Private ACF(1 To 5)      As Double
Private VSHIFT(1 To 5)   As Double
Private OmegaA(1 To 5)   As Double
Private OmegaB(1 To 5)   As Double

'--- Refitted dimensionless Cp-polynomial coefficients (1=CO2 … 5=Gas) ---
Private Cp_a_(1 To 5)    As Double
Private Cp_b_(1 To 5)    As Double
Private Cp_c_(1 To 5)    As Double
Private Cp_d_(1 To 5)    As Double
Private Cp_e_(1 To 5)    As Double

'--- VCVIS array for LBC (must be filled with the same five values used in Python) ---
Private VCVIS(1 To 5)    As Double


'===================================================================
'   Initialize all module-level arrays (run once per session)
'===================================================================
Private Sub InitConstants()

    ' 1) Molar masses (lb-mol basis)
    mws(1) = 44.01:   mws(2) = 34.082:  mws(3) = 28.014
    mws(4) = 2.016:   mws(5) = 0#

    ' 2) Critical temperatures (°R) and pressures (psia)
    tcs(1) = 547.416: tcs(2) = 672.12: tcs(3) = 227.16
    tcs(4) = 47.43:   tcs(5) = 1#       ' placeholder for "Gas"
    pcs(1) = 1069.51: pcs(2) = 1299.97: pcs(3) = 492.84
    pcs(4) = 187.53:  pcs(5) = 1#       ' placeholder for "Gas"

    ' 3) Acentric factors
    ACF(1) = 0.12256: ACF(2) = 0.04909: ACF(3) = 0.037
    ACF(4) = -0.217:  ACF(5) = -0.03899

    ' 4) Volume shift parameters
    VSHIFT(1) = -0.27607: VSHIFT(2) = -0.22901
    VSHIFT(3) = -0.21066: VSHIFT(4) = -0.3627
    VSHIFT(5) = -0.19076

    ' 5) OmegaA and OmegaB for Peng–Robinson
    OmegaA(1) = 0.427671: OmegaA(2) = 0.436725
    OmegaA(3) = 0.457236: OmegaA(4) = 0.457236
    OmegaA(5) = 0.457236
    OmegaB(1) = 0.0696397: OmegaB(2) = 0.0724345
    OmegaB(3) = 0.0777961: OmegaB(4) = 0.0777961
    OmegaB(5) = 0.0777961

    ' 6) Refitted Cp polynomial coefficients (dimensionless)
    Cp_a_(1) = 2.725473196:  Cp_b_(1) = 0.004103751
    Cp_c_(1) = 0.000015602:  Cp_d_(1) = -0.0000000419321
    Cp_e_(1) = 3.10542E-11

    Cp_a_(2) = 4.446031265:  Cp_b_(2) = -0.005296052
    Cp_c_(2) = 0.000020533:  Cp_d_(2) = -0.0000000258993
    Cp_e_(2) = 1.25555E-11

    Cp_a_(3) = 3.423811591:  Cp_b_(3) = 0.001007461
    Cp_c_(3) = -0.00000458491: Cp_d_(3) = 0.0000000084252
    Cp_e_(3) = -4.38083E-12

    Cp_a_(4) = 1.421468418:  Cp_b_(4) = 0.018192108
    Cp_c_(4) = -0.0000604285: Cp_d_(4) = 0.0000000908033
    Cp_e_(4) = -5.18972E-11

    Cp_a_(5) = 5.369051342:  Cp_b_(5) = -0.014851371
    Cp_c_(5) = 0.0000486358: Cp_d_(5) = -0.0000000370187
    Cp_e_(5) = 1.80641E-12

    ' 7) VCVIS array
    VCVIS(1) = 1.46352
    VCVIS(2) = 1.46808
    VCVIS(3) = 1.35526
    VCVIS(4) = 0.68473
    VCVIS(5) = 0#
    
End Sub

'===================================================================
'   1) tc_pc: Given sg_hc and AG (True for Associated Gas),
'      return tpc_hc (°R) and ppc_hc (psia)
'===================================================================
Private Sub tc_pc(ByVal sg_hc As Double, ByVal AG As Boolean, ByRef tpc_hc As Double, ByRef ppc_hc As Double)
    Dim a As Double, b As Double, c As Double
    Dim SlopeVc As Double
    Dim mw_gas As Double, X As Double, vc_on_zc As Double

    Const MW_AIR As Double = 28.9625 ' (replace with your value if different)
    'Const R_field As Double = 10.731577089016 ' (psia·ft³)/(lb·mol·°R)

    mw_gas = MW_AIR * sg_hc
    X = mw_gas - 16.0425

    If AG Then
        ' Associated Gas
        a = 2695.14765: b = 274.341701: c = 343.008: SlopeVc = 0.177497835
    Else
        ' Gas Condensate
        a = 1098.10948: b = 101.529237: c = 343.008: SlopeVc = 0.170931432
    End If

    ' Critical Temperature
    tpc_hc = a * X / (b + X) + c
    ' Linear Vc/Zc
    vc_on_zc = SlopeVc * X + 5.51852587149056
    ' Critical Pressure
    ppc_hc = tpc_hc * R_field / vc_on_zc

End Sub


'===================================================================
'   2) calc_bips: Build the 5×5 BIP matrix (+ derivatives)
'===================================================================
Private Sub calc_bips( _
        ByVal degR As Double, _
        ByVal tpc_hc As Double, _
        ByVal Thermo As Boolean, _
        ByRef kij() As Double, _
        ByRef dkij_dT() As Double, _
        ByRef d2kij_dT2() As Double)

    Dim components(1 To 5) As String
    Dim i As Long, j As Long
    Dim bip_constant As Double, bip_Tr_slope As Double, bip_tc As Double
    Dim Tr As Double, rawKey As String, Key As String, revKey As String

    components(1) = "CO2": components(2) = "H2S"
    components(3) = "N2":  components(4) = "H2"
    components(5) = "Gas"

    ReDim kij(1 To 5, 1 To 5)
    If Thermo Then
        ReDim dkij_dT(1 To 5, 1 To 5)
        ReDim d2kij_dT2(1 To 5, 1 To 5)
    End If

    ' Initialize everything to zero
    For i = 1 To 5
        For j = 1 To 5
            kij(i, j) = 0#
            If Thermo Then
                dkij_dT(i, j) = 0#: d2kij_dT2(i, j) = 0#
            End If
        Next j
    Next i

    ' Loop only over i < j; fill both [i,j] and [j,i]
    For i = 1 To 5
        For j = i + 1 To 5

            rawKey = components(i) & "_" & components(j)
            revKey = components(j) & "_" & components(i)

            ' Decide which case matches: try rawKey first, then revKey
            Select Case True
                Case (rawKey = "Gas_CO2" Or revKey = "Gas_CO2")
                    bip_constant = -0.145561: bip_Tr_slope = 0.276572: bip_tc = tpc_hc
                Case (rawKey = "Gas_H2S" Or revKey = "Gas_H2S")
                    bip_constant = 0.16852: bip_Tr_slope = -0.122378: bip_tc = tpc_hc
                Case (rawKey = "Gas_N2" Or revKey = "Gas_N2")
                    bip_constant = -0.108: bip_Tr_slope = 0.0605506: bip_tc = tpc_hc
                Case (rawKey = "Gas_H2" Or revKey = "Gas_H2")
                    bip_constant = -0.0620119:  bip_Tr_slope = 0.0427873: bip_tc = tpc_hc

                Case (rawKey = "CO2_H2S" Or revKey = "CO2_H2S")
                    bip_constant = 0.248638: bip_Tr_slope = -0.138185: bip_tc = 547.416
                Case (rawKey = "CO2_N2" Or revKey = "CO2_N2")
                    bip_constant = -0.25:      bip_Tr_slope = 0.11602: bip_tc = 547.416
                Case (rawKey = "CO2_H2" Or revKey = "CO2_H2")
                    bip_constant = -0.247153: bip_Tr_slope = 0.16377: bip_tc = 547.416

                Case (rawKey = "H2S_N2" Or revKey = "H2S_N2")
                    bip_constant = -0.204414:  bip_Tr_slope = 0.234417: bip_tc = 672.12
                Case (rawKey = "H2S_H2" Or revKey = "H2S_H2")
                    bip_constant = 0#:        bip_Tr_slope = 0#:        bip_tc = 672.12

                Case (rawKey = "N2_H2" Or revKey = "N2_H2")
                    bip_constant = -0.166253: bip_Tr_slope = 0.0788129: bip_tc = 227.16

                Case Else
                    bip_constant = 0#: bip_Tr_slope = 0#: bip_tc = 0#
            End Select

            ' If bip_tc is still zero, we skip dividing
            If bip_tc <> 0# Then
                Tr = degR / bip_tc
                kij(i, j) = bip_constant + bip_Tr_slope / Tr
            Else
                kij(i, j) = 0#
            End If
            kij(j, i) = kij(i, j)

            If Thermo Then
                If bip_tc <> 0# Then
                    dkij_dT(i, j) = -bip_Tr_slope * bip_tc / (degR ^ 2)
                    d2kij_dT2(i, j) = 2# * bip_Tr_slope * bip_tc / (degR ^ 3)
                Else
                    dkij_dT(i, j) = 0#
                    d2kij_dT2(i, j) = 0#
                End If
                dkij_dT(j, i) = dkij_dT(i, j)
                d2kij_dT2(j, i) = d2kij_dT2(i, j)
            End If

        Next j
    Next i
End Sub


' -----------------------------------------------------------------------
' 1) CUBIC ROOT SOLVER
' -----------------------------------------------------------------------
'   Solves a cubic polynomial a(0)*Z^3 + a(1)*Z^2 + a(2)*Z + a(3) = 0.
'   - If flag= 1, returns maximum real root
'   - If flag=-1, returns minimum real root
'   - If flag= 0, returns all real roots (Variant array)
'
' Though not used directly by the user, you can call it in a sheet:
'   =CubicRoot(A1:A4, 1)
'   (where A1:A4 hold the polynomial coefficients) -> max real root
' -----------------------------------------------------------------------
Public Function CubicRoot(aCoeffs As Range, Optional ByVal flag As Long = 0) As Variant
    Dim i As Long
    If aCoeffs.Cells.Count <> 4 Then
        CubicRoot = CVErr(xlErrRef)
        Exit Function
    End If
    
    ' Copy input from the range
    Dim a(0 To 3) As Double
    For i = 0 To 3
        a(i) = aCoeffs.Cells(i + 1, 1).Value
    Next i
    
    ' Solve
    CubicRoot = CubicRootArrayDirect(a, flag)
End Function

' -----------------------------------------------------------------------
' Helper for CubicRoot: does the actual math on a 4-element Double array
' -----------------------------------------------------------------------
Private Function CubicRootArrayDirect(a() As Double, ByVal flag As Long) As Variant
    ' Ensure leading coefficient is non-zero
    If a(0) = 0# Then
        CubicRootArrayDirect = CVErr(xlErrDiv0)
        Exit Function
    End If
    
    ' Normalize leading coefficient to 1 if needed
    If a(0) <> 1# Then
        Dim lead As Double
        lead = a(0)
        Dim j As Long
        For j = 1 To 3
            a(j) = a(j) / lead
        Next j
        a(0) = 1#
    End If
    
    Dim a0 As Double, a1 As Double, a2 As Double, a3 As Double
    a0 = a(0) ' should be 1
    a1 = a(1)
    a2 = a(2)
    a3 = a(3)
    
    ' Depressed cubic: p, q
    Dim P As Double, Q As Double
    P = (3# * a2 - a1 ^ 2) / 3#
    Q = (2# * a1 ^ 3 - 9# * a1 * a2 + 27# * a3) / 27#
    
    ' Discriminant
    Dim disc As Double
    disc = (Q ^ 2 / 4#) + (P ^ 3 / 27#)
    
    Dim roots(1 To 3) As Double
    Dim nRoots As Long
    
    If disc < 0# Then
        ' Three real roots
        Dim m As Double
        m = 2# * Sqr(-P / 3#)
        
        Dim tmp As Double
        tmp = P * m
        Dim qpm As Double
        If Abs(tmp) < 0.00000000000001 Then
            qpm = 0#
        Else
            qpm = 3# * Q / tmp
        End If
        If qpm > 1# Then qpm = 1#
        If qpm < -1# Then qpm = -1#
        
        Dim theta As Double
        theta = Application.WorksheetFunction.Acos(qpm) / 3#
        
        roots(1) = m * Cos(theta)
        roots(2) = m * Cos(theta + 4# * Atn(1#) * 2# / 3#)
        roots(3) = m * Cos(theta + 2# * Atn(1#) * 2# / 3#)
        nRoots = 3
        
        Dim shift As Double
        shift = a1 / 3#
        Dim k As Long
        For k = 1 To 3
            roots(k) = roots(k) - shift
        Next k
    Else
        ' One real root
        Dim temp1 As Double, temp2 As Double
        temp1 = -Q / 2# + Sqr(disc)
        temp2 = -Q / 2# - Sqr(disc)
        
        P = cubicSignCbrt(temp1)
        Q = cubicSignCbrt(temp2)
        
        roots(1) = P + Q - (a1 / 3#)
        nRoots = 1
    End If
    
    ' Decide which root(s) to return
    Dim retVal As Variant
    Select Case flag
        Case -1 ' min
            Dim mn As Double
            mn = roots(1)
            Dim h As Long
            For h = 2 To nRoots
                If roots(h) < mn Then mn = roots(h)
            Next h
            retVal = mn
        Case 1  ' max
            Dim mx As Double
            mx = roots(1)
            Dim hh As Long
            For hh = 2 To nRoots
                If roots(hh) > mx Then mx = roots(hh)
            Next hh
            retVal = mx
        Case Else
            If nRoots = 1 Then
                retVal = roots(1)
            Else
                Dim arrA() As Variant
                ReDim arrA(1 To nRoots, 1 To 1)
                Dim hh2 As Long
                For hh2 = 1 To nRoots
                    arrA(hh2, 1) = roots(hh2)
                Next hh2
                retVal = arrA
            End If
    End Select
    
    CubicRootArrayDirect = retVal
End Function

' -----------------------------------------------------------------------
' sign-based cubic root helper
' -----------------------------------------------------------------------
Private Function cubicSignCbrt(X As Double) As Double
    If X >= 0# Then
        cubicSignCbrt = X ^ (1# / 3#)
    Else
        cubicSignCbrt = -(Abs(X) ^ (1# / 3#))
    End If
End Function

'===================================================================
'   4) Compute_da_mix_dT: returns a_i(), da_i_dT(), da_mix_dT
'===================================================================
Private Sub Compute_da_mix_dT( _
        ByRef zf As Variant, _
        ByRef a_c As Variant, _
        ByRef m_arr As Variant, _
        ByRef Tc_arr As Variant, _
        ByRef Tr_arr As Variant, _
        ByRef kij As Variant, _
        ByRef dkij_dT As Variant, _
        ByRef a_i As Variant, _
        ByRef da_i_dT As Variant, _
        ByRef da_mix_dT As Double)

    Dim i As Long, j As Long
    Dim sqrt_tr(1 To 5)        As Double
    Dim u_arr(1 To 5)          As Double
    Dim d_alpha_dTr(1 To 5)    As Double
    Dim d_alpha_dT_elem(1 To 5) As Double
    Dim N_mat(1 To 25)         As Double
    Dim sqrt_ai_aj             As Double
    Dim d_aij_dT               As Double
    Dim tempSum                As Double

    ' --- Step 1: build a_i and da_i_dT arrays ---
    For i = 1 To 5
        sqrt_tr(i) = Sqr(Tr_arr(i))
        ' u_i = 1 + m_i * (1 - sqrt_tr)
        u_arr(i) = 1# + m_arr(i) * (1# - sqrt_tr(i))
        ' da/dT_r = -m_i * u_i / sqrt_tr
        If sqrt_tr(i) = 0 Then
            d_alpha_dTr(i) = 0#
        Else
            d_alpha_dTr(i) = -m_arr(i) * u_arr(i) / sqrt_tr(i)
        End If
        ' chain-rule into da/dT
        d_alpha_dT_elem(i) = d_alpha_dTr(i) / Tc_arr(i)
        ' a_i = u_i^2, so a_i = a_c * a_i
        a_i(i) = a_c(i) * (u_arr(i) ^ 2)
        da_i_dT(i) = a_c(i) * d_alpha_dT_elem(i)
    Next i

    ' --- Step 2: populate N_mat = da_i_dT(i)*a_i(j) + a_i(i)*da_i_dT(j) ---
    For i = 1 To 5
        For j = 1 To 5
            N_mat((i - 1) * 5 + j) = da_i_dT(i) * a_i(j) + a_i(i) * da_i_dT(j)
        Next j
    Next i

    ' --- Step 3: assemble da_mix_dT using the exact same formula as Python ---
    tempSum = 0#
    For i = 1 To 5
        For j = 1 To 5
            sqrt_ai_aj = Sqr(a_i(i) * a_i(j))
            If sqrt_ai_aj = 0 Then
                d_aij_dT = 0#
            Else
                d_aij_dT = -dkij_dT(i, j) * sqrt_ai_aj _
                           + (1# - kij(i, j)) * 0.5 * (N_mat((i - 1) * 5 + j) / sqrt_ai_aj)
            End If
            tempSum = tempSum + zf(i) * zf(j) * d_aij_dT
        Next j
    Next i

    da_mix_dT = tempSum

End Sub




'===================================================================
'   5) Compute_d2a_mix_dT2: returns d2a_mix_dT2, d2a_i_dT2()
'===================================================================
Private Sub Compute_d2a_mix_dT2( _
    ByRef zf As Variant, _
    ByRef a_c As Variant, _
    ByRef m_arr As Variant, _
    ByRef Tc_arr As Variant, _
    ByRef Tr_arr As Variant, _
    ByRef kij As Variant, _
    ByRef dkij_dT As Variant, _
    ByRef d2kij_dT2 As Variant, _
    ByRef d2a_mix_dT2 As Double, _
    ByRef d2a_i_dT2 As Variant _
)

    Dim i As Long, j As Long
    Dim sqrt_tr(1 To 5)     As Double
    Dim u_arr(1 To 5)       As Double
    Dim alpha_arr(1 To 5)   As Double
    Dim da_dTr(1 To 5)      As Double
    Dim du_dTr(1 To 5)      As Double
    Dim d2u_dTr2(1 To 5)    As Double
    Dim d2alpha_dTr2(1 To 5) As Double
    Dim d_alpha_dT(1 To 5)  As Double
    Dim d2alpha_dT2(1 To 5) As Double

    Dim a_i_local(1 To 5)        As Double
    Dim da_i_dT_local(1 To 5)    As Double
    Dim d2a_i_dT2_local(1 To 5)  As Double

    Dim N_mat(1 To 25)        As Double
    Dim cross_2_da(1 To 25)    As Double
    Dim d2a_term(1 To 25)     As Double

    Dim sqrt_ai_aj As Double
    Dim term1 As Double, term2 As Double, term3 As Double
    Dim tempSum As Double

    ' First, compute each component's a, da_dT, d²a_dT²
    For i = 1 To 5
        sqrt_tr(i) = Sqr(Tr_arr(i))            ' Tr_arr is now a Variant containing [1..5] array
        u_arr(i) = 1# + m_arr(i) * (1# - sqrt_tr(i))
        alpha_arr(i) = u_arr(i) ^ 2

        If sqrt_tr(i) = 0 Then
            da_dTr(i) = 0#
            du_dTr(i) = 0#
            d2u_dTr2(i) = 0#
        Else
            da_dTr(i) = -m_arr(i) * u_arr(i) / sqrt_tr(i)
            du_dTr(i) = -m_arr(i) / (2# * sqrt_tr(i))
            d2u_dTr2(i) = m_arr(i) / (4# * (Tr_arr(i) ^ 1.5))
        End If

        d2alpha_dTr2(i) = 2# * (du_dTr(i) ^ 2) + 2# * u_arr(i) * d2u_dTr2(i)

        If Tc_arr(i) = 0 Then
            d2alpha_dT2(i) = 0#
            d_alpha_dT(i) = 0#
        Else
            d2alpha_dT2(i) = d2alpha_dTr2(i) / (Tc_arr(i) ^ 2)
            d_alpha_dT(i) = da_dTr(i) / Tc_arr(i)
        End If

        a_i_local(i) = a_c(i) * alpha_arr(i)
        da_i_dT_local(i) = a_c(i) * d_alpha_dT(i)
        d2a_i_dT2_local(i) = a_c(i) * d2alpha_dT2(i)
    Next i

    ' Build the N_mat, cross_2_da, and d2a_term matrices
    For i = 1 To 5
        For j = 1 To 5
            N_mat((i - 1) * 5 + j) = da_i_dT_local(i) * a_i_local(j) + a_i_local(i) * da_i_dT_local(j)
            cross_2_da((i - 1) * 5 + j) = 2# * da_i_dT_local(i) * da_i_dT_local(j)
            d2a_term((i - 1) * 5 + j) = _
                 d2a_i_dT2_local(i) * a_i_local(j) _
               + cross_2_da((i - 1) * 5 + j) _
               + a_i_local(i) * d2a_i_dT2_local(j)
        Next j
    Next i

    ' Sum up d²a_mix/dT² as S_i S_j zf(i)*zf(j)*(term1+term2+term3)
    tempSum = 0#
    For i = 1 To 5
        For j = 1 To 5
            sqrt_ai_aj = Sqr(a_i_local(i) * a_i_local(j))
            If sqrt_ai_aj = 0# Then
                term1 = 0#: term2 = 0#: term3 = 0#
            Else
                term1 = -(N_mat((i - 1) * 5 + j) / sqrt_ai_aj) * dkij_dT(i, j)
                term2 = (1# - kij(i, j)) * 0.5 * ( _
                             (d2a_term((i - 1) * 5 + j) / sqrt_ai_aj) _
                           - ((N_mat((i - 1) * 5 + j) ^ 2) / (2# * (sqrt_ai_aj ^ 3))) _
                         )
                term3 = -sqrt_ai_aj * d2kij_dT2(i, j)
            End If

            tempSum = tempSum + zf(i) * zf(j) * (term1 + term2 + term3)
        Next j
    Next i

    d2a_mix_dT2 = tempSum

    ' Populate d2a_i_dT2(1 To 5) from the local second-derivative array
    ReDim d2a_i_dT2(1 To 5)

    For i = 1 To 5
        d2a_i_dT2(i) = d2a_i_dT2_local(i)
    Next i
End Sub




'===================================================================
'   6) Compute_dz_dT_constP
'===================================================================
Private Function Compute_dz_dT_constP( _
        ByVal T As Double, _
        ByVal P As Double, _
        ByVal a_mix_val As Double, _
        ByVal da_mix_dT_val As Double, _
        ByVal b_mix_val As Double, _
        ByVal z_val As Double) As Double

    Dim A_ As Double, B_ As Double
    Dim dB_dT_ As Double, dA_dT_ As Double
    Dim dF_dz As Double, dF_dA As Double, dF_dB As Double
    Dim dF_dT_ As Double

    A_ = a_mix_val * P / (R_field * T) ^ 2
    B_ = b_mix_val * P / (R_field * T)

    dB_dT_ = -B_ / T
    dA_dT_ = (P / (R_field * T) ^ 2) * da_mix_dT_val - 2# * a_mix_val * P / (R_field ^ 2 * T ^ 3)

    dF_dz = 3# * z_val ^ 2 - 2# * (1# - B_) * z_val + (A_ - 3# * B_ ^ 2 - 2# * B_)
    dF_dA = z_val - B_
    dF_dB = z_val ^ 2 - (6# * B_ + 2#) * z_val - A_ + 2# * B_ + 3# * B_ ^ 2

    dF_dT_ = dF_dA * dA_dT_ + dF_dB * dB_dT_
    Compute_dz_dT_constP = -dF_dT_ / dF_dz
End Function

'===================================================================
'   7) Main UDF: BNS_Full
'      Inputs: Temperature, Pressure, sg, co2, h2s, n2, h2, AG, Vis, Den, Thermo
'              degF or degC, psia or MPa
'      Returns: Scripting.Dictionary with keys
'    {'Z':        Z-Factor                                     (Dimensionless)
'    'H':         Enthalpy relative to 60 degF and 14.606 psia (Btu/(lb-mol.R) or kJ/(kmol K))
'    'Cp':        Isobaric heat capacity                       (Btu/(lb-mol.R) or kJ/(kmol K))
'    'Cv':        Isochoric heat capacity                      (Btu/(lb-mol·R) or kJ/(kmol K))
'    'mu_JT':        Joule-Thomson COefficient                    (Btu/(lb-mol.R) / degC/MPa)
'    'viscosity': Gas viscosity                                (cP or mPa.s)}
'    'density':   Gas Density                                  (lbm/cuft or kg/m3)
'
'===================================================================


Public Function BNS_Full( _
    ByVal temp As Double, _
    ByVal pres As Double, _
    ByVal sg As Double, _
    Optional ByVal co2 As Double = 0#, _
    Optional ByVal h2s As Double = 0#, _
    Optional ByVal n2 As Double = 0#, _
    Optional ByVal h2 As Double = 0#, _
    Optional ByVal AG As Boolean = False, _
    Optional ByVal Vis As Boolean = True, _
    Optional ByVal Den As Boolean = False, _
    Optional ByVal Thermo As Boolean = True, _
    Optional ByVal Metric As Boolean = False _
) As Variant

    Dim degF As Double
    Dim psia As Double
    If Metric = True Then
        degF = temp * 9 / 5 + 32
        psia = pres * 145.0377
    Else
        degF = temp
        psia = pres
    End If
    

    InitConstants

    Dim dict As Dictionary
    Set dict = New Dictionary

    ' 1) Convert degF -> °R
    Dim degR As Double
    degR = degF + DEG_F_TO_R

    ' 2) Build composition array zf(1 To 5)
    Dim zf(1 To 5) As Double
    Dim sumNonHC As Double
    sumNonHC = co2 + h2s + n2 + h2
    If sumNonHC > 1# Then sumNonHC = 1#
    zf(1) = co2: zf(2) = h2s: zf(3) = n2: zf(4) = h2
    zf(5) = 1# - sumNonHC

    ' 3) Hydrocarbon pseudo-critical if hc_fraction > 0
    Dim hc_fraction As Double, sg_hc As Double
    hc_fraction = zf(5)
    If hc_fraction > 0# Then
        Dim sumNonhwm As Double
        sumNonhwm = co2 * mws(1) + h2s * mws(2) + n2 * mws(3) + h2 * mws(4)
        sg_hc = (sg - sumNonhwm / MW_AIR) / hc_fraction
    Else
        sg_hc = 0.75
    End If
    sg_hc = WorksheetFunction.Max(sg_hc, 16.0425 / MW_AIR)
    Dim hc_mw As Double
    hc_mw = sg_hc * MW_AIR
    
    mws(5) = hc_mw

    ' 4) Override "Gas" critical properties locally
    Dim tcs_local(1 To 5) As Double, pcs_local(1 To 5) As Double
    Dim tpc_hc As Double, ppc_hc As Double
    tc_pc sg_hc, AG, tpc_hc, ppc_hc
    Dim i As Long, j As Long
    For i = 1 To 5
        tcs_local(i) = tcs(i)
        pcs_local(i) = pcs(i)
    Next i
    tcs_local(5) = tpc_hc
    pcs_local(5) = ppc_hc


 
    ' 5) Compute reduced temperatures & pressures
    Dim trs(1 To 5) As Double, prs(1 To 5) As Double
    For i = 1 To 5
        trs(i) = degR / tcs_local(i)
        prs(i) = psia / pcs_local(i)
    Next i
 
    ' 6) Compute BIPs + derivatives
    Dim kijArr() As Double, dkij_dTArr() As Double, d2kij_dT2Arr() As Double
    If Thermo Then
        calc_bips degR, tpc_hc, True, kijArr, dkij_dTArr, d2kij_dT2Arr
    Else
        calc_bips degR, tpc_hc, False, kijArr, dkij_dTArr, d2kij_dT2Arr
    End If

    ' 7) a_i, a_c_i, a_i, b_i, a_mix, b_mix
    Dim m_i(1 To 5) As Double
    Dim alpha_i(1 To 5) As Double
    Dim a_c_i(1 To 5) As Double
    Dim a_i(1 To 5) As Double
    Dim b_i(1 To 5) As Double
    Dim aij(1 To 5, 1 To 5) As Double
    Dim a_mix As Double, b_mix As Double

    For i = 1 To 5
        m_i(i) = 0.37464 + 1.54226 * ACF(i) - 0.26992 * (ACF(i) ^ 2)
        alpha_i(i) = (1# + m_i(i) * (1# - Sqr(trs(i)))) ^ 2
        a_c_i(i) = OmegaA(i) * R_field ^ 2 * tcs_local(i) ^ 2 / pcs_local(i)
        a_i(i) = a_c_i(i) * alpha_i(i)
        b_i(i) = OmegaB(i) * R_field * tcs_local(i) / pcs_local(i)
    Next i

    a_mix = 0#: b_mix = 0#
    For i = 1 To 5
        For j = 1 To 5
            aij(i, j) = Sqr(a_i(i) * a_i(j)) * (1# - kijArr(i, j))
            a_mix = a_mix + zf(i) * zf(j) * aij(i, j)
        Next j
        b_mix = b_mix + zf(i) * b_i(i)
    Next i

    ' 8) Solve cubic for z_eos
    Dim RT As Double, A_param As Double, B_param As Double
    Dim coeffs(0 To 3) As Double
    Dim z_eos As Double
    RT = R_field * degR
    A_param = a_mix * psia / (RT * RT)
    B_param = b_mix * psia / (RT)

    coeffs(0) = 1#
    coeffs(1) = -(1# - B_param)
    coeffs(2) = A_param - 3# * B_param ^ 2 - 2# * B_param
    coeffs(3) = -(A_param * B_param - B_param ^ 2 - B_param ^ 3)
    
    z_eos = CubicRootArrayDirect(coeffs, 1)

    ' 9) Volume shift
    Dim Bi(1 To 5) As Double
    For i = 1 To 5
        Bi(i) = OmegaB(i) * (prs(i) / trs(i))
    Next i
    Dim sum_zVBi As Double
    sum_zVBi = 0#
    For i = 1 To 5
        sum_zVBi = sum_zVBi + zf(i) * VSHIFT(i) * Bi(i)
    Next i
    Dim z_vshift As Double
    z_vshift = z_eos - sum_zVBi

    ' Peneloux translation as a constant molar volume offset (ft3/lb-mol):
    '   V_shifted = V_eos - c_mix,  c_mix = SUM(zf * VSHIFT * b_i)
    ' The p/T dependence in Bi cancels against RT/p exactly, so c_mix is constant.
    ' A constant translation leaves Cp, Cv and entropy unchanged, but shifts
    ' enthalpy by -c_mix*p and the JT coefficient by +c_mix/Cp.
    Dim c_mix As Double
    c_mix = 0#
    For i = 1 To 5
        c_mix = c_mix + zf(i) * VSHIFT(i) * b_i(i)
    Next i
    
    'Debug.Print "Adding Z", z_vshift
    dict.Add "Z", z_vshift
    'Debug.Print "Done"
    
    '--- Calculate density if requested ---
    If Den Then
        Dim mwt As Double
        mwt = 0#
        For i = 1 To 5
            mwt = mwt + zf(i) * mws(i)
        Next i
        Dim density As Double
        ' lbm/ft³
        density = mwt * psia / (z_vshift * R_field * degR)

        dict.Add "Dens", density
    End If



    '===================================================================
    '   If Thermo, compute H, Cp, Cv, mu_JT
    '===================================================================
    If Thermo Then
        ' (a) Recompute a_i, da_i_dT, da_mix_dT
        Dim da_i_dT(1 To 5) As Double
        Dim da_mix_dT_val As Double
        
        Call Compute_da_mix_dT(zf, a_c_i, m_i, tcs_local, trs, kijArr, dkij_dTArr, a_i, da_i_dT, da_mix_dT_val)


        ' (b) Compute d2a_mix_dT2
        Dim d2a_mix_dT2_val As Double
        Dim d2a_i_dT2(1 To 5) As Double
        
        Dim d2a_i(1 To 5) As Double
        Dim v_d2a   As Variant

        
        'Call Compute_d2a_mix_dT2(zf, a_c_i, m_i, tcs_local, trs, kijArr, dkij_dTArr, d2kij_dT2Arr, d2a_mix_dT2_val, d2a_i_dT2)
        v_d2a = d2a_i      ' assign the fixed array into a Variant
        Call Compute_d2a_mix_dT2( _
        zf, a_c_i, m_i, tcs_local, trs, _
        kijArr, dkij_dTArr, d2kij_dT2Arr, _
        d2a_mix_dT2_val, v_d2a _
        )

        ' (c) Compute dz/dT_constP
        Dim dz_dT As Double
        
        dz_dT = Compute_dz_dT_constP(degR, psia, a_mix, da_mix_dT_val, b_mix, z_eos)
    

        ' (d) Bdim
        Dim sqrt2 As Double: sqrt2 = Sqr(2#)
        Dim Bdim As Double: Bdim = (b_mix * psia) / (R_field * degR)

        ' (e) Adjust Cp coefficients for heavy hydrocarbon
    Dim CpA(1 To 5) As Double, CpB(1 To 5) As Double
    Dim CpC(1 To 5) As Double, CpD(1 To 5) As Double, CpE(1 To 5) As Double
    
    For i = 1 To 5
        CpA(i) = Cp_a_(i)
        CpB(i) = Cp_b_(i)
        CpC(i) = Cp_c_(i)
        CpD(i) = Cp_d_(i)
        CpE(i) = Cp_e_(i)
    Next i
    
    ' --- adjust *only* the 5th (pseudo-comp) entry locally ---
    Dim a0(1 To 5) As Double, a1(1 To 5) As Double, X As Double
    a0(1) = 0.0007857:   a1(1) = -0.0081649
    a0(2) = 0.0013123:   a1(2) = 0.0055485
    a0(3) = 0.00098133:  a1(3) = 0.083258
    a0(4) = 0.0016463:   a1(4) = 0.20635
    a0(5) = 0.017306:    a1(5) = 2.5551
    
    X = hc_mw - 16.0425
    CpA(5) = CpA(5) * (a0(1) * X ^ 2 + a1(1) * X + 1#)
    CpB(5) = CpB(5) * (a0(2) * X ^ 2 + a1(2) * X + 1#)
    CpC(5) = CpC(5) * (a0(3) * X ^ 2 + a1(3) * X + 1#)
    CpD(5) = CpD(5) * (a0(4) * X ^ 2 + a1(4) * X + 1#)
    CpE(5) = CpE(5) * (a0(5) * X ^ 2 + a1(5) * X + 1#)
    
    ' --- now compute the ideal-gas heat capacity ---
    Dim T_K As Double
    T_K = degR * 5# / 9#  ' °R -> K
    
    Dim Cp_poly As Double, Cp_comp As Double
    Cp_poly = 0#
    For i = 1 To 5
        Cp_comp = (CpA(i) + CpB(i) * T_K + CpC(i) * T_K ^ 2 + CpD(i) * T_K ^ 3 + CpE(i) * T_K ^ 4)
        Cp_poly = Cp_poly + zf(i) * Cp_comp
    Next i
    
    Dim Cp_IG_Btu As Double
    Cp_IG_Btu = Cp_poly * R_THERMO  ' Btu/(lb-mol·°R)
    
        ' (g) Departure enthalpy
        Dim H_dep_ft3psia As Double
        H_dep_ft3psia = R_field * degR * (z_eos - 1#) _
                        + (degR * da_mix_dT_val - a_mix) / (2# * sqrt2 * b_mix) * _
                          Log((z_eos + (sqrt2 + 1#) * Bdim) / (z_eos - (sqrt2 - 1#) * Bdim))
        Dim H_dep_Btu As Double
        H_dep_Btu = H_dep_ft3psia / FT3_PSIA_TO_BTU


        ' (h) Ideal-gas enthalpy from reference T_ref = 60°F
        Dim T_ref_F As Double: T_ref_F = 60#
        Dim T_ref_K As Double: T_ref_K = (T_ref_F + DEG_F_TO_R) * 5# / 9#
        Dim H_IG_Btu As Double: H_IG_Btu = 0#

        For i = 1 To 5
            H_IG_Btu = H_IG_Btu + zf(i) * ( _
                CpA(i) * (T_K - T_ref_K) _
            + CpB(i) / 2# * (T_K ^ 2 - T_ref_K ^ 2) _
            + CpC(i) / 3# * (T_K ^ 3 - T_ref_K ^ 3) _
            + CpD(i) / 4# * (T_K ^ 4 - T_ref_K ^ 4) _
            + CpE(i) / 5# * (T_K ^ 5 - T_ref_K ^ 5) _
            )
        Next i

        ' now convert from Btu/(lb-mol·K) back to Btu/(lb-mol·°R)
        H_IG_Btu = H_IG_Btu * R_THERMO * 9# / 5#

        ' (i) Total enthalpy
        ' Absolute enthalpy at 60 degF and 14.696psia
        Dim H0_hc As Double, H0 As Double
        H0_hc = -0.015774 * (hc_mw - 16.0425) ^ 2 - 0.646645 * (hc_mw - 16.0425) - 8.2551915 ' Poyfit of enthalpy of pseudocomponent at SC
        H0 = -16.6022 * zf(1) + -21.5512 * zf(2) + -3.57757 * zf(3) + 0.008054 * zf(4) + H0_hc * zf(5) ' Sumproduct of Pure component Enthalpy at 60 degF and 14.696 psia
        'H0 = 0
        
        ' Volume-shift contribution to enthalpy: H_shifted = H_eos - c_mix*p.
        ' H is reported relative to 60 degF and 14.696 psia, so the reference
        ' pressure term cancels and the H0 constants above stay valid unchanged.
        Const P_REF_PSIA As Double = 14.696
        Dim H_vshift_Btu As Double
        H_vshift_Btu = -c_mix * (psia - P_REF_PSIA) / FT3_PSIA_TO_BTU

        Dim H_total_Btu As Double
        H_total_Btu = H_IG_Btu + H_dep_Btu + H_vshift_Btu - H0


        ' (j) Cp_total = Cp_IG + dH_dep/dT
        Dim dH_dep_dT As Double
        dH_dep_dT = HDepDeriv(degR, z_eos, dz_dT, a_mix, da_mix_dT_val, b_mix, sqrt2, Bdim, d2a_mix_dT2_val)
        Dim Cp_total_Btu As Double
        Cp_total_Btu = Cp_IG_Btu + dH_dep_dT


        ' (k) Molar volume V = z_eos·R·T / p, untranslated, plus the translated
        '     volume. Cp and Cv are evaluated on the untranslated EOS deliberately:
        '     a constant translation leaves U(T,V) and S(T,p) unchanged, so both
        '     heat capacities are invariant under it. Only V itself is translated.
        Dim V As Double
        V = z_eos * R_field * degR / psia
        Dim V_shifted As Double
        V_shifted = V - c_mix

        ' (l) dP/dT at constant V
        Dim dP_dT_constV As Double
        dP_dT_constV = R_field / (V - b_mix) - da_mix_dT_val / (V ^ 2 + 2# * b_mix * V - b_mix ^ 2)


        ' (m) dV/dT = (R/psia)*(z_eos + T·dz/dT)
        Dim dV_dT As Double
        dV_dT = (R_field / psia) * (z_eos + degR * dz_dT)

        ' analytic dP/dV at constant T from PR
        Dim Dens As Double
        Dens = (V ^ 2 + 2 * b_mix * V - b_mix ^ 2) ^ 2
        Dim dP_dV_constT As Double
        dP_dV_constT = (-R_field * degR / (V - b_mix) ^ 2 + 2 * a_mix * (V + b_mix) / Dens)

        Dim dV_dP_constT As Double
        dV_dP_constT = -1 / dP_dV_constT
        

        ' (n) Cv_total = Cp_total – T·(dV/dT)·(dP/dT)|V
        ' Everything in ft³·psia/(lb-mol·°R) for the volume derivatives, so you must divide by FT3_PSIA_TO_BTU at the end to get Btu/(lb-mol·°R)
        Dim Cv_total_Btu As Double
        Cv_total_Btu = Cp_total_Btu + (degR * dV_dT ^ 2 * dP_dV_constT) / FT3_PSIA_TO_BTU


        ' (o) mu_JT = (T·dV/dT – V_shifted)/(Cp_total * FT3_PSIA_TO_BTU)
        '     Uses the translated volume, consistent with the density and the
        '     enthalpy the model reports.
        Dim mu_JT As Double
        mu_JT = (degR * dV_dT - V_shifted) / (Cp_total_Btu * FT3_PSIA_TO_BTU)

        'Debug.Print "Adding H", H_total_Btu
        dict.Add "H", H_total_Btu
        'Debug.Print "Adding Cp", Cp_total_Btu
        dict.Add "Cp", Cp_total_Btu
        'Debug.Print "Adding Cv", Cv_total_Btu
        dict.Add "Cv", Cv_total_Btu
        'Debug.Print "Adding mu_JT", mu_JT
        dict.Add "mu_JT", mu_JT
    End If
    

    
    '===================================================================
    '   8) If Vis, call LBC_Scalar
    '===================================================================
    If Vis Then
        Dim mu_cP As Double
        mu_cP = LBC_Scalar(z_vshift, degF, psia, sg, co2, h2s, n2, h2, AG)
        'Debug.Print "Adding viscosity", mu_cP
        dict.Add "viscosity", mu_cP
    End If


    '------------------------------------------------------------------------
    ' 2) Build a 1×N Variant array whose size depends on Vis/Thermo
    '------------------------------------------------------------------------

    Dim colCount As Long
    colCount = 1                ' always have at least 1 column (Z)

    If Vis Then colCount = colCount + 1
    If Den Then colCount = colCount + 1
    If Thermo Then colCount = colCount + 4     ' H, Cp, Cv, mu_JT adds 4 more columns

    Dim outArr() As Variant
    ReDim outArr(1 To 1, 1 To colCount)

    ' Fill in outArr(1,1) = Z-factor
    'Debug.Print "Adding dict Z"
    outArr(1, 1) = dict.Item("Z")
    'Debug.Print "Done"

    Dim idx As Long
    idx = 1

    ' If viscosity was requested, put it in the next column:
    If Vis Then
        idx = idx + 1
        'Debug.Print "Adding dict viscosity"
        outArr(1, idx) = dict.Item("viscosity")
        'Debug.Print "Done"
    End If
    
    ' If density was requested, put it in the next column:
    If Den Then
        idx = idx + 1
        outArr(1, idx) = dict.Item("Dens")
    End If

    ' If thermo was requested, put H, Cp, Cv, mu_JT in the following columns (in that order):
    If Thermo Then
        idx = idx + 1
        outArr(1, idx) = dict.Item("H")

        idx = idx + 1
        outArr(1, idx) = dict.Item("Cp")

        idx = idx + 1
        outArr(1, idx) = dict.Item("Cv")

        idx = idx + 1
        outArr(1, idx) = dict.Item("mu_JT")
    End If
    
    '******* Units conversion
    If Metric Then
        ' Convert outputs to metric units
        idx = 1
        If Vis Then idx = idx + 1 ' No need to convert, but do need to increment index
        
        If Den Then
            idx = idx + 1
            outArr(1, idx) = outArr(1, idx) * 16.01846337396 ' kg/m³
        End If
    
        If Thermo Then
            idx = idx + 1
            outArr(1, idx) = dict.Item("H") * 2.326 ' kJ/(kmol)
            idx = idx + 1
            outArr(1, idx) = dict.Item("Cp") * 4.186800585 ' kJ/(kmol·K)
            idx = idx + 1
            outArr(1, idx) = dict.Item("Cv") * 4.186800585 ' kJ/(kmol·K)
            idx = idx + 1
            outArr(1, idx) = dict.Item("mu_JT") * 80.5765 ' degC/MPa
        End If
    End If
    ' ******* End units conversion

    BNS_Full = outArr
    Exit Function

    
End Function

'===================================================================
'   8a) HDepDeriv: computes dH_dep/dT exactly as in Python's Hdep_deriv
'===================================================================
Private Function HDepDeriv( _
        ByVal T As Double, _
        ByVal z As Double, _
        ByVal dzdT As Double, _
        ByVal a_mix As Double, _
        ByVal da_mix_dT As Double, _
        ByVal b_mix As Double, _
        ByVal sqrt2 As Double, _
        ByVal Bdim As Double, _
        ByVal d2a_mix_dT2 As Double _
    ) As Double

    Dim dF1dT As Double, N As Double, D As Double
    Dim dBdim_dT As Double, dN_dT As Double, dD_dT As Double
    Dim dln_dT As Double
    Dim X As Double, dX_dT As Double
    Dim term_dF2 As Double

    dF1dT = R_field * (z - 1#) + R_field * T * dzdT
    N = z + (sqrt2 + 1#) * Bdim
    D = z - (sqrt2 - 1#) * Bdim
    dBdim_dT = -Bdim / T
    dN_dT = dzdT + (sqrt2 + 1#) * dBdim_dT
    dD_dT = dzdT - (sqrt2 - 1#) * dBdim_dT
    dln_dT = (1# / N) * dN_dT - (1# / D) * dD_dT

    X = (T * da_mix_dT - a_mix) / (2# * sqrt2 * b_mix)
    dX_dT = (da_mix_dT + T * d2a_mix_dT2 - da_mix_dT) / (2# * sqrt2 * b_mix)
    term_dF2 = dX_dT * Log(N / D) + X * dln_dT

    HDepDeriv = (dF1dT + term_dF2) / FT3_PSIA_TO_BTU
End Function

'===================================================================
'   9) LBC_Scalar
'===================================================================
Private Function LBC_Scalar( _
        ByVal Z_scalar As Double, _
        ByVal degF_scalar As Double, _
        ByVal psia_scalar As Double, _
        ByVal sg As Double, _
        ByVal co2 As Double, _
        ByVal h2s As Double, _
        ByVal n2 As Double, _
        ByVal h2 As Double, _
        ByVal AG As Boolean _
    ) As Double

    Dim zi(1 To 5) As Double
    Dim sumNonHC As Double
    sumNonHC = co2 + h2s + n2 + h2

    ' (a) Validate composition
    If sumNonHC > 1# Or co2 < 0# Or h2s < 0# Or n2 < 0# Or h2 < 0# Then
        LBC_Scalar = Empty
        Exit Function
    End If

    ' (b) Build composition array
    zi(1) = co2: zi(2) = h2s: zi(3) = n2: zi(4) = h2: zi(5) = 1# - sumNonHC

    ' (c) Weighted SG for HC, with lower bound 0.553779772
    Dim sg_hc As Double
    If sumNonHC > 1# Then
        sg_hc = 0.75
    ElseIf zi(5) > 0# Then
        Dim sum_nonhwm As Double
        sum_nonhwm = co2 * mws(1) + h2s * mws(2) + n2 * mws(3) + h2 * mws(4)
        sg_hc = (sg - sum_nonhwm / MW_AIR) / zi(5)
    Else
        sg_hc = 0.75
    End If
    sg_hc = WorksheetFunction.Max(sg_hc, 0.553779772)
    Dim hc_gas_mw As Double
    hc_gas_mw = sg_hc * MW_AIR

    ' (d) Local copies of tcs and pcs, override "Gas"
    Dim tcs_l(1 To 5) As Double, pcs_l(1 To 5) As Double
    Dim tpc_hc As Double, ppc_hc As Double
    Dim i As Long
    For i = 1 To 5
        tcs_l(i) = tcs(i)
        pcs_l(i) = pcs(i)
    Next i
    tc_pc sg_hc, AG, tpc_hc, ppc_hc
    tcs_l(5) = tpc_hc
    pcs_l(5) = ppc_hc

    ' (e) Local VCVIS and update "Gas"
    Dim VCVIS_l(1 To 5) As Double
    For i = 1 To 5: VCVIS_l(i) = VCVIS(i): Next i
    VCVIS_l(5) = 0.057671 * (hc_gas_mw - 16.0425) + 1.44383

    ' (f) Local mws and update "Gas"
    Dim mws_l(1 To 5) As Double
    For i = 1 To 5: mws_l(i) = mws(i): Next i
    mws_l(5) = hc_gas_mw

    ' (g) Temperature in °R
    Dim degR As Double
    degR = degF_scalar + DEG_F_TO_R

    ' (h) Stiel-Thodos for each component
    Dim mu_components() As Double
    mu_components = StielThodosViscosity(degR, mws_l, tcs_l, pcs_l)

    ' (i) Herning-Zippener mixing for dilute gas
    Dim sqrt_mw(1 To 5) As Double
    Dim numer As Double, denom As Double
    numer = 0#: denom = 0#
    For i = 1 To 5
        sqrt_mw(i) = Sqr(mws_l(i))
        numer = numer + zi(i) * mu_components(i) * sqrt_mw(i)
        denom = denom + zi(i) * sqrt_mw(i)
    Next i
    Dim mu_dilute As Double
    If denom = 0 Then
        mu_dilute = 0
    Else
        mu_dilute = numer / denom
    End If

    ' (j) LBC correlation part
    Dim rhoc As Double
    Dim sumVC As Double
    sumVC = 0#
    For i = 1 To 5
        sumVC = sumVC + VCVIS_l(i) * zi(i)
    Next i
    If sumVC = 0 Then
        LBC_Scalar = mu_dilute
        Exit Function
    End If
    rhoc = 1# / sumVC

    Dim gas_density As Double
    gas_density = psia_scalar / (Z_scalar * R_field * degR)  ' lb-mol/ft³
    Dim rhor As Double
    rhor = gas_density / rhoc

    ' (k) Lorenz–Bray–Clark polynomial coefficients
    Dim a(1 To 5) As Double
    a(1) = 0.1023:       a(2) = 0.023364
    a(3) = 0.058533:     a(4) = -0.0392852
    a(5) = 0.00926279

    Dim lhs As Double
    lhs = a(1) + a(2) * rhor + a(3) * (rhor ^ 2) + a(4) * (rhor ^ 3) + a(5) * (rhor ^ 4)

    ' (l) Mixture pseudo-criticals for "eta" group
    Dim Tc_mix As Double, Pc_mix As Double, Mw_mix As Double
    Tc_mix = 0#: Pc_mix = 0#: Mw_mix = 0#
    For i = 1 To 5
        Tc_mix = Tc_mix + zi(i) * tcs_l(i)
        Pc_mix = Pc_mix + zi(i) * pcs_l(i)
        Mw_mix = Mw_mix + zi(i) * mws_l(i)
    Next i
    Tc_mix = Tc_mix * (5# / 9#)   ' °R
    Pc_mix = Pc_mix / 14.696     ' psia

    Dim eta_mix As Double
    If Mw_mix <= 0# Or Pc_mix <= 0# Then
        eta_mix = 1#
    Else
        eta_mix = (Tc_mix ^ (1# / 6#)) / (Sqr(Mw_mix) * (Pc_mix ^ (2# / 3#)))
    End If

    ' (m) Final viscosity (cP)
    Dim viscosity As Double
    viscosity = ((lhs ^ 4) - 0.0001) / eta_mix + mu_dilute

    LBC_Scalar = viscosity
End Function

'===================================================================
'   9a) StielThodosViscosity: exact port of Python's stiel_thodos_viscosity
'      Returns a 1D array(1 To 5) of µ_component (cP)
'===================================================================
Private Function StielThodosViscosity( _
        ByVal degR As Double, _
        ByRef mwsArr() As Double, _
        ByRef tcsArr() As Double, _
        ByRef pcsArr() As Double _
    ) As Variant

    Dim i As Long
    Dim Tr(1 To 5) As Double
    Dim muC(1 To 5) As Double
    Dim Tc_k(1 To 5) As Double
    Dim Pc_atm(1 To 5) As Double
    Dim eta_factor(1 To 5) As Double
    Dim Tred As Double

    For i = 1 To 5
        Tr(i) = degR / tcsArr(i)
        Tc_k(i) = tcsArr(i) * (5# / 9#)      ' °R
        Pc_atm(i) = pcsArr(i) / 14.696       ' psia
        If mwsArr(i) <= 0 Or Pc_atm(i) <= 0 Then
            eta_factor(i) = 1#
        Else
            eta_factor(i) = Tc_k(i) ^ (1# / 6#) / (Sqr(mwsArr(i)) * (Pc_atm(i) ^ (2# / 3#)))
        End If
    Next i

    For i = 1 To 5
        Tred = Tr(i)
        If Tred <= 1.5 Then
            muC(i) = 0.00034 * (Tred ^ 0.94) / eta_factor(i)
        Else
            muC(i) = 0.0001778 * ((4.58 * Tred) - 1.67) ^ (5# / 8#) / eta_factor(i)
        End If
    Next i

    StielThodosViscosity = muC
End Function

Public Function BNS_Get( _
        ByVal whichProp As String, _
        ByVal temp As Double, _
        ByVal pres As Double, _
        ByVal sg As Double, _
        Optional ByVal co2 As Double = 0#, _
        Optional ByVal h2s As Double = 0#, _
        Optional ByVal n2 As Double = 0#, _
        Optional ByVal h2 As Double = 0#, _
        Optional ByVal AG As Boolean = False, _
        Optional ByVal Metric As Boolean = False _
    ) As Variant

    Dim arr As Variant
    Dim idx As Long
    Dim propName As String

    propName = Trim(UCase$(whichProp))

    Select Case propName
        Case "Z"
            ' only Z
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, False, False, Metric)
            idx = 1

        Case "Vis", "VIS"
            ' ask exactly for viscosity (and Z)
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    True, False, False, Metric)
            idx = 2
            
        Case "Den", "DEN"
            ' ask exactly for density (and Z)
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, True, False, Metric)
            idx = 2

        Case "H"
            ' ask for thermo only: {Z, H, Cp, Cv, mu_JT}
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, False, True, Metric)
            idx = 2

        Case "Cp", "CP"
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, False, True, Metric)
            idx = 3

        Case "Cv", "CV"
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, False, True, Metric)
            idx = 4

        Case "JT", "MU-JT"
            arr = BNS_Full(temp, pres, sg, co2, h2s, n2, h2, AG, _
                                    False, False, True, Metric)
            idx = 5

        Case Else
            ' Unknown property -> return #N/A
            BNS_Get = CVErr(xlErrNA)
            Exit Function
    End Select

    ' Extract the single value from the 1×N array:
    ' arr is dimensioned (1 To 1, 1 To colCount)
    ' so arr(1, idx) is the desired scalar.
    On Error GoTo OutputError
    BNS_Get = arr(1, idx)
    Exit Function

OutputError:
    BNS_Get = CVErr(xlErrNA)
End Function



