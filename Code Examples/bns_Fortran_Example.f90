! bns_pr_fixed2.f90 — Fortran Peng–Robinson + Thermo + Viscosity
! NOTE: compile as FREE-FORM (.f90) or add compiler flag -ffree-form
 
module bns
  use iso_fortran_env,  only: real64
  use ieee_arithmetic,  only: ieee_is_finite
  implicit none
  private
  public :: pr_properties

  integer, parameter :: ncomp = 5
  real(real64), parameter :: SQRT2 = 1.4142135623730951_real64
  real(real64), parameter :: PI    = 3.1415926535897932_real64

  real(real64), parameter :: Rgas       = 10.731577089016_real64
  real(real64), parameter :: Rthermo    = 1.98588_real64
  real(real64), parameter :: MW_AIR     = 28.97_real64
  real(real64), parameter :: DEG_F_TO_R = 459.67_real64
  real(real64), parameter :: FT3PSIA_TO_BTU = 5.403_real64
  real(real64), parameter :: mwCH4      = 16.0425_real64
  real(real64), parameter :: TcCH4      = 343.008_real64
  real(real64), parameter :: PcCH4      = 667.029_real64
  real(real64), parameter :: VcZcCH4    = (Rgas*TcCH4)/PcCH4

contains

function clamp_min(x, lo) result(y)
  real(real64), intent(in) :: x, lo
  real(real64) :: y
  if (x >= lo) then; y = x; else; y = lo; end if
end function

function dot5(a, b) result(s)
  real(real64), intent(in) :: a(ncomp), b(ncomp)
  real(real64) :: s
  s = sum(a*b)
end function

function poly5_eval(a, x) result(p)
  real(real64), intent(in) :: a(5), x
  real(real64) :: p, xn
  integer :: k
  p = 0.0_real64; xn = 1.0_real64
  do k = 1, 5
    p = p + a(k)*xn
    xn = xn * x
  end do
end function

function cube_root(x) result(y)
  real(real64), intent(in) :: x
  real(real64) :: y
  if (x == 0.0_real64) then
    y = 0.0_real64
  else
    y = sign(1.0_real64, x) * abs(x)**(1.0_real64/3.0_real64)
  end if
end function

subroutine base_props(mws, tcs, pcs, acf, vshift, omega_a, omega_b, vcvis, cp_poly)
  real(real64), intent(out) :: mws(ncomp), tcs(ncomp), pcs(ncomp), acf(ncomp), vshift(ncomp)
  real(real64), intent(out) :: omega_a(ncomp), omega_b(ncomp), vcvis(ncomp)
  real(real64), intent(out) :: cp_poly(ncomp,5)

  mws  = [44.01_real64, 34.082_real64, 28.014_real64,  2.016_real64,  0.0_real64]
  tcs  = [547.416_real64, 672.120_real64, 227.160_real64, 47.430_real64, 1.0_real64]
  pcs  = [1069.51_real64, 1299.97_real64, 492.84_real64,  187.53_real64, 1.0_real64]
  acf  = [0.12253_real64, 0.04909_real64, 0.037_real64,  -0.217_real64, -0.03899_real64]
  vshift = [-0.27607_real64, -0.22901_real64, -0.21066_real64, -0.36270_real64, -0.19076_real64]
  omega_a = [0.427671_real64, 0.436725_real64, 0.457236_real64, 0.457236_real64, 0.457236_real64]
  omega_b = [0.0696397_real64,0.0724345_real64,0.0777961_real64,0.0777961_real64,0.0777961_real64]
  vcvis = [1.46352_real64, 1.46808_real64, 1.35526_real64, 0.68473_real64, 1.44383_real64]

  cp_poly(1,:) = [2.725473196_real64,  0.004103751_real64,  1.5602e-05_real64, &
                   -4.19321e-08_real64,  3.10542e-11_real64]
  cp_poly(2,:) = [4.446031265_real64, -0.005296052_real64,  2.0533e-05_real64, &
                   -2.58993e-08_real64,  1.25555e-11_real64]
  cp_poly(3,:) = [3.423811591_real64,  0.001007461_real64, -4.58491e-06_real64, &
                    8.4252e-09_real64, -4.38083e-12_real64]
  cp_poly(4,:) = [1.421468418_real64,  0.018192108_real64, -6.04285e-05_real64, &
                    9.08033e-08_real64, -5.18972e-11_real64]
  cp_poly(5,:) = [5.369051342_real64, -0.014851371_real64,  4.86358e-05_real64, &
                   -3.70187e-08_real64,  1.80641e-12_real64]
end subroutine

function tc_ag(x) result(tc)
  real(real64), intent(in) :: x
  real(real64) :: tc
  tc = 2695.14765_real64 * x / (274.341701_real64 + x) + TcCH4
end function
function tc_gc(x) result(tc)
  real(real64), intent(in) :: x
  real(real64) :: tc
  tc = 1098.10948_real64 * x / (101.529237_real64 + x) + TcCH4
end function
function pc_from_vc_zc(x, vc_slope, tc) result(pc)
  real(real64), intent(in) :: x, vc_slope, tc
  real(real64) :: pc, vc_over_zc
  vc_over_zc = vc_slope * x + VcZcCH4
  pc = Rgas * tc / vc_over_zc
end function

subroutine pseudo_critical(sg_hc, ag, tpc, ppc)
  real(real64), intent(in)  :: sg_hc
  logical,      intent(in)  :: ag
  real(real64), intent(out) :: tpc, ppc
  real(real64) :: x, slope
  x = max(0.0_real64, MW_AIR*sg_hc - mwCH4)
  if (ag) then
    tpc = tc_ag(x); slope = 0.177497835_real64
  else
    tpc = tc_gc(x); slope = 0.170931432_real64
  end if
  ppc = pc_from_vc_zc(x, slope, tpc)
end subroutine

subroutine update_hydrocarbon(mws, tcs, pcs, vcvis, sg_hc, ag)
  real(real64), intent(inout) :: mws(ncomp), tcs(ncomp), pcs(ncomp), vcvis(ncomp)
  real(real64), intent(in)    :: sg_hc
  logical,      intent(in)    :: ag
  real(real64) :: tpc, ppc, hc_mw
  call pseudo_critical(sg_hc, ag, tpc, ppc)
  hc_mw = sg_hc * MW_AIR
  tcs(5) = tpc; pcs(5) = ppc; mws(5) = hc_mw
  vcvis(5) = 0.0576710_real64 * (hc_mw - mwCH4) + 1.44383_real64
end subroutine

subroutine get_z_fractions(co2, h2s, n2, h2, zf)
  real(real64), intent(in)  :: co2, h2s, n2, h2
  real(real64), intent(out) :: zf(ncomp)
  zf = [co2, h2s, n2, h2, 1.0_real64 - (co2+h2s+n2+h2)]
end subroutine

function hydrocarbon_sg(sg_bulk, zf, mws) result(sg_hc)
  real(real64), intent(in) :: sg_bulk, zf(ncomp), mws(ncomp)
  real(real64) :: sg_hc, frac_hc, sum_nonhc
  frac_hc = zf(5)
  sum_nonhc = zf(1)*mws(1) + zf(2)*mws(2) + zf(3)*mws(3) + zf(4)*mws(4)
  if (frac_hc > 0.0_real64) then
    sg_hc = (sg_bulk - sum_nonhc/MW_AIR)/frac_hc
  else
    sg_hc = 0.75_real64
  end if
end function

function pr_m_param(acf) result(m)
  real(real64), intent(in) :: acf
  real(real64) :: m
  m = 0.37464_real64 + 1.54226_real64*acf - 0.26992_real64*acf*acf
end function

function fugacity_coeff(z, a, b) result(phi)
  real(real64), intent(in) :: z, a, b
  real(real64) :: phi, arg, lnphi
  if (z <= b) then; phi = huge(1.0_real64); return; end if
  arg = (z + (1.0_real64 + SQRT2)*b) / (z + (1.0_real64 - SQRT2)*b)
  if (arg <= 0.0_real64) then; phi = huge(1.0_real64); return; end if
  lnphi = (z - 1.0_real64) - log(z - b) - (a / (2.0_real64*SQRT2*b))*log(arg)
  phi = exp(lnphi)
end function

subroutine solve_pr_cubic_select(c2, c1, c0, a, b, zroot, ideal_fallback)
  real(real64), intent(in)  :: c2, c1, c0, a, b
  real(real64), intent(out) :: zroot
  logical,      intent(out) :: ideal_fallback
  real(real64) :: p, q, disc, a_div3, t, z, r, phiA, mA
  real(real64) :: roots(3), bestphi, phi_i
  integer :: n, i

  a_div3 = c2/3.0_real64
  p = c1 - c2*c2/3.0_real64
  q = (2.0_real64*c2*c2*c2)/27.0_real64 - (c2*c1)/3.0_real64 + c0
  disc = (q*q)/4.0_real64 + (p*p*p)/27.0_real64

  roots = 0.0_real64; n = 0
  if (disc > 0.0_real64) then
    t = cube_root(-q/2.0_real64 + sqrt(disc)) + cube_root(-q/2.0_real64 - sqrt(disc))
    z = t - a_div3
    if (z > 0.0_real64 .and. ieee_is_finite(z)) then; n=1; roots(1)=z; end if
  else if (abs(disc) <= 1.0e-14_real64) then
    t = cube_root(-q/2.0_real64)
    z =  2.0_real64*t - a_div3; if (z > 0.0_real64 .and. ieee_is_finite(z)) then; n=n+1; roots(n)=z; end if
    z = -1.0_real64*t - a_div3; if (z > 0.0_real64 .and. ieee_is_finite(z)) then; n=n+1; roots(n)=z; end if
  else
    r   = sqrt(-p*p*p/27.0_real64)
    phiA = acos( (-q/2.0_real64)/r )
    mA  = 2.0_real64*sqrt(-p/3.0_real64)
    do i=0,2
      t = mA * cos( (phiA + 2.0_real64*PI*real(i,real64))/3.0_real64 )
      z = t - a_div3
      if (z > 0.0_real64 .and. ieee_is_finite(z)) then; n=n+1; roots(n)=z; end if
    end do
  end if

  if (n == 0) then; zroot = 1.0_real64; ideal_fallback = .true.; return; end if
  if (n == 1 .or. .not.(a>0.0_real64 .and. b>0.0_real64)) then
    zroot = roots(1); ideal_fallback = .false.; return
  end if

  zroot = roots(1); ideal_fallback = .false.
  bestphi = fugacity_coeff(zroot, a, b)
  do i=2,n
    phi_i = fugacity_coeff(roots(i), a, b)
    if (phi_i < bestphi) then; bestphi = phi_i; zroot = roots(i); end if
  end do
end subroutine

subroutine calc_bips(degR, tpc_hc, kij, dkij_dT, d2kij_dT2)
  real(real64), intent(in)  :: degR, tpc_hc
  real(real64), intent(out) :: kij(ncomp,ncomp), dkij_dT(ncomp,ncomp), d2kij_dT2(ncomp,ncomp)
  integer :: i, j
  real(real64) :: c, s, tc, trm

  kij = 0.0_real64; dkij_dT = 0.0_real64; d2kij_dT2 = 0.0_real64
  do i=1,ncomp
    do j=1,ncomp
      if (i == j) cycle
      call pair_coeff(i,j,tpc_hc, c, s, tc)
      trm = degR / tc
      kij(i,j)      = c + s / trm
      dkij_dT(i,j)  = -s * tc / (degR*degR)
      d2kij_dT2(i,j)=  2.0_real64 * s * tc / (degR*degR*degR)
    end do
  end do
contains
  subroutine pair_coeff(i,j,tpc, c, s, tc)
    integer, intent(in) :: i,j
    real(real64), intent(in) :: tpc
    real(real64), intent(out) :: c, s, tc
    c=0.0_real64; s=0.0_real64; tc=1.0_real64
    if ((i==5 .and. j==1) .or. (i==1 .and. j==5)) then
      c=-0.145561_real64; s=0.276572_real64; tc=tpc
    else if ((i==5 .and. j==2) .or. (i==2 .and. j==5)) then
      c= 0.16852_real64;  s=-0.122378_real64; tc=tpc
    else if ((i==5 .and. j==3) .or. (i==3 .and. j==5)) then
      c=-0.108_real64;    s= 0.0605506_real64; tc=tpc
    else if ((i==5 .and. j==4) .or. (i==4 .and. j==5)) then
      c=-0.0620119_real64;s= 0.0427873_real64; tc=tpc
    else if ((i==1 .and. j==2) .or. (i==2 .and. j==1)) then
      c= 0.248638_real64; s=-0.138185_real64; tc=547.416_real64
    else if ((i==1 .and. j==3) .or. (i==3 .and. j==1)) then
      c=-0.25_real64;     s= 0.11602_real64;  tc=547.416_real64
    else if ((i==1 .and. j==4) .or. (i==4 .and. j==1)) then
      c=-0.247153_real64; s= 0.16377_real64;  tc=547.416_real64
    else if ((i==2 .and. j==3) .or. (i==3 .and. j==2)) then
      c=-0.204414_real64; s= 0.234417_real64; tc=672.12_real64
    else if ((i==2 .and. j==4) .or. (i==4 .and. j==2)) then
      c= 0.0_real64;      s= 0.0_real64;      tc=672.12_real64
    else if ((i==3 .and. j==4) .or. (i==4 .and. j==3)) then
      c=-0.166253_real64; s= 0.0788129_real64;tc=227.16_real64
    end if
  end subroutine
end subroutine

subroutine stiel_thodos_viscosity(tR, mws, tcs, pcs, mu)
  real(real64), intent(in)  :: tR, mws(ncomp), tcs(ncomp), pcs(ncomp)
  real(real64), intent(out) :: mu(ncomp)
  integer :: i
  ! NOTE: local named 'trm' (not 'tr') to avoid case-insensitive clash with dummy tR
  real(real64) :: trm, tc_k, pc_atm, eta
  do i=1,ncomp
    trm = tR / tcs(i)
    tc_k = tcs(i) * (5.0_real64/9.0_real64)
    pc_atm = pcs(i) / 14.696_real64
    eta = tc_k**(1.0_real64/6.0_real64) / (sqrt(mws(i)) * pc_atm**(2.0_real64/3.0_real64))
    if (trm <= 1.5_real64) then
      mu(i) = 34.0e-5_real64 * trm**0.94_real64 / eta
    else
      mu(i) = 17.78e-5_real64 * ((4.58_real64*trm) - 1.67_real64)**(5.0_real64/8.0_real64) / eta
    end if
  end do
end subroutine

subroutine methane_adjust(cp_poly, hc_mw, out_poly)
  real(real64), intent(in)  :: cp_poly(ncomp,5), hc_mw
  real(real64), intent(out) :: out_poly(ncomp,5)
  real(real64), parameter :: a0(5) = [7.8570e-04_real64, 1.3123e-03_real64, 9.8133e-04_real64, &
                                      1.6463e-03_real64, 1.7306e-02_real64]
  real(real64), parameter :: a1(5) = [-8.1649e-03_real64, 5.5485e-03_real64, 8.3258e-02_real64, &
                                      2.0635e-01_real64, 2.5551e+00_real64]
  real(real64) :: x, scale
  integer :: k
  out_poly = cp_poly
  x = hc_mw - mwCH4
  do k=1,5
    scale = a0(k)*x*x + a1(k)*x + 1.0_real64
    out_poly(5,k) = cp_poly(5,k) * scale
  end do
end subroutine

function lbc_viscosity(z, degF, psia, sg, co2, h2s, n2, h2, ag) result(mu_mix)
  real(real64), intent(in) :: z, degF, psia, sg, co2, h2s, n2, h2
  logical,      intent(in) :: ag
  real(real64) :: mu_mix
  real(real64) :: zf(ncomp), mws(ncomp), tcs(ncomp), pcs(ncomp), acf(ncomp), vshift(ncomp)
  real(real64) :: omega_a(ncomp), omega_b(ncomp), vcvis(ncomp), cp_poly(ncomp,5)
  real(real64) :: tR, mu_i(ncomp), sqrt_mw(ncomp), numerator, denom, rhoc, gasdens, rhor
  real(real64) :: a(5), lhs, tc_mix_k, pc_mix_atm, mw_mix, eta_mix
  real(real64) :: sg_hc, hc_mw
  integer :: i

  call base_props(mws,tcs,pcs,acf,vshift,omega_a,omega_b,vcvis,cp_poly)
  call get_z_fractions(co2,h2s,n2,h2, zf)

  sg_hc = hydrocarbon_sg(sg, zf, mws)
  sg_hc = clamp_min(sg_hc, mwCH4/MW_AIR)
  hc_mw = sg_hc * MW_AIR
  call update_hydrocarbon(mws,tcs,pcs,vcvis, sg_hc, ag)

  tR = degF + DEG_F_TO_R
  call stiel_thodos_viscosity(tR, mws, tcs, pcs, mu_i)

  do i=1,ncomp
    sqrt_mw(i) = sqrt(mws(i))
  end do
  numerator = sum( zf * mu_i * sqrt_mw )
  denom     = sum( zf * sqrt_mw )
  mu_mix    = numerator / denom

  rhoc = 1.0_real64 / dot5(vcvis, zf)
  gasdens = psia / (z * Rgas * tR)
  rhor = gasdens / rhoc

  a = [0.1023_real64, 0.023364_real64, 0.058533_real64, -0.0392852_real64, 0.00926279_real64]
  lhs = a(1) + a(2)*rhor + a(3)*rhor**2 + a(4)*rhor**3 + a(5)*rhor**4

  tc_mix_k = dot5(tcs, zf) * (5.0_real64/9.0_real64)
  pc_mix_atm = dot5(pcs, zf) / 14.696_real64
  mw_mix = dot5(mws, zf)
  eta_mix = tc_mix_k**(1.0_real64/6.0_real64) / (sqrt(mw_mix) * pc_mix_atm**(2.0_real64/3.0_real64))

  mu_mix = (lhs**4 - 1.0e-4_real64)/eta_mix + mu_mix
end function

subroutine pr_properties( &
    temp, pres, sg, co2, h2s, n2, h2, &
    ag, viscosity, density, thermo, metric, verbose, &
    z, rho_out, h_out, cp_out, cv_out, jt_out, mu_out)

  real(real64), intent(in)  :: temp, pres, sg, co2, h2s, n2, h2
  logical,      intent(in)  :: ag, viscosity, density, thermo, metric, verbose
  real(real64), intent(out) :: z, rho_out, h_out, cp_out, cv_out, jt_out, mu_out

  real(real64) :: zf(ncomp), mws(ncomp), tcs(ncomp), pcs(ncomp), acf(ncomp), vshift(ncomp)
  real(real64) :: omega_a(ncomp), omega_b(ncomp), vcvis(ncomp), cp_poly0(ncomp,5), cp_poly(ncomp,5)
  real(real64) :: degF, psia, degR
  real(real64) :: sg_hc, hc_mw, trs(ncomp), prs_r(ncomp), kij(ncomp,ncomp), dkij(ncomp,ncomp), d2kij(ncomp,ncomp)
  real(real64) :: m_i(ncomp), alpha(ncomp), a_c_i(ncomp), a_i(ncomp), b_i(ncomp)
  real(real64) :: a_mix, b_mix, rt, a_dim, b_dim, c2, c1, c0, z_eos
  logical :: ideal_fallback
  real(real64) :: bi(ncomp), shift, mwt

  real(real64) :: d_alpha_dT(ncomp), d2_alpha_dT2(ncomp), da_i_dT(ncomp), d2a_i_dT2(ncomp)
  real(real64) :: da_mix_dT, d2a_mix_dT2, sqrt_ai_aj, nterm, d2term, d, daij, d2aij
  real(real64) :: aA, bB, dB_dT, dA_dT, dF_dz, dF_dA, dF_dB, dF_dT, dz_dT
  real(real64) :: t_k, cp_vals(ncomp), cp_ig, bdim, h_dep_ft3psia, h_dep_btu
  real(real64) :: h0(ncomp), h0_mix, t_ref_k, term, h_ig_btu, h_total_btu
  real(real64) :: dBdim_dT, num_n, den_d, dln_dT, x, dx_dT, dHdep_dT_btu, cp_total_btu
  real(real64) :: v, dV_dT, den4, dP_dV_T, cv_total_btu, mu_jt
  integer :: i, j

  if (co2 + h2s + n2 + h2 > 1.0_real64 + 1.0e-12_real64) then
    z = 1.0_real64; rho_out=0.0_real64; h_out=0.0_real64
    cp_out=0.0_real64; cv_out=0.0_real64; jt_out=0.0_real64; mu_out=0.0_real64
    return
  end if

  if (metric) then
    degF = (temp * 9.0_real64/5.0_real64) + 32.0_real64
    psia = pres * 145.0377_real64
  else
    degF = temp
    psia = pres
  end if
  degR = degF + DEG_F_TO_R

  call base_props(mws,tcs,pcs,acf,vshift,omega_a,omega_b,vcvis,cp_poly0)
  call get_z_fractions(co2,h2s,n2,h2, zf)

  sg_hc = hydrocarbon_sg(sg, zf, mws)
  sg_hc = clamp_min(sg_hc, mwCH4/MW_AIR)
  hc_mw = sg_hc * MW_AIR
  call update_hydrocarbon(mws,tcs,pcs,vcvis, sg_hc, ag)

  do i=1,ncomp
    trs(i) = degR / tcs(i)
    prs_r(i) = psia / pcs(i)
    m_i(i) = pr_m_param(acf(i))
    alpha(i) = (1.0_real64 + m_i(i)*(1.0_real64 - sqrt(trs(i))))**2
    a_c_i(i) = omega_a(i) * Rgas*Rgas * tcs(i)*tcs(i) / pcs(i)
    a_i(i)   = a_c_i(i) * alpha(i)
    b_i(i)   = omega_b(i) * Rgas * tcs(i) / pcs(i)
  end do

  call calc_bips(degR, TcCH4, kij, dkij, d2kij)

  a_mix = 0.0_real64
  do i=1,ncomp
    do j=1,ncomp
      a_mix = a_mix + zf(i)*zf(j)*sqrt(a_i(i)*a_i(j))*(1.0_real64 - kij(i,j))
    end do
  end do
  b_mix = dot5(b_i, zf)

  rt   = Rgas * degR
  a_dim = a_mix * psia / (rt*rt)
  b_dim = b_mix * psia / rt

  c2 = -(1.0_real64 - b_dim)
  c1 = a_dim - 3.0_real64*b_dim*b_dim - 2.0_real64*b_dim
  c0 = -(a_dim*b_dim - b_dim*b_dim - b_dim*b_dim*b_dim)
  call solve_pr_cubic_select(c2,c1,c0, a_dim, b_dim, z_eos, ideal_fallback)

  shift = 0.0_real64
  do i=1,ncomp
    bi(i) = omega_b(i) * (prs_r(i) / trs(i))
    shift = shift + zf(i) * vshift(i) * bi(i)
  end do
  z = z_eos - shift

  mwt = dot5(mws, zf)
  rho_out = 0.0_real64
  if (density) rho_out = mwt * psia / (z * Rgas * degR)

  h_out = 0.0_real64; cp_out = 0.0_real64; cv_out = 0.0_real64; jt_out = 0.0_real64
  if (thermo) then
    do i=1,ncomp
      d_alpha_dT(i)  = (-m_i(i) * (1.0_real64 + m_i(i)*(1.0_real64 - sqrt(trs(i)))) / &
                        sqrt(trs(i))) / tcs(i)
      d2_alpha_dT2(i)= ( 2.0_real64 * (-m_i(i)/(2.0_real64*sqrt(trs(i))))**2 + &
                         2.0_real64 * (1.0_real64 + m_i(i)*(1.0_real64 - sqrt(trs(i)))) * &
                         (m_i(i)/(4.0_real64*trs(i)**1.5_real64)) ) / (tcs(i)*tcs(i))
      da_i_dT(i)   = a_c_i(i) * d_alpha_dT(i)
      d2a_i_dT2(i) = a_c_i(i) * d2_alpha_dT2(i)
    end do

    da_mix_dT = 0.0_real64; d2a_mix_dT2 = 0.0_real64
    do i=1,ncomp
      do j=1,ncomp
        sqrt_ai_aj = sqrt(a_i(i)*a_i(j))
        nterm = da_i_dT(i)*a_i(j) + a_i(i)*da_i_dT(j)
        daij = -dkij(i,j)*sqrt_ai_aj + (1.0_real64 - kij(i,j))*0.5_real64*(nterm/sqrt_ai_aj)
        da_mix_dT = da_mix_dT + zf(i)*zf(j)*daij

        d2term = d2a_i_dT2(i)*a_i(j) + 2.0_real64*da_i_dT(i)*da_i_dT(j) + a_i(i)*d2a_i_dT2(j)
        d = sqrt_ai_aj
        d2aij = -(nterm/d)*dkij(i,j) + (1.0_real64 - kij(i,j))*0.5_real64 * &
                ((d2term/d) - (nterm*nterm)/(2.0_real64*d**3)) - d*d2kij(i,j)
        d2a_mix_dT2 = d2a_mix_dT2 + zf(i)*zf(j)*d2aij
      end do
    end do

    aA = a_mix * psia / (Rgas*degR)**2
    bB = b_mix * psia / (Rgas*degR)
    dB_dT = -bB / degR
    dA_dT = (psia/(Rgas*degR)**2) * da_mix_dT - 2.0_real64 * a_mix*psia / (Rgas*Rgas*degR**3)

    dF_dz = 3.0_real64*z_eos*z_eos - 2.0_real64*(1.0_real64 - bB)*z_eos + &
            (aA - 3.0_real64*bB*bB - 2.0_real64*bB)
    dF_dA = z_eos - bB
    dF_dB = z_eos*z_eos - (6.0_real64*bB + 2.0_real64)*z_eos - aA + 2.0_real64*bB + 3.0_real64*bB*bB
    dF_dT = dF_dA*dA_dT + dF_dB*dB_dT
    dz_dT = - dF_dT / dF_dz

    call methane_adjust(cp_poly0, hc_mw, cp_poly)
    t_k = degR * (5.0_real64/9.0_real64)
    do i=1,ncomp; cp_vals(i) = poly5_eval(cp_poly(i,:), t_k); end do
    cp_ig = dot5(cp_vals, zf) * Rthermo

    bdim = (b_mix * psia) / (Rgas * degR)
    h_dep_ft3psia = Rgas*degR*(z_eos - 1.0_real64) + &
                    (degR*da_mix_dT - a_mix) / (2.0_real64*SQRT2*b_mix) * &
                    log( (z_eos + (SQRT2+1.0_real64)*bdim) / (z_eos - (SQRT2-1.0_real64)*bdim) )
    h_dep_btu = h_dep_ft3psia / FT3PSIA_TO_BTU

    h0 = [-16.6022_real64, -21.5512_real64, -3.57757_real64, 0.008054_real64, 0.0_real64]
    h0(5) = -0.015774_real64*(hc_mw - mwCH4)**2 - 0.646645_real64*(hc_mw - mwCH4) - 8.2551915_real64
    h0_mix = dot5(h0, zf)

    t_ref_k = (60.0_real64 + DEG_F_TO_R) * (5.0_real64/9.0_real64)
    h_ig_btu = 0.0_real64
    do i=1,ncomp
      term = cp_poly(i,1)*(t_k - t_ref_k) + &
             cp_poly(i,2)/2.0_real64*(t_k*t_k - t_ref_k*t_ref_k) + &
             cp_poly(i,3)/3.0_real64*(t_k**3 - t_ref_k**3) + &
             cp_poly(i,4)/4.0_real64*(t_k**4 - t_ref_k**4) + &
             cp_poly(i,5)/5.0_real64*(t_k**5 - t_ref_k**5)
      h_ig_btu = h_ig_btu + zf(i)*Rthermo*(9.0_real64/5.0_real64)*term
    end do

    h_total_btu = h_ig_btu + h_dep_btu - h0_mix

    dBdim_dT = - bdim / degR
    num_n = z_eos + (SQRT2 + 1.0_real64)*bdim
    den_d = z_eos - (SQRT2 - 1.0_real64)*bdim
    dln_dT = (dz_dT + (SQRT2+1.0_real64)*dBdim_dT)/num_n - &
             (dz_dT - (SQRT2-1.0_real64)*dBdim_dT)/den_d
    x = (degR*da_mix_dT - a_mix) / (2.0_real64*SQRT2*b_mix)
    dx_dT = (degR * d2a_mix_dT2) / (2.0_real64*SQRT2*b_mix)
    dHdep_dT_btu = ( Rgas*(z_eos - 1.0_real64) + Rgas*degR*dz_dT + &
                     dx_dT*log(num_n/den_d) + x*dln_dT ) / FT3PSIA_TO_BTU

    cp_total_btu = cp_ig + dHdep_dT_btu

    v = z_eos * Rgas * degR / psia
    dV_dT = (Rgas/psia) * (z_eos + degR*dz_dT)
    den4 = (v*v + 2.0_real64*b_mix*v - b_mix*b_mix)**2
    dP_dV_T = - Rgas*degR / (v - b_mix)**2 + 2.0_real64*a_mix*(v + b_mix)/den4
    cv_total_btu = cp_total_btu + (degR * dV_dT*dV_dT * dP_dV_T) / FT3PSIA_TO_BTU

    mu_jt = (degR * dV_dT - v) / (cp_total_btu * FT3PSIA_TO_BTU)

    h_out  = h_total_btu
    cp_out = cp_total_btu
    cv_out = cv_total_btu
    jt_out = mu_jt
  end if

  mu_out = 0.0_real64
  if (viscosity) mu_out = lbc_viscosity(z, degF, psia, sg, co2, h2s, n2, h2, ag)

  if (metric) then
    if (density) rho_out = rho_out * 16.01846337396_real64
    if (thermo) then
      h_out  = h_out  * 2.326_real64
      cp_out = cp_out * 4.186800585_real64
      cv_out = cv_out * 4.186800585_real64
      jt_out = jt_out * 80.576521_real64
    end if
  end if
end subroutine

end module bns

program demo
  use iso_fortran_env, only: real64
  use bns
  implicit none
  real(real64) :: z, rho, h, cp, cv, jt, mu

  call pr_properties( 100.0_real64, 1500.0_real64, 0.65_real64, &
                      0.05_real64, 0.0_real64, 0.02_real64, 0.0_real64, &
                      .false., &
                      .true., .true., .true., .false., .false., &
                      z, rho, h, cp, cv, jt, mu)

  print '(A,F18.10)',    'Z   = ', z
  print '(A,F18.10, A)', 'rho = ', rho, ' (lbm/ft^3)'
  print '(A,F18.10, A)', 'mu  = ', mu,  ' (cP)'
  print '(A,F18.10, A)', 'H  = ', cp,  ' (Btu/(lb-mol))'
  print '(A,F18.10, A)', 'Cp  = ', cp,  ' (Btu/(lb-mol·R))'
  print '(A,F18.10, A)', 'Cv  = ', cv,  ' (Btu/(lb-mol·R))'
  print '(A,F18.10, A)', 'JT  = ', jt,  ' (F/psi)'
end program
