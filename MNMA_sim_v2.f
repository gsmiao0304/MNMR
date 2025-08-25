program MNMA_Model
! Parameters
integer, parameter :: ns = 103, nsc = 79, K = 40, nx = 10, nT = 10, &
     nTT = 11, nrep = 10000, J = 3, nthin = 5, nwarm = 5000, &
     nbeta = 10, nb = 4, nphi = 6, nbbeta = nbeta + nb + nT, ngrp = 1
real(8), parameter :: epi = 3.141592741012573d0, &
     tol = 2.220446049250313d-14
! Integer arrays
integer :: nldlc(ns), nhdlc(ns), ntg(ns), grp(ns), bk(ns)
integer :: ids(ns), iarm(ns), iT(ns), narms(K), iobs(K)
! Real arrays
real(8) :: x(nx, ns), xc(nx, ns), ny(ns), xx(nbbeta, ns), xpre(nx, ns)
real(8) :: xmean(nx), xsd(nx), x2mean(nx)
real(8) :: w(ns, nphi), xk(K, ns), xt(nT, ns), xb(nb, ns)
real(8) :: yldlc(ns), yhdlc(ns), ytg(ns), sdldlc(ns), sdhdlc(ns), &
     sdtg(ns)
real(8) :: y(J, ns), grp_n(ngrp)
real(8) :: beta(nbeta, J), bbeta(nbbeta, J), seqbeta(J, nbbeta, nrep), beta_true(J,nbbeta)
real(8) :: phi(nphi, J), seqphi(J, nphi, nrep) , phi_true(J,nphi)
real(8) :: alpha(J), seqalpha(J, nrep)
real(8) :: tauk(ns, J, J), Rhok(ns, J, J)
real(8) :: Rrho(J), Erho(ns, 3), Srho(ns, 3)
real(8) :: seqRrho(J, nrep), seqErho(ns, 3, nrep), seqSrho(ns, 3, nrep)
real(8) :: Rgam(ns, J), seqRgam(ns, J, nrep)
real(8) :: lam(K, J), seqlam(J, K, nrep)
real(8) :: Omega(J, nTT, nTT), kOmega(J, K, nTT, nTT), pOmega(J, nTT, nTT)
real(8) :: Rho(nTT, nTT), pRho(nTT, nTT), seqrho(nrep, J, 55)
real(8) :: eta(J), seqeta(J, nrep)
real(8) :: VoS(ns, J, J), VRV(ns, J, J), seqVRV(nrep, ns, J, J), pVRV(ns, J, J)
real(8) :: pRoS(ns, J, J), RoS(ns, J, J), zRoS(ns, J)
real(8) :: SigInv(ns, J, J), sigk(ns, J, J)
real(8) :: sig(ngrp, J, J), sig0(J, J)
real(8) :: sigma(ns, 6), seqsigma(ns, 6, nrep)
real(8) :: delta(ngrp, J), gamma(ngrp, J)
real(8) :: seqdelta(ngrp, J, nrep), seqgamma(ngrp, J, nrep)
real(8) :: seqzRoS(J, ns, nrep)
real(8) :: E(ns, nTT), EO(J, 12), EO2(J, 6), EO3(J, 9), EO4(J, 12)
real(8) :: v, seqv(nrep),ve
integer :: narmsk, icount, temp, icount4, idim, tv2, dim1, dim2,jj
real(8) :: tmp, c, phitmp, a, b, det
real(8) :: ss1, ss2, ss3, ss4, ss5, ss6, ss7
real(8) :: lg1, lg2, lg3
real(8) :: m1, m2, m3, m4, m5, m6, m7
real(8) :: sigtmp(J, J), sigitmp(J, J), ssig(J, J)
real(8) :: LR(J, J), LRvec(J), tv(J)
real(8) :: t3(1, 1), t4(1, 1), t5(1, 1), t6(1, 1), t7(1, 1)
real(8) :: yDIC(ns), ypD(ns), ybarD(ns), yDbar(ns)
real(8) :: sDIC(ns), spD(ns), sbarD(ns), sDbar(ns)
real(8) :: y1DIC(ns), y1pD(ns), y1barD(ns), y1Dbar(ns)
real(8) :: s1DIC(ns), s1pD(ns), s1barD(ns), s1Dbar(ns)
real(8) :: sumDIC(ns), sumpD(ns), sumbarD(ns), sumDbar(ns)
real(8) :: lDIC(ns), lpD(ns), lbarD(ns), lDbar(ns)
real(8) :: hDIC(ns), hpD(ns), hbarD(ns), hDbar(ns)
real(8) :: tDIC(ns), tpD(ns), tbarD(ns), tDbar(ns)
real(8) :: DICall, pDall, DICall1, pDall1
real(8) :: ipd_dic_s, ipd_pd_s, ipd_dic_y, ipd_pd_y,ipd_dic,ipd_pd
real(8) :: lDICall, hDICall, tDICall, lpDall, hpDall, tpDall
real(8) :: seqphi_gam(ns, J, nrep), seqsig(nrep, ns, J, J)
real(8) :: pbeta(nbbeta, J), pphi_gam(ns, J), psig(ns, J, J)
real(8) :: pphi(nphi, J), ptauk(ns, J, J)
real(8) :: beta_out(nbbeta, J, 7), beta_tmp(nrep)
real(8) :: phi_out(nphi, J, 7), phi_tmp(nrep)
real(8) :: Sigk2(K, 6, 6), Sigk3(K, 9, 9), Sigk4(K, 12, 12)
real(8) :: WOW(12, 12)
real(8) :: WOW2(6, 6), WOW3(9, 9), WOW4(12, 12)
real(8) :: WOWInv2(6, 6), WOWInv3(9, 9), WOWInv4(12, 12)
real(8) :: Sigpre2(6, 6), Sigpre3(9, 9), Sigpre4(12, 12)
real(8) :: SigGam2(6, 6), SigGam3(9, 9), SigGam4(12, 12)
real(8) :: LR2(6, 6), LR3(9, 9), LR4(12, 12)
real(8) :: LRvec2(6), LRvec3(9), LRvec4(12)
real(8) :: tv72(6), tv73(9), tv74(12), tv1(12, 12), tv7(3)
real(8) :: Ome(J, nTT, nTT), p(K, nrep), logp(K, nrep)
real(8) :: InvSig_t(J, J), Sig_t(J, J), tau(J, J)
real(8) :: s12(1, 2), s22(2, 2)
real(8) :: sig12(1, 2), sig22(2, 2), sig22i(2, 2)
real(8) :: ldenom(nrep, ns), lpdenom(ns)
real(8) :: lprod1(nrep, ns), lpprod1(ns)
real(8) :: lprod2(nrep, ns), lpprod2(ns)
real(8) :: lprod3(nrep, ns), lpprod3(ns)
real(8) :: lprod4(nrep, ns), lpprod4(ns)
real(8) :: hdenom(nrep, ns), hpdenom(ns)
real(8) :: hprod1(nrep, ns), hpprod1(ns)
real(8) :: hprod2(nrep, ns), hpprod2(ns)
real(8) :: hprod3(nrep, ns), hpprod3(ns)
real(8) :: hprod4(nrep, ns), hpprod4(ns)
real(8) :: tdenom(nrep, ns), tpdenom(ns)
real(8) :: tprod1(nrep, ns), tpprod1(ns)
real(8) :: tprod2(nrep, ns), tpprod2(ns)
real(8) :: tprod3(nrep, ns), tpprod3(ns)
real(8) :: tprod4(nrep, ns), tpprod4(ns)
real(8) :: lmu1(nrep, ns), lpmu1(ns)
real(8) :: hmu1(nrep, ns), hpmu1(ns)
real(8) :: tmu1(nrep, ns), tpmu1(ns)
real(8) :: y2(2, 1), mu2(2, 1)
real(8) :: smat1(2, 2), smat2(2, 2), smat3(2, 2)
real(8) :: smat4(1, 2), smat5(1, 2), smat6(1, 2)
real(8) :: pa1(2, 2), pa2(2, 2), pa3(2, 2), pa4(1, 2), pa5(1, 2)
real(8) :: lhdenom(nrep, ns, 2, 2), lhpdenom(ns, 2, 2)
real(8) :: lhprod1(nrep, ns, 2, 2), lhpprod1(ns, 2, 2)
real(8) :: lhprod2(nrep, ns, 2, 2), lhpprod2(ns, 2, 2)
real(8) :: lhprod3(nrep, ns, 1, 2), lhpprod3(ns, 1, 2)
real(8) :: lhprod4(nrep, ns, 1, 2), lhpprod4(ns, 1, 2)
real(8) :: lhmu(nrep, ns, 1, 2), lhpmu(ns, 1, 2)
real(8) :: lhDIC(ns), lhpD(ns), lhbarD(ns), lhDbar(ns),lhDICall,lhpDall
real(8) :: htdenom(nrep, ns, 2, 2), htpdenom(ns, 2, 2)
real(8) :: htprod1(nrep, ns, 2, 2), htpprod1(ns, 2, 2)
real(8) :: htprod2(nrep, ns, 2, 2), htpprod2(ns, 2, 2)
real(8) :: htprod3(nrep, ns, 1, 2), htpprod3(ns, 1, 2)
real(8) :: htprod4(nrep, ns, 1, 2), htpprod4(ns, 1, 2)
real(8) :: htmu(nrep, ns, 1, 2), htpmu(ns, 1, 2)
real(8) :: htDIC(ns), htpD(ns), htbarD(ns), htDbar(ns),htDICall,htpDall
real(8) :: ltdenom(nrep, ns, 2, 2), ltpdenom(ns, 2, 2)
real(8) :: ltprod1(nrep, ns, 2, 2), ltpprod1(ns, 2, 2)
real(8) :: ltprod2(nrep, ns, 2, 2), ltpprod2(ns, 2, 2)
real(8) :: ltprod3(nrep, ns, 1, 2), ltpprod3(ns, 1, 2)
real(8) :: ltprod4(nrep, ns, 1, 2), ltpprod4(ns, 1, 2)
real(8) :: ltmu(nrep, ns, 1, 2), ltpmu(ns, 1, 2)
real(8) :: ltDIC(ns), ltpD(ns), ltbarD(ns), ltDbar(ns),ltDICall,ltpDall
real(8) :: DIC_out(24)
real(8) :: aupp(2), alow(2)
real(8) :: sum, sumsq, var, mse
character(100) :: datafile, infile, outfile
! COMMON blocks for arrays
common /vecy/ y
common /vecny/ ny
common /vecx/ x
common /vecxx/ xx
common /vecw/ w
common /vecids/ ids
common /veciarm/ iarm
common /vecnarms/ narms
common /vecobs/ iobs
common /vecgrp/ grp
common /vecgrpn/ grp_n
common /vecbbeta/ bbeta
common /vecRgam/ Rgam
common /veclam/ lam
common /vecphi/ phi
common /vecalpha/ alpha
common /vecdelta/ delta
common /vecgamma/ gamma
common /vectauk/ tauk
common /vecrhok/ Rhok
common /vecOmega/ Omega
common /vecpOmega/ pOmega
common /veckOmega/ kOmega
common /vecSigInv/ SigInv
common /vecsigk/ sigk
common /vecSigma/ sigma
common /vecSig/ sig
common /vecSigo/ sig0
common /veceta/ eta
common /vecVoS/ VoS
common /vecpRoS/ pRoS
common /veczRoS/ zRoS
common /vecRoS/ RoS
common /vecVRV/ VRV
common /vecRrho/ Rrho
common /vecErho/ Erho
common /vecSrho/ Srho
common /vecE/ E
common /numv/ v
common /numve/ ve

! External functions
real(8), external :: alngam, FindDet

! Time variables
integer :: now(3), nowb(3), nowe(3), ntoday(3), ntodayb(3)
real :: etime, elapsed(2), total

! Character variables
character(18) :: xname(10)
character(15) :: tname(nT)

! Initialize time variables
call idate(ntoday)
ntodayb(1) = ntoday(1)
ntodayb(2) = ntoday(2)
ntodayb(3) = ntoday(3)

call itime(now)
nowb(1) = now(1)
nowb(2) = now(2)
nowb(3) = now(3)

! Set initial random number seed
iseed = 2333

! Read the data
infile = '~/MMNAC.txt'
!infile = './Data/' // trim(datafile) 

open(unit=12, file=infile, status='old', iostat=ios)
if (ios /= 0) then
  print *, 'Error opening file: ', infile
  stop
end if

! Read data from file
do i = 1, ns
  read(12, *) ids(i), iarm(i),iT(i), bk(i), nldlc(i), yldlc(i), &
              sdldlc(i), nhdlc(i), yhdlc(i), sdhdlc(i), ntg(i), ytg(i), sdtg(i), &
              (xpre(jj, i), jj=1, nx)

  ! Write formatted output for verification
  write(*, 700) ids(i), iarm(i),iT(i), bk(i), nldlc(i), yldlc(i), &
                sdldlc(i), nhdlc(i), yhdlc(i), sdhdlc(i), ntg(i), ytg(i), sdtg(i), &
                (xpre(jj, i), jj=1, nx)
end do

close(12)
write(*, *) 'end data reading'

700  format(I3,1x,I3,1x,I3,1x,I3,1x,I5,1x,f8.2,1x,f8.2,1x,I5,1x,&
         f8.2,1x,f8.2,1x,I5,1x,f8.2,1x,f8.2,10f8.2)
           



! Calculate number of arms for each trial
do kk = 1, K
  narms(kk) = 0
  iobs(kk) = 0
end do

do i = 1, ns
  narms(ids(i)) = narms(ids(i)) + 1
end do

temp = 0
iobs(1) = 1
do kk = 2, K
  temp = narms(kk - 1)
  iobs(kk) = iobs(kk - 1) + temp
end do

do kk = 1, K
  write(*, 701) kk, narms(kk), iobs(kk)
end do

701 format('k=', I3, ', # of arms=', I3, ', 1st_record', I3)

! Prepare response variables
do i = 1, ns
  y(1, i) = yldlc(i)
  y(2, i) = yhdlc(i)
  y(3, i) = ytg(i)
  write(*, *) 'obs=', i, ' y=', y(:, i)
end do

! Prepare covariate matrices
xx(:, :) = 0.0d0
do jj = 1, nx
  do i = 1, ns
    ! Standardization code commented out
    ! x(jj, i) = (xpre(jj, i) - xmean(jj)) / xsd(jj)
    x(jj, i) = xpre(jj, i) !xpre is pre-scaled
  end do
end do

write(*, *) 'xpre obs=8, ', xpre(:, 8)
write(*, *) 'standardize obs=8', x(:, 8)

! Prepare design matrices
xb(:, :) = 0.0d0
!xb(1, :) = 1.0d0

do l1 = 1, ns
  select case (bk(l1))
    case (1)
      xb(1, l1) = 1.0d0    
    case (5)
      xb(2, l1) = 1.0d0
    case (6)
      xb(3, l1) = 1.0d0
    case (7)
      xb(4, l1) = 1.0d0
  end select
  
  dim1 = bk(l1)
  dim2 = iT(l1)
  
  xt(:, l1) = 0.0d0
  if (dim1 /= dim2) then
    xt(dim2 - 1, l1) = 1.0d0
    if (dim1 /= 1) then
      xt(dim1 - 1, l1) = -1.0d0
    end if
  end if
  !write(*,*) 'i=',l1, xt(:,l1)
end do

! Combine covariates into design matrix (xx)
do i = 1, ns
  do j2 = 1, nb
    xx(j2, i) = xb(j2, i)
  end do
  
  do j3 = 1, nx
    xx(nb + j3, i) = x(j3, i)
  end do
  
  do j1 = 1, nT
    xx(nx + nb + j1, i) = xt(j1, i)
  end do
  write(*,'(24F10.2)') (xx(ll, i),ll=1,nbbeta) 
end do

! Prepare covariates Z_kt
w(:, :) = 0.0d0
do i=1,ns
    w(i,1)=1.0d0
    w(i,2)=x(1,i)
    w(i,3)=x(3,i)
    if (iT(i) == 1) w(i,4)=1.0d0
    if (iT(i) == 7) w(i,5)=1.0d0
    if ((iT(i) == 8) .or. (iT(i) == 9) .or. (iT(i) == 10) .or. (iT(i) == 11)) then
        w(i,6)=1.0d0
    end if
end do

! Initialize covariance-related matrices
vos(:, :, :) = 0.0d0
pRoS(:, :, :) = 0.0d0
RoS(:, :, :) = 0.0d0
VRV(:, :, :) = 0.0d0

zRoS(:, 1) = -1.0d0
zRoS(:, 2) = 1.0d0
zRoS(:, 3) = 0.0d0

do i = 1, ns
  vos(i, 1, 1) = sdldlc(i)
  vos(i, 2, 2) = sdhdlc(i)
  vos(i, 3, 3) = sdtg(i)
end do

do i = 1, 87
  do j1 = 1, J
    pRoS(i, j1, j1) = 1.0d0
    RoS(i, j1, j1) = 1.0d0
    VRV(i, j1, j1) = vos(i, j1, j1)**2.0d0
  end do
  
  pRoS(i, 1, 2) = -0.3d0
  pRoS(i, 1, 3) = 0.6d0
  pRoS(i, 2, 3) = -0.6d0
  pRoS(i, 2, 1) = pRoS(i, 1, 2)   ! Maintain symmetry
  pRoS(i, 3, 1) = pRoS(i, 1, 3)
  pRoS(i, 3, 2) = pRoS(i, 2, 3)
end do
! Set initial values
v = 2.5d0
ve = exp(v) + 4.0d0

! Initialize phi matrix
phi(:, :) = 0.0d0 

! Initialize other parameters
alpha(:) = 0.0d0
lam(:, :) = 1.0d0 
Rgam(:, :) = 0.0d0   
Rhok(:, :, :) = 0.0d0
SigInv(:, :, :) = 0.0d0
sigk(:, :, :) = 0.0d0
sig(:, :, :) = 0.0d0
sig0(:, :) = 0.0d0
delta(:, :) = 200.0d0
gamma(:, :) = 0.0d0
Srho(:, :) = 0.0d0

! Initialize sig matrices
do i1 = 1, ngrp
  sig(i1, 1, 1) = 1.0d0
  sig(i1, 2, 2) = 1.0d0
  sig(i1, 3, 3) = 1.0d0
end do

! Initialize Rhok and SigInv matrices
do i1 = 1, ns
  do j2 = 1, J
    Rhok(i1, j2, j2) = 1.0d0
    SigInv(i1, j2, j2) = 1.0d0
  end do
  sigk(i1, 1, 1) = 200.0d0
  sigk(i1, 2, 2) = 150.0d0
  sigk(i1, 3, 3) = 400.0d0
end do

! Initialize sig0 matrix
do jj = 1, J
  sig0(jj, jj) = 1.0d0 / 10.0d0
end do

! Initialize bbeta
bbeta(:, :) = 0.0d0

! Initialize tauk matrix
do ii = 1, ns
  tauk(ii, :, :) = 0.0d0
  do l1 = 1, J
    ! Calculate tauk using dot product
    tauk(ii, l1, l1) = exp(dot_product(w(ii, :), phi(:, l1)))
  end do
end do

! Initialize correlation parameters
seqRrho(:, :) = 0.0d0
Erho(:, :) = 0.0d0
seqErho(:, :, :) = 0.0d0
seqSrho(:, :, :) = 0.0d0

! Calculate ny values
do i = 1, ns
  ny(i) = real(min(nldlc(i), nhdlc(i), ntg(i)))
end do
!write(*, *) 'ny', ny(:)

! Set up Omega matrix
kOmega(:, :, :, :) = 1.0d0
Omega(:, :, :) = 1.0d0

do l1 = 1, J
  do l2 = 1, nTT - 1
    do l3 = l2 + 1, nTT
      Omega(l1, l2, l3) = 0.1d0
      Omega(l1, l3, l2) = Omega(l1, l2, l3)
    end do
  end do
  
  do kk = 1, K
    kOmega(l1, kk, :, :) = Omega(l1, :, :) / lam(kk, l1)
  end do
end do

pOmega(:, :, :) = Omega(:, :, :)

! Initialize entire E matrix to zero
E = 0.0d0

do i = 1, ns
  E(i, iT(i)) = 1.0d0
!  print *, 'i=', i,' E(i,)=', E(i,:)
end do

! Call subroutine
call vars1()
! Warm up Gibbs sampling
do i1 = 1, nwarm
   write(*, *) 'iwarm=', i1
!  print *, 'Main: iseed=', iseed
  call gibbs(iseed)
end do

! Write Gibbs samples
do i1 = 1, nrep
  write(*, *) 'irep=', i1
  
  ! Thinning loop
  do ithin = 1, nthin
    call gibbs(iseed)
  end do
  
  ! Store beta parameters
  do i2 = 1, nbbeta
    do jj = 1, J
      seqbeta(jj, i2, i1) = bbeta(i2, jj)
    end do
  end do
  write(*,*) 'LDLC Gamma-fixed:'
  write(*, '(10F10.4)') (bbeta(ll, 1),ll=15,nbbeta) 
  
  ! Store alpha and Rrho parameters
  do jj = 1, J
    seqalpha(jj, i1) = alpha(jj)
    seqRrho(jj, i1) = Rrho(jj)
  end do
  
  ! Store delta and gamma parameters
  do ii = 1, ngrp
    do jj = 1, J
      seqdelta(ii, jj, i1) = delta(ii, jj)
      seqgamma(ii, jj, i1) = gamma(ii, jj)
    end do
  end do
  write(*,*) 'delta'
  write(*, '(3F10.4)') (delta(1, ll),ll=1,J) 
  
  ! Store Erho, Srho, and Rgam parameters
  do i2 = 1, ns
    do jj = 1, J
      seqErho(i2, jj, i1) = Erho(i2, jj)
      seqSrho(i2, jj, i1) = Srho(i2, jj)
      seqRgam(i2, jj, i1) = Rgam(i2, jj)
    end do
  end do
  
  ! Store sigma parameters
  do i2 = 1, ns
    do jj = 1, 6
      seqsigma(i2, jj, i1) = sigma(i2, jj)
    end do
    !if (i2 > 79) then
    !  write(*, '(6F10.4)') (sigma(i2, ll),ll=1,6) 
    !end if 
  end do
  !write(*,*) 'sigma obs=1'
  
  ! Store phi parameters
  do i2 = 1, nphi
    do jj = 1, J
      seqphi(jj, i2, i1) = phi(i2, jj)
    end do
  end do
  
  ! Store v and lam parameters
  seqv(i1) = v
  write(*, *) 'df=', v
  
  do i2 = 1, K
    do jj = 1, J
      seqlam(jj, i2, i1) = lam(i2, jj)
    end do
  end do
  
  ! Store VRV matrices
  do i2 = 1, ns
    seqVRV(i1, i2, :, :) = VRV(i2, :, :)
  end do
  
  ! store Rho parameters
  do jj =1, J
    icount = 1
    do l1 = 1, (nTT - 1)
      do l2 = (l1+1), nTT
        seqrho(i1,jj,icount) = Omega(jj,l1,l2)
        icount = icount + 1
      end do 
    end do
  end do
end do
! Calculate posterior statistics for beta parameters
pbeta(:, :) = 0.0d0  ! Initialize posterior mean matrix for beta
do i1 = 1, nbbeta    ! Loop over all beta parameters
  do i2 = 1, J       ! Loop over outcomes (J dimensions)
    tmp = 0.0d0       ! Temporary sum for mean calculation
    sumsq = 0.0d0     ! Temporary sum for sum of squares
    mse = 0.0d0
    
    ! Accumulate samples for mean and variance calculation
    do i3 = 1, nrep  ! Loop through all MCMC samples
      tmp = tmp + seqbeta(i2, i1, i3)           ! Sum of samples
      sumsq = sumsq + seqbeta(i2, i1, i3)**2.0d0    ! Sum of squares
     ! mse = mse + (seqbeta(i2, i1, i3)-beta_true(i2,i1))**2.0d0
    end do
    mse = mse / real(nrep, 8)
    beta_out(i1, i2, 3) = mse
    
    ! Calculate posterior mean
    pbeta(i1, i2) = tmp / real(nrep, 8)
    
    ! Calculate sample variance (unbiased estimator)
    var = (sumsq - real(nrep, 8) * pbeta(i1, i2)**2.0d0) / real(nrep - 1, 8)
    
    ! Store posterior mean and standard deviation
    beta_out(i1, i2, 1) = pbeta(i1, i2)   ! Posterior mean
    beta_out(i1, i2, 2) = sqrt(var)        ! Posterior standard deviation
    
    ! Prepare array for HPD interval calculation
    beta_tmp = seqbeta(i2, i1, :)  ! Extract all samples for this parameter
    
    ! Calculate Highest Posterior Density (HPD) interval
    call hpd_interval(nrep, 0.05d0, beta_tmp, alow, aupp)
    
    ! Store HPD interval results (two methods)
    beta_out(i1, i2, 4) = alow(1)  ! Lower bound (method 1)
    beta_out(i1, i2, 5) = aupp(1)  ! Upper bound (method 1)
    beta_out(i1, i2, 6) = alow(2)  ! Lower bound (method 2)
    beta_out(i1, i2, 7) = aupp(2)  ! Upper bound (method 2)
  end do
end do
!print *, 'pbeta=', pbeta

! Calculate posterior statistics for phi parameters
do i1 = 1, nphi     ! Loop over all phi parameters
  do i2 = 1, J      ! Loop over outcomes (J dimensions)
    tmp = 0.0d0     ! Temporary sum for mean calculation
    sumsq = 0.0d0   ! Temporary sum for sum of squares
    mse = 0.0d0
    
    ! Accumulate samples for mean and variance calculation
    do i3 = 1, nrep  ! Loop through all MCMC samples
      tmp = tmp + seqphi(i2, i1, i3)           ! Sum of samples
      sumsq = sumsq + seqphi(i2, i1, i3)**2.0d0    ! Sum of squares
     ! mse = mse + (seqphi(i2, i1, i3)-phi_true(i2,i1))**2.0d0
    end do
    mse = mse / real(nrep, 8)
    phi_out(i1, i2, 3) = mse
    
    ! Calculate posterior mean
    phi_out(i1, i2, 1) = tmp / real(nrep, 8)
    
    ! Calculate sample variance (unbiased estimator)
    var = (sumsq - real(nrep, 8) * phi_out(i1, i2, 1)**2.0d0) / real(nrep - 1, 8)
    
    ! Store posterior mean and standard deviation
    phi_out(i1, i2, 2) = sqrt(var)  ! Posterior standard deviation
    
    ! Prepare array for HPD interval calculation
    phi_tmp = seqphi(i2, i1, :)  ! Extract all samples for this parameter
    
    ! Calculate Highest Posterior Density (HPD) interval
    call hpd_interval(nrep, 0.05d0, phi_tmp, alow, aupp)
    
    ! Store HPD interval results (two methods)
    phi_out(i1, i2, 4) = alow(1)  ! Lower bound (method 1)
    phi_out(i1, i2, 5) = aupp(1)  ! Upper bound (method 1)
    phi_out(i1, i2, 6) = alow(2)  ! Lower bound (method 2)
    phi_out(i1, i2, 7) = aupp(2)  ! Upper bound (method 2)
  end do
end do
!print *, 'pphi=', phi_out(:,:,1)

! Calculate posterior mean for phi_gamma parameters
pphi_gam(:, :) = 0.0d0  ! Initialize posterior mean matrix
do i1 = 1, ns      ! Loop over all obs
  do i2 = 1, J     ! Loop over outcomes (J dimensions)
    tmp = 0.0d0    ! Temporary sum for mean calculation
    
    do i3 = 1, nrep  ! Loop through all MCMC samples
      phitmp = 0.0d0  ! Temporary variable for linear combination
      
      ! Calculate linear combination: w * phi
      do i4 = 1, nphi  ! Loop over phi dimensions
        phitmp = phitmp + w(i1, i4) * seqphi(i2, i4, i3)
      end do
      
      ! Apply exponential transformation and multiply by Rgam
      seqphi_gam(i1, i2, i3) = exp(phitmp) * seqRgam(i1, i2, i3)
      
      ! Accumulate for posterior mean
      tmp = tmp + seqphi_gam(i1, i2, i3)
    end do
    
    ! Calculate posterior mean
    pphi_gam(i1, i2) = tmp / real(nrep, 8)
  end do
end do
!print *, 'pphi_gam=', pphi_gam

! Calculate posterior mean for covariance matrices
psig(:, :, :) = 0.0d0       ! Initialize posterior mean array
seqsig(:, :, :, :) = 0.0d0  ! Initialize sample storage array
do i1 = 1, ns  ! Loop over all obs
  ssig(:, :) = 0.0d0  ! Initialize sum matrix for this subject
  
  do i2 = 1, nrep  ! Loop through all MCMC samples
    sigtmp(:, :) = 0.0d0  ! Initialize temporary covariance matrix
    
    ! Set diagonal elements of covariance matrix
    do l1 = 1, J  ! Loop over dimensions
      sigtmp(l1, l1) = seqsigma(i1, l1, i2)
    end do
    
    ! Set off-diagonal elements of covariance matrix
    do l1 = 1, J - 1        ! Loop over rows
      do l2 = l1 + 1, J     ! Loop over columns (upper triangle)
        ! Calculate storage index: l1 + l2 + 1
        sigtmp(l1, l2) = seqsigma(i1, l1 + l2 + 1, i2)
        ! Mirror to lower triangle for symmetry
        sigtmp(l2, l1) = sigtmp(l1, l2)
      end do
    end do
    
    ! Store reconstructed covariance matrix
    seqsig(i2, i1, :, :) = sigtmp
    
    ! Accumulate for posterior mean
    ssig = ssig + sigtmp
  end do
  
  ! Calculate posterior mean covariance matrix
  psig(i1, :, :) = ssig / real(nrep, 8)
!  print *, 'i=',i1, ' psig=', psig(i1, :, :)
end do

! Calculate posterior mean for VRV matrices
pVRV(:, :, :) = 0.0d0  ! Initialize posterior mean array
do i1 = 1, ns  ! Loop over all obs
  ssig(:, :) = 0.0d0  ! Initialize sum matrix for this subject
  
  do i2 = 1, nrep  ! Loop through all MCMC samples
    ! Extract stored VRV matrix
    sigtmp(:, :) = seqVRV(i2, i1, :, :)
    
    ! Accumulate for posterior mean
    ssig = ssig + sigtmp
  end do
  
  ! Calculate posterior mean VRV matrix
  pVRV(i1, :, :) = ssig / real(nrep, 8)
!  print *, 'i=',i1, ' pVRV=', pVRV(i1, :, :)
end do
! Calculate posterior statistics for LDL-C
do i = 1, nsc
  ss1 = 0.0d0; ss2 = 0.0d0; ss3 = 0.0d0
  ss4 = 0.0d0; ss5 = 0.0d0; ss6 = 0.0d0
  
  do j1 = 1, nrep
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Extract components for LDL-C conditional distribution
    sig11 = sigtmp(1, 1)                        ! LDL-C variance
    sig12(1, 1) = sigtmp(1, 2); sig12(1, 2) = sigtmp(1, 3)  ! Covariance vector
    sig22(1, 1) = sigtmp(2, 2); sig22(2, 2) = sigtmp(3, 3)  ! HDL-C and TG variances
    sig22(1, 2) = sigtmp(2, 3); sig22(2, 1) = sig22(1, 2)   ! HDL-TG covariance
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(1, 1)                           ! LDL-C VRV component
    s12(1, 1) = ssig(1, 2); s12(1, 2) = ssig(1, 3) ! Covariance vector
    s22(1, 1) = ssig(2, 2); s22(2, 2) = ssig(3, 3) ! HDL-C and TG VRV components
    s22(1, 2) = ssig(2, 3); s22(2, 1) = s22(1, 2)  ! HDL-TG VRV covariance
    
    ! Invert sig22 matrix (covariance of HDL-C and TG)
    call inverse(sig22, sig22i, 2)       ! Matrix inversion         
    
    ! Calculate conditional variance components
    t7 = matmul(sig12, matmul(sig22i, transpose(sig12)))
    ldenom(j1, i) = sig11 - t7(1, 1)            ! Conditional variance
    
    ! Extract observed values
    y1 = y(1, i)                                ! Observed LDL-C
    y2(1, 1) = y(2, i); y2(2, 1) = y(3, i)     ! Observed HDL-C and TG
    
    ! Calculate conditional mean components
    lmu1(j1, i) = dot_product(xx(:, i), seqbeta(1, :, j1)) + seqphi_gam(i, 1, j1)
    mu2(1, 1) = dot_product(xx(:, i), seqbeta(2, :, j1)) + seqphi_gam(i, 2, j1)
    mu2(2, 1) = dot_product(xx(:, i), seqbeta(3, :, j1)) + seqphi_gam(i, 3, j1)
    
    ! Calculate intermediate products
    t3 = matmul(sig12, matmul(sig22i, transpose(s12)))
    t4 = matmul(sig12, matmul(sig22i, matmul(s22, matmul(sig22i, transpose(sig12)))))
    t5 = matmul(sig12, matmul(sig22i, y2))
    t6 = matmul(sig12, matmul(sig22i, mu2))
    
    ! Store intermediate products
    lprod1(j1, i) = t3(1, 1)
    lprod2(j1, i) = t4(1, 1)
    lprod3(j1, i) = t5(1, 1)
    lprod4(j1, i) = t6(1, 1)
    
    ! Accumulate for posterior means
    ss1 = ss1 + ldenom(j1, i)
    ss2 = ss2 + lprod1(j1, i)
    ss3 = ss3 + lprod2(j1, i)
    ss4 = ss4 + lprod3(j1, i)
    ss5 = ss5 + lprod4(j1, i)
    ss6 = ss6 + lmu1(j1, i)
  end do
  
  ! Calculate posterior means
  lpdenom(i) = ss1 / real(nrep, 8)
  lpprod1(i) = ss2 / real(nrep, 8)
  lpprod2(i) = ss3 / real(nrep, 8)
  lpprod3(i) = ss4 / real(nrep, 8)
  lpprod4(i) = ss5 / real(nrep, 8)
  lpmu1(i) = ss6 / real(nrep, 8)
end do

! Calculate posterior statistics for HDL-C
do i = 1, nsc
  ss1 = 0.0d0; ss2 = 0.0d0; ss3 = 0.0d0
  ss4 = 0.0d0; ss5 = 0.0d0; ss6 = 0.0d0
  
  do j1 = 1, nrep
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Extract components for HDL-C conditional distribution
    sig11 = sigtmp(2, 2)                        ! HDL-C variance
    sig12(1, 1) = sigtmp(2, 1); sig12(1, 2) = sigtmp(2, 3)  ! Covariance vector
    sig22(1, 1) = sigtmp(1, 1); sig22(2, 2) = sigtmp(3, 3)  ! LDL-C and TG variances
    sig22(1, 2) = sigtmp(1, 3); sig22(2, 1) = sigtmp(3, 1)   ! LDL-TG covariance
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(2, 2)                           ! HDL-C VRV component
    s12(1, 1) = ssig(2, 1); s12(1, 2) = ssig(2, 3) ! Covariance vector
    s22(1, 1) = ssig(1, 1); s22(2, 2) = ssig(3, 3) ! LDL-C and TG VRV components
    s22(1, 2) = ssig(1, 3); s22(2, 1) = ssig(3, 1) ! LDL-TG VRV covariance
    
    ! Invert sig22 matrix (covariance of LDL-C and TG)
    call inverse(sig22, sig22i, 2)               ! Matrix inversion
    
    ! Calculate conditional variance components
    t7 = matmul(sig12, matmul(sig22i, transpose(sig12)))
    hdenom(j1, i) = sig11 - t7(1, 1)            ! Conditional variance
    
    ! Extract observed values
    y1 = y(2, i)                                ! Observed HDL-C
    y2(1, 1) = y(1, i); y2(2, 1) = y(3, i)     ! Observed LDL-C and TG
    
    ! Calculate conditional mean components
    hmu1(j1, i) = dot_product(xx(:, i), seqbeta(2, :, j1)) + seqphi_gam(i, 2, j1)
    mu2(1, 1) = dot_product(xx(:, i), seqbeta(1, :, j1)) + seqphi_gam(i, 1, j1)
    mu2(2, 1) = dot_product(xx(:, i), seqbeta(3, :, j1)) + seqphi_gam(i, 3, j1)
    
    ! Calculate intermediate products
    t3 = matmul(sig12, matmul(sig22i, transpose(s12)))
    t4 = matmul(sig12, matmul(sig22i, matmul(s22, matmul(sig22i, transpose(sig12)))))
    t5 = matmul(sig12, matmul(sig22i, y2))
    t6 = matmul(sig12, matmul(sig22i, mu2))
    
    ! Store intermediate products
    hprod1(j1, i) = t3(1, 1)
    hprod2(j1, i) = t4(1, 1)
    hprod3(j1, i) = t5(1, 1)
    hprod4(j1, i) = t6(1, 1)
    
    ! Accumulate for posterior means
    ss1 = ss1 + hdenom(j1, i)
    ss2 = ss2 + hprod1(j1, i)
    ss3 = ss3 + hprod2(j1, i)
    ss4 = ss4 + hprod3(j1, i)
    ss5 = ss5 + hprod4(j1, i)
    ss6 = ss6 + hmu1(j1, i)
  end do
  
  ! Calculate posterior means
  hpdenom(i) = ss1 / real(nrep, 8)
  hpprod1(i) = ss2 / real(nrep, 8)
  hpprod2(i) = ss3 / real(nrep, 8)
  hpprod3(i) = ss4 / real(nrep, 8)
  hpprod4(i) = ss5 / real(nrep, 8)
  hpmu1(i) = ss6 / real(nrep, 8)
end do

! Calculate posterior statistics for Triglycerides (TG)
do i = 1, nsc
  ss1 = 0.0d0; ss2 = 0.0d0; ss3 = 0.0d0
  ss4 = 0.0d0; ss5 = 0.0d0; ss6 = 0.0d0
  
  do j1 = 1, nrep
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Extract components for TG conditional distribution
    sig11 = sigtmp(3, 3)                        ! TG variance
    sig12(1, 1) = sigtmp(1, 3); sig12(1, 2) = sigtmp(2, 3)  ! Covariance vector
    sig22 = sigtmp(1:2, 1:2)                    ! LDL-C and HDL-C covariance matrix
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(3, 3)                           ! TG VRV component
    s12(1, 1) = ssig(1, 3); s12(1, 2) = ssig(2, 3) ! Covariance vector
    s22 = ssig(1:2, 1:2)
    
    ! Invert sig22 matrix (covariance of LDL-C and HDL-C)
    call inverse(sig22, sig22i, 2)               ! Matrix inversion
    
    ! Calculate conditional variance components
    t7 = matmul(sig12, matmul(sig22i, transpose(sig12)))
    tdenom(j1, i) = sig11 - t7(1, 1)            ! Conditional variance
    
    ! Extract observed values
    y1 = y(3, i)                                ! Observed TG
    y2(1, 1) = y(1, i); y2(2, 1) = y(2, i)     ! Observed LDL-C and HDL-C
    
    ! Calculate conditional mean components
    tmu1(j1, i) = dot_product(xx(:, i), seqbeta(3, :, j1)) + seqphi_gam(i, 3, j1)
    mu2(1, 1) = dot_product(xx(:, i), seqbeta(1, :, j1)) + seqphi_gam(i, 1, j1)
    mu2(2, 1) = dot_product(xx(:, i), seqbeta(2, :, j1)) + seqphi_gam(i, 2, j1)
    
    ! Calculate intermediate products
    t3 = matmul(sig12, matmul(sig22i, transpose(s12)))
    t4 = matmul(sig12, matmul(sig22i, matmul(s22, matmul(sig22i, transpose(sig12)))))
    t5 = matmul(sig12, matmul(sig22i, y2))
    t6 = matmul(sig12, matmul(sig22i, mu2))
    
    ! Store intermediate products
    tprod1(j1, i) = t3(1, 1)
    tprod2(j1, i) = t4(1, 1)
    tprod3(j1, i) = t5(1, 1)
    tprod4(j1, i) = t6(1, 1)
    
    ! Accumulate for posterior means
    ss1 = ss1 + tdenom(j1, i)
    ss2 = ss2 + tprod1(j1, i)
    ss3 = ss3 + tprod2(j1, i)
    ss4 = ss4 + tprod3(j1, i)
    ss5 = ss5 + tprod4(j1, i)
    ss6 = ss6 + tmu1(j1, i)
  end do
  
  ! Calculate posterior means
  tpdenom(i) = ss1 / real(nrep, 8)
  tpprod1(i) = ss2 / real(nrep, 8)
  tpprod2(i) = ss3 / real(nrep, 8)
  tpprod3(i) = ss4 / real(nrep, 8)
  tpprod4(i) = ss5 / real(nrep, 8)
  tpmu1(i) = ss6 / real(nrep, 8)
end do

! Calculate DIC (Deviance Information Criterion) components
do icount = 1, nsc
    sumDbar(icount) = 0.0d0
    t1 = ny(icount) - 1.0d0  ! Degrees of freedom
    
    ! Get posterior mean covariance matrix
    sigtmp = psig(icount, :, :)
!    print *, 'i=',icount, ' psig=', sigtmp
    
    ! Invert covariance matrix
    call inverse(sigtmp, sigitmp, J)
    
    ! Calculate determinant of covariance matrix
    det = FindDet(tol, J, sigtmp)
!    print *, 'i=', icount, 'det=', det
    
    ! Calculate residuals
    do i1 = 1, J 
      tv(i1) = y(i1, icount) - &
               dot_product(xx(:, icount), pbeta(:, i1)) - &
               pphi_gam(icount, i1)
    end do
    
    ! Calculate log-likelihood components for yDbar
    tmp = -real(J, 8)/2.0d0 * log(2.0d0*epi) - &
          0.5d0 * log(det) + &
          1.5d0 * log(ny(icount)) - &
          0.5d0 * ny(icount) * dot_product(matmul(sigitmp, tv), tv)
    
    yDbar(icount) = -2.0d0 * tmp  ! Deviance for yDbar
    
    ! Calculate log-likelihood components for y1Dbar
    tmp = -real(J*ny(icount), 8)/2.0d0 * log(2.0d0*epi) - &
          0.5d0 * log(det) - &
          0.5d0 * ny(icount) * dot_product(matmul(sigitmp, tv), tv)
    
    y1Dbar(icount) = -2.0d0 * tmp  ! Deviance for y1Dbar
    sumDbar(icount) = sumDbar(icount) -2.0d0 * tmp
    
    ! for S part
      sigtmp = pVRV(icount, :, :) * t1  ! Scale VRV matrix
      det1 = FindDet(tol, J, sigtmp)       ! Determinant of scaled VRV
      ssig = matmul(sigitmp, sigtmp)    ! Matrix product
      ss1 = ssig(1,1) + ssig(2,2) + ssig(3,3)  ! Trace of product
      
      ! Calculate gamma function terms
      lg1 = alngam(t1/2.0d0, ifault)
      lg2 = alngam((t1-1.0d0)/2.0d0, ifault)
      lg3 = alngam((t1-2.0d0)/2.0d0, ifault)
      
      ! Calculate log-likelihood components for sDbar
      tmp = real(J*(J+1), 8)/2.0d0 * log(t1) - &
            t1*real(J, 8)/2.0d0 * log(2.0d0) - &
            real(J*(J-1), 8)/4.0d0 * log(epi) - &
            lg1 - lg2 - lg3 - &
            (t1/2.0d0) * log(det) + &
            real(ny(icount)-J-2, 8)/2.0d0 * log(det1) - &
            0.5d0 * ss1
      
      sDbar(icount) = -2.0d0 * tmp  ! Deviance for sDbar
      
      ! Simplified calculation for s1Dbar
      tmp = - (t1/2.0d0) * log(det) - 0.5d0 * ss1
      s1Dbar(icount) = -2.0d0 * tmp  ! Deviance for s1Dbar
      sumDbar(icount) = sumDbar(icount) -2.0d0 * tmp
end do 


! Calculate posterior mean deviances and DIC components
do icount = 1, nsc
    t1 = ny(icount) - 1.0d0  ! Degrees of freedom

    ! Initialize accumulators
    ss1 = 0.0d0; ss2 = 0.0d0; ss5 = 0.0d0; ss6 = 0.0d0; ss7 = 0.0d0
    
    do j1 = 1, nrep  ! Loop over MCMC samples
      ! Get covariance matrix for this sample
      sigtmp = seqsig(j1, icount, :, :)
      
      ! Invert covariance matrix
      call inverse(sigtmp, sigitmp, J)
      
      ! Calculate determinant
      det = FindDet(tol, J, sigtmp)
      
      ! Calculate residuals
      do i1 = 1, J 
        ss3 = 0.0d0
        do i2 = 1, nbbeta
          ss3 = ss3 + xx(i2, icount) * seqbeta(i1, i2, j1)
        end do
        tv(i1) = y(i1, icount) - ss3 - seqphi_gam(icount, i1, j1)
      end do

      ! Calculate log-likelihood for ybarD
      tmp = -real(J, 8)/2.0d0 * log(2.0d0*epi) - &
            0.5d0 * log(det) + &
            1.5d0 * log(ny(icount)) - &
            0.5d0 * ny(icount) * dot_product(matmul(sigitmp, tv), tv)
      
      ss1 = ss1 - 2.0d0 * tmp  ! Accumulate deviance
      
      ! Calculate log-likelihood for y1barD
      tmp = -real(J*ny(icount), 8)/2.0d0 * log(2.0d0*epi) - &
            0.5d0 * log(det) - &
            0.5d0 * ny(icount) * dot_product(matmul(sigitmp, tv), tv)
      
      ss6 = ss6 - 2.0d0 * tmp  ! Accumulate deviance
      ss7 = ss7 - 2.0d0 * tmp
      
      ! Additional calculations for obs with sufficient data

        sigtmp = seqVRV(j1, icount, :, :) * t1  ! Scale VRV matrix
        det1 = FindDet(tol, J, sigtmp)             ! Determinant
        ssig = matmul(sigitmp, sigtmp)          ! Matrix product
        ss4 = ssig(1,1) + ssig(2,2) + ssig(3,3) ! Trace of product
        
        ! Calculate gamma function terms
        lg1 = alngam(t1/2.0d0, ifault)
        lg2 = alngam((t1-1.0d0)/2.0d0, ifault)
        lg3 = alngam((t1-2.0d0)/2.0d0, ifault)
        
        ! Calculate log-likelihood for sbarD
        tmp = real(J*(J+1), 8)/2.0d0 * log(t1) - &
              t1*real(J, 8)/2.0d0 * log(2.0d0) - &
              real(J*(J-1), 8)/4.0d0 * log(epi) - &
              lg1 - lg2 - lg3 - &
              (t1/2.0d0) * log(det) + &
              real(ny(icount)-J-2, 8)/2.0d0 * log(det1) - &
              0.5d0 * ss4
        
        ss2 = ss2 - 2.0d0 * tmp  ! Accumulate deviance
        
        ! Simplified calculation for s1barD
        tmp = - (t1/2.0d0) * log(det) - 0.5d0 * ss4
        ss5 = ss5 - 2.0d0 * tmp  ! Accumulate deviance
        ss7 = ss7 - 2.0d0 * tmp
      
    end do
    
    ! Calculate posterior mean deviances
    ybarD(icount) = ss1 / real(nrep, 8)
    y1barD(icount) = ss6 / real(nrep, 8)
    sumbarD(icount) = ss7 / real(nrep, 8)
    
    ! Calculate effective number of parameters (pD)
    ypD(icount) = ybarD(icount) - yDbar(icount)
    y1pD(icount) = y1barD(icount) - y1Dbar(icount)
    
    
    ! Calculate DIC = Dbar + pD
    yDIC(icount) = yDbar(icount) + 2.0d0 * ypD(icount)
    y1DIC(icount) = y1Dbar(icount) + 2.0d0 * y1pD(icount)
    
    ! Additional DIC calculations for obs with sufficient data
      sbarD(icount) = ss2 / real(nrep, 8)
      s1barD(icount) = ss5 / real(nrep, 8)
      
      spD(icount) = sbarD(icount) - sDbar(icount)
      s1pD(icount) = s1barD(icount) - s1Dbar(icount)
      sumpD(icount) = sumbarD(icount) - sumDbar(icount)
      
      sDIC(icount) = sDbar(icount) + 2.0d0 * spD(icount)
      s1DIC(icount) = s1Dbar(icount) + 2.0d0 * s1pD(icount)
      sumDIC(icount) = sumDbar(icount) + 2.0d0 * sumpD(icount)
end do

! Calculate conditional DIC for LDL-C
! -----------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 1, 1)  ! Posterior mean VRV component for LDL-C
  y1 = y(1, i)         ! Observed LDL-C value
  
  ! Calculate log-likelihood components
  tmp = -t1/2.0d0 * log(2.0d0*epi) - &
        t1/2.0d0 * log(lpdenom(i)) - &
        0.5d0 * t2 * (s11 - 2.0d0*lpprod1(i) + lpprod2(i)) / lpdenom(i) - &
        0.5d0 * t1 / lpdenom(i) * &
        ((y1 - lpprod3(i) - lpmu1(i) + lpprod4(i))**2)
  
  ! Store deviance at posterior mean
  lDbar(i) = -2.0d0 * tmp
end do

! Calculate conditional DIC for HDL-C
! -----------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 2, 2)  ! Posterior mean VRV component for HDL-C
  y1 = y(2, i)         ! Observed HDL-C value
  
  ! Calculate log-likelihood components
  tmp = -t1/2.0d0 * log(2.0d0*epi) - &
        t1/2.0d0 * log(hpdenom(i)) - &
        0.5d0 * t2 * (s11 - 2.0d0*hpprod1(i) + hpprod2(i)) / hpdenom(i) - &
        0.5d0 * t1 / hpdenom(i) * &
        ((y1 - hpprod3(i) - hpmu1(i) + hpprod4(i))**2)
  
  ! Store deviance at posterior mean
  hDbar(i) = -2.0d0 * tmp
end do

! Calculate conditional DIC for Triglycerides (TG)
! ----------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 3, 3)  ! Posterior mean VRV component for TG
  y1 = y(3, i)         ! Observed TG value
  
  ! Calculate log-likelihood components
  tmp = -t1/2.0d0 * log(2.0d0*epi) - &
        t1/2.0d0 * log(tpdenom(i)) - &
        0.5d0 * t2 * (s11 - 2.0d0*tpprod1(i) + tpprod2(i)) / tpdenom(i) - &
        0.5d0 * t1 / tpdenom(i) * &
        ((y1 - tpprod3(i) - tpmu1(i) + tpprod4(i))**2)
  
  ! Store deviance at posterior mean
  tDbar(i) = -2.0d0 * tmp
end do

! Calculate posterior mean deviance for LDL-C
! -------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 1, 1)  ! Posterior mean VRV component for LDL-C
  y1 = y(1, i)         ! Observed LDL-C value
  ss1 = 0.0d0          ! Accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Calculate log-likelihood for this sample
    tmp = -t1/2.0d0 * log(2.0d0*epi) - &
          t1/2.0d0 * log(ldenom(j1, i)) - &
          0.5d0 * t2 * (s11 - 2.0d0*lprod1(j1, i) + lprod2(j1, i)) / ldenom(j1, i) - &
          0.5d0 * t1 / ldenom(j1, i) * &
          ((y1 - lprod3(j1, i) - lmu1(j1, i) + lprod4(j1, i))**2)
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  lbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  lpD(i) = lbarD(i) - lDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  lDIC(i) = lDbar(i) + 2.0d0 * lpD(i)
end do

! Calculate posterior mean deviance for HDL-C
! -------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 2, 2)  ! Posterior mean VRV component for HDL-C
  y1 = y(2, i)         ! Observed HDL-C value
  ss1 = 0.0d0          ! Accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Calculate log-likelihood for this sample
    tmp = -t1/2.0d0 * log(2.0d0*epi) - &
          t1/2.0d0 * log(hdenom(j1, i)) - &
          0.5d0 * t2 * (s11 - 2.0d0*hprod1(j1, i) + hprod2(j1, i)) / hdenom(j1, i) - &
          0.5d0 * t1 / hdenom(j1, i) * &
          ((y1 - hprod3(j1, i) - hmu1(j1, i) + hprod4(j1, i))**2)
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  hbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  hpD(i) = hbarD(i) - hDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  hDIC(i) = hDbar(i) + 2.0d0 * hpD(i)
end do

! Calculate posterior mean deviance for Triglycerides (TG)
! ------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  s11 = pVRV(i, 3, 3)  ! Posterior mean VRV component for TG
  y1 = y(3, i)         ! Observed TG value
  ss1 = 0.0d0          ! Accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Calculate log-likelihood for this sample
    tmp = -t1/2.0d0 * log(2.0d0*epi) - &
          t1/2.0d0 * log(tdenom(j1, i)) - &
          0.5d0 * t2 * (s11 - 2.0d0*tprod1(j1, i) + tprod2(j1, i)) / tdenom(j1, i) - &
          0.5d0 * t1 / tdenom(j1, i) * &
          ((y1 - tprod3(j1, i) - tmu1(j1, i) + tprod4(j1, i))**2)
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  tbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  tpD(i) = tbarD(i) - tDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  tDIC(i) = tDbar(i) + 2.0d0 * tpD(i)
end do
! Calculate posterior statistics for LDL-C and HDL-C conditional on TG
! ----------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  ! Initialize accumulation matrices
  smat1 = 0.0d0; smat2 = 0.0d0; smat3 = 0.0d0
  smat4 = 0.0d0; smat5 = 0.0d0; smat6 = 0.0d0
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Decompose covariance matrix into blocks
    sig11 = sigtmp(3, 3)  ! TG variance
    sig12(1, 1) = sigtmp(1, 3)  ! Cov(LDL-C, TG)
    sig12(1, 2) = sigtmp(2, 3)  ! Cov(HDL-C, TG)
    sig22(1, 1) = sigtmp(1, 1)  ! LDL-C variance
    sig22(2, 2) = sigtmp(2, 2)  ! HDL-C variance
    sig22(1, 2) = sigtmp(1, 2)  ! Cov(LDL-C, HDL-C)
    sig22(2, 1) = sig22(1, 2)   ! Ensure symmetry
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(3, 3)           ! TG VRV component
    s12(1, 1) = ssig(1, 3)    ! VRV Cov(LDL-C, TG)
    s12(1, 2) = ssig(2, 3)    ! VRV Cov(HDL-C, TG)
    s22(1, 1) = ssig(1, 1)    ! LDL-C VRV component
    s22(2, 2) = ssig(2, 2)    ! HDL-C VRV component
    s22(1, 2) = ssig(1, 2)    ! VRV Cov(LDL-C, HDL-C)
    s22(2, 1) = s22(1, 2)     ! Ensure symmetry
    
    ! Calculate conditional covariance matrix
    pa1 = matmul(transpose(sig12), sig12) / sig11
    lhdenom(j1, i, :, :) = sig22 - pa1
    
    ! Extract observed values
    y1 = y(3, i)          ! TG
    y2(1, 1) = y(1, i)    ! LDL-C
    y2(2, 1) = y(2, i)    ! HDL-C
    
    ! Calculate conditional means
    lhmu(j1, i, 1, 1) = dot_product(xx(:, i), seqbeta(1, :, j1)) + &
                         seqphi_gam(i, 1, j1)  ! LDL-C mean
    lhmu(j1, i, 1, 2) = dot_product(xx(:, i), seqbeta(2, :, j1)) + &
                         seqphi_gam(i, 2, j1)  ! HDL-C mean
    mu1 = dot_product(xx(:, i), seqbeta(3, :, j1)) + &
          seqphi_gam(i, 3, j1)  ! TG mean
    
    ! Calculate intermediate products
    pa2 = matmul(transpose(sig12), s12) / sig11
    pa3 = matmul(transpose(sig12), sig12) * s11 / (sig11 * sig11)
    pa4 = sig12 * y1 / sig11
    pa5 = sig12 * mu1 / sig11
    
    ! Store intermediate products
    lhprod1(j1, i, :, :) = pa2
    lhprod2(j1, i, :, :) = pa3
    lhprod3(j1, i, :, :) = pa4
    lhprod4(j1, i, :, :) = pa5
    
    ! Accumulate for posterior means
    smat1 = smat1 + lhdenom(j1, i, :, :)
    smat2 = smat2 + lhprod1(j1, i, :, :)
    smat3 = smat3 + lhprod2(j1, i, :, :)
    smat4 = smat4 + lhprod3(j1, i, :, :)
    smat5 = smat5 + lhprod4(j1, i, :, :)
    smat6 = smat6 + lhmu(j1, i, :, :)
  end do
  
  ! Calculate posterior means
  lhpdenom(i, :, :) = smat1 / real(nrep, 8)  ! Mean conditional covariance
  lhpprod1(i, :, :) = smat2 / real(nrep, 8)   ! Mean product 1
  lhpprod2(i, :, :) = smat3 / real(nrep, 8)   ! Mean product 2
  lhpprod3(i, :, :) = smat4 / real(nrep, 8)   ! Mean product 3
  lhpprod4(i, :, :) = smat5 / real(nrep, 8)   ! Mean product 4
  lhpmu(i, :, :) = smat6 / real(nrep, 8)      ! Mean conditional mean
end do

! Calculate Dbar (deviance at posterior mean) for LDL-C and HDL-C | TG
! ------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract posterior mean VRV matrix for LDL-C and HDL-C
  s22 = pVRV(i, 1:2, 1:2)
  
  ! Extract observed values
  y2(1, 1) = y(1, i)  ! LDL-C
  y2(2, 1) = y(2, i)  ! HDL-C
  
  ! Get posterior mean conditional covariance matrix
  smat2 = lhpdenom(i, :, :)
  
  ! Calculate determinant of conditional covariance matrix
  det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
  
  ! Invert conditional covariance matrix
  call inverse(smat2, sig22i, 2)
  
  ! Calculate trace component
  smat1 = s22 + lhpprod2(i, :, :) - 2.0d0 * lhpprod1(i, :, :)
  smat3 = matmul(sig22i, smat1)
  a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
  
  ! Calculate quadratic form component
  smat4 = transpose(y2) - lhpprod3(i, :, :) - lhpmu(i, :, :) + lhpprod4(i, :, :)
  t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
  
  ! Calculate log-likelihood
  tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
        t1 / 2.0d0 * log(det) - &    ! Determinant term
        0.5d0 * t2 * a - &            ! Trace term
        0.5d0 * t1 * t7(1, 1)        ! Quadratic form term
  
  ! Store deviance at posterior mean
  lhDbar(i) = -2.0d0 * tmp
end do

! Calculate posterior mean deviance (barD) for LDL-C and HDL-C | TG
! ----------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract observed values
  y2(1, 1) = y(1, i)  ! LDL-C
  y2(2, 1) = y(2, i)  ! HDL-C
  
  ss1 = 0.0d0  ! Initialize accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract VRV matrix for this sample
    s22 = seqVRV(j1, i, 1:2, 1:2)
    
    ! Get conditional covariance matrix for this sample
    smat2 = lhdenom(j1, i, :, :)
    
    ! Calculate determinant of conditional covariance matrix
    det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
    
    ! Invert conditional covariance matrix
    call inverse(smat2, sig22i, 2)
    
    ! Calculate trace component
    smat1 = s22 + lhprod2(j1, i, :, :) - 2.0d0 * lhprod1(j1, i, :, :)
    smat3 = matmul(sig22i, smat1)
    a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
    
    ! Calculate quadratic form component
    smat4 = transpose(y2) - lhprod3(j1, i, :, :) - lhmu(j1, i, :, :) + lhprod4(j1, i, :, :)
    t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
    
    ! Calculate log-likelihood for this sample
    tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
          t1 / 2.0d0 * log(det) - &    ! Determinant term
          0.5d0 * t2 * a - &           ! Trace term
          0.5d0 * t1 * t7(1, 1)       ! Quadratic form term
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  lhbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  lhpD(i) = lhbarD(i) - lhDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  lhDIC(i) = lhDbar(i) + 2.0d0 * lhpD(i)
end do
! Calculate posterior statistics for HDL-C and TG conditional on LDL-C
! ----------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  ! Initialize accumulation matrices
  smat1 = 0.0d0; smat2 = 0.0d0; smat3 = 0.0d0
  smat4 = 0.0d0; smat5 = 0.0d0; smat6 = 0.0d0
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Decompose covariance matrix into blocks
    sig11 = sigtmp(1, 1)  ! LDL-C variance
    sig12(1, 1) = sigtmp(1, 2)  ! Cov(LDL-C, HDL-C)
    sig12(1, 2) = sigtmp(1, 3)  ! Cov(LDL-C, TG)
    sig22(1, 1) = sigtmp(2, 2)  ! HDL-C variance
    sig22(2, 2) = sigtmp(3, 3)  ! TG variance
    sig22(1, 2) = sigtmp(2, 3)  ! Cov(HDL-C, TG)
    sig22(2, 1) = sig22(1, 2)   ! Ensure symmetry
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(1, 1)           ! LDL-C VRV component
    s12(1, 1) = ssig(1, 2)    ! VRV Cov(LDL-C, HDL-C)
    s12(1, 2) = ssig(1, 3)    ! VRV Cov(LDL-C, TG)
    s22(1, 1) = ssig(2, 2)    ! HDL-C VRV component
    s22(2, 2) = ssig(3, 3)    ! TG VRV component
    s22(1, 2) = ssig(2, 3)    ! VRV Cov(HDL-C, TG)
    s22(2, 1) = s22(1, 2)     ! Ensure symmetry
    
    ! Calculate conditional covariance matrix
    pa1 = matmul(transpose(sig12), sig12) / sig11
    htdenom(j1, i, :, :) = sig22 - pa1
    
    ! Extract observed values
    y1 = y(1, i)          ! LDL-C
    y2(1, 1) = y(2, i)    ! HDL-C
    y2(2, 1) = y(3, i)    ! TG
    
    ! Calculate conditional means
    htmu(j1, i, 1, 1) = dot_product(xx(:, i), seqbeta(2, :, j1)) + &
                         seqphi_gam(i, 2, j1)  ! HDL-C mean
    htmu(j1, i, 1, 2) = dot_product(xx(:, i), seqbeta(3, :, j1)) + &
                         seqphi_gam(i, 3, j1)  ! TG mean
    mu1 = dot_product(xx(:, i), seqbeta(1, :, j1)) + &
          seqphi_gam(i, 1, j1)  ! LDL-C mean
    
    ! Calculate intermediate products
    pa2 = matmul(transpose(sig12), s12) / sig11
    pa3 = matmul(transpose(sig12), sig12) * s11 / (sig11 * sig11)
    pa4 = sig12 * y1 / sig11
    pa5 = sig12 * mu1 / sig11
    
    ! Store intermediate products
    htprod1(j1, i, :, :) = pa2
    htprod2(j1, i, :, :) = pa3
    htprod3(j1, i, :, :) = pa4
    htprod4(j1, i, :, :) = pa5
    
    ! Accumulate for posterior means
    smat1 = smat1 + htdenom(j1, i, :, :)
    smat2 = smat2 + htprod1(j1, i, :, :)
    smat3 = smat3 + htprod2(j1, i, :, :)
    smat4 = smat4 + htprod3(j1, i, :, :)
    smat5 = smat5 + htprod4(j1, i, :, :)
    smat6 = smat6 + htmu(j1, i, :, :)
  end do
  
  ! Calculate posterior means
  htpdenom(i, :, :) = smat1 / real(nrep, 8)  ! Mean conditional covariance
  htpprod1(i, :, :) = smat2 / real(nrep, 8)   ! Mean product 1
  htpprod2(i, :, :) = smat3 / real(nrep, 8)   ! Mean product 2
  htpprod3(i, :, :) = smat4 / real(nrep, 8)   ! Mean product 3
  htpprod4(i, :, :) = smat5 / real(nrep, 8)   ! Mean product 4
  htpmu(i, :, :) = smat6 / real(nrep, 8)      ! Mean conditional mean
end do

! Calculate Dbar (deviance at posterior mean) for HDL-C and TG | LDL-C
! ------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract posterior mean VRV matrix for HDL-C and TG
  s22 = pVRV(i, 2:3, 2:3)
  
  ! Extract observed values
  y2(1, 1) = y(2, i)  ! HDL-C
  y2(2, 1) = y(3, i)  ! TG
  
  ! Get posterior mean conditional covariance matrix
  smat2 = htpdenom(i, :, :)
  
  ! Calculate determinant of conditional covariance matrix
  det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
  
  ! Invert conditional covariance matrix
  call inverse(smat2, sig22i, 2)
  
  ! Calculate trace component
  smat1 = s22 + htpprod2(i, :, :) - 2.0d0 * htpprod1(i, :, :)
  smat3 = matmul(sig22i, smat1)
  a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
  
  ! Calculate quadratic form component
  smat4 = transpose(y2) - htpprod3(i, :, :) - htpmu(i, :, :) + htpprod4(i, :, :)
  t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
  
  ! Calculate log-likelihood
  tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
        t1 / 2.0d0 * log(det) - &    ! Determinant term
        0.5d0 * t2 * a - &           ! Trace term
        0.5d0 * t1 * t7(1, 1)        ! Quadratic form term
  
  ! Store deviance at posterior mean
  htDbar(i) = -2.0d0 * tmp
end do

! Calculate posterior mean deviance (barD) for HDL-C and TG | LDL-C
! ----------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract observed values
  y2(1, 1) = y(2, i)  ! HDL-C
  y2(2, 1) = y(3, i)  ! TG
  
  ss1 = 0.0d0  ! Initialize accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract VRV matrix for this sample
    s22 = seqVRV(j1, i, 2:3, 2:3)
    
    ! Get conditional covariance matrix for this sample
    smat2 = htdenom(j1, i, :, :)
    
    ! Calculate determinant of conditional covariance matrix
    det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
    
    ! Invert conditional covariance matrix
    call inverse(smat2, sig22i, 2)
    
    ! Calculate trace component
    smat1 = s22 + htprod2(j1, i, :, :) - 2.0d0 * htprod1(j1, i, :, :)
    smat3 = matmul(sig22i, smat1)
    a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
    
    ! Calculate quadratic form component
    smat4 = transpose(y2) - htprod3(j1, i, :, :) - htmu(j1, i, :, :) + htprod4(j1, i, :, :)
    t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
    
    ! Calculate log-likelihood for this sample
    tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
          t1 / 2.0d0 * log(det) - &   ! Determinant term
          0.5d0 * t2 * a - &           ! Trace term
          0.5d0 * t1 * t7(1, 1)       ! Quadratic form term
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  htbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  htpD(i) = htbarD(i) - htDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  htDIC(i) = htDbar(i) + 2.0d0 * htpD(i)
end do
! Calculate posterior statistics for LDL-C and TG conditional on HDL-C
! ----------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  ! Initialize accumulation matrices
  smat1 = 0.0d0; smat2 = 0.0d0; smat3 = 0.0d0
  smat4 = 0.0d0; smat5 = 0.0d0; smat6 = 0.0d0
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract covariance matrix for this sample
    sigtmp = seqsig(j1, i, :, :)
    
    ! Decompose covariance matrix into blocks
    sig11 = sigtmp(2, 2)  ! HDL-C variance
    sig12(1, 1) = sigtmp(1, 2)  ! Cov(LDL-C, HDL-C)
    sig12(1, 2) = sigtmp(2, 3)  ! Cov(HDL-C, TG)
    sig22(1, 1) = sigtmp(1, 1)  ! LDL-C variance
    sig22(2, 2) = sigtmp(3, 3)  ! TG variance
    sig22(1, 2) = sigtmp(1, 3)  ! Cov(LDL-C, TG)
    sig22(2, 1) = sig22(1, 2)   ! Ensure symmetry
    
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s11 = ssig(2, 2)           ! HDL-C VRV component
    s12(1, 1) = ssig(1, 2)    ! VRV Cov(LDL-C, HDL-C)
    s12(1, 2) = ssig(2, 3)    ! VRV Cov(HDL-C, TG)
    s22(1, 1) = ssig(1, 1)    ! LDL-C VRV component
    s22(2, 2) = ssig(3, 3)    ! TG VRV component
    s22(1, 2) = ssig(1, 3)    ! VRV Cov(LDL-C, TG)
    s22(2, 1) = s22(1, 2)     ! Ensure symmetry
    
    ! Calculate conditional covariance matrix
    pa1 = matmul(transpose(sig12), sig12) / sig11
    ltdenom(j1, i, :, :) = sig22 - pa1
    
    ! Extract observed values
    y1 = y(2, i)          ! HDL-C
    y2(1, 1) = y(1, i)    ! LDL-C
    y2(2, 1) = y(3, i)    ! TG
    
    ! Calculate conditional means
    ltmu(j1, i, 1, 1) = dot_product(xx(:, i), seqbeta(1, :, j1)) + &
                         seqphi_gam(i, 1, j1)  ! LDL-C mean
    ltmu(j1, i, 1, 2) = dot_product(xx(:, i), seqbeta(3, :, j1)) + &
                         seqphi_gam(i, 3, j1)  ! TG mean
    mu1 = dot_product(xx(:, i), seqbeta(2, :, j1)) + &
          seqphi_gam(i, 2, j1)  ! HDL-C mean
    
    ! Calculate intermediate products
    pa2 = matmul(transpose(sig12), s12) / sig11
    pa3 = matmul(transpose(sig12), sig12) * s11 / (sig11 * sig11)
    pa4 = sig12 * y1 / sig11
    pa5 = sig12 * mu1 / sig11
    
    ! Store intermediate products
    ltprod1(j1, i, :, :) = pa2
    ltprod2(j1, i, :, :) = pa3
    ltprod3(j1, i, :, :) = pa4
    ltprod4(j1, i, :, :) = pa5
    
    ! Accumulate for posterior means
    smat1 = smat1 + ltdenom(j1, i, :, :)
    smat2 = smat2 + ltprod1(j1, i, :, :)
    smat3 = smat3 + ltprod2(j1, i, :, :)
    smat4 = smat4 + ltprod3(j1, i, :, :)
    smat5 = smat5 + ltprod4(j1, i, :, :)
    smat6 = smat6 + ltmu(j1, i, :, :)
  end do
  
  ! Calculate posterior means
  ltpdenom(i, :, :) = smat1 / real(nrep, 8)  ! Mean conditional covariance
  ltpprod1(i, :, :) = smat2 / real(nrep, 8)   ! Mean product 1
  ltpprod2(i, :, :) = smat3 / real(nrep, 8)   ! Mean product 2
  ltpprod3(i, :, :) = smat4 / real(nrep, 8)   ! Mean product 3
  ltpprod4(i, :, :) = smat5 / real(nrep, 8)   ! Mean product 4
  ltpmu(i, :, :) = smat6 / real(nrep, 8)      ! Mean conditional mean
end do

! Calculate Dbar (deviance at posterior mean) for LDL-C and TG | HDL-C
! ------------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract posterior mean VRV matrix for LDL-C and TG
  ssig = pVRV(i, :, :)
  s22(1, 1) = ssig(1, 1)  ! LDL-C VRV component
  s22(2, 2) = ssig(3, 3)  ! TG VRV component
  s22(1, 2) = ssig(1, 3)  ! VRV Cov(LDL-C, TG)
  s22(2, 1) = s22(1, 2)   ! Ensure symmetry
  
  ! Extract observed values
  y2(1, 1) = y(1, i)  ! LDL-C
  y2(2, 1) = y(3, i)  ! TG
  
  ! Get posterior mean conditional covariance matrix
  smat2 = ltpdenom(i, :, :)
  
  ! Calculate determinant of conditional covariance matrix
  det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
  
  ! Invert conditional covariance matrix
  call inverse(smat2, sig22i, 2)
  
  ! Calculate trace component
  smat1 = s22 + ltpprod2(i, :, :) - 2.0d0 * ltpprod1(i, :, :)
  smat3 = matmul(sig22i, smat1)
  a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
  
  ! Calculate quadratic form component
  smat4 = transpose(y2) - ltpprod3(i, :, :) - ltpmu(i, :, :) + ltpprod4(i, :, :)
  t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
  
  ! Calculate log-likelihood
  tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
        t1 / 2.0d0 * log(det) - &    ! Determinant term
        0.5d0 * t2 * a - &           ! Trace term
        0.5d0 * t1 * t7(1, 1)       ! Quadratic form term
  
  ! Store deviance at posterior mean
  ltDbar(i) = -2.0d0 * tmp
end do

! Calculate posterior mean deviance (barD) for LDL-C and TG | HDL-C
! ----------------------------------------------------------------
do i = 1, nsc  ! Loop over obs
  t1 = ny(i)      ! Number of measurements
  t2 = ny(i) - 1.0d0  ! Degrees of freedom
  
  ! Extract observed values
  y2(1, 1) = y(1, i)  ! LDL-C
  y2(2, 1) = y(3, i)  ! TG
  
  ss1 = 0.0d0  ! Initialize accumulator for posterior mean deviance
  
  do j1 = 1, nrep  ! Loop over MCMC samples
    ! Extract VRV matrix for this sample
    ssig = seqVRV(j1, i, :, :)
    s22(1, 1) = ssig(1, 1)  ! LDL-C VRV component
    s22(2, 2) = ssig(3, 3)  ! TG VRV component
    s22(1, 2) = ssig(1, 3)  ! VRV Cov(LDL-C, TG)
    s22(2, 1) = s22(1, 2)   ! Ensure symmetry
    
    ! Get conditional covariance matrix for this sample
    smat2 = ltdenom(j1, i, :, :)
    
    ! Calculate determinant of conditional covariance matrix
    det = smat2(1, 1) * smat2(2, 2) - smat2(1, 2)**2.0d0
    
    ! Invert conditional covariance matrix
    call inverse(smat2, sig22i, 2)
    
    ! Calculate trace component
    smat1 = s22 + ltprod2(j1, i, :, :) - 2.0d0 * ltprod1(j1, i, :, :)
    smat3 = matmul(sig22i, smat1)
    a = smat3(1, 1) + smat3(2, 2)  ! Trace of matrix product
    
    ! Calculate quadratic form component
    smat4 = transpose(y2) - ltprod3(j1, i, :, :) - ltmu(j1, i, :, :) + ltprod4(j1, i, :, :)
    t7 = matmul(smat4, matmul(sig22i, transpose(smat4)))
    
    ! Calculate log-likelihood for this sample
    tmp = -t1 * log(2.0d0 * epi) - &  ! Constant term
          t1 / 2.0d0 * log(det) - &    ! Determinant term
          0.5d0 * t2 * a - &           ! Trace term
          0.5d0 * t1 * t7(1, 1)       ! Quadratic form term
    
    ! Accumulate deviance
    ss1 = ss1 - 2.0d0 * tmp
  end do
  
  ! Calculate posterior mean deviance
  ltbarD(i) = ss1 / real(nrep, 8)
  
  ! Calculate effective number of parameters (pD)
  ltpD(i) = ltbarD(i) - ltDbar(i)
  
  ! Calculate DIC = Dbar + 2*pD
  ltDIC(i) = ltDbar(i) + 2.0d0 * ltpD(i)
end do
! Aggregate DIC results across all obs
! -----------------------------------------

! Initialize all DIC and pD accumulators
DICall = 0.0d0; pDall = 0.0d0
DICall1 = 0.0d0; pDall1 = 0.0d0
ipd_dic_s = 0.0d0; ipd_pd_s = 0.0d0
ipd_dic_y = 0.0d0; ipd_pd_y = 0.0d0
lDICall = 0.0d0; hDICall = 0.0d0; tDICall = 0.0d0
lpDall = 0.0d0; hpDall = 0.0d0; tpDall = 0.0d0
lhDICall = 0.0d0; lhpDall = 0.0d0
htDICall = 0.0d0; htpDall = 0.0d0
ltDICall = 0.0d0; ltpDall = 0.0d0
ipd_dic = 0.0d0; ipd_pd = 0.0d0

! Loop over all obs
do kk = 1, nsc
  ! Only process obs with complete data (nsc = obs with complete data)
    ! Aggregate main DIC models
    DICall = DICall + yDIC(kk)      ! Total DIC for primary model
    DICall1 = DICall1 + sDIC(kk)    ! Total DIC for alternative model
    
    ! Aggregate IPD DIC models
    ipd_dic_s = ipd_dic_s + s1DIC(kk)  ! IPD DIC for alternative model
    ipd_pd_s = ipd_pd_s + s1pD(kk)     ! IPD pD for alternative model
    ipd_dic_y = ipd_dic_y + y1DIC(kk)  ! IPD DIC for primary model
    ipd_pd_y = ipd_pd_y + y1pD(kk)     ! IPD pD for primary model
    ipd_dic = ipd_dic + sumDIC(kk)
    ipd_pd = ipd_pd + sumpD(kk)
    
    ! Aggregate univariate DIC models
    lDICall = lDICall + lDIC(kk)    ! LDL-C DIC
    hDICall = hDICall + hDIC(kk)    ! HDL-C DIC
    tDICall = tDICall + tDIC(kk)    ! TG DIC
    
    ! Aggregate bivariate DIC models
    lhDICall = lhDICall + lhDIC(kk)  ! LDL-C & HDL-C DIC
    htDICall = htDICall + htDIC(kk)  ! HDL-C & TG DIC
    ltDICall = ltDICall + ltDIC(kk)  ! LDL-C & TG DIC
    
    ! Aggregate pD values
    pDall = pDall + ypD(kk)         ! Primary model pD
    pDall1 = pDall1 + spD(kk)       ! Alternative model pD
    lpDall = lpDall + lpD(kk)       ! LDL-C pD
    hpDall = hpDall + hpD(kk)       ! HDL-C pD
    tpDall = tpDall + tpD(kk)       ! TG pD
    lhpDall = lhpDall + lhpD(kk)    ! LDL-C & HDL-C pD
    htpDall = htpDall + htpD(kk)    ! HDL-C & TG pD
    ltpDall = ltpDall + ltpD(kk)    ! LDL-C & TG pD
end do

! Store aggregated results in output array
DIC_out(1) = DICall        ! Total DIC (primary model)
DIC_out(2) = pDall         ! Total pD (primary model)
DIC_out(3) = DICall1       ! Total DIC (alternative model)
DIC_out(4) = pDall1        ! Total pD (alternative model)
DIC_out(5) = DICall + DICall1  ! Combined DIC (both models)
DIC_out(6) = pDall + pDall1    ! Combined pD (both models)
DIC_out(7) = ipd_dic_y     ! IPD DIC (primary model)
DIC_out(8) = ipd_pd_y      ! IPD pD (primary model)
DIC_out(9) = ipd_dic_s     ! IPD DIC (alternative model)
DIC_out(10) = ipd_pd_s     ! IPD pD (alternative model)
DIC_out(11) = ipd_dic  ! Double IPD primary DIC (potential error)
DIC_out(12) = ipd_pd   ! IPD alternative DIC + pD (potential error)
DIC_out(13) = lDICall      ! LDL-C DIC
DIC_out(14) = lpDall       ! LDL-C pD
DIC_out(19) = lhDICall     ! LDL-C & HDL-C DIC
DIC_out(20) = lhpDall      ! LDL-C & HDL-C pD
DIC_out(15) = hDICall      ! HDL-C DIC
DIC_out(16) = hpDall       ! HDL-C pD
DIC_out(21) = htDICall    ! HDL-C & TG DIC
DIC_out(22) = htpDall     ! HDL-C & TG pD
DIC_out(17) = tDICall      ! TG DIC
DIC_out(18) = tpDall       ! TG pD
DIC_out(23) = ltDICall    ! LDL-C & TG DIC
DIC_out(24) = ltpDall     ! LDL-C & TG pD

! Output aggregated DIC results
write(*, *) 'DIC:', DIC_out
! Write simulated data for convergence diagnostics (beta parameters)
outfile = './Results/' // 'NMAbeta.out' 
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=20, file=outfile, access='sequential', status='unknown')

do i2 = 1, nbbeta  ! Loop over output statistics (mean, sd, etc.)
  do jj = 1, J  ! Loop over outcomes
    write(20, 3001) (beta_out(i2, jj, i1), i1 = 1, 7)
  end do
end do
3001 format(7f25.4)  ! Format for beta parameters
close(20)

! Write phi parameters to file
outfile = './Results/' // 'NMAphi.out' 
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=21, file=outfile, access='sequential', status='unknown')
do i2 = 1, nphi  ! Loop over output statistics
  do jj = 1, J  ! Loop over outcomes
    write(21, 3002) (phi_out(i2, jj, i1), i1 = 1, 7)
  end do
end do
3002 format(7f25.4)  ! Format for phi parameters
close(21)

! beta
outfile = './Results/' // 'NMAseqbeta.out' 
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=22,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
  do jj=1,J
  write(22,3003) (seqbeta(jj,i2,i1),i2=1,nbbeta)
  enddo
enddo
3003 format(24f14.4)
close(22)

! phi
outfile = './Results/NMAseqphi.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=23,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
  do jj=1,J
  write(23,3004) (seqphi(jj,i2,i1),i2=1,nphi)
  enddo
enddo
3004 format(6f14.4)
close(23)

! gamma
outfile = './Results/NMAseqgamma.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=24,file=outfile, access='sequential',status='unknown')
do i1=1,nrep
    do i2=1,ngrp
        write(24,3005) (seqgamma(i2,jj,i1),jj=1,J)
    enddo
enddo
3005  format(3f25.4)
close(24)

! delta
outfile = './Results/NMAseqdelta.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=25,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
    do i2=1,ngrp
        write(25,3006) (seqdelta(i2,jj,i1),jj=1,J)
     enddo
3006  format(3f25.4)
enddo
close(25)

! sigma
outfile = './Results/NMAseqsigma.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=26,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
    do i2=1,ns
    write(26,3007) (seqsigma(i2,jj,i1),jj=1,6)
    enddo
3007  format(6f25.4)
enddo
close(26)

! degree of freedom
outfile = './Results/NMAseqv.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=27,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
    write(27,3008) seqv(i1)
enddo
3008  format(1f25.4)
close(27)

! R correlation
outfile = './Results/NMAseqErho.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=28,file=outfile, access='sequential',status='unknown')
do i1=1,nrep
    do i3=1,ns
        write(28,3009) (seqErho(i3,i2,i1),i2=1,3)
    enddo
enddo
3009  format(3f25.4)
close(28)

! Rho correlation

outfile = './Results/NMAseqRho.out'
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=29,file=outfile,access='sequential',status='unknown')
do i1=1,nrep
  do jj=1,J
  write(29,3010) (seqrho(i1,jj,i2),i2=1,55)
  enddo
enddo
3010 format(55f14.4)
close(29)

! Write DIC results to file
outfile = './Results/' // 'DIC.out' 
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=30, file=outfile, access='sequential', status='unknown')
do i1 = 1, 24  ! Loop over DIC output components
  write(30, 3014) DIC_out(i1)
end do
3014 format(1f25.4)  ! Format for DIC values
close(30)

! Write comprehensive diagnostics to main output file
outfile = './Results/' // 'MNMA.out' 
print *, 'outfile name is: [', trim(outfile), ']'
open(unit=15, file=outfile, access='sequential', status='unknown')

! Write run information
write(15, *) 'Gibbs sample size=', nrep
write(15, *) '--------------------------------------------------'

! Get current date and time
call idate(ntoday)
call itime(now)

! Write timing information
write(15, *) 'date: month/day/year and time: hour, minute, and second'
write(15, *) 'Beginning at date=', ntodayb(2), '/', ntodayb(1), '/', ntodayb(3)
write(15, *) 'time at', nowb(1), ':', nowb(2), ':', nowb(3)
write(15, *) 'Ending at date=', ntoday(2), '/', ntoday(1), '/', ntoday(3)
write(15, *) 'time at', now(1), ':', now(2), ':', now(3)

! Calculate and write elapsed time
total = etime(elapsed)
write(15, *) 'elapsed time in minutes'
write(15, 2216) total/60.0d0, elapsed(1)/60.0d0, elapsed(2)/60.0d0
2216 format('end: total=', f12.4, ' user=', f12.4, ' system=', f12.4)

! Write aggregated DIC results
write(15, *) '--------------------------------------------------'
write(15, *) 'AD yDIC=', DICall, ' pD', pDall
write(15, *) 'AD sDIC=', DICall1, ' pD', pDall1
write(15, *) 'IPD yDIC=', ipd_dic_y, ' pD', ipd_pd_y
write(15, *) 'IPD sDIC=', ipd_dic_s, ' pD', ipd_pd_s
write(15, *) 'L|HT LDLC=', lDICall, ' pD', lpDall
write(15, *) 'LH|T: DIC=', lhDICall, ' pD', lhpDall
write(15, *) 'H|LT HDLC=', hDICall, ' pD', hpDall
write(15, *) 'HT|L: DIC=', htDICall, ' pD', htpDall
write(15, *) 'T|LH TG=', tDICall, ' pD', tpDall
write(15, *) 'LT|H: DIC=', ltDICall, ' pD', ltpDall
write(15, *) 'AD DIC=', DICall + DICall1, ' pD', pDall + pDall1
write(15, *) 'IPD DIC=', ipd_dic_y + ipd_dic_s, ' pD', ipd_pd_y + ipd_pd_s
write(15, *) '--------------------------------------------------'

! Write subject-specific diagnostics for each model
! Primary model (y)
do kk = 1, nsc
  write(15, *) 'i=', kk, ' ypD=', ypD(kk), ' yDIC=', yDIC(kk), &
               ' barD', ybarD(kk), ' Dbar', yDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! Alternative model (s)
do kk = 1, nsc
  write(15, *) 'i=', kk, ' spD=', spD(kk), ' sDIC=', sDIC(kk), &
               ' barD', sbarD(kk), ' Dbar', sDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! IPD primary model (y1)
do kk = 1, nsc
  write(15, *) 'i=', kk, ' y1pD=', y1pD(kk), ' y1DIC=', y1DIC(kk)
end do
write(15, *) '--------------------------------------------------'

! IPD alternative model (s1)
do kk = 1, nsc
  write(15, *) 'i=', kk, ' s1pD=', s1pD(kk), ' s1DIC=', s1DIC(kk)
end do
write(15, *) '--------------------------------------------------'

! LDL-C model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' lpD=', lpD(kk), ' lDIC=', lDIC(kk), &
               ' barD', lbarD(kk), ' Dbar', lDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! LDL-C & HDL-C model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' lhpD=', lhpD(kk), ' lhDIC=', lhDIC(kk), &
               ' barD', lhbarD(kk), ' Dbar', lhDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! HDL-C model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' hpD=', hpD(kk), ' hDIC=', hDIC(kk), &
               ' barD', hbarD(kk), ' Dbar', hDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! HDL-C & TG model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' htpD=', htpD(kk), ' htDIC=', htDIC(kk), &
               ' barD', htbarD(kk), ' Dbar', htDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! TG model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' tpD=', tpD(kk), ' tDIC=', tDIC(kk), &
               ' barD', tbarD(kk), ' Dbar', tDbar(kk)
end do
write(15, *) '--------------------------------------------------'

! LDL-C & TG model
do kk = 1, nsc
  write(15, *) 'i=', kk, ' ltpD=', ltpD(kk), ' ltDIC=', ltDIC(kk), &
               ' barD', ltbarD(kk), ' Dbar', ltDbar(kk)
end do

! Close output file 
close(15)
stop

end program
