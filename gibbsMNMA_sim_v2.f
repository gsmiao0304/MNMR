subroutine gibbs(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed
  real(8) :: r

!  print *, 'Gibbs using seed=', iseed
  call r8_uniform_01_sample(r, iseed)
  
  print *, 'Ggamma'
  ! Generate gamma parameters
  call Ggamma(iseed)
  write(*,*) 'Exiting Ggamma'
  flush(6)
  
  print *, 'Gbeta'
  ! Generate beta parameters (treatment effects)
  call Gbeta(iseed)
  write(*,*) 'Exiting Gbeta'
  flush(6)
  
  
  ! Update derived variables after beta and gamma updates
  call vars1()
  write(*,*) 'Exiting vars1'
  flush(6)
  
  print *, 'GR'
  ! Generate R matrices (random error correlation)
  call GR(iseed)
  write(*,*) 'Exiting GR'
  flush(6)

  print *, 'GRho'
  ! Generate rho parameters (correlations)
  call GRho(iseed)
  write(*,*) 'Exiting GRho'
  flush(6)

  print *, 'Gphi'
  ! Generate phi parameters (covariate effects)
  call Gphi(iseed)
  write(*,*) 'Exiting Gphi'
  flush(6)

  print *, 'Gsigma'
  ! Generate Sigma matrices (residual covariance)
  call Gsigma(iseed)
  write(*,*) 'Exiting Gsigma'
  flush(6)

  print *, 'Gsiggrp'
  ! Generate group-level sigma parameters
  call Gsiggrp(iseed)
  write(*,*) 'Exiting Gsiggrp'
  flush(6)

  print *, 'Gv'
  ! Generate degrees of freedom parameter (v)
  call Gv(iseed)
  write(*,*) 'Exiting Gv'
  flush(6)

end subroutine gibbs

subroutine Ggamma(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3
  integer, parameter :: nbeta = 10, nb = 4, nbbeta = nx + nb + nT
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K)
  integer :: icount, icount4, idim, tv2, ifault, ij, nn, NULLTY, irank, lda, ldr, nr
  real(8) :: x(nx, ns), xx(nbbeta, ns), G(J*nbbeta, J)
  real(8) :: y(J, ns), ny(ns)
  real(8) :: beta(nbeta, J), bbeta(nbbeta, J)
  real(8) :: lam(K, J), Rgam(ns, J)
  real(8) :: kOmega(J, K, nTT, nTT)
  real(8) :: SigInv(ns, J, J), Rho(ns, J, J), Q(ns, J, J)
  real(8) :: sigma_temp(J, J), Rho_temp(J, J), Q_temp(J, J)
  real(8) :: InvSig_t(J, J)
  real(8) :: E(ns, nTT)
  real(8) :: tau(J, J), tauk(ns, J, J)
  real(8) :: tv12(6, 6), tv13(9, 9), tv14(12, 12)
  real(8) :: tv62(6), tv63(9), tv64(12)
  real(8) :: tv1_all(K, 12, 12), tv7_234(K, 12), tv7_all(ns, J), tv7(J)
  real(8) :: tv72(6, J), tv73(9, J), tv74(12, J)
  real(8) :: EO2(J, 6), EO3(J, 9), EO4(J, 12)
  real(8) :: WOW(12, 12), WOW2(6, 6), WOW3(9, 9), WOW4(12, 12)
  real(8) :: WOWInv2(6, 6), WOWInv3(9, 9), WOWInv4(12, 12)
  real(8) :: Sigpre2(6, 6), Sigpre3(9, 9), Sigpre4(12, 12)
  real(8) :: SigGam2(6, 6), SigGam3(9, 9), SigGam4(12, 12)
  real(8) :: MuGam2(6), MuGam3(9), MuGam4(12)
  real(8) :: SigGk2(K, 6, 6), SigGk3(K, 9, 9), SigGk4(K, 12, 12)
  real(8) :: tol, r2(6), r3(9), r4(12)
  real(8) :: Sig2chol(6, 6), Sig3chol(9, 9), Sig4chol(12, 12)
  real(8) :: RV2(6), RV3(9), RV4(12)
  real(8) :: AI2vec(21), AI3vec(45), AI4vec(78)
  real(8) :: upper2(21), upper3(45), upper4(78)
  real(8) :: tmp, bn, ss
  
  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecx/ x
  common /vecxx/ xx
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecbbeta/ bbeta
  common /vecRgam/ Rgam
  common /veclam/ lam
  common /vectauk/ tauk
  common /veckOmega/ kOmega
  common /vecSigInv/ SigInv
  common /vecE/ E
  common /vecSigGam2/ SigGk2
  common /vecSigGam3/ SigGk3
  common /vecSigGam4/ SigGk4
  common /vectv/ tv1_all, tv7_234, tv7_all

  ! Initialize variables
  icount4 = 0
  icount = 0
  SigGk2 = 0.0d0
  SigGk3 = 0.0d0
  SigGk4 = 0.0d0

  ! Loop over trials
  do kk = 1, K
    ! Initialize temporary variables
    tv62 = 0.0d0; tv63 = 0.0d0; tv64 = 0.0d0
    tv12 = 0.0d0; tv13 = 0.0d0; tv14 = 0.0d0
    tv72 = 0.0d0; tv73 = 0.0d0; tv74 = 0.0d0
    WOW = 0.0d0

    ! Get number of arms for this trial
    idim = narms(kk)
    tv2 = 3 * idim  ! Total dimension for gamma vector

    ! Construct WOW matrix = E' * Omega * E
    ! For LDL-C component
    WOW(1:idim, 1:idim) = matmul(matmul( &
      E((icount4+1):(icount4+idim), :), &
      kOmega(1, kk, :, :)), &
      transpose(E((icount4+1):(icount4+idim), :)))
    
    ! For HDL-C component
    WOW((idim+1):(2*idim), (idim+1):(2*idim)) = matmul(matmul( &
      E((icount4+1):(icount4+idim), :), &
      kOmega(2, kk, :, :)), &
      transpose(E((icount4+1):(icount4+idim), :)))
    
    ! For TG component
    WOW((2*idim+1):(3*idim), (2*idim+1):(3*idim)) = matmul(matmul( &
      E((icount4+1):(icount4+idim), :), &
      kOmega(3, kk, :, :)), &
      transpose(E((icount4+1):(icount4+idim), :)))

    ! Compute inverse based on dimension
    if (idim == 2) then
      WOW2 = WOW(1:tv2, 1:tv2)
      call inverse(WOW2, WOWInv2, tv2)
    else if (idim == 3) then
      WOW3 = WOW(1:tv2, 1:tv2)
      call inverse(WOW3, WOWInv3, tv2)
    else if (idim == 4) then
      WOW4 = WOW(1:tv2, 1:tv2)
      call inverse(WOW4, WOWInv4, tv2)
    end if

    ! Loop over arms within trial
    do jj = 1, idim
      icount = icount + 1
      tau = tauk(icount, :, :)  ! Get tau matrix for this subject
      InvSig_t = SigInv(icount, :, :)  ! Inverse Sigma for this subject
      tv7 = tv7_all(icount, :)  ! Residuals y - Xβ

      ! Construct EO matrix based on dimension
      if (idim == 2) then
        EO2 = 0.0d0
        EO2(1, jj) = 1.0d0
        EO2(2, jj + idim) = 1.0d0
        EO2(3, jj + 2*idim) = 1.0d0
        
        ! Compute tv72 = EO2' * tau * InvSig_t
        tv72 = matmul(matmul(transpose(EO2), tau), InvSig_t)
        
        ! Accumulate precision matrix component
        tv12 = tv12 + ny(icount) * matmul(tv72, matmul(tau, EO2))
        
        ! Accumulate mean component
        do i1 = 1, tv2
          ss = 0.0d0
          do i2 = 1, J
            ss = ss + tv72(i1, i2) * tv7(i2)
          end do
          tv62(i1) = tv62(i1) + ny(icount) * ss
        end do
        
      else if (idim == 3) then
        ! Similar calculations for idim=3
        EO3 = 0.0d0
        EO3(1, jj) = 1.0d0
        EO3(2, jj + idim) = 1.0d0
        EO3(3, jj + 2*idim) = 1.0d0
        tv73 = matmul(matmul(transpose(EO3), tau), InvSig_t)
        tv13 = tv13 + ny(icount) * matmul(tv73, matmul(tau, EO3))
        do i1 = 1, tv2
          ss = 0.0d0
          do i2 = 1, J
            ss = ss + tv73(i1, i2) * tv7(i2)
          end do
          tv63(i1) = tv63(i1) + ny(icount) * ss
        end do
        
      else if (idim == 4) then
        ! Similar calculations for idim=4
        EO4 = 0.0d0
        EO4(1, jj) = 1.0d0
        EO4(2, jj + idim) = 1.0d0
        EO4(3, jj + 2*idim) = 1.0d0
        tv74 = matmul(matmul(transpose(EO4), tau), InvSig_t)
        tv14 = tv14 + ny(icount) * matmul(tv74, matmul(tau, EO4))
        do i1 = 1, tv2
          ss = 0.0d0
          do i2 = 1, J
            ss = ss + tv74(i1, i2) * tv7(i2)
          end do
          tv64(i1) = tv64(i1) + ny(icount) * ss
        end do
      end if
    end do  ! End arm loop

    ! Compute posterior covariance and mean for gamma
    if (idim == 2) then
      Sigpre2 = tv12 + WOWInv2
      call inverse(Sigpre2, SigGam2, tv2)
      SigGk2(kk, :, :) = SigGam2
      
      ! Compute MuGam2 = SigGam2 * tv62
      do i1 = 1, tv2
        ss = 0.0d0
        do i2 = 1, tv2
          ss = ss + SigGam2(i1, i2) * tv62(i2)
        end do
        MuGam2(i1) = ss
      end do
      
      ! Sample from multivariate normal
      call cholesky(tv2, SigGam2, Sig2chol, ifault)
      call multinormal(tv2, MuGam2, Sig2chol, RV2, iseed)
      
      ! Store sampled gamma values
      do i2 = 1, idim
        Rgam(i2 + icount4, 1) = RV2(i2)
        Rgam(i2 + icount4, 2) = RV2(i2 + idim)
        Rgam(i2 + icount4, 3) = RV2(i2 + 2*idim)
      end do
      
    else if (idim == 3) then
      ! Similar process for idim=3
      Sigpre3 = tv13 + WOWInv3
      call inverse(Sigpre3, SigGam3, tv2)
      SigGk3(kk, :, :) = SigGam3
      do i1 = 1, tv2
        ss = 0.0d0
        do i2 = 1, tv2
          ss = ss + SigGam3(i1, i2) * tv63(i2)
        end do
        MuGam3(i1) = ss
      end do
      call cholesky(tv2, SigGam3, Sig3chol, ifault)
      call multinormal(tv2, MuGam3, Sig3chol, RV3, iseed)
      do i2 = 1, idim
        Rgam(i2 + icount4, 1) = RV3(i2)
        Rgam(i2 + icount4, 2) = RV3(i2 + idim)
        Rgam(i2 + icount4, 3) = RV3(i2 + 2*idim)
      end do
      
    else if (idim == 4) then
      ! Similar process for idim=4
      Sigpre4 = tv14 + WOWInv4
      call inverse(Sigpre4, SigGam4, tv2)
      SigGk4(kk, :, :) = SigGam4
      do i1 = 1, tv2
        ss = 0.0d0
        do i2 = 1, tv2
          ss = ss + SigGam4(i1, i2) * tv64(i2)
        end do
        MuGam4(i1) = ss
      end do
      call cholesky(tv2, SigGam4, Sig4chol, ifault)
      call multinormal(tv2, MuGam4, Sig4chol, RV4, iseed)
      do i2 = 1, idim
        Rgam(i2 + icount4, 1) = RV4(i2)
        Rgam(i2 + icount4, 2) = RV4(i2 + idim)
        Rgam(i2 + icount4, 3) = RV4(i2 + 2*idim)
      end do
    end if
    
    icount4 = icount4 + idim  ! Update subject counter
  end do  ! End trial loop

end subroutine Ggamma

subroutine Gbeta(iseed)
  implicit real(8) (a-h, o-z)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3
  integer, parameter :: nb = 4, nbbeta = nx + nb + nT
  real(8), parameter :: c01 = 100000.0d0
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K)
  integer :: icount, idim, ndim, ifault, i1, i2, j1, j2, kk, jj, tv2
  real(8) :: xx(nbbeta, ns), G(J*nbbeta, J)
  real(8) :: y(J, ns), ny(ns)
  real(8) :: bbeta(nbbeta, J)
  real(8) :: tauk(ns, J, J), tau(J, J)
  real(8) :: SigInv(ns, J, J)
  real(8) :: Sigbeta(J*nbbeta, J*nbbeta), Mubeta(J*nbbeta)
  real(8) :: Sigbeta_pre(J*nbbeta, J*nbbeta)
  real(8) :: tv1(J*nbbeta, J*nbbeta), tv3(J*nbbeta, J*nbbeta)
  real(8) :: tv4inv(J, J), tv5(J*nbbeta), tv8(J*nbbeta)
  real(8) :: temp1(J*nbbeta, J*nbbeta), temp2(J*nbbeta), temp3(J*nbbeta, J)
  
  ! Dimension-specific variables
  real(8) :: tv62(J*nbbeta, 6), tv63(J*nbbeta, 9), tv64(J*nbbeta, 12)
  real(8) :: tv72(6), tv73(9), tv74(12)
  real(8) :: tv92(6, J), tv93(9, J), tv94(12, J)
  real(8) :: tv82(J*nbbeta, 6), tv83(J*nbbeta, 9), tv84(J*nbbeta, 12)
  real(8) :: EO2(J, 6), EO3(J, 9), EO4(J, 12)
  real(8) :: SigGam2(6, 6), SigGam3(9, 9), SigGam4(12, 12)
  real(8) :: SigGk2(K, 6, 6), SigGk3(K, 9, 9), SigGk4(K, 12, 12)
  
  ! Sampling variables
  real(8) :: Sigchol(J*nbbeta, J*nbbeta), RV(J*nbbeta)
  integer, intent(inout) :: iseed

  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecxx/ xx
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vectauk/ tauk
  common /vecSigInv/ SigInv
  common /vecbbeta/ bbeta
  common /vecSigGam2/ SigGk2
  common /vecSigGam3/ SigGk3
  common /vecSigGam4/ SigGk4

  ! Initialize accumulators
  ndim = J * nbbeta
  tv1 = 0.0d0
  tv8 = 0.0d0
  icount = 0

  ! Loop over trials
  do kk = 1, K
    idim = narms(kk)  ! Number of arms in current trial
    tv2 = J * idim
    
    ! Initialize trial-specific accumulators
    tv3 = 0.0d0
    tv5 = 0.0d0
    tv62 = 0.0d0
    tv63 = 0.0d0
    tv64 = 0.0d0
    tv72 = 0.0d0
    tv73 = 0.0d0
    tv74 = 0.0d0

    ! Loop over arms within trial
    do jj = 1, idim
      icount = icount + 1
      
      ! Construct design matrix G for current subject
      G = 0.0d0
      do j1 = 1, J
        do j2 = 1, nbbeta
          G((j1-1)*nbbeta + j2, j1) = xx(j2, icount)
        end do
      end do
      
      ! Compute tau matrix and inverse Sigma for current subject
      tau = 0.0d0
      do i1 = 1, J
        do i2 = 1, J
          tv4inv(i1, i2) = SigInv(icount, i1, i2) * ny(icount)
          tau(i1, i2) = tauk(icount, i1, i2)
        end do
      end do
      
      ! Compute intermediate matrices
      temp1 = 0.0d0
      temp2 = 0.0d0
      temp3 = 0.0d0
      temp3 = matmul(G, tv4inv)                     
      temp1 = matmul(temp3, transpose(G))            
      
      ! Compute mean component
      do i1 = 1, ndim
        do i2 = 1, J
          temp2(i1) = temp2(i1) + temp3(i1, i2) * y(i2, icount)
        end do
      end do
      
      ! Accumulate precision and mean components
      tv3 = tv3 + temp1
      tv5 = tv5 + temp2
      
      ! Dimension-specific calculations
      select case(idim)
        case(2)  ! 2-arm trial
          EO2 = 0.0d0
          EO2(1, jj) = 1.0d0
          EO2(2, jj + idim) = 1.0d0
          EO2(3, jj + 2*idim) = 1.0d0
          tv62 = tv62 + matmul(matmul(temp3, tau), EO2)
          tv92 = matmul(matmul(transpose(EO2), tau), tv4inv)
          do i1 = 1, tv2
            do i2 = 1, J
              tv72(i1) = tv72(i1) + tv92(i1, i2) * y(i2, icount)
            end do
          end do
          
        case(3)  ! 3-arm trial
          EO3 = 0.0d0
          EO3(1, jj) = 1.0d0
          EO3(2, jj + idim) = 1.0d0
          EO3(3, jj + 2*idim) = 1.0d0
          tv63 = tv63 + matmul(matmul(temp3, tau), EO3)
          tv93 = matmul(matmul(transpose(EO3), tau), tv4inv)
          do i1 = 1, tv2
            do i2 = 1, J
              tv73(i1) = tv73(i1) + tv93(i1, i2) * y(i2, icount)
            end do
          end do
          
        case(4)  ! 4-arm trial
          EO4 = 0.0d0
          EO4(1, jj) = 1.0d0
          EO4(2, jj + idim) = 1.0d0
          EO4(3, jj + 2*idim) = 1.0d0
          tv64 = tv64 + matmul(matmul(temp3, tau), EO4)
          tv94 = matmul(matmul(transpose(EO4), tau), tv4inv)
          do i1 = 1, tv2
            do i2 = 1, J
              tv74(i1) = tv74(i1) + tv94(i1, i2) * y(i2, icount)
            end do
          end do
      end select
    end do  ! End arm loop
    
    ! Apply gamma adjustments based on number of arms
    select case(idim)
      case(2)  ! 2-arm adjustments
        SigGam2 = SigGk2(kk, :, :)
        tv82 = matmul(tv62, SigGam2)
        tv1 = tv1 + tv3 - matmul(tv82, transpose(tv62))
        do i1 = 1, ndim
          do i2 = 1, tv2
            tv8(i1) = tv8(i1) - tv82(i1, i2) * tv72(i2)
          end do
          tv8(i1) = tv8(i1) + tv5(i1)
        end do
        
      case(3)  ! 3-arm adjustments
        SigGam3 = SigGk3(kk, :, :)
        tv83 = matmul(tv63, SigGam3)
        tv1 = tv1 + tv3 - matmul(tv83, transpose(tv63))
        do i1 = 1, ndim
          do i2 = 1, tv2
            tv8(i1) = tv8(i1) - tv83(i1, i2) * tv73(i2)
          end do
          tv8(i1) = tv8(i1) + tv5(i1)
        end do
        
      case(4)  ! 4-arm adjustments
        SigGam4 = SigGk4(kk, :, :)
        tv84 = matmul(tv64, SigGam4)
        tv1 = tv1 + tv3 - matmul(tv84, transpose(tv64))
        do i1 = 1, ndim
          do i2 = 1, tv2
            tv8(i1) = tv8(i1) - tv84(i1, i2) * tv74(i2)
          end do
          tv8(i1) = tv8(i1) + tv5(i1)
        end do
    end select
  end do  ! End trial loop
  
  ! Add prior precision to diagonal
  do i1 = 1, ndim
    do i2 = 1, ndim
      Sigbeta_pre(i1, i2) = tv1(i1, i2)
    end do
    Sigbeta_pre(i1, i1) = tv1(i1, i1) + 1.0d0 / c01
  end do
  
  ! Compute posterior covariance matrix
  call inverse(Sigbeta_pre, Sigbeta, ndim)
  
  ! Compute posterior mean
  Mubeta = 0.0d0
  do i1 = 1, ndim
    do i2 = 1, ndim
      Mubeta(i1) = Mubeta(i1) + Sigbeta(i1, i2) * tv8(i2)
    end do
  end do
  
  ! Sample from multivariate normal distribution
  call cholesky(ndim, Sigbeta, Sigchol, ifault)
  call multinormal(ndim, Mubeta, Sigchol, RV, iseed)
  
  ! Store sampled beta values
  do i1 = 1, nbbeta
    do i2 = 1, J
      bbeta(i1, i2) = RV((i2-1)*nbbeta + i1)
    end do
  end do

end subroutine Gbeta

subroutine GR(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, J = 3
  real(8), parameter :: reqmin = 1.0d-10  ! Optimization tolerance
  integer, parameter :: konvge = 5       ! Convergence check interval
  integer, parameter :: kcount = 1000    ! Max optimization iterations
  real(8), parameter :: tol = 2.220446049250313d-14 
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K), iobs(K), icount5
  integer :: kk, jj, iR, iC, i, ni, ifault1, icount1, numres, nopt
  real(8) :: ny(ns), Erho(ns, 3), Srho(ns, 3)
  real(8) :: VoS(ns, J, J), VRV(ns, J, J), VRVtemp(J, J)
  real(8) :: RoS(ns, J, J), pRoS(ns, J, J), pRho(J, J), Rho(J, J)
  real(8) :: SigInv(ns, J, J), sigk(ns, J, J), sigtemp(J,J)
  real(8) :: tempv(J,J),temps(J,J)
  real(8) :: tempmat1(2,2),tempmat2(2,2),tv(1,1),tv2(2,2)
  real(8) :: psi111,psi112(1,2),psi122(2,2),psi122inv(2,2),sig1(1,1)
  real(8) :: s111,s112(1,2),s122(2,2),s122inv(2,2),mu1(1,2),mu1_v(2)
  real(8) :: rmvn(2),cholvec(2,2),mvnsig(2,2)
  real(8) :: psi211(2,2),psi221(1,2),psi222,psi211inv(2,2),sig2(2,2)
  real(8) :: s211(2,2),s221(1,2),s222,s211inv(2,2),mu2(1,2),mu2_v(2), rn(1,1),wis1(2,2)
  real(8) :: start, zprhomax, plmax, step, e1, sigmaa, zprhostar
  real(8) :: r1, r2, rat1, step1, current_corr, proposed_corr
  real(8) :: cl(5, 3), dl(5, 1), al(3, 3), alinv(3, 3), fl(3, 1), el(3, 1)
  integer :: index1, index2
  
  ! Common blocks
  common /vecny/ ny
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecobs/ iobs
  common /vecSigInv/ SigInv
  common /vecsigk/ sigk
  common /d3/ icount5
  common /vecVoS/ VoS
  common /vecpRoS/ pRoS
  common /vecRoS/ RoS
  common /vecpRho/ pRho
  common /vecVRV/ VRV
  common /vecErho/ Erho
  common /vecSrho/ Srho
  common  /vecindex1/index1
  common  /vecindex2/index2

  ! External functions and subroutines
  real(8), external :: bNloglike, r8_normal_01_sample,r8_chi_sample,FindDet

  ! Loop over trials
  do kk = 1, 34
    icount5 = iobs(kk) - 1  ! Starting index for current trial
    if (kk == 32) then
      do jj = 1, narms(kk)    ! Loop over arms in trial
        icount5 = icount5 + 1
        sigtemp=sigk(icount5,:,:)
        pRho=RoS(icount5,:,:)
        VRVtemp=0.0d0
        iR=2
        iC=3
        index1=iR
        index2=iC
        nopt = 1
        step = 0.02d0 
        current_corr = pRho(iR, iC)
        start = 0.5d0 * log((1.0d0 + current_corr) / (1.0d0 - current_corr))
        call nelmin(bNloglike, nopt, start, zprhomax, plmax, reqmin, &
                      step, konvge, kcount, icount1, numres, ifault1)
        step1 = 0.2d0
        eps   = 1.0d-10
        iter  = 0
        refine: do
          iter = iter + 1
            ! Evaluate log-likelihood at points around mode
            do i = 1, 5
              e1 = real(i - 3, 8)
              cl(i, 1) = (zprhomax + e1 * step1)**2
              cl(i, 2) = zprhomax + e1 * step1
              cl(i, 3) = 1.0d0
              dl(i, 1) = bNloglike(zprhomax + e1 * step1)
            end do
            
            ! Check if refinement is needed
            if (any( (/ dl(1,1), dl(2,1), dl(4,1), dl(5,1) /) <= plmax + eps )) then
              step1 = max(step1 * 0.80d0, 1.0d-4)
              if (iter < 30) cycle refine
            end if
          exit refine
        end do refine
          
          ! Fit quadratic model to log-likelihood
          al = matmul(transpose(cl), cl)
          el = matmul(transpose(cl), dl)
          call inverse(al, alinv, 3)
          fl = matmul(alinv, el)
          
          ! Compute proposal distribution parameters
          sigmaa = 1.0d0 / sqrt(2.0d0 * fl(1, 1))
          
          ! Generate proposal on Fisher z-scale
          r1 = r8_normal_01_sample(iseed)
          zprhostar = zprhomax + sigmaa * r1
          
          ! Compute Metropolis acceptance ratio
          rat1 = -bNloglike(zprhostar) + bNloglike(start) &
                 - 0.5d0 * (start - zprhomax)**2 / sigmaa**2 &
                 + 0.5d0 * (zprhostar - zprhomax)**2 / sigmaa**2
          
          ! Transform proposal back to correlation scale
          proposed_corr = (exp(2.0d0 * zprhostar) - 1.0d0) / (exp(2.0d0 * zprhostar) + 1.0d0)
          
          ! Generate uniform random number for Metropolis acceptance
          call r8_uniform_01_sample(r2, iseed)  ! Subroutine call
          
          ! Metropolis acceptance step
          if (rat1 >= 0.0d0 .or. log(r2) <= rat1) then
            pRho(iR, iC) = proposed_corr
          else
            pRho(iR, iC) = current_corr
          end if
          
          ! Ensure symmetry and store result
          pRho(iC, iR) = pRho(iR, iC) 
          pRoS(icount5,iC, iR)=pRho(iR, iC); pRoS(icount5,iR, iC)=pRoS(icount5,iC, iR)
          RoS(icount5,iC, iR)=pRho(iR, iC); RoS(icount5,iR, iC) = RoS(icount5,iC, iR)
          tempmat1=VoS(icount5,2:3,2:3)
          tempmat2=pRho(2:3,2:3)
          s122=matmul(tempmat1,matmul(tempmat2,tempmat1))
          psi122=sigtemp(2:3,2:3)
          psi112(1,1:2)=sigtemp(1,2:3)
          psi111=sigtemp(1,1)
          call inverse(s122, s122inv, 2)
          call inverse(psi122,psi122inv,2)
          tv=matmul(psi112,matmul(psi122inv,transpose(psi112)))
          sig1(1,1)=psi111-tv(1,1)
          mu1=transpose(matmul(s122,matmul(psi122inv,transpose(psi112))))
          mu1_v = mu1(1,:)
          mvnsig=sig1(1,1)*s122 / (ny(icount5) - 1.0d0)
          df=ny(icount5)-3.0d0
          call cholesky(2, mvnsig, cholvec, ifault)
          call multinormal(2, mu1_v, cholvec, rmvn, iseed)   
          ! call wishart_sample(1, int(df), sig1, rn, iseed)
          rn(1,1) = sig1(1,1) * r8_chi_sample ( df,iseed ) / (ny(icount5) - 1.0d0)
          s112(1,:) = rmvn
          tv=matmul(s112,matmul(s122inv,transpose(s112)))
          s111=rn(1,1)+tv(1,1)
          VRVtemp(2:3,2:3)=s122
          VRVtemp(1,1)=s111
          do ll=2,3
            VRVtemp(1,ll)=s112(1,ll-1)
            VRVtemp(ll,1)=s112(1,ll-1)
          end do    
          VRV(icount5,:,:)=VRVtemp 
          Srho(icount5,1)=dsqrt(VRV(icount5,1,1))
          Erho(icount5,1)=VRV(icount5,1,2)/Srho(icount5,1)/VoS(icount5,2,2)
          Erho(icount5,2)=VRV(icount5,1,3)/Srho(icount5,1)/VoS(icount5,3,3)
          Erho(icount5,3)=pRho(2,3)
      end do 
    else if (kk == 33 .or. kk == 34) then
      do jj = 1, narms(kk)    ! Loop over arms in trial
        icount5 = icount5 + 1
        VRVtemp = 0.0d0
        sigtemp=sigk(icount5,:,:)
        !write(*,*) 'i=', icount5, ' sigma=', sigtemp
        psi221(1,1:2)=sigtemp(3,1:2)
        !write(*,*) 'psi221',psi221
        psi222=sigtemp(3,3)
        !write(*,*) 'psi222',psi222
        psi211=sigtemp(1:2,1:2)
        !write(*,*) 'psi211',psi211
        s222=VoS(icount5,3,3)**2.0d0
        !write(*,*) 's222',s222
        sig2=psi211-matmul(transpose(psi221),psi221) / psi222
        ! write(*,*) 'sig2',sig2
        mvnsig=sig2*s222 / (ny(icount5)-1.0d0)
        ! write(*,*) 'mvnsig',mvnsig
        mu2=s222/psi222 *psi221
        mu2_v=mu2(1,:)
        !write(*,*) 'mu2',mu2_v
        df=ny(icount5)-2.0d0
        !write(*,*) 'df',df
        call cholesky(2, mvnsig, cholvec, ifault)
        call multinormal(2, mu2_v, cholvec, rmvn, iseed) 
        s221(1,:) = rmvn
        !write(*,*) 's221',s221
        call wishart_sample(2, int(df), sig2, wis1, iseed)
        !write(*,*) 'wis1',wis1
        tv2=matmul(transpose(s221),s221)/s222
        VRVtemp(1:2,1:2)=wis1 / (ny(icount5)-1.0d0)+tv2
        do ll=1,2
          VRVtemp(3,ll)=s221(1,ll)
          VRVtemp(ll,3)=s221(1,ll)
        end do
        VRVtemp(3,3)=s222
        VRV(icount5,:,:)=VRVtemp 
        do l1=1,2
          Srho(icount5,l1)=sqrt(VRV(icount5,l1,l1))
        end do
        Erho(icount5,1)=VRV(icount5,1,2)/Srho(icount5,1)/Srho(icount5,2)
        Erho(icount5,2)=VRV(icount5,1,3)/Srho(icount5,1)/VoS(icount5,3,3)     
        Erho(icount5,3)=VRV(icount5,2,3)/Srho(icount5,2)/VoS(icount5,3,3)
        ! write(*,*) 'Srho=',  Srho(icount5,:)
      det1 = FindDet(tol, J, VRVtemp)
      if (det1 <= 0.0d0) then
        print *, 'Warning: Invalid VRV at icount=', icount5
      end if
      end do
    else
      do jj = 1, narms(kk)    ! Loop over arms in trial
        icount5 = icount5 + 1  ! Current subject index
      
      ! Initialize correlation matrix for current subject
        pRho = pRoS(icount5, :, :)

      
      ! Sample correlations for each pair of response variables
        do iR = 1, 2
          do iC = iR + 1, 3
          index1=iR
          index2=iC
          ! Transform current correlation to Fisher z-scale
          current_corr = pRho(iR, iC)
          start = 0.5d0 * log((1.0d0 + current_corr) / (1.0d0 - current_corr))
          
          ! Find mode of log-likelihood using Nelder-Mead optimization
          nopt = 1  ! Number of optimization variables
          step = 0.02d0  ! Initial step size
          call nelmin(bNloglike, nopt, start, zprhomax, plmax, reqmin, &
                      step, konvge, kcount, icount1, numres, ifault1)
          
          ! Refine optimization with quadratic approximation
          step1 = 0.2d0
          eps   = 1.0d-10
          iter  = 0
          ref_loop: do
            iter = iter + 1
            ! Evaluate log-likelihood at points around mode
            do i = 1, 5
              e1 = real(i - 3, 8)
              cl(i, 1) = (zprhomax + e1 * step1)**2
              cl(i, 2) = zprhomax + e1 * step1
              cl(i, 3) = 1.0d0
              dl(i, 1) = bNloglike(zprhomax + e1 * step1)
            end do
            
            ! Check if refinement is needed
          if (any( (/ dl(1,1), dl(2,1), dl(4,1), dl(5,1) /) <= plmax + eps )) then
            step1 = max(step1 * 0.80d0, 1.0d-4)
            if (iter < 30) cycle ref_loop
          end if
          exit ref_loop
          end do ref_loop
          
          ! Fit quadratic model to log-likelihood
          al = matmul(transpose(cl), cl)
          el = matmul(transpose(cl), dl)
          call inverse(al, alinv, 3)
          fl = matmul(alinv, el)
          
          ! Compute proposal distribution parameters
          sigmaa = 1.0d0 / sqrt(2.0d0 * fl(1, 1))
          
          ! Generate proposal on Fisher z-scale
          r1 = r8_normal_01_sample(iseed)
          zprhostar = zprhomax + sigmaa * r1
          
          ! Compute Metropolis acceptance ratio
          rat1 = -bNloglike(zprhostar) + bNloglike(start) &
                 - 0.5d0 * (start - zprhomax)**2 / sigmaa**2 &
                 + 0.5d0 * (zprhostar - zprhomax)**2 / sigmaa**2
          
          ! Transform proposal back to correlation scale
          proposed_corr = (exp(2.0d0 * zprhostar) - 1.0d0) / (exp(2.0d0 * zprhostar) + 1.0d0)
          
          ! Generate uniform random number for Metropolis acceptance
          call r8_uniform_01_sample(r2, iseed)  ! Subroutine call
          
          ! Metropolis acceptance step
          if (rat1 >= 0.0d0 .or. log(r2) <= rat1) then
            pRho(iR, iC) = proposed_corr
          else
            pRho(iR, iC) = current_corr
          end if
          
          ! Ensure symmetry and store result
          pRho(iC, iR) = pRho(iR, iC)
          Erho(icount5, iC + iR - 2) = pRho(iR, iC)
        end do
      end do
      
      ! Convert partial correlations to full correlation matrix
      ! DIRECT CALL TO EXTERNAL SUBROUTINE AS IN ORIGINAL CODE
      call prhoTorho(J, pRho, Rho)
      
      ! Store results for current subject
      RoS(icount5, :, :) = Rho
      pRoS(icount5, :, :) = pRho
      VRV(icount5, :, :) = matmul(VoS(icount5, :, :), matmul(Rho, VoS(icount5, :, :)))
    
    end do
    end if
  end do

end subroutine GR


real(8) function bNloglike(zprho)
  implicit real(8) (a-h, o-z)
  real(8), intent(in) :: zprho  ! Fisher z-transformed correlation
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, J = 3
  real(8), parameter :: c10 = 4.0d0  ! Prior precision constant
  real(8), parameter :: tol = 2.220446049250313d-14  ! Tolerance for determinant calculation
  
  ! Declare variables
  integer :: ids(ns), icount5, index1, index2, k1, j1, j2, i1, i2
  real(8) :: ny(ns), SigInv(ns, J, J), VoS(ns, J, J), pRho(J, J)
  real(8) :: temppRho(J, J), tempRho(J, J), siginv1(J, J), VRV1(J, J), temppRho2(2, 2), tempRho2(2, 2)
  real(8) :: temp1(J, J), det2, tr, s, zprho_corr, log_jacobian
  real(8) :: siginv2(2,2), VRV2(2,2),temp2(2,2)
  real(8) :: E0, E1(1,3),E2(2,3),E3(3,3)
  
  ! Common blocks
  common /d3/ icount5
  common /vecids/ ids
  common /vecny/ ny
  common /vecindex1/ index1
  common /vecindex2/ index2
  common /vecVoS/ VoS
  common /vecpRho/ pRho
  common /vecSigInv/ SigInv

  ! External functions
  real(8), external :: FindDet  ! Determinant calculation function
  external :: prhoTorho      ! Partial to full correlation conversion

  ! Get subject ID and current correlation matrix
  k1 = ids(icount5)
  do j1=1,J
    do j2=1,J
      E3(j1,j2)=0.0d0
    end do
    E3(j1,j1)=1.0d0
  end do
  E0=0.0d0
  E1=E3(3:3,:)
  E2(1,:)=E3(2,:)
  E2(2,:)=E3(3,:)
! Compute log Jacobian of transformation
  log_jacobian = 2.0d0 * zprho - 2.0d0 * log(exp(2.0d0 * zprho) + 1.0d0)
  
  temppRho = pRho
  
  ! Update the specified correlation with proposed value
  zprho_corr = (exp(2.0d0 * zprho) - 1.0d0) / (exp(2.0d0 * zprho) + 1.0d0)
  temppRho(index1, index2) = zprho_corr
  temppRho(index2, index1) = zprho_corr
  
  ! Convert partial correlations to full correlation matrix
  call prhoTorho(J, temppRho, tempRho)
  ! missing one sd
  if (k1 == 32) then
    siginv2=SigInv(icount5,2:3,2:3)
    s=0.5 * (ny(icount5)-2.0d0-real(2))  * log(1.0d0 - zprho_corr**2) &
      - (ny(icount5)-1.0d0) * siginv2(1,2) * VoS(icount5,2,2) * VoS(icount5,3,3) * zprho_corr &
      + log_jacobian
  end if
  if (k1 < 32) then
  ! Get inverse covariance and VoS matrix for current subject
    siginv1 = SigInv(icount5, :, :)
  
  ! Compute VRV1 = VoS * Rho * VoS 
    VRV1 = matmul(VoS(icount5, :, :), matmul(tempRho, transpose(VoS(icount5, :, :))))
  
  ! Compute determinant of Rho matrix
    det2 = FindDet(tol, J, tempRho)
  
  ! Compute trace term: tr(SigInv * VRV1)
    temp1 = matmul(siginv1, VRV1)
    tr = temp1(1, 1) + temp1(2, 2) + temp1(3, 3)
  
  
  ! Compute log-likelihood components
    s = 0.5d0 * (ny(icount5) - 2.0d0 - real(J, 8)) * log(det2) &
      - 0.5d0 * (ny(icount5) - 1.0d0) * tr &
      + 0.5d0 * (real(J-1) - real(index2 - index1, 8)) * log(1.0d0 - zprho_corr**2) &
      + log_jacobian
  end if
  ! Return negative log-likelihood (for minimization)
  bNloglike = -s
end function bNloglike


subroutine Gphi(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, J = 3, nphi = 6
  integer, parameter :: nb = 4, nbbeta = nx + nb + nT
  real(8), parameter :: reqmin = 1.0d-10  ! Optimization tolerance
  integer, parameter :: konvge = 5       ! Convergence check interval
  integer, parameter :: kcount = 1000    ! Max optimization iterations
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K), inphi1, inphi2
  integer :: j1, j2, i, ni, ifault1, icount1, numres, nopt, ial
  real(8) :: x(nx, ns), ny(ns), xx(nbbeta, ns), w(ns, nphi)
  real(8) :: y(J, ns), phi(nphi, J), tauk(ns, J, J)
  real(8) :: start, ommax, plmax, step, e1, sigmaa, omkstar
  real(8) :: r1, r2, rat1, step1, current_phi
  real(8) :: cl(5, 3), dl(5, 1), al(3, 3), alinv(3, 3), fl(3, 1), el(3, 1)
  
  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecxx/ xx
  common /vecw/ w
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecphi/ phi
  common /vectauk/ tauk
  common /dummy3/ inphi1
  common /dummy4/ inphi2

  ! External functions
  real(8), external :: cNloglike, r8_normal_01_sample
  external :: r8_uniform_01_sample  ! Uniform random number generator

  ! Loop over phi parameters (covariate effects)
  do j1 = 1, nphi    ! Loop over covariates
    do j2 = 1, J     ! Loop over response variables
      ! Set current parameter indices
      inphi1 = j1
      inphi2 = j2
      current_phi = phi(j1, j2)
      
      ! Initialize optimization parameters
      nopt = 1       ! Number of optimization variables
      step = 0.05d0  ! Initial step size
      
      ! Find mode of log-likelihood using Nelder-Mead optimization
      call nelmin(cNloglike, nopt, current_phi, ommax, plmax, reqmin, &
                  step, konvge, kcount, icount1, numres, ifault1)
      
      ! Refine optimization with quadratic approximation
      step1 = 0.20d0
      iter  = 0
      eps   = 1.0d-10
      refine: do
        iter = iter + 1
        ! Evaluate log-likelihood at points around mode
        do i = 1, 5
          e1 = real(i - 3, 8)
          cl(i, 1) = (ommax + e1 * step1)**2
          cl(i, 2) = ommax + e1 * step1
          cl(i, 3) = 1.0d0
          dl(i, 1) = cNloglike(ommax + e1 * step1)
        end do
        
        ! Check if refinement is needed
        if (any( (/ dl(1,1), dl(2,1), dl(4,1), dl(5,1) /) <= plmax + eps)) then
          step1 = max(step1 * 0.80d0, 1.0d-4)
          if (iter < 30) cycle refine
        end if
        exit refine
      end do refine
      
      ! Fit quadratic model to log-likelihood
      al = matmul(transpose(cl), cl)
      el = matmul(transpose(cl), dl)
      call inverse(al, alinv, 3)
      fl = matmul(alinv, el)
      
      ! Compute proposal distribution parameters
      sigmaa = 1.0d0 / sqrt(2.0d0 * fl(1, 1))
      
      ! Generate proposal from normal distribution
!      print *, 'normal in iseed=', iseed
      r1 = r8_normal_01_sample(iseed)
!      print *, 'normal out iseed=', iseed
      omkstar = ommax + sigmaa * r1
      
      ! Compute Metropolis acceptance ratio
      rat1 = -cNloglike(omkstar) + cNloglike(current_phi) &
             - 0.5d0 * (current_phi - ommax)**2.0d0 / sigmaa**2.0d0 &
             + 0.5d0 * (omkstar - ommax)**2.0d0 / sigmaa**2.0d0
      
      ! Generate uniform random number for Metropolis acceptance
!      print *, 'uniform in iseed=', iseed
      call r8_uniform_01_sample(r2, iseed)
!      print *, 'uniform out iseed=', iseed
      
      ! Metropolis acceptance step
      if (rat1 >= 0.0d0 .or. log(r2) <= rat1) then
        phi(j1, j2) = omkstar
      else
        phi(j1, j2) = current_phi
      end if
    end do
  end do
  
  ! Update tau_k matrices with new phi values
  do i = 1, ns
    do l1 = 1, J
      tauk(i, l1, l1) = exp(dot_product(w(i, :), phi(:, l1)))
    end do
  end do

end subroutine Gphi

real(8) function cNloglike(phiij)
  implicit real(8) (a-h, o-z)
  real(8), intent(in) :: phiij  ! Current phi parameter value
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3, nphi = 6
  integer, parameter :: nb = 4, nbbeta = nx + nb + nT, narms_max = 4
  real(8), parameter :: c02 = 4.0d0  ! Prior precision constant
  real(8), parameter :: tol = 2.220446049250313d-14  ! Tolerance for determinant calculation
  
real(8) :: x(nx, ns), xx(nbbeta, ns), w(ns, nphi)
real(8) :: y(J, ns), ny(ns)
real(8) :: bbeta(nbbeta, J)
real(8) :: tOmega(nTT, nTT), kOmega(J, K, nTT, nTT)
real(8) :: SigInv(ns, J, J), InvSig_t(J, J)
real(8) :: E(ns, nTT), EO(J, 12), EO2(J, 6), EO3(J, 9), EO4(J, 12)
real(8) :: tau(J, J)
real(8) :: WOW(12, 12)
real(8) :: WOW2(6, 6), WOW3(9, 9), WOW4(12, 12)
real(8) :: WOWInv2(6, 6), WOWInv3(9, 9), WOWInv4(12, 12)
real(8) :: Sigpre2(6, 6), Sigpre3(9, 9), Sigpre4(12, 12)
real(8) :: SigGam2(6, 6), SigGam3(9, 9), SigGam4(12, 12)
real(8) :: tv72(6), tv73(9), tv74(12)
real(8) :: LR2(6, 6), LR3(9, 9), LR4(12, 12)
real(8) :: LRvec2(6), LRvec3(9), LRvec4(12)
real(8) :: s1, s2, det2
real(8) :: tv1(12, 12), tv4(J, J), tv4inv(J, J)
real(8) :: tv1_all(K, 12, 12), tv7_234(K, 12), tv7_all(ns, J), tv7(J)
real(8) :: phi(nphi, J), tphi(nphi, J)
integer :: ids(ns), iarm(ns), narms(K)
integer :: icount, icount4, tv2, idim, i1, i2, kk, jj, l1, l2
real(8), external :: FindDet
external :: inverse

! Common blocks
common /vecy/ y
common /vecny/ ny
common /vecxx/ xx
common /vecw/ w
common /vecids/ ids
common /veciarm/ iarm
common /vecnarms/ narms
common /vecphi/ phi
common /veckOmega/ kOmega
common /vecSigInv/ SigInv
common /vecE/ E
common /vecbbeta/ bbeta
common /dummy3/ inphi1
common /dummy4/ inphi2
common /vectv/ tv1_all, tv7_234, tv7_all

! Copy phi to tphi
do i1 = 1, nphi
    do i2 = 1, J
        tphi(i1, i2) = phi(i1, i2)
    end do
end do

! Update tphi with proposed value
tphi(inphi1, inphi2) = phiij

s1 = 0.0d0
s2 = 0.5d0 * phiij**2.0d0 / c02
icount = 0
icount4 = 0

do kk = 1, K
    tv1(:, :) = 0.0d0
    tv72(:) = 0.0d0
    tv73(:) = 0.0d0
    tv74(:) = 0.0d0
    idim = narms(kk)
    tv2 = 3 * idim
    WOW(:, :) = 0.0d0
    WOW(1:idim, 1:idim) = matmul(matmul( &
        E((icount4+1):(icount4+idim), :), &
        kOmega(1, kk, :, :)), &
        transpose(E((icount4+1):(icount4+idim), :)))
    
    WOW((idim+1):(2*idim), (idim+1):(2*idim)) = matmul(matmul( &
        E((icount4+1):(icount4+idim), :), &
        kOmega(2, kk, :, :)), &
        transpose(E((icount4+1):(icount4+idim), :)))
    
    WOW((2*idim+1):(3*idim), (2*idim+1):(3*idim)) = matmul(matmul( &
        E((icount4+1):(icount4+idim), :), &
        kOmega(3, kk, :, :)), &
        transpose(E((icount4+1):(icount4+idim), :)))
    
    if (idim == 2) then
        WOW2 = WOW(1:tv2, 1:tv2)
        call inverse(WOW2, WOWInv2, tv2)
    end if
    if (idim == 3) then
        WOW3 = WOW(1:tv2, 1:tv2)
        call inverse(WOW3, WOWInv3, tv2)
    end if
    if (idim == 4) then
        WOW4 = WOW(1:tv2, 1:tv2)
        call inverse(WOW4, WOWInv4, tv2)
    end if
    
    do jj = 1, idim
        EO2(:, :) = 0.0d0
        EO3(:, :) = 0.0d0
        EO4(:, :) = 0.0d0
        icount = icount + 1
        do l1 = 1, J
            do l2 = 1, J
                tau(l1, l2) = 0.0d0
            end do
            tau(l1, l1) = exp(dot_product(w(icount, :), tphi(:, l1)))
        end do
        
        ! Set up EO matrix
        do l1 = 1, J
            do l2 = 1, 12
                EO(l1, l2) = 0.0d0
            end do
        end do
        EO(1, jj) = 1.0d0
        EO(2, jj + idim) = 1.0d0
        EO(3, jj + 2*idim) = 1.0d0
        
        tv7 = tv7_all(icount, :)
        
        InvSig_t(:, :) = 0.0d0
        InvSig_t = SigInv(icount, :, :)
        
        tv1(1:tv2, 1:tv2) = tv1(1:tv2, 1:tv2) + ny(icount) * &
            matmul(matmul(matmul(transpose(EO(:, 1:tv2)), tau), &
            InvSig_t), matmul(tau, EO(:, 1:tv2)))
        
        if (idim == 2) then
            EO2(1, jj) = 1.0d0
            EO2(2, jj + idim) = 1.0d0
            EO2(3, jj + 2*idim) = 1.0d0
            tv72 = tv72 + ny(icount) * &
                matmul(matmul(matmul(transpose(EO2), tau), InvSig_t), tv7)
        end if
        if (idim == 3) then
            EO3(1, jj) = 1.0d0
            EO3(2, jj + idim) = 1.0d0
            EO3(3, jj + 2*idim) = 1.0d0
            tv73 = tv73 + ny(icount) * &
                matmul(matmul(matmul(transpose(EO3), tau), InvSig_t), tv7)
        end if
        if (idim == 4) then
            EO4(1, jj) = 1.0d0
            EO4(2, jj + idim) = 1.0d0
            EO4(3, jj + 2*idim) = 1.0d0
            tv74 = tv74 + ny(icount) * &
                matmul(matmul(matmul(transpose(EO4), tau), InvSig_t), tv7)
        end if
    end do
    
    if (idim == 2) then
        Sigpre2 = tv1(1:tv2, 1:tv2) + WOWInv2
        call inverse(Sigpre2, SigGam2, tv2)
        det2 = FindDet(tol, tv2, SigGam2)
        s1 = s1 + 0.5d0 * log(det2) + 0.5d0 * dot_product( &
            matmul(SigGam2, tv72), tv72)
    else if (idim == 3) then
        Sigpre3 = tv1(1:tv2, 1:tv2) + WOWInv3
        call inverse(Sigpre3, SigGam3, tv2)
        det2 = FindDet(tol, tv2, SigGam3)
        s1 = s1 + 0.5d0 * log(det2) + 0.5d0 * dot_product( &
            matmul(SigGam3, tv73), tv73)
    else
        Sigpre4 = tv1(1:tv2, 1:tv2) + WOWInv4
        call inverse(Sigpre4, SigGam4, tv2)
        det2 = FindDet(tol, tv2, SigGam4)
        s1 = s1 + 0.5d0 * log(det2) + 0.5d0 * dot_product( &
            matmul(SigGam4, tv74), tv74)
    end if
    
    icount4 = icount4 + idim
end do

  
  ! Final log-likelihood
  cNloglike = -s1 + s2

end function cNloglike


subroutine GRho(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3
  real(8), parameter :: reqmin = 1.0d-10  ! Optimization tolerance
  integer, parameter :: konvge = 5        ! Convergence check interval
  integer, parameter :: kcount = 1000     ! Max optimization iterations
  
  ! Declare variables
  integer :: jj, iR, iC, i, ni, ifault1, icount, numres, nopt, ial
  real(8) :: Omega(J, nTT, nTT), pOmega(J, nTT, nTT)
  real(8) :: Rho(nTT, nTT), pRho(nTT, nTT)
  real(8) :: kOmega(J, K, nTT, nTT), lam(K, J)
  real(8) :: start, zprhomax, plmax, step, e1, sigmaa, zprhostar
  real(8) :: r1, r2, rat1, step1, current_corr, proposed_corr
  real(8) :: cl(5, 3), dl(5, 1), al(3, 3), alinv(3, 3), fl(3, 1), el(3, 1)
  integer :: index3, index4, index5
  
  ! Common blocks
  common /vecOmega/ Omega
  common /vecpOmega/ pOmega
  common /veckOmega/ kOmega
  common /vecindex3/ index3
  common /vecindex4/ index4
  common /vecindex5/ index5
  common /Rho/ Rho
  common /pRho/ pRho
  common /veclam/ lam

  ! External functions
  real(8), external :: gNloglike, r8_normal_01_sample
  external :: r8_uniform_01_sample, prhoTorho

  ! Loop over response variables
  do jj = 1, J
    index5 = jj  ! Current response index
    pRho = pOmega(jj, :, :)  ! Get current partial correlation matrix
    
    ! Loop over pairs of time points
    do iR = 1, nTT - 1
      do iC = iR + 1, nTT
!        print *, 'j=', jj,' iR=',iR, ' iC=',iC
        index3 = iR  ! Row index
        index4 = iC  ! Column index
        
        ! Transform current correlation to Fisher z-scale
        current_corr = pRho(iR, iC)
        start = 0.5d0 * log((1.0d0 + current_corr) / (1.0d0 - current_corr))
        
        ! Find mode of log-likelihood using Nelder-Mead optimization
        nopt = 1  ! Number of optimization variables
        step = 0.02d0  ! Initial step size
        call nelmin(gNloglike, nopt, start, zprhomax, plmax, reqmin, &
                    step, konvge, kcount, icount, numres, ifault1)
        
        ! Refine optimization with quadratic approximation
        step1 = 0.2d0
        eps   = 1.0d-10
        iter  = 0
        refine: do
          iter = iter + 1
          ! Evaluate log-likelihood at points around mode
          do i = 1, 5
            e1 = real(i - 3, 8)
            cl(i, 1) = (zprhomax + e1 * step1)**2
            cl(i, 2) = zprhomax + e1 * step1
            cl(i, 3) = 1.0d0
            dl(i, 1) = gNloglike(zprhomax + e1 * step1)
          end do
          
          ! Check if refinement is needed
          if (any( (/ dl(1,1), dl(2,1), dl(4,1), dl(5,1) /) <= plmax + eps )) then
            step1 = max(step1 * 0.80d0, 1.0d-4)
            if (iter < 30) cycle refine
          end if
          exit refine
        end do refine
        
        ! Fit quadratic model to log-likelihood
        al = matmul(transpose(cl), cl)
        el = matmul(transpose(cl), dl)
        call inverse(al, alinv, 3)
        fl = matmul(alinv, el)
        
        ! Compute proposal distribution parameters
        sigmaa = 1.0d0 / sqrt(2.0d0 * fl(1, 1))
        
        ! Generate proposal on Fisher z-scale
        r1 = r8_normal_01_sample(iseed)
        zprhostar = zprhomax + sigmaa * r1
        
        ! Compute Metropolis acceptance ratio
        rat1 = -gNloglike(zprhostar) + gNloglike(start) &
               - 0.5d0 * (start - zprhomax)**2 / sigmaa**2 &
               + 0.5d0 * (zprhostar - zprhomax)**2 / sigmaa**2
        
        ! Transform proposal back to correlation scale
        proposed_corr = (exp(2.0d0 * zprhostar) - 1.0d0) / (exp(2.0d0 * zprhostar) + 1.0d0)
        
        ! Generate uniform random number for Metropolis acceptance
        call r8_uniform_01_sample(r2, iseed)
        
        ! Metropolis acceptance step
        if (rat1 >= 0.0d0 .or. log(r2) <= rat1) then
          pRho(iR, iC) = proposed_corr
          pRho(iC, iR) = proposed_corr
        end if
!        print *, ' pRho(iR, iC)=',  pRho(iR, iC)
        ! Convert partial correlations to full correlation matrix
        call prhoTorho(nTT, pRho, Rho)
      end do
    end do
    
    ! Store updated matrices
    pOmega(jj, :, :) = pRho
    Omega(jj, :, :) = Rho
  end do
  
  ! Update kOmega matrices
  do l1 = 1, J
    do kk = 1, K
      kOmega(l1, kk, :, :) = Omega(l1, :, :) / lam(kk, l1)
    end do
  end do

end subroutine GRho

real(8) function gNloglike(zprho)
  implicit real(8) (a-h, o-z)
  real(8), intent(in) :: zprho  ! Input parameter (Fisher z-transformed correlation)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3
  integer, parameter :: nb = 4, nbbeta = nx + nb + nT
  real(8), parameter :: c03 = 4.0d0, tol = 2.220446049250313d-14
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K), iobs(K), icount, icount4, idim, tv2
  integer :: j1, j2, j3, kk, iL, index3, index4, index5
  real(8) :: x(nx, ns), ny(ns), xx(nbbeta, ns), y(J, ns)
  real(8) :: bbeta(nbbeta, J), Omega(J, nTT, nTT), Ome(J, nTT, nTT), kOme(J, nTT, nTT)
  real(8) :: rho(nTT, nTT), prho(nTT, nTT), temppRho(nTT, nTT), tempRho(nTT, nTT)
  real(8) :: SigInv(ns, J, J), E(ns, nTT), tau(J, J), tauk(ns, J, J), lam(K, J)
  real(8), allocatable :: WOW(:,:), E_sub(:,:)
  real(8), allocatable :: WOW_now(:,:), WOWInv(:,:)
  real(8), allocatable :: tv7(:), Sigpre(:,:), SigGam(:,:)
  real(8) :: tv1_all(K, 12, 12), tv7_234(K, 12), tv7_all(ns, J)
  real(8) :: s1, s2, s3, det1, det2, fisher_z, temp_corr, quad_form
  real(8) :: kOme_j(nTT,nTT)
  
  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecxx/ xx
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecobs/ iobs
  common /veclam/ lam
  common /vectauk/ tauk
  common /vecSigInv/ SigInv
  common /vecE/ E
  common /vecbbeta/ bbeta
  common /vecindex3/ index3
  common /vecindex4/ index4
  common /vecindex5/ index5
  common /Rho/ Rho
  common /pRho/ pRho
  common /vecOmega/ Omega
  common /vectv/ tv1_all, tv7_234, tv7_all

  ! External functions
  real(8), external :: FindDet  ! Determinant calculation function
  external :: inverse, prhoTorho  ! Matrix inversion and correlation conversion

  ! Update partial correlation matrix with proposed value
  temppRho = pRho
  fisher_z = exp(2.0d0 * zprho)
  temp_corr = (fisher_z - 1.0d0) / (fisher_z + 1.0d0)
  temppRho(index3, index4) = temp_corr
  temppRho(index4, index3) = temp_corr

  ! Convert partial correlations to full correlation matrix
  call prhoTorho(nTT, temppRho, tempRho)
  
  ! Update group-level covariance matrix
  Ome = Omega
  Ome(index5, :, :) = tempRho
  
  ! Compute log Jacobian term
  iL = index4 - index3
  s2 = 0.50d0 * (real(nT) - real(iL, 8)) * log(1.0d0 - temp_corr**2.0d0) &
       + 2.0d0 * zprho - 2.0d0 * log(fisher_z + 1.0d0)
  
  ! Initialize log-likelihood accumulator
  s3 = 0.0d0
  icount4 = 0

  ! Loop over trials
  do kk = 1, K
    idim = narms(kk)  ! Number of arms in current trial
!    print *, 'kk=', kk,' dim=',idim
    tv2 = 3 * idim    ! Total dimension for gamma vector
    if (allocated(E_sub)) deallocate(E_sub)
    allocate(E_sub(idim,nTT))

    E_sub = E(icount4+1:icount4+idim, :)    
    ! Scale Omega by lambda for this trial
    do j1 = 1, J
      kOme(j1, :, :) = Ome(j1, :, :) / lam(kk, j1)
    end do
  
    
    if (allocated(WOW)) deallocate(WOW)
    allocate(WOW(tv2, tv2))

    if (allocated(WOWInv)) deallocate(WOWInv)
    allocate(WOWInv(tv2, tv2))

    if (allocated(Sigpre)) deallocate(Sigpre)
    allocate(Sigpre(tv2, tv2))

    if (allocated(SigGam)) deallocate(SigGam)
    allocate(SigGam(tv2, tv2))

    if (allocated(tv7)) deallocate(tv7)
    allocate(tv7(tv2))
    ! Construct WOW matrix = E' * kOme * E
    WOW = 0.0d0
    do j1 = 1, J
       kOme_j = kOme(j1, :, :)
       idx_start = (j1 - 1) * idim + 1
       idx_end = j1 * idim
       WOW(idx_start:idx_end, idx_start:idx_end) = matmul( &
         matmul(E_sub, kOme_j), transpose(E_sub))
    end do

    ! Compute determinant and inverse based on number of arms
    WOW_now = WOW

    det1 = FindDet(tol, tv2, WOW_now)
    call inverse(WOW_now, WOWInv, tv2)

    ! Compute log-likelihood components based on number of arms
    Sigpre = tv1_all(kk, 1:tv2, 1:tv2) + WOWInv
    call inverse(Sigpre, SigGam, tv2)
    det2 = FindDet(tol, tv2, Sigpre)
    tv7 = tv7_234(kk, 1:tv2)
    quad_form = dot_product(matmul(SigGam, tv7), tv7)
    s1 = -0.5d0 * (log(det1) + log(det2)) + 0.5d0 * quad_form
    
    ! Accumulate log-likelihood
    s3 = s3 + s1
    icount4 = icount4 + idim
  end do
  
  ! Combine all components
  gNloglike = -s3 - s2

end function gNloglike

subroutine Gsigma(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, J = 3, ngrp = 1
  integer, parameter :: nb = 4, nbbeta = nx + nb + nT
  real(8), parameter :: tol = 2.220446049250313d-14  ! Tolerance for determinant calculation
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K), grp(ns), icount, kk, jj, j1, j2, l1, l2, idim
  real(8) :: ny(ns), xx(nbbeta, ns), y(J, ns), Rgam(ns, J), df, ve
  real(8) :: bbeta(nbbeta, J), SigInv(ns, J, J), sigk(ns, J, J), sigma(ns, 6)
  real(8) :: tauk(ns, J, J), tau(J, J), VRV(ns, J, J), tempvrv(J, J)
  real(8) :: sig(ngrp, J, J), ptk(J, J), ptk_inv(J, J), wis(J, J), sigtest(J, J)
  real(8) :: tv(J), tv1(J, J), det1, residual
  
  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecxx/ xx
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecRgam/ Rgam
  common /vecbbeta/ bbeta
  common /vecSigInv/ SigInv
  common /vecsigk/ sigk
  common /vecSigma/ sigma
  common /vectauk/ tauk
  common /vecVRV/ VRV
  common /vecSig/ sig
  common /vecgrp/ grp
  common /numve/ ve

  ! External subroutines
  external :: inverse, wishart_sample
  real(8), external :: FindDet

  ! Initialize observation counter
  icount = 0

  ! Loop over trials
  do kk = 1, K
    idim = narms(kk)  ! Number of arms in current trial
    
    ! Loop over arms within trial
    do jj = 1, idim
      icount = icount + 1
      
      ! Get tau matrix and VRV for current subject
      tau = tauk(icount, :, :)
      tempvrv = VRV(icount, :, :)
      
      
      ! Set degrees of freedom for Wishart distribution
      if (kk > 34) then
        df = 1.0d0 + ve  ! Reduced DF for trials 34+
      else
        df = ny(icount) + ve  ! Full DF for other trials
      end if
      
      
      ! Compute residuals for each response variable
      do j1 = 1, J
        residual = y(j1, icount) - dot_product(xx(:, icount), bbeta(:, j1)) - tau(j1, j1) * Rgam(icount, j1)
        tv(j1) = residual
      end do
      
      ! Compute outer product of residuals
      do j1 = 1, J
        do j2 = 1, J
          tv1(j1, j2) = tv(j1) * tv(j2)
        end do
      end do
      
      
      ! Construct precision matrix for Wishart distribution
      ptk = ny(icount) * tv1 + (ny(icount) - 1.0d0) * tempvrv + (ve - real(J + 1, 8)) * sig(1, :, :)
      
      
      det1 = FindDet(tol, J, ptk)
      if (det1 <= 0) then
        print *, 'Warning: Invalid ptk at icount=', icount
        print *, 'kk=', kk, 'icount=', icount, 'tempvrv=', tempvrv
        print *, 'ny=', ny(icount), 've=', ve, 'df=', df
        print *, 'tv1=', tv1
        print *, 'ptk=', ptk
      end if
      
      
      ! Invert to get scale matrix
      call inverse(ptk, ptk_inv, J)
      
      ! Sample from Wishart distribution
      call wishart_sample(J, int(df), ptk_inv, wis, iseed)
      
      ! Compute covariance matrix (inverse of Wishart sample)
      call inverse(wis, sigtest, J)
      
      ! Check determinant for positive definiteness
      det1 = FindDet(tol, J, wis)
      if (det1 > 0.0d0) then
        sigk(icount,:,:) = sigtest
        SigInv(icount,:,:) = wis
        sigma(icount,1) = sigtest(1,1)
        sigma(icount,2) = sigtest(2,2)
        sigma(icount,3) = sigtest(3,3)
        sigma(icount,4) = sigtest(1,2)
        sigma(icount,5) = sigtest(1,3)
        sigma(icount,6) = sigtest(2,3)
      else
        print *, 'Warning: Invalid Wishart sample at icount=', icount
      end if
    end do
  end do

end subroutine Gsigma
        
subroutine Gsiggrp(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, nx = 10, nT = 10, J = 3, ngrp = 1
  real(8), parameter :: v0 = 5.0d0  ! Prior degrees of freedom
  
  ! Declare variables
  integer :: grp(ns), narms(K), icount, i1, kk, idim, l, l2, l3, l4
  real(8) :: ve, df, sig(ngrp, J, J), SigInv(ns, J, J), grp_n(ngrp)
  real(8) :: delta(ngrp, J), gamma(ngrp, J), sig0(J, J)
  real(8) :: sigtemp(J, J), sigtemp_inv(J, J), wis(J, J), sigtemp1(J, J)
  
  ! Common blocks
  common /vecnarms/ narms
  common /vecdelta/ delta
  common /vecgamma/ gamma
  common /vecSig/ sig
  common /vecSigInv/ SigInv
  common /vecgrp/ grp
  common /vecgrpn/ grp_n
  common /vecSigo/ sig0
  common /numve/ ve

  ! External subroutines
  external :: inverse, wishart_sample

  ! Loop over groups (only 1 group in this implementation)
  do i1 = 1, ngrp
    icount = 0
    sigtemp = 0.0d0
    
    ! Calculate degrees of freedom for Wishart distribution
    df = ns * ve + v0
    
    ! Accumulate precision matrices across all obs
    do kk = 1, K
      idim = narms(kk)  ! Number of arms in current trial
      do l = 1, idim
        icount = icount + 1
        sigtemp = sigtemp + SigInv(icount, :, :)
      end do
    end do
    
    ! Construct scale matrix for Wishart distribution
    sigtemp1 = (ve - real(J + 1, 8)) * sigtemp + sig0
    
    ! Compute inverse of scale matrix
    call inverse(sigtemp1, sigtemp_inv, J)
    
    ! Sample from Wishart distribution
    call wishart_sample(J, int(df), sigtemp_inv, wis, iseed)
    
    ! Store sampled covariance matrix for the group
    sig(i1, :, :) = wis
  end do

  ! Extract variances and covariances for storage
  do i1 = 1, ngrp
    ! Store variances (diagonal elements)
    do l2 = 1, J
      delta(i1, l2) = sig(i1, l2, l2)
    end do
    
    ! Store covariances (off-diagonal elements)
    do l3 = 1, J - 1
      do l4 = l3 + 1, J
        gamma(i1, l3 + l4 - 2) = sig(i1, l3, l4)
      end do
    end do
  end do

end subroutine Gsiggrp

subroutine Gv(iseed)
  implicit real(8) (a-h, o-z)
  integer, intent(inout) :: iseed  ! Random number seed (input/output)
  integer, parameter :: J = 3
  
  ! Declare variables
  real(8) :: v, ve
  real(8) :: start, ommax, plmax, step, step1, e1, sigmaa
  real(8) :: omkstar, rat1, r1, r2
  real(8) :: cl(5, 3), dl(5, 1), al(3, 3), alinv(3, 3), fl(3, 1), el(3, 1)
  integer :: nopt, konvge, kcount, icount1, numres, ifault1, ial, i, ni
  
  ! Common blocks
  common /numv/ v
  common /numve/ ve

  ! External functions
  real(8), external :: aNloglike, r8_normal_01_sample
  external :: inverse, r8_uniform_01_sample  ! Matrix inversion subroutine

  ! Initialize optimization parameters
  nopt = 1             ! Number of optimization variables
  reqmin = 1.0d-10      ! Optimization tolerance
  konvge = 5           ! Convergence check interval
  kcount = 1000        ! Max optimization iterations
  step = 0.1d0         ! Initial step size
  start = v            ! Starting value for optimization

  ! Find mode of log-likelihood using Nelder-Mead optimization
  call nelmin(aNloglike, nopt, start, ommax, plmax, reqmin, &
              step, konvge, kcount, icount1, numres, ifault1)

  ! Refine optimization with quadratic approximation
  step1 = 0.2d0
  eps   = 1.0d-10
  iter  = 0
  refine: do
    iter = iter + 1
    ! Evaluate log-likelihood at points around mode
    do i = 1, 5
      e1 = real(i - 3, 8)
      cl(i, 1) = (ommax + e1 * step1)**2
      cl(i, 2) = ommax + e1 * step1
      cl(i, 3) = 1.0d0
      dl(i, 1) = aNloglike(ommax + e1 * step1)
    end do
    
    ! Check if refinement is needed (any point has lower likelihood than mode)
    if (any( (/ dl(1,1), dl(2,1), dl(4,1), dl(5,1) /) <= plmax + eps)) then
      step1 = max(step1 * 0.80d0, 1.0d-4)
      if (iter < 30) cycle refine
    end if
    exit refine
  end do refine

  ! Fit quadratic model to log-likelihood
  al = matmul(transpose(cl), cl)
  el = matmul(transpose(cl), dl)
  ial = 3
  call inverse(al, alinv, ial)
  fl = matmul(alinv, el)
  
  ! Compute proposal distribution parameters
  sigmaa = 1.0d0 / sqrt(2.0d0 * fl(1, 1))
  
  ! Generate proposal from normal distribution
  r1 = r8_normal_01_sample(iseed)
  omkstar = ommax + sigmaa * r1
  
  ! Compute Metropolis acceptance ratio
  rat1 = -aNloglike(omkstar) + aNloglike(start) &
         - 0.5d0 * (start - ommax)**2 / sigmaa**2 &
         + 0.5d0 * (omkstar - ommax)**2 / sigmaa**2
  
  ! Generate uniform random number for Metropolis acceptance
  call r8_uniform_01_sample(r2, iseed)
  
  ! Metropolis acceptance step
  if (rat1 >= 0.0d0 .or. log(r2) <= rat1) then
    v = omkstar
  else
    v = start
  end if
  
  ! Transform v to ve (degrees of freedom parameter)
  ve = exp(v) + (J+1)

end subroutine Gv

real(8) function aNloglike(iv)
  implicit real(8) (a-h, o-z)
  real(8), intent(in) :: iv  ! Input parameter (log-scale)
  
  ! Model dimensions and parameters
  integer, parameter :: ns = 103, K = 40, ngrp = 1, J = 3
  real(8), parameter :: c04 = 10.0d0, epi = 3.141592741012573d0
  real(8), parameter :: tol = 2.220446049250313d-14  ! Tolerance for determinant calculation
  
  ! Declare variables
  integer :: ids(ns), iarm(ns), narms(K), icount, kk, idim, l, ifault
  real(8) :: SigInv(ns, J, J), sig(ngrp, J, J), sig1(J, J), sig2(J, J)
  real(8) :: sigtemp(J, J), sigtemp1(J, J), z, tz, det1, det2
  real(8) :: s1, s2, s4, ll, lg1, lg2, lg3, trace_term
  real(8) :: log_gamma1, log_gamma2, log_gamma3
  
  ! Common blocks
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecSigInv/ SigInv
  common /vecSig/ sig

  ! External functions
  real(8), external :: alngam  ! Logarithm of gamma function
  real(8), external :: FindDet    ! Determinant calculation function

  ! Transform input parameter
  z = exp(iv)
  tz = z + (J+1)  ! Degrees of freedom parameter
  
  ! Compute prior term
  s4 = -0.5d0 * iv**2 / c04
  
  ! Get group-level covariance matrix
  sig1 = sig(1, :, :)
  
  ! Compute determinant of group-level covariance
  det1 = FindDet(tol, J, sig1)
  
  ! Initialize accumulators
  icount = 0
  sigtemp = 0.0d0
  s1 = 0.0d0
  
  ! Precompute log-gamma terms
  log_gamma1 = alngam(tz/2.0d0, ifault)
  log_gamma2 = alngam((tz - 1.0d0)/2.0d0, ifault)
  log_gamma3 = alngam((tz - 2.0d0)/2.0d0, ifault)
  
  ! Loop over trials and arms
  do kk = 1, K
    idim = narms(kk)  ! Number of arms in current trial
    do l = 1, idim
      icount = icount + 1
      
      ! Get subject-level precision matrix
      sig2 = SigInv(icount, :, :)
      
      ! Compute determinant of precision matrix
      det2 = FindDet(tol, J, sig2)
      
      ! Accumulate precision matrices
      sigtemp = sigtemp + sig2
      
      ! Accumulate log-determinant terms
      s1 = s1 + 0.5d0 * (tz + real(J + 1, 8)) * log(det2)
    end do
  end do
  
  ! Compute trace term: tr(sig1 * sigtemp)
  sigtemp1 = matmul(sig1, sigtemp)
  trace_term = sigtemp1(1, 1) + sigtemp1(2, 2) + sigtemp1(3, 3)
  
  ! Compute log-likelihood components
  s2 = 0.5d0 * tz * ns * log(det1) &  ! Group covariance term
       + s1 &                          ! Subject precision terms
       - 0.5d0 * (tz - real(J + 1, 8)) * trace_term &     ! Trace term
       + 0.5d0 * tz * ns * real(J, 8) * log(tz - real(J + 1, 8)) &  ! DF scaling
       - 0.5d0 * tz * real(J, 8) * ns * log(2.0d0) &  ! Constant factor
       - ns * (log_gamma1 + log_gamma2 + log_gamma3) &  ! Gamma function terms
       + s4  ! Prior term
  
  ! Return negative log-likelihood
  aNloglike = -s2

end function aNloglike

subroutine vars1()
  implicit real*8 (a-h,o-z)
  parameter (ns = 103, K = 40, nx = 10, nT = 10, nTT = 11, J = 3)
  parameter (nbeta = 10, nb = 4, nbbeta=nbeta+nb+nT)
  
  ! Declare variables
  integer ids(ns), iarm(ns), narms(K), iobs(K)
  integer icount, idim, tv2, kk, jj, l1, l2
  real*8 ny(ns), xx(nbbeta, ns), y(J, ns)
  real*8 bbeta(nbbeta, J), SigInv(ns, J, J), InvSig_t(J, J)
  real*8 E(ns, nTT), EO(J, 12), EO2(J, 6), EO3(J, 9), EO4(J, 12)
  real*8 tau(J, J), tauk(ns, J, J)
  real*8 tv1_all(K, 12, 12), tv7_234(K, 12), tv7_all(ns, J), tv7(J)
  
  ! Common blocks
  common /vecy/ y
  common /vecny/ ny
  common /vecxx/ xx
  common /vecids/ ids
  common /veciarm/ iarm
  common /vecnarms/ narms
  common /vecobs/ iobs
  common /vectauk/ tauk
  common /vecSigInv/ SigInv
  common /vecE/ E
  common /vecbbeta/ bbeta
  common /vectv/ tv1_all, tv7_234, tv7_all

  ! Initialize accumulators
  icount = 0
  tv1_all = 0.0d0
  tv7_234 = 0.0d0

  ! Loop over trials
  do kk = 1, K
    idim = narms(kk)
    tv2 = 3 * idim
    
    ! Loop over arms within trial
    do jj = 1, idim
      icount = icount + 1
      
      ! Initialize EO matrices
      EO2 = 0.0d0
      EO3 = 0.0d0
      EO4 = 0.0d0
      
      ! Get tau matrix for current subject
      tau = tauk(icount, :, :)
      
      ! Initialize EO matrix
      do l1 = 1, J
        do l2 = 1, 12
          EO(l1, l2) = 0.0d0
        end do
      end do
      
      ! Set up EO matrix for current arm
      EO(1, jj) = 1.0d0
      EO(2, jj + idim) = 1.0d0
      EO(3, jj + 2 * idim) = 1.0d0
      
      ! Compute residuals: y - Xβ
      do ii = 1, J 
        tv7_all(icount, ii) = y(ii, icount) - dot_product(xx(:, icount), bbeta(:, ii))
      end do
      tv7 = tv7_all(icount, :)
      
      ! Get inverse Sigma for current subject
      InvSig_t = SigInv(icount, :, :)
      
      ! Update precision matrix component
      tv1_all(kk, 1:tv2, 1:tv2) = tv1_all(kk, 1:tv2, 1:tv2) + ny(icount) * &
        matmul(matmul(matmul(transpose(EO(:, 1:tv2)), tau), InvSig_t), &
              matmul(tau, EO(:, 1:tv2)))
      
      ! Update mean vector component based on number of arms
      if (idim == 2) then
        EO2(1, jj) = 1.0d0
        EO2(2, jj + idim) = 1.0d0
        EO2(3, jj + 2 * idim) = 1.0d0
        tv7_234(kk, 1:tv2) = tv7_234(kk, 1:tv2) + ny(icount) * &
          matmul(matmul(matmul(transpose(EO2), tau), InvSig_t), tv7)
      else if (idim == 3) then
        EO3(1, jj) = 1.0d0
        EO3(2, jj + idim) = 1.0d0
        EO3(3, jj + 2 * idim) = 1.0d0
        tv7_234(kk, 1:tv2) = tv7_234(kk, 1:tv2) + ny(icount) * &
          matmul(matmul(matmul(transpose(EO3), tau), InvSig_t), tv7)
      else if (idim == 4) then
        EO4(1, jj) = 1.0d0
        EO4(2, jj + idim) = 1.0d0
        EO4(3, jj + 2 * idim) = 1.0d0
        tv7_234(kk, 1:tv2) = tv7_234(kk, 1:tv2) + ny(icount) * &
          matmul(matmul(matmul(transpose(EO4), tau), InvSig_t), tv7)
      end if
    end do
  end do
end subroutine vars1