       subroutine prhoTorho(nT, pRho, Rho)
c      Hao Li, Uconn
c      March 28, 2016
       implicit double precision (a-h,o-z)
       double precision  Rho(nT,nT), pRho(nT,nT)
       double precision rvec1(nT-2), rvec3(nT-2), rr11, rr13, rr33
       double precision rmat1(1,1), rmatinv1(1,1)
       double precision rmat2(2,2), rmatinv2(2,2)
       double precision rmat3(3,3), rmatinv3(3,3)
       double precision rmat4(4,4), rmatinv4(4,4)
       double precision rmat5(5,5), rmatinv5(5,5)
       double precision rmat6(6,6), rmatinv6(6,6)
       double precision rmat7(7,7), rmatinv7(7,7)
       double precision rmat8(8,8), rmatinv8(8,8)
       double precision rmat9(9,9), rmatinv9(9,9)
c
        do j=1,nT
         Rho(j,j)=pRho(j,j)
        enddo
c
        do j=1,nT-1
         Rho(j,j+1)=pRho(j,j+1)
         Rho(j+1,j)=Rho(j,j+1)
        enddo
c
        do 555 iL=2,nT-1
         do 556 iR=1,nT-iL
          do i1=1,nT-2
           rvec1(i1)=0.0d0
           rvec3(i1)=0.0d0
          enddo
          do i1=1,iL-1
           rvec1(i1)=Rho(iR,iR+i1)
           rvec3(i1)=Rho(iR+i1,iR+iL)
          enddo
          rr11=0.0d0
          rr13=0.0d0
          rr33=0.0d0
c.......iL=2         
         if (iL .eq. 2) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat1(i1,j1)=0.0d0
           rmatinv1(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat1(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat1, rmatinv1, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv1(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv1(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv1(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=3
         if (iL .eq. 3) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat2(i1,j1)=0.0d0
           rmatinv2(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat2(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat2, rmatinv2, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv2(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv2(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv2(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=4
         if (iL .eq. 4) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat3(i1,j1)=0.0d0
           rmatinv3(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat3(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat3, rmatinv3, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv3(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv3(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv3(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=5
         if (iL .eq. 5) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat4(i1,j1)=0.0d0
           rmatinv4(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat4(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat4, rmatinv4, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv4(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv4(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv4(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=6
         if (iL .eq. 6) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat5(i1,j1)=0.0d0
           rmatinv5(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat5(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat5, rmatinv5, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv5(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv5(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv5(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=7
         if (iL .eq. 7) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat6(i1,j1)=0.0d0
           rmatinv6(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat6(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat6, rmatinv6, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv6(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv6(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv6(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=8
         if (iL .eq. 8) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat7(i1,j1)=0.0d0
           rmatinv7(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat7(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat7, rmatinv7, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv7(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv7(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv7(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif         
c.......iL=9
         if (iL .eq. 9) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat8(i1,j1)=0.0d0
           rmatinv8(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat8(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat8, rmatinv8, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv8(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv8(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv8(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c.......iL=10
         if (iL .eq. 10) then           
          do i1=1,iL-1
           do j1=1,iL-1
           rmat9(i1,j1)=0.0d0
           rmatinv9(i1,j1)=0.0d0
           enddo
          enddo
          do i1=1,iL-1
           do j1=1,iL-1
            rmat9(i1,j1)=Rho(iR+i1,iR+j1)
           enddo
          enddo
          call inverse(rmat9, rmatinv9, iL-1)
          do i1=1,iL-1
           do j1=1,iL-1
            rr11=rr11 + rvec1(i1)*rmatinv9(i1,j1)*rvec1(j1)
            rr13=rr13 + rvec1(i1)*rmatinv9(i1,j1)*rvec3(j1)
            rr33=rr33 + rvec3(i1)*rmatinv9(i1,j1)*rvec3(j1)
           enddo
          enddo
         endif
c
        Rho(iR,iR+iL) = rr13 +
     1                pRho(iR,iR+iL)*dsqrt((1.0d0-rr11)*(1.0d0-rr33))
        Rho(iR+iL,iR) = Rho(iR,iR+iL)

 556    continue
 555   continue

        return
        end
c
      subroutine wishart_sample ( m, df, sigma, a ,iseed)

c*********************************************************************72
c
cc WISHART_SAMPLE samples the Wishart distribution.
c
c  Discussion:
c
c    This function requires functions from the PDFLIB and RNGLIB libraries.
c
c    The "initialize()" function from RNGLIB must be called before using
c    this function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Patrick Odell, Alan Feiveson,
c    A numerical procedure to generate a sample covariance matrix,
c    Journal of the American Statistical Association,
c    Volume 61, Number 313, March 1966, pages 199-203.
c
c    Stanley Sawyer,
c    Wishart Distributions and Inverse-Wishart Sampling,
c    Washington University,
c    30 April 2007, 12 pages.
c
c  Parameters:
c
c    Input, integer M, the order of the matrix.
c
c    Input, integer DF, the number of degrees of freedom.
c    M <= DF.
c
c    Input, double precision SIGMA(M,M), the covariance matrix, which should be 
c    a symmetric positive definite matrix.
c
c    Output, double precision A(M,M), the sample matrix from 
c    the Wishart distribution.
c
      implicit none

      integer m

      double precision a(m,m)
      double precision au(m,m)
      double precision aur(m,m)
      integer df
      integer flag
      double precision r(m,m)
      double precision sigma(m,m)
      integer iseed

      if ( df .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WISHART_SAMPLE - Fatal error!'
        write ( *, '(a,i6)' ) '  DF = ', df
        write ( *, '(a,i4)' ) ' < M = ', m
        stop 1
      end if
c
c  Get R, the upper triangular Cholesky factor of SIGMA.
c
      call cholesky ( m, sigma, r, flag )

      if ( flag .ne. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'WISHART_SAMPLE - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Unexpected error return from cholesky.'
        write ( *, '(a,i4)' ) '  FLAG = ', flag
        stop 1
      end if
c
c  Get AU, a sample from the unit Wishart distribution.
c
      call wishart_unit_sample ( m, df, au,iseed )
c
c  Construct the matrix A = R' * AU * R.
c
      call r8mat_mm ( m, m, m, au, r, aur )

      call r8mat_mtm ( m, m, m, r, aur, a )

      return
      end
c
      subroutine wishart_unit_sample ( m, df, a ,iseed)

c*********************************************************************72
c
cc WISHART_UNIT_SAMPLE samples the unit Wishart distribution.
c
c  Discussion:
c
c    This function requires functions from the PDFLIB and RNGLIB libraries.
c
c    The "initialize()" function from RNGLIB must be called before using
c    this function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 October 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Patrick Odell, Alan Feiveson,
c    A numerical procedure to generate a sample covariance matrix,
c    Journal of the American Statistical Association,
c    Volume 61, Number 313, March 1966, pages 199-203.
c
c    Stanley Sawyer,
c    Wishart Distributions and Inverse-Wishart Sampling,
c    Washington University,
c    30 April 2007, 12 pages.
c
c  Parameters:
c
c    Input, integer M, the order of the matrix.
c
c    Input, integer DF, the number of degrees of freedom.
c    M <= DF.
c
c    Output, double A(M,M), the sample matrix from the 
c    unit Wishart distribution.
c
      implicit none

      integer m

      double precision a(m,m)
      double precision c(m,m)
      integer df
      double precision df_chi
      integer i
      integer j
      double precision r8_chi_sample
      double precision r8_normal_01_sample
      integer iseed

      if ( df .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WISHART_UNIT_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
        stop 1
      end if

      do i = 1, m

        do j = 1, i - 1
          c(i,j) = 0.0D+00
        end do

        df_chi = dble ( df + 1 - i )
        c(i,i) = sqrt ( r8_chi_sample ( df_chi,iseed ) )

        do j = i + 1, m
          c(i,j) = r8_normal_01_sample (iseed )
        end do

      end do

      call r8mat_mtm ( m, m, m, c, c, a )

      return
      end
c
      subroutine r8mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MM multiplies two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 July 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n2)
      double precision b(n2,n3)
      double precision c(n1,n3)
      double precision c1(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c1(i,j) = 0.0D+00
          do k = 1, n2
            c1(i,j) = c1(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      do j = 1, n3
        do i = 1, n1
          c(i,j) = c1(i,j)
        end do
      end do

      return
      end
c
      subroutine r8mat_mtm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MTM multiplies computes C = A' * B for two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = matmul ( transpose ( A(1:N2,1:N1) ), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N2,N1), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A' * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n2,n1)
      double precision b(n2,n3)
      double precision c(n1,n3)
      double precision c1(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c1(i,j) = 0.0D+00
          do k = 1, n2
            c1(i,j) = c1(i,j) + a(k,i) * b(k,j)
          end do
        end do
      end do

      call r8mat_copy ( n1, n3, c1, c )

      return
      end
c
      function r8_normal_01_sample (iseed )

c*********************************************************************72
c
cc R8_NORMAL_01_SAMPLE returns a unit pseudonormal R8.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    The Box-Muller method is used, which is efficient, but
c    generates two values at a time.
c
c    Typically, we would use one value and save the other for the next call.
c    However, the fact that this function has saved memory makes it difficult
c    to correctly handle cases where we want to re-initialize the code,
c    or to run in parallel.  Therefore, we will instead use the first value
c    and DISCARD the second.
c
c    EFFICIENCY must defer to SIMPLICITY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_NORMAL_01_SAMPLE, a sample of the standard
c    normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01_sample
      double precision x
      integer iseed


      call r8_uniform_01_sample(r1, iseed)
      call r8_uniform_01_sample(r2, iseed)

      x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )

      r8_normal_01_sample = x

      return
      end
c
      function r8_chi_sample ( df ,iseed)

c*********************************************************************72
c
cc R8_CHI_SAMPLE generates a Chi-Square random deviate.
c
c  Discussion:
c
c    This procedure generates a random deviate from the chi square distribution
c    with DF degrees of freedom random variable.
c
c    The algorithm exploits the relation between chisquare and gamma.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    This version by John Burkardt.
c
c  Parameters:
c
c    Input, double precision DF, the degrees of freedom.
c    0.0 .lt. DF.
c
c    Output, double precision R8_CHI_SAMPLE, a random deviate
c    from the distribution.
c
      implicit none

      double precision arg1
      double precision arg2
      double precision df
      double precision r8_chi_sample
      double precision r8_gamma_sample
      integer iseed

      if ( df .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CHI_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  DF .le. 0.'
        write ( *, '(a,g14.6)' ) '  Value of DF: ', df
        stop
      end if

      arg1 = 1.0D+00
      arg2 = df / 2.0D+00

      r8_chi_sample = 2.0D+00 * r8_gamma_sample ( arg1, arg2 ,iseed)

      return
      end
c
      function r8_gamma_sample ( a, r ,iseed)

c*********************************************************************72
c
cc R8_GAMMA_SAMPLE generates a Gamma random deviate.
c
c  Discussion:
c
c    This procedure generates random deviates from the gamma distribution whose
c    density is (A^R)/Gamma(R) * X^(R-1) * Exp(-A*X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Generating Gamma Variates by a Modified Rejection Technique,
c    Communications of the ACM,
c    Volume 25, Number 1, January 1982, pages 47-54.
c
c    Joachim Ahrens, Ulrich Dieter,
c    Computer Methods for Sampling from Gamma, Beta, Poisson and
c    Binomial Distributions,
c    Computing,
c    Volume 12, Number 3, September 1974, pages 223-246.
c
c  Parameters:
c
c    Input, double precision A, the rate parameter.
c
c    Input, double precision R, the shape parameter.
c
c    Output, double precision R8_GAMMA_SAMPLE, a random deviate
c    from the distribution.
c
      implicit none

      double precision a
      double precision r
      double precision r8_gamma_sample
      double precision r8_gamma_01_sample
      integer iseed

      r8_gamma_sample = r8_gamma_01_sample ( r ,iseed) / a

      return
      end
c
      subroutine r8mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc R8MAT_COPY copies an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision A1(M,N), the matrix to be copied.
c
c    Output, double precision A2(M,N), a copy of the matrix.
c
      implicit none

      integer m
      integer n

      double precision a1(m,n)
      double precision a2(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
        end do
      end do

      return
      end
c
      subroutine cholesky ( n, a, c, flag )

!*********************************************************************72
!
!! cholesky: upper Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is an upper triangular matrix R such that:
!
!      A = R * R'
!
!    The lower Cholesky factor is a lower triangular matrix L such that
!
!      A = L * L'
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of
!    the matrix A.
!
!    Input, double precision A(N,N), the N by N matrix.
!
!    Output, double precision C(N,N), the N by N upper triangular
!    Cholesky factor.
!
!    Output, integer FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!    2, the matrix is not nonnegative definite.
!
      implicit none

      integer n

      double precision a(n,n)
      double precision c(n,n)
      integer flag
      integer i
      integer j
      integer k
      double precision sum2
      double precision tol

      flag = 0

      do j = 1, n
        do i = 1, n
          c(i,j) = a(i,j)
        end do
      end do

      do j = 1, n

        do i = 1, j - 1
          c(j,i) = 0.0D+00
       end do

        do i = j, n

          sum2 = c(i,j)
          do k = 1, j - 1
            sum2 = sum2 - c(k,j) * c(k,i)
          end do

          if ( i .eq. j ) then
            if ( sum2 .le. 0.0D+00 ) then
              flag = 1
              return
            else
              c(j,i) = sqrt ( sum2 )
            end if
          else
            if ( c(j,j) .ne. 0.0D+00 ) then
              c(j,i) = sum2 / c(j,j)
            else
              c(j,i) = 0.0D+00
            end if
          end if

        end do

      end do

      return
      end
c
      function r8_gamma_01_sample ( a ,iseed)

c*********************************************************************72
c
cc R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.
c
c  Discussion:
c
c    This procedure corresponds to algorithm GD in the reference.
c
c    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 April 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Generating Gamma Variates by a Modified Rejection Technique,
c    Communications of the ACM,
c    Volume 25, Number 1, January 1982, pages 47-54.
c
c  Parameters:
c
c    Input, double precision A, the shape parameter of the standard gamma
c    distribution. 0 .lt. A.
c
c    Output, double precision R8_GAMMA_01_SAMPLE, a random deviate
c    from the distribution.
c
      implicit none

      double precision a
      double precision a1
      parameter ( a1 =  0.3333333D+00 )
      double precision a2
      parameter ( a2 = -0.2500030D+00 )
      double precision a3
      parameter ( a3 =  0.2000062D+00 )
      double precision a4
      parameter ( a4 = -0.1662921D+00 )
      double precision a5
      parameter ( a5 =  0.1423657D+00 )
      double precision a6
      parameter ( a6 = -0.1367177D+00 )
      double precision a7
      parameter ( a7 =  0.1233795D+00 )
      double precision b
      double precision c
      double precision d
      double precision e
      double precision, parameter :: e1 = 1.0D+00
      double precision, parameter :: e2 = 0.4999897D+00
      double precision, parameter :: e3 = 0.1668290D+00
      double precision, parameter :: e4 = 0.0407753D+00
      double precision, parameter :: e5 = 0.0102930D+00
      double precision p
      double precision q
      double precision q0
      double precision, parameter :: q1 =  0.04166669D+00
      double precision, parameter :: q2 =  0.02083148D+00
      double precision, parameter :: q3 =  0.00801191D+00
      double precision, parameter :: q4 =  0.00144121D+00
      double precision, parameter :: q5 = -0.00007388D+00
      double precision, parameter :: q6 =  0.00024511D+00
      double precision, parameter :: q7 =  0.00024240D+00
      double precision r
      double precision r8_exponential_01_sample
      double precision r8_gamma_01_sample
      double precision r8_normal_01_sample
      double precision s
      double precision s2
      double precision si
      double precision, parameter :: sqrt32 = 5.656854D+00
      double precision t
      double precision u
      double precision v
      double precision w
      double precision x
      double precision ur
      integer iseed

      if ( 1.0D+00 .le. a ) then

        s2 = a - 0.5D+00
        s = sqrt ( s2 )
        d = sqrt32 - 12.0D+00 * s
c
c  Immediate acceptance.
c
        t = r8_normal_01_sample (iseed )
        x = s + 0.5D+00 * t
        r8_gamma_01_sample = x * x

        if ( 0.0D+00 .le. t ) then
          return
        end if
c
c  Squeeze acceptance.
c
        call r8_uniform_01_sample (u,iseed)
        if ( d * u .le. t * t * t ) then
          return
        end if

        r = 1.0D+00 / a
        q0 = (((((( q7 
     &    * r + q6 ) 
     &    * r + q5 ) 
     &    * r + q4 ) 
     &    * r + q3 ) 
     &    * r + q2 ) 
     &    * r + q1 ) 
     &    * r
c
c  Approximation depending on size of parameter A.
c
        if ( 13.022D+00 .lt. a ) then
          b = 1.77D+00
          si = 0.75D+00
          c = 0.1515D+00 / s
        else if ( 3.686D+00 .lt. a ) then
          b = 1.654D+00 + 0.0076D+00 * s2
          si = 1.68D+00 / s + 0.275D+00
          c = 0.062D+00 / s + 0.024D+00
        else
          b = 0.463D+00 + s + 0.178D+00 * s2
          si = 1.235D+00
          c = 0.195D+00 / s - 0.079D+00 + 0.16D+00 * s
        end if
c
c  Quotient test.
c
        if ( 0.0D+00 .lt. x ) then

          v = 0.5D+00 * t / s

          if ( 0.25D+00 .lt. abs ( v ) ) then
            q = q0 - s * t + 0.25D+00 * t * t 
     &        + 2.0D+00 * s2 * log ( 1.0D+00 + v )
          else
            q = q0 + 0.5D+00 * t * t * (((((( a7 
     &        * v + a6 ) 
     &        * v + a5 ) 
     &        * v + a4 ) 
     &        * v + a3 ) 
     &        * v + a2 ) 
     &        * v + a1 ) 
     &        * v
          end if

          if ( log ( 1.0D+00 - u ) .le. q ) then
            return
          end if

        end if

10      continue

          e = r8_exponential_01_sample (iseed )
          call r8_uniform_01_sample (ur,iseed)
          u = 2.0D+00 * ur - 1.0D+00

          if ( 0.0D+00 .le. u ) then
            t = b + abs ( si * e )
          else
            t = b - abs ( si * e )
          end if
c
c  Possible rejection.
c
          if ( t .lt. -0.7187449D+00 ) then
            go to 10
          end if
c
c  Calculate V and quotient Q.
c
          v = 0.5D+00 * t / s

          if ( 0.25D+00 .lt. abs ( v ) ) then
            q = q0 - s * t + 0.25D+00 * t * t 
     &        + 2.0D+00 * s2 * log ( 1.0D+00 + v )
          else
            q = q0 + 0.5D+00 * t * t * (((((( a7 
     &        * v + a6 ) 
     &        * v + a5 ) 
     &        * v + a4 ) 
     &        * v + a3 ) 
     &        * v + a2 ) 
     &        * v + a1 ) 
     &        * v
          end if
c
c  Hat acceptance.
c
          if ( q .le. 0.0D+00 ) then
            go to 10
          end if

          if ( 0.5D+00 .lt. q ) then
            w = exp ( q ) - 1.0D+00
          else
            w = (((( e5 * q + e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q
          end if
c
c  May have to sample again.
c
          if ( c * abs ( u ) .le. w * exp ( e - 0.5D+00 * t * t ) ) then
            go to 20
          end if

        go to 10

20      continue

        x = s + 0.5D+00 * t
        r8_gamma_01_sample = x * x

        return
c
c  Method for A .lt. 1.
c
      else

        b = 1.0D+00 + 0.3678794D+00 * a

30      continue
          call r8_uniform_01_sample (ur,iseed)
          p = b * ur

          if ( p .lt. 1.0D+00 ) then

            r8_gamma_01_sample = exp ( log ( p ) / a )

            if ( r8_gamma_01_sample .le. 
     &        r8_exponential_01_sample ( iseed) ) then
              return
            end if

            go to 30

          end if

          r8_gamma_01_sample = - log ( ( b - p ) / a )

          if ( ( 1.0D+00 - a ) * log ( r8_gamma_01_sample ) .le. 
     &      r8_exponential_01_sample ( iseed) ) then
            go to 40
          end if

        go to 30

40      continue

      end if

      return
      end
c
      function r8_exponential_01_sample (seed )

c*********************************************************************72
c
cc R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 April 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
c
      implicit none

      double precision r
      double precision r8_exponential_01_sample
      integer seed
      external r8_uniform_01_sample

      call r8_uniform_01_sample ( r,seed)

      r8_exponential_01_sample = - log ( r )

      return
      END FUNCTION r8_exponential_01_sample
c=======================================================================
c  Below is the uniform random number generator subroutine
c=======================================================================
      subroutine r8_uniform_01_sample( r,seed)
c     seed: input/output
c     r   : output
      integer seed
      double precision r
      integer k, i4_huge
      parameter (i4_huge = 2147483647)
      
c      print *, 'uniform received seed=', seed
      if (seed .eq. 0) then
         write(*,*) 'Warning: seed was zero. Using default 123456789'
         seed = 123456789
      end if

      k = seed / 127773
      seed = 16807 * (seed - k * 127773) - k * 2836

      if (seed .lt. 0) seed = seed + i4_huge
c      print *, 'uniform output seed=', seed

      r = dble(seed) * 4.656612875D-10

      return
      end
c
      subroutine multinormal ( n, mu, r, x ,iseed)

c*********************************************************************72
c
cc multinormal samples a multivariate normal PDF.
c
c  Discussion:
c
c    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) =
c      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / det ( C )
c      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )
c
c    Here,
c
c      X is the argument vector of length N,
c      MU is the mean vector of length N,
c      C is an N by N positive definite symmetric covariance matrix.
c
c    The properties of C guarantee that it has an upper triangular
c    matrix R, the Cholesky factor, such that C = R' * R.  It is the
c    matrix R that is required by this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 June 2013
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision MU(N), the mean vector.
c
c    Input, double precision R(N,N), the upper triangular Cholesky
c    factor of the covariance matrix C.
c
c    Output, double precision X(N), a sample of the distribution.
c
      implicit none
      integer iseed
      integer n

      integer i
      integer j
      double precision mu(n)
      double precision r(n,n)
      double precision r8_normal_01_sample
      double precision x(n)
      double precision z(n)
c
c  Compute X = MU + R' * Z
c  where Z is a vector of standard normal variates.
c
      do j = 1, n
        z(j) = r8_normal_01_sample (iseed )
      end do

      do i = 1, n
        x(i) = mu(i)
        do j = 1, i
          x(i) = x(i) + r(j,i) * z(j)
        end do
      end do

      return
      end
c
          subroutine inverse(w,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
        implicit double precision (a-h,o-z)
        integer n
        double precision a(n,n), c(n,n), w(n,n)
        double precision L(n,n), U(n,n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k
        
        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0
        do i1=1,n
         do j1=1,n
          a(i1,j1)=w(i1,j1)
         enddo
        enddo
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do
        
        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do
        
        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        end subroutine inverse   

        Subroutine TSRGT(eps, n, A, it, C, Kp, Lp)
        implicit double precision (a-h,o-z)        
          double precision eps
          integer n,it
          double precision A(n,n), C(n,n)
          integer Kp(n),Lp(n) 
          double precision  po,t0
          C=A; it=1; k=1
          do while (it==1.and.k<n)
            po=C(k,k); lo=k; ko=k
            do i=k, n
              do j=k, n
                if (dabs(C(i,j))>dabs(po)) then
                  po=C(i,j); lo=i; ko=j
                end if
              end do
            end do
            Lp(k)=lo; Kp(k)=ko
            if (dabs(po)<eps) then
              it=0
            else
              if (lo.ne.k) then
                do j=k, n
                  t0=C(k,j); C(k,j)=C(lo,j); C(lo,j)=t0
                end do
              end if
              if (ko.ne.k) then
                do i=1, n
                  t0=C(i,k); C(i,k)=C(i,ko); C(i,ko)=t0
                end do
              end if 
              do i=k+1, n
                C(i,k)=C(i,k)/po
                do j=k+1, n
                  C(i,j)=C(i,j)-C(i,k)*C(k,j)
                end do 
              end do
              k=k+1
            end if
          end do
          if (it==1.and.dabs(C(n,n))<eps)  it=0
          return
        End !TSRGT
        
        !The function DMGT returns the determinant of a real square matrix
        !A(n,n) by Gauss method with full pivoting.
        !----------------------------------------------------------------------------
        !  Input parameters:
        !  eps        precision (double precision)
        !  n          size of A matrix (integer)
        !  A          pointer to input real square matrix
        !  Output parameters:
        !  None
        !-----------------------------------------------------------------------------
        !The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
        !Output variables are it(integer), C(n,n), Kp(n) and Lp(n).
        !If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
        !at location i,j (j>=i) the corresponding element of the upper triangular matrix.
        !Tables Lp and Kp contain informations relative to exchanges of line or column
        !that occured during the process. For instance, the element number k of Lp is
        !an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
        !The number of exchanges of lines and columns is stored in integer L. the
        !determinant of A matrix is stored in d0 (double precision).
        !-----------------------------------------------------------------------------
        Function DMGT(eps, n, A)
        implicit double precision (a-h,o-z)        
        integer n
        double precision eps, A(n,n)
        double precision d0
c        
        double precision C(n,n)
        integer Kp(n), Lp(n)
        call TSRGT(eps,n,A,it,C,Kp,Lp)  !call triangularization subroutine
          if (it==0) then
            d0=0.d0  !matrix singular, det=0
          else       !matrix regular, det<>0
            d0=1.d0
            do k=1, n
              d0=d0*C(k,k)
            end do
            l=0
            do k=1, n-1
              if (Lp(k).ne.k)  l=l+1
              if (Kp(k).ne.k)  l=l+1
            end do
            if (MOD(l,2).ne.0) d0=-d0  !l is odd
          end if
          DMGT=d0   !return determinant
          return
        End
c
      function alngam ( xvalue, ifault )

c*********************************************************************72
c
cc ALNGAM computes the logarithm of the gamma function.
c
c  Modified:
c
c    30 March 1999
c
c  Author:
c
c    Allan Macleod
c    Modifications by John Burkardt
c
c  Reference:
c
c    Allan Macleod,
c    Algorithm AS 245,
c    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
c    Applied Statistics,
c    Volume 38, Number 2, 1989, pages 397-402.
c
c  Parameters:
c
c    Input, double precision XVALUE, the argument of the Gamma function.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c    1, XVALUE is less than or equal to 0.
c    2, XVALUE is too big.
c
c    Output, double precision ALNGAM, the logarithm of the gamma function of X.
c
      implicit none

      double precision alngam
      double precision alr2pi
      parameter ( alr2pi = 0.918938533204673D+00 )
      integer ifault
      double precision r1(9)
      double precision r2(9)
      double precision r3(9)
      double precision r4(5)
      double precision x
      double precision x1
      double precision x2
      double precision xlge
      parameter ( xlge = 5.10D+06 )
      double precision xlgst
      parameter ( xlgst = 1.0D+30 )
      double precision xvalue
      double precision y

      data r1 /
     &  -2.66685511495D+00,
     &  -24.4387534237D+00,
     &  -21.9698958928D+00,
     &   11.1667541262D+00,
     &   3.13060547623D+00,
     &   0.607771387771D+00,
     &   11.9400905721D+00,
     &   31.4690115749D+00,
     &   15.2346874070D+00 /

      data r2 /
     &  -78.3359299449D+00,
     &  -142.046296688D+00,
     &   137.519416416D+00,
     &   78.6994924154D+00,
     &   4.16438922228D+00,
     &   47.0668766060D+00,
     &   313.399215894D+00,
     &   263.505074721D+00,
     &   43.3400022514D+00 /

      data r3 /
     &  -2.12159572323D+05,
     &   2.30661510616D+05,
     &   2.74647644705D+04,
     &  -4.02621119975D+04,
     &  -2.29660729780D+03,
     &  -1.16328495004D+05,
     &  -1.46025937511D+05,
     &  -2.42357409629D+04,
     &  -5.70691009324D+02 /

      data r4 / 
     &   0.279195317918525D+00, 
     &   0.4917317610505968D+00,
     &   0.0692910599291889D+00, 
     &   3.350343815022304D+00,
     &   6.012459259764103D+00 /

      x = xvalue
      alngam = 0.0D+00
c
c  Check the input.
c
      if ( xlgst .le. x ) then
        ifault = 2
        return
      end if

      if ( x .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0
c
c  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
c
      if ( x .lt. 1.5D+00 ) then

        if ( x .lt. 0.5D+00 ) then

          alngam = - dlog ( x )
          y = x + 1.0D+00
c
c  Test whether X < machine epsilon.
c
          if ( y .eq. 1.0D+00 ) then
            return
          end if

        else

          alngam = 0.0D+00
          y = x
          x = ( x - 0.5D+00 ) - 0.5D+00

        end if

        alngam = alngam + x * ((((
     &      r1(5)   * y 
     &    + r1(4) ) * y 
     &    + r1(3) ) * y
     &    + r1(2) ) * y 
     &    + r1(1) ) / ((((
     &                y 
     &    + r1(9) ) * y 
     &    + r1(8) ) * y
     &    + r1(7) ) * y 
     &    + r1(6) )

        return

      end if
c
c  Calculation for 1.5 <= X < 4.0.
c
      if ( x .lt. 4.0D+00 ) then

        y = ( x - 1.0D+00 ) - 1.0D+00

        alngam = y * ((((
     &      r2(5)   * x 
     &    + r2(4) ) * x 
     &    + r2(3) ) * x 
     &    + r2(2) ) * x
     &    + r2(1) ) / ((((
     &                x 
     &    + r2(9) ) * x 
     &    + r2(8) ) * x 
     &    + r2(7) ) * x
     &    + r2(6) )
c
c  Calculation for 4.0 <= X < 12.0.
c
      else if ( x .lt. 12.0D+00 ) then

        alngam = ((((
     &      r3(5)   * x 
     &    + r3(4) ) * x 
     &    + r3(3) ) * x 
     &    + r3(2) ) * x 
     &    + r3(1) ) / (((( 
     &                x 
     &    + r3(9) ) * x 
     &    + r3(8) ) * x 
     &    + r3(7) ) * x 
     &    + r3(6) )
c
c  Calculation for X >= 12.0.
c
      else

        y = dlog ( x )
        alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

        if ( x .le. xlge ) then

          x1 = 1.0D+00 / x
          x2 = x1 * x1

          alngam = alngam + x1 * ( ( 
     &           r4(3)   * 
     &      x2 + r4(2) ) * 
     &      x2 + r4(1) ) / ( ( 
     &      x2 + r4(5) ) * 
     &      x2 + r4(4) )

        end if

      end if

      return
      end
c
      REAL*8 FUNCTION FindDet(tol,N, A)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N, I, J, K
      REAL*8 A(N,N), M(N,N), TEMP, MULT,tol
      INTEGER ISWAP
      LOGICAL SINGULAR

      FindDet = 1.0D0
      SINGULAR = .FALSE.
      ISWAP = 0

C     Make a local copy of A to preserve input
      DO I = 1, N
         DO J = 1, N
            M(I,J) = A(I,J)
         END DO
      END DO

C     Gaussian elimination to upper triangular form
      DO K = 1, N-1
         IF (ABS(M(K,K)) .LT. tol) THEN
            SINGULAR = .TRUE.
            DO I = K+1, N
               IF (ABS(M(I,K)) .GT. tol) THEN
C                 Swap rows K and I
                  DO J = 1, N
                     TEMP = M(K,J)
                     M(K,J) = M(I,J)
                     M(I,J) = TEMP
                  END DO
                  ISWAP = ISWAP + 1
                  SINGULAR = .FALSE.
                  EXIT
               END IF
            END DO
            IF (SINGULAR) THEN
               FindDet = 0.0D0
               RETURN
            END IF
         END IF

C        Eliminate entries below pivot
         DO I = K+1, N
            MULT = M(I,K) / M(K,K)
            DO J = K+1, N
               M(I,J) = M(I,J) - MULT * M(K,J)
            END DO
         END DO
      END DO

C     Multiply diagonals
      DO I = 1, N
         FindDet = FindDet * M(I,I)
      END DO

C     Adjust sign for number of swaps
      IF (MOD(ISWAP, 2) .EQ. 1) THEN
         FindDet = -FindDet
      END IF

      RETURN
      END