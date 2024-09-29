subroutine changepoints_(n,hA,hT,p,A,T,X,W,Xb,AT,ATi,ATs,AT1,AT1i,AT1s,&
critical,criticali,criticals,nboot,kbin,kernel,nh, umatrixA, umatrixT)

!!DEC$ ATTRIBUTES DLLEXPORT:: change.points
!!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'changepoints_' :: change.points
use iso_fortran_env
implicit none
integer(kind=int32) i,p,iboot,j,icont,imin,imax,Xmin(200),&
Xmax(200),nh,kernel,n,kbin,nboot

double precision A(n),T(n),X(n),W(n),hA,hT,Signi(5000),&
critical(kbin),vmax,criticalboot(kbin,nboot), &
criticali(kbin),criticals(kbin),Ati(kbin),Ats(kbin),At1i(kbin),&
At1s(kbin),xb(kbin),At(kbin),At1(kbin), umatrixA(n, nboot),&
umatrixT(n, nboot), VauxInt(n), VauxBoot(nboot), VauxKbin(kbin)
double precision,allocatable::PbA(:,:),PbT(:,:),pred(:),predA(:),PredT(:),&
Aboot(:),Tboot(:),Atboot(:,:),At1boot(:,:)
double precision,external::Bernoulli



allocate (pred(n),predA(n),PredT(n),Aboot(n),Tboot(n))


allocate(PbA(kbin,3),PbT(kbin,3),Atboot(kbin,nboot),At1boot(kbin,nboot))



call GRID(X,W,n,Xb,kbin)
call RFAST_H(X,A,W,n,hA,p,Xb,PbA,kbin,kernel,nh)
!hT=hA
call RFAST_H(X,T,W,n,hT,p,Xb,PbT,kbin,kernel,nh)






do i=1,kbin
 AT(i)=(PbA(i,1)-PbT(i,1))/(PbA(i,1)+PbT(i,1))
end do

call Interpola (Xb,PbA(1,1),kbin,X,PredA,VauxInt,n) ! Vaux
call Interpola (Xb,PbT(1,1),kbin,X,PredT,VauxInt,n) ! Vaux




do iboot=1,nboot
 do i=1,n
  Aboot(i)=Bernoulli(predA(i), umatrixA(i, iboot))
  Tboot(i)=Bernoulli(predT(i), umatrixT(i, iboot))
 end do
 call RFAST_H(X,Aboot,W,n,hA,p,Xb,PbA,kbin,kernel,nh)
 call RFAST_H(X,Tboot,W,n,hT,p,Xb,PbT,kbin,kernel,nh)
 do i=1, kbin
  ATboot(i,iboot)=(PbA(i,1)-PbT(i,1))/(PbA(i,1)+PbT(i,1))
 end do
end do


do i=1,kbin
 do j=1,nboot
  VauxBoot(j)=Atboot(i,j) ! Vaux
 end do
 call IC(At(i),VauxBoot,nboot,ATi(i),Ats(i))
end do




call Interpola (Xb,At,kbin,Xb,VauxKbin,At1,kbin) !Vaux
do iboot=1,nboot
 do i=1,kbin
  VauxKbin(i)=Atboot(i,iboot)
 end do
 call Interpola (Xb,VauxKbin,kbin,Xb,VauxKbin,At1boot(1,iboot),kbin)
end do

do i=1,kbin
 do j=1,nboot
  VauxBoot(j)=At1boot(i,j)
 end do
 call IC(At1(i),VauxBoot,nboot,AT1i(i),At1s(i))
end do

do i=1,kbin
 if (At1i(i).gt.0.or.At1s(i).lt.0) then
  Signi(i)=1
 else
  Signi(i)=0
 end if
end do

icont=0
do i=2,kbin
 if (Signi(i).eq.1.and.Signi(i-1).eq.0) then
  icont=icont+1
  xmin(icont)=1
  do j=i-1,1,-1
   if (Signi(j).eq.0) then
    xmin(icont)=j
    goto 11
   end if
  end do
11      continue
xmax(icont)=kbin
do j=i+1,kbin
 if (Signi(j).eq.0) then
  xmax(icont)=j
  goto 12
 end if
end do
12 continue
end if
end do

do i=1,icont
 vmax=-999
 do j=xmin(i),xmax(i)
  if (abs(at1(j)).gt.vmax) then
   critical(i)=Xb(j)
   vmax=abs(at1(j))
  end if
 end do
end do

do i=1,icont
 do iboot=1,nboot
  vmax=-999
   do j=xmin(i),xmax(i)
    if (abs(at1boot(j,iboot)).gt.vmax) then
     criticalboot(i,iboot)=Xb(j)
     vmax=abs(at1boot(j,iboot))
    end if
   end do
 end do
end do

do i=1,icont
 do j=1,nboot
  VauxBoot(j)=criticalboot(i,j)
 end do
 call IC(critical(i),VauxBoot,nboot,criticali(i),criticals(i))
 imin=xmin(i)
 imax=xmax(i)
 criticali(i)=max(criticali(i),xb(imin))
 criticals(i)=min(criticals(i),xb(imax))
end do




deallocate(PbT,At1boot,pred,PredT,Tboot) !atboot,aboot,predA,PbA,
end



















double precision function Bernoulli(p, u)
use iso_fortran_env
implicit none
double precision p, u
!REAL*8 u

!u=RAND()

Bernoulli=0
if (u.le.p) Bernoulli=1
end


integer function which_min(X,n)
use iso_fortran_env
implicit none
integer(kind=int32) n,i
double precision X(n),aux

aux=X(1)
which_min=1
do i=2,n
if (X(i).le.aux) then
aux=X(i)
which_min=i
end if
end do

end

subroutine IC(X0,X,nboot,li,ls)
use iso_fortran_env
implicit none
integer(kind=int32) nboot,nalfa,i
double precision X0,X(nboot),Dif(nboot),li,ls,alfa(3),Q(3)
do i=1,nboot
Dif(i)=X(i)-X0
end do

alfa(1)=0.025
alfa(2)=0.5
alfa(3)=0.975
nalfa=3
call quantile (Dif,nboot,alfa,nalfa,Q)


li=x0-Q(3)+Q(2)
ls=x0-Q(1)+Q(2)
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE lsq

!  Module for unconstrained linear least-squares calculations.
!  The algorithm is suitable for updating LS calculations as more
!  data are added.   This is sometimes called recursive estimation.
!  Only one dependent variable is allowed.
!  Based upon Applied Statistics algorithm AS 274.
!  Translation from Fortran 77 to Fortran 90 by Alan Miller.
!  A function, VARPRD, has been added for calculating the variances
!  of predicted values, and this uses a subroutine BKSUB2.

!  Version 1.14, 19 August 2002 - ELF90 compatible version
!  Author: Alan Miller
!  e-mail : amiller @ bigpond.net.au
!  WWW-pages: http://www.ozemail.com.au/~milleraj
!             http://users.bigpond.net.au/amiller/

!  Bug fixes:
!  1. In REGCF a call to TOLSET has been added in case the user had
!     not set tolerances.
!  2. In SING, each time a singularity is detected, unless it is in the
!     variables in the last position, INCLUD is called.   INCLUD assumes
!     that a new observation is being added and increments the number of
!     cases, NOBS.   The line:  nobs = nobs - 1 has been added.
!  3. row_ptr was left out of the DEALLOCATE statement in routine startup
!     in version 1.07.
!  4. In COV, now calls SS if rss_set = .FALSE.  29 August 1997
!  5. In TOLSET, correction to accomodate negative values of D.  19 August 2002

!  Other changes:
!  1. Array row_ptr added 18 July 1997.   This points to the first element
!     stored in each row thus saving a small amount of time needed to
!     calculate its position.
!  2. Optional parameter, EPS, added to routine TOLSET, so that the user
!     can specify the accuracy of the input data.
!  3. Cosmetic change of lsq_kind to dp (`Double precision')
!  4. Change to routine SING to use row_ptr rather than calculate the position
!     of first elements in each row.

!  The PUBLIC variables are:
!  dp       = a KIND parameter for the floating-point quantities calculated
!             in this module.   See the more detailed explanation below.
!             This KIND parameter should be used for all floating-point
!             arguments passed to routines in this module.

!  nobs    = the number of observations processed to date.
!  ncol    = the total number of variables, including one for the constant,
!            if a constant is being fitted.
!  r_dim   = the dimension of array r = ncol*(ncol-1)/2
!  vorder  = an integer vector storing the current order of the variables
!            in the QR-factorization.   The initial order is 0, 1, 2, ...
!            if a constant is being fitted, or 1, 2, ... otherwise.
!  initialized = a logical variable which indicates whether space has
!                been allocated for various arrays.
!  tol_set = a logical variable which is set when subroutine TOLSET has
!            been called to calculate tolerances for use in testing for
!            singularities.
!  rss_set = a logical variable indicating whether residual sums of squares
!            are available and usable.
!  d()     = array of row multipliers for the Cholesky factorization.
!            The factorization is X = Q.sqrt(D).R where Q is an ortho-
!            normal matrix which is NOT stored, D is a diagonal matrix
!            whose diagonal elements are stored in array d, and R is an
!            upper-triangular matrix with 1's as its diagonal elements.
!  rhs()   = vector of RHS projections (after scaling by sqrt(D)).
!            Thus Q'y = sqrt(D).rhs
!  r()     = the upper-triangular matrix R.   The upper triangle only,
!            excluding the implicit 1's on the diagonal, are stored by
!            rows.
!  tol()   = array of tolerances used in testing for singularities.
!  rss()   = array of residual sums of squares.   rss(i) is the residual
!            sum of squares with the first i variables in the model.
!            By changing the order of variables, the residual sums of
!            squares can be found for all possible subsets of the variables.
!            The residual sum of squares with NO variables in the model,
!            that is the total sum of squares of the y-values, can be
!            calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
!            is a constant, then rss(1) is the sum of squares of
!            (y - ybar) where ybar is the average value of y.
!  sserr   = residual sum of squares with all of the variables included.
!  row_ptr() = array of indices of first elements in each row of R.
!
!--------------------------------------------------------------------------

!     General declarations
use iso_fortran_env
IMPLICIT NONE

INTEGER(kind=int32), SAVE                :: nobs, ncol, r_dim
INTEGER(kind=int32), ALLOCATABLE, SAVE   :: vorder(:), row_ptr(:)
LOGICAL, SAVE                :: initialized = .false.,                  &
                                tol_set = .false., rss_set = .false.

! Note. dp is being set to give at least 12 decimal digit
!       representation of floating point numbers.   This should be adequate
!       for most problems except the fitting of polynomials.   dp is
!       being set so that the same code can be run on PCs and Unix systems,
!       which will usually represent floating-point numbers in `double
!       precision', and other systems with larger word lengths which will
!       give similar accuracy in `single precision'.

INTEGER(kind=int32), PARAMETER           :: dp = SELECTED_REAL_KIND(12,60)
double precision, ALLOCATABLE, SAVE :: d(:), rhs(:), r(:), tol(:), rss(:)
double precision, SAVE              :: zero = 0.0_dp, one = 1.0_dp, vsmall
double precision, SAVE              :: sserr, toly

PUBLIC                       :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                initialized, tol_set, rss_set,          &
                                d, rhs, r, tol, rss, sserr
PRIVATE                      :: zero, one, vsmall


CONTAINS

SUBROUTINE startup(nvar, fit_const)

!     Allocates dimensions for arrays and initializes to zero
!     The calling program must set nvar = the number of variables, and
!     fit_const = .true. if a constant is to be included in the model,
!     otherwise fit_const = .false.
!
!--------------------------------------------------------------------------
use iso_fortran_env
IMPLICIT NONE
INTEGER(kind=int32), INTENT(IN)  :: nvar
LOGICAL, INTENT(IN)  :: fit_const

!     Local variable
INTEGER(kind=int32)   :: i

vsmall = 10. * TINY(zero)

nobs = 0
IF (fit_const) THEN
  ncol = nvar + 1
ELSE
  ncol = nvar
END IF

IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder, row_ptr)
r_dim = ncol * (ncol - 1)/2
ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol),  &
          row_ptr(ncol) )

d = zero
rhs = zero
r = zero
sserr = zero

IF (fit_const) THEN
  DO i = 1, ncol
    vorder(i) = i-1
  END DO
ELSE
  DO i = 1, ncol
    vorder(i) = i
  END DO
END IF ! (fit_const)

! row_ptr(i) is the position of element R(i,i+1) in array r().

row_ptr(1) = 1
DO i = 2, ncol-1
  row_ptr(i) = row_ptr(i-1) + ncol - i + 1
END DO
row_ptr(ncol) = 0

initialized = .true.
tol_set = .false.
rss_set = .false.

RETURN
END SUBROUTINE startup




SUBROUTINE includ(weight, xrow, yelem)

!     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

!     Calling this routine updates D, R, RHS and SSERR by the
!     inclusion of xrow, yelem with the specified weight.

!     *** WARNING  Array XROW is overwritten.

!     N.B. As this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------

use iso_fortran_env
IMPLICIT NONE
double precision,INTENT(IN)                    :: weight, yelem
double precision, DIMENSION(:), INTENT(IN OUT) :: xrow

!     Local variables

INTEGER(kind=int32)     :: i, k, nextr
double precision   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

nobs = nobs + 1
w = weight
y = yelem
rss_set = .false.
nextr = 1
DO i = 1, ncol

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (ABS(w) < vsmall) RETURN
  xi = xrow(i)
  IF (ABS(xi) < vsmall) THEN
    nextr = nextr + ncol - i
  ELSE
    di = d(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    d(i) = dpi
    DO k = i+1, ncol
      xk = xrow(k)
      xrow(k) = xk - xi * r(nextr)
      r(nextr) = cbar * r(nextr) + sbar * xk
      nextr = nextr + 1
    END DO
    xk = y
    y = xk - xi * rhs(i)
    rhs(i) = cbar * rhs(i) + sbar * xk
  END IF
END DO ! i = 1, ncol

!     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
!     residual.

sserr = sserr + w * y * y

RETURN
END SUBROUTINE includ



SUBROUTINE regcf(beta, nreq, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------
use iso_fortran_env
IMPLICIT NONE
INTEGER(kind=int32), INTENT(IN)                  :: nreq
INTEGER(kind=int32), INTENT(OUT)                 :: ifault
double precision, DIMENSION(:), INTENT(OUT) :: beta

!     Local variables

INTEGER(kind=int32)   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
  IF (SQRT(d(i)) < tol(i)) THEN
    beta(i) = zero
    d(i) = zero
    ifault = -i
  ELSE
    beta(i) = rhs(i)
    nextr = row_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
  END IF
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE regcf



SUBROUTINE tolset(eps)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Sets up array TOL for testing for zeroes in an orthogonal
!     reduction formed using AS75.1.

double precision, INTENT(IN), OPTIONAL :: eps

!     Unless the argument eps is set, it is assumed that the input data are
!     recorded to full machine accuracy.   This is often not the case.
!     If, for instance, the data are recorded to `single precision' of about
!     6-7 significant decimal digits, then singularities will not be detected.
!     It is suggested that in this case eps should be set equal to
!     10.0 * EPSILON(1.0)
!     If the data are recorded to say 4 significant decimals, then eps should
!     be set to 1.0E-03
!     The above comments apply to the predictor variables, not to the
!     dependent variable.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Local variables.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER(kind=int32)    :: col, row, pos
double precision  :: eps1, ten = 10.0, total, work(ncol)

!     EPS is a machine-dependent constant.

IF (PRESENT(eps)) THEN
  eps1 = MAX(ABS(eps), ten * EPSILON(ten))
ELSE
  eps1 = ten * EPSILON(ten)
END IF

!     Set tol(i) = sum of absolute values in column I of R after
!     scaling each element by the square root of its row multiplier,
!     multiplied by EPS1.

work = SQRT(ABS(d))
DO col = 1, ncol
  pos = col - 1
  total = work(col)
  DO row = 1, col-1
    total = total + ABS(r(pos)) * work(row)
    pos = pos + ncol - row - 1
  END DO
  tol(col) = eps1 * total
END DO

tol_set = .TRUE.
RETURN
END SUBROUTINE tolset




SUBROUTINE sing(lindep, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Checks for singularities, reports, and adjusts orthogonal
!     reductions produced by AS75.1.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Auxiliary routines called: INCLUD, TOLSET
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(OUT)                :: ifault
LOGICAL, DIMENSION(:), INTENT(OUT)  :: lindep

!     Local variables

double precision  :: temp, x(ncol), work(ncol), y, weight
INTEGER(kind=int32)    :: pos, row, pos2

ifault = 0

work = SQRT(ABS(d))
IF (.NOT. tol_set) CALL tolset()

DO row = 1, ncol
  temp = tol(row)
  pos = row_ptr(row)         ! pos = location of first element in row

!     If diagonal element is near zero, set it to zero, set appropriate
!     element of LINDEP, and use INCLUD to augment the projections in
!     the lower rows of the orthogonalization.

  lindep(row) = .FALSE.
  IF (work(row) <= temp) THEN
    lindep(row) = .TRUE.
    ifault = ifault - 1
    IF (row < ncol) THEN
      pos2 = pos + ncol - row - 1
      x = zero
      x(row+1:ncol) = r(pos:pos2)
      y = rhs(row)
      weight = d(row)
      r(pos:pos2) = zero
      d(row) = zero
      rhs(row) = zero
      CALL includ(weight, x, y)
                             ! INCLUD automatically increases the number
                             ! of cases each time it is called.
      nobs = nobs - 1
    ELSE
      sserr = sserr + d(row) * rhs(row)**2
    END IF ! (row < ncol)
  END IF ! (work(row) <= temp)
END DO ! row = 1, ncol

RETURN
END SUBROUTINE sing



SUBROUTINE ss()

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculates partial residual sums of squares from an orthogonal
!     reduction from AS75.1.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER(kind=int32)    :: i
double precision  :: total

total = sserr
rss(ncol) = sserr
DO i = ncol, 2, -1
  total = total + d(i) * rhs(i)**2
  rss(i-1) = total
END DO

rss_set = .TRUE.
RETURN
END SUBROUTINE ss



SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate covariance matrix for regression coefficients for the
!     first nreq variables, from an orthogonal reduction produced from
!     AS75.1.

!     Auxiliary routine called: INV
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                   :: nreq, dimcov
INTEGER(kind=int32), INTENT(OUT)                  :: ifault
double precision, INTENT(OUT)                :: var
double precision, DIMENSION(:), INTENT(OUT)  :: covmat, sterr

!     Local variables.

INTEGER(kind=int32)                :: dim_rinv, pos, row, start, pos2, col, pos1, k
double precision              :: total
double precision, ALLOCATABLE :: rinv(:)

!     Check that dimension of array covmat is adequate.

IF (dimcov < nreq*(nreq+1)/2) THEN
  ifault = 1
  RETURN
END IF

!     Check for small or zero multipliers on the diagonal.

ifault = 0
DO row = 1, nreq
  IF (ABS(d(row)) < vsmall) ifault = -row
END DO
IF (ifault /= 0) RETURN

!     Calculate estimate of the residual variance.

IF (nobs > nreq) THEN
  IF (.NOT. rss_set) CALL ss()
  var = rss(nreq) / (nobs - nreq)
ELSE
  ifault = 2
  RETURN
END IF

dim_rinv = nreq*(nreq-1)/2
ALLOCATE ( rinv(dim_rinv) )

CALL INV(nreq, rinv)
pos = 1
start = 1
DO row = 1, nreq
  pos2 = start
  DO col = row, nreq
    pos1 = start + col - row
    IF (row == col) THEN
      total = one / d(col)
    ELSE
      total = rinv(pos1-1) / d(col)
    END IF
    DO K = col+1, nreq
      total = total + rinv(pos1) * rinv(pos2) / d(k)
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    END DO ! K = col+1, nreq
    covmat(pos) = total * var
    IF (row == col) sterr(row) = SQRT(covmat(pos))
    pos = pos + 1
  END DO ! col = row, nreq
  start = start + nreq - row
END DO ! row = 1, nreq

DEALLOCATE(rinv)
RETURN
END SUBROUTINE cov



SUBROUTINE inv(nreq, rinv)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Invert first nreq rows and columns of Cholesky factorization
!     produced by AS 75.1.
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                  :: nreq
double precision, DIMENSION(:), INTENT(OUT) :: rinv

!     Local variables.

INTEGER(kind=int32)    :: pos, row, col, start, k, pos1, pos2
double precision  :: total

!     Invert R ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
DO row = nreq-1, 1, -1
  start = row_ptr(row)
  DO col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    total = zero
    DO k = row+1, col-1
      pos2 = pos2 + nreq - k
      total = total - r(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    END DO ! k = row+1, col-1
    rinv(pos) = total - r(pos1)
    pos = pos - 1
  END DO ! col = nreq, row+1, -1
END DO ! row = nreq-1, 1, -1

RETURN
END SUBROUTINE inv



SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

!     Replaces subroutines PCORR and COR of:
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate partial correlations after the variables in rows
!     1, 2, ..., IN have been forced into the regression.
!     If IN = 1, and the first row of R represents a constant in the
!     model, then the usual simple correlations are returned.

!     If IN = 0, the value returned in array CORMAT for the correlation
!     of variables Xi & Xj is:
!       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

!     On return, array CORMAT contains the upper triangle of the matrix of
!     partial correlations stored by rows, excluding the 1's on the diagonal.
!     e.g. if IN = 2, the consecutive elements returned are:
!     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
!     Array YCORR stores the partial correlations with the Y-variable
!     starting with YCORR(IN+1) = partial correlation with the variable in
!     position (IN+1).
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                  :: in, dimc
INTEGER(kind=int32), INTENT(OUT)                 :: ifault
double precision, DIMENSION(:), INTENT(OUT) :: cormat, ycorr

!     Local variables.

INTEGER(kind=int32)    :: base_pos, pos, row, col, col1, col2, pos1, pos2
double precision  :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

!     Some checks.

ifault = 0
IF (in < 0 .OR. in > ncol-1) ifault = ifault + 4
IF (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
IF (ifault /= 0) RETURN

!     Base position for calculating positions of elements in row (IN+1) of R.

base_pos = in*ncol - (in+1)*(in+2)/2

!     Calculate 1/RMS of elements in columns from IN to (ncol-1).

IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
DO col = in+2, ncol
  pos = base_pos + col
  sumxx = d(col)
  DO row = in+1, col-1
    sumxx = sumxx + d(row) * r(pos)**2
    pos = pos + ncol - row - 1
  END DO ! row = in+1, col-1
  IF (sumxx > zero) THEN
    rms(col) = one / SQRT(sumxx)
  ELSE
    rms(col) = zero
    ifault = -col
  END IF ! (sumxx > zero)
END DO ! col = in+1, ncol-1

!     Calculate 1/RMS for the Y-variable

sumyy = sserr
DO row = in+1, ncol
  sumyy = sumyy + d(row) * rhs(row)**2
END DO ! row = in+1, ncol
IF (sumyy > zero) sumyy = one / SQRT(sumyy)

!     Calculate sums of cross-products.
!     These are obtained by taking dot products of pairs of columns of R,
!     but with the product for each row multiplied by the row multiplier
!     in array D.

pos = 1
DO col1 = in+1, ncol
  sumxy = zero
  work(col1+1:ncol) = zero
  pos1 = base_pos + col1
  DO row = in+1, col1-1
    pos2 = pos1 + 1
    DO col2 = col1+1, ncol
      work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
      pos2 = pos2 + 1
    END DO ! col2 = col1+1, ncol
    sumxy = sumxy + d(row) * r(pos1) * rhs(row)
    pos1 = pos1 + ncol - row - 1
  END DO ! row = in+1, col1-1

!     Row COL1 has an implicit 1 as its first element (in column COL1)

  pos2 = pos1 + 1
  DO col2 = col1+1, ncol
    work(col2) = work(col2) + d(col1) * r(pos2)
    pos2 = pos2 + 1
    cormat(pos) = work(col2) * rms(col1) * rms(col2)
    pos = pos + 1
  END DO ! col2 = col1+1, ncol
  sumxy = sumxy + d(col1) * rhs(col1)
  ycorr(col1) = sumxy * rms(col1) * sumyy
END DO ! col1 = in+1, ncol-1

ycorr(1:in) = zero

RETURN
END SUBROUTINE partial_corr




SUBROUTINE vmove(from, to, ifault)

!     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

!     Move variable from position FROM to position TO in an
!     orthogonal reduction produced by AS75.1.
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)    :: from, to
INTEGER(kind=int32), INTENT(OUT)   :: ifault

!     Local variables

double precision  :: d1, d2, x, d1new, d2new, cbar, sbar, y
INTEGER(kind=int32)    :: m, first, last, inc, m1, m2, mp1, col, pos, row

!     Check input parameters

ifault = 0
IF (from < 1 .OR. from > ncol) ifault = ifault + 4
IF (to < 1 .OR. to > ncol) ifault = ifault + 8
IF (ifault /= 0) RETURN

IF (from == to) RETURN

IF (.NOT. rss_set) CALL ss()

IF (from < to) THEN
  first = from
  last = to - 1
  inc = 1
ELSE
  first = from - 1
  last = to
  inc = -1
END IF

DO m = first, last, inc

!     Find addresses of first elements of R in rows M and (M+1).

  m1 = row_ptr(m)
  m2 = row_ptr(m+1)
  mp1 = m + 1
  d1 = d(m)
  d2 = d(mp1)

!     Special cases.

  IF (d1 < vsmall .AND. d2 < vsmall) GO TO 40
  x = r(m1)
  IF (ABS(x) * SQRT(d1) < tol(mp1)) THEN
    x = zero
  END IF
  IF (d1 < vsmall .OR. ABS(x) < vsmall) THEN
    d(m) = d2
    d(mp1) = d1
    r(m1) = zero
    DO col = m+2, ncol
      m1 = m1 + 1
      x = r(m1)
      r(m1) = r(m2)
      r(m2) = x
      m2 = m2 + 1
    END DO ! col = m+2, ncol
    x = rhs(m)
    rhs(m) = rhs(mp1)
    rhs(mp1) = x
    GO TO 40
  ELSE IF (d2 < vsmall) THEN
    d(m) = d1 * x**2
    r(m1) = one / x
    r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
    rhs(m) = rhs(m) / x
    GO TO 40
  END IF

!     Planar rotation in regular case.

  d1new = d2 + d1*x**2
  cbar = d2 / d1new
  sbar = x * d1 / d1new
  d2new = d1 * cbar
  d(m) = d1new
  d(mp1) = d2new
  r(m1) = sbar
  DO col = m+2, ncol
    m1 = m1 + 1
    y = r(m1)
    r(m1) = cbar*r(m2) + sbar*y
    r(m2) = y - x*r(m2)
    m2 = m2 + 1
  END DO ! col = m+2, ncol
  y = rhs(m)
  rhs(m) = cbar*rhs(mp1) + sbar*y
  rhs(mp1) = y - x*rhs(mp1)

!     Swap columns M and (M+1) down to row (M-1).

  40 pos = m
  DO row = 1, m-1
    x = r(pos)
    r(pos) = r(pos-1)
    r(pos-1) = x
    pos = pos + ncol - row - 1
  END DO ! row = 1, m-1

!     Adjust variable order (VORDER), the tolerances (TOL) and
!     the vector of residual sums of squares (RSS).

  m1 = vorder(m)
  vorder(m) = vorder(mp1)
  vorder(mp1) = m1
  x = tol(m)
  tol(m) = tol(mp1)
  tol(mp1) = x
  rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
END DO

RETURN
END SUBROUTINE vmove



SUBROUTINE reordr(list, n, pos1, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Re-order the variables in an orthogonal reduction produced by
!     AS75.1 so that the N variables in LIST start at position POS1,
!     though will not necessarily be in the same order as in LIST.
!     Any variables in VORDER before position POS1 are not moved.

!     Auxiliary routine called: VMOVE
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)               :: n, pos1
INTEGER(kind=int32), DIMENSION(:), INTENT(IN) :: list
INTEGER(kind=int32), INTENT(OUT)              :: ifault

!     Local variables.

INTEGER(kind=int32)    :: next, i, l, j

!     Check N.

ifault = 0
IF (n < 1 .OR. n > ncol+1-pos1) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     Work through VORDER finding variables which are in LIST.

next = pos1
i = pos1
10 l = vorder(i)
DO j = 1, n
  IF (l == list(j)) GO TO 40
END DO
30 i = i + 1
IF (i <= ncol) GO TO 10

!     If this point is reached, one or more variables in LIST has not
!     been found.

ifault = 8
RETURN

!     Variable L is in LIST; move it up to position NEXT if it is not
!     already there.

40 IF (i > next) CALL vmove(i, next, ifault)
next = next + 1
IF (next < n+pos1) GO TO 30

RETURN
END SUBROUTINE reordr



SUBROUTINE hdiag(xrow, nreq, hii, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
!
!                         -1           -1
! The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

!              -1
! where z = x'R

! Here we only calculate the diagonal element hii corresponding to one
! row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                  :: nreq
INTEGER(kind=int32), INTENT(OUT)                 :: ifault
double precision, DIMENSION(:), INTENT(IN)  :: xrow
double precision, INTENT(OUT)               :: hii

!     Local variables

INTEGER(kind=int32)    :: col, row, pos
double precision  :: total, wk(ncol)

!     Some checks

ifault = 0
IF (nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

hii = zero
DO col = 1, nreq
  IF (SQRT(d(col)) <= tol(col)) THEN
    wk(col) = zero
  ELSE
    pos = col - 1
    total = xrow(col)
    DO row = 1, col-1
      total = total - wk(row)*r(pos)
      pos = pos + ncol - row - 1
    END DO ! row = 1, col-1
    wk(col) = total
    hii = hii + total**2 / d(col)
  END IF
END DO ! col = 1, nreq

RETURN
END SUBROUTINE hdiag



FUNCTION varprd(x, nreq) RESULT(fn_val)

!     Calculate the variance of x'b where b consists of the first nreq
!     least-squares regression coefficients.
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                  :: nreq
double precision, DIMENSION(:), INTENT(IN)  :: x
double precision                            :: fn_val

!     Local variables

INTEGER(kind=int32)    :: ifault, row
double precision  :: var, wk(nreq)

!     Check input parameter values

fn_val = zero
ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (nobs <= nreq) ifault = ifault + 8
IF (ifault /= 0) THEN
 ! WRITE(*, '(1x, a, i4)') 'Error in function VARPRD: ifault =', ifault
  RETURN
END IF

!     Calculate the residual variance estimate.

var = sserr / (nobs - nreq)

!     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
!     First call BKSUB2 to calculate (inv R')x by back-substitution.

CALL BKSUB2(x, wk, nreq)
DO row = 1, nreq
  IF(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
END DO

fn_val = fn_val * var

RETURN
END FUNCTION varprd



SUBROUTINE bksub2(x, b, nreq)

!     Solve x = R'b for b given x, using only the first nreq rows and
!     columns of R, and only the first nreq elements of R.
!
!--------------------------------------------------------------------------

INTEGER(kind=int32), INTENT(IN)                  :: nreq
double precision, DIMENSION(:), INTENT(IN)  :: x
double precision, DIMENSION(:), INTENT(OUT) :: b

!     Local variables

INTEGER(kind=int32)    :: pos, row, col
double precision  :: temp

!     Solve by back-substitution, starting from the top.

DO row = 1, nreq
  pos = row - 1
  temp = x(row)
  DO col = 1, row-1
    temp = temp - r(pos)*b(col)
    pos = pos + ncol - col - 1
  END DO
  b(row) = temp
END DO

RETURN
END SUBROUTINE bksub2


END MODULE lsq






!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	WREGRESION
!*******************************************************
subroutine WRegresion(X,Y,W,n,nvar,beta,sterr,se,r2,iopt)
USE lsq
use iso_fortran_env
IMPLICIT NONE
INTEGER(kind=int32)             :: i, ier, j, m, n,nvar,iopt
double precision          :: x(n,nvar), y(n),W(n), xrow(0:nvar+1), beta(0:nvar+1),var, &
covmat(231), sterr(0:nvar+1), totalSS,se,r2
LOGICAL             :: fit_const = .TRUE., lindep(0:20)



! Least-squares calculations
m=nvar
CALL startup(m, fit_const)
DO i = 1, n
  xrow(0) = 1.0_dp
  DO j = 1, m
    xrow(j) = x(i,j)
  END DO
  CALL includ(W(i), xrow, y(i))
END DO


if (iopt.gt.0) then
CALL sing(lindep, ier)
IF (ier /= 0) THEN
  DO i = 0, m
   ! IF (lindep(i)) WRITE(*, '(a, i3)') ' Singularity detected for power: ', i
    !IF (lindep(i)) WRITE(9, '(a, i3)') ' Singularity detected for power: ', i
  END DO
END IF
end if

! Calculate progressive residual sums of squares
CALL ss()
var = rss(m+1) / (n - m - 1)

! Calculate least-squares regn. coeffs.
CALL regcf(beta, m+1, ier)

if (iopt.gt.0) then
! Calculate covariance matrix, and hence std. errors of coeffs.
CALL cov(m+1, var, covmat, 231, sterr, ier)

!WRITE(*, *) 'Least-squares coefficients & std. errors'
!WRITE(9, *) 'Least-squares coefficients & std. errors'
!WRITE(*, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
!WRITE(9, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
!DO i = 0, m
! WRITE(*, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
!      i, beta(i), sterr(i), rss(i+1)
!  WRITE(9, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
!  !      i, beta(i), sterr(i), rss(i+1)
!END DO

!WRITE(*, *)

!WRITE(*, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)

se=SQRT(var)
totalSS = rss(1)
!WRITE(*, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS

r2=(totalSS - rss(m+1))/totalSS

end if
END


subroutine quantile (X,n,alfa,nalfa,Q)
use iso_fortran_env
implicit none
INTEGER(kind=int32) n,nalfa,ip,j,ind(n)
double precision X(n),alfa(nalfa),Q(nalfa),R,xest
call qsortd(x,ind,n)

do j=1,nalfa
IP=int(alfa(j)*(n+1.))
XEST=alfa(j)*(n+1.)
IF(ip .lt. 1) then
Q(j)=X(ind(1))
elseif (ip.ge.n) then
Q(j)=X(ind(n))
else
  R=alfa(j)*(n+1.)-IP
Q(j)=(1.-R)*X(ind(IP)) + R*X(ind(IP+1))
end if
end do
end

double precision function Cuant (X,n,alfa)
use iso_fortran_env
implicit none
INTEGER(kind=int32) n,ip,ind(n)
double precision X(n),alfa,Q,R,xest
call qsortd(x,ind,n)
IP=int(alfa*(n+1.))
XEST=alfa*(n+1.)
IF(ip .lt. 1) then
Q=X(ind(1))
elseif (ip.ge.n) then
Q=X(ind(n))
else
  R=alfa*(n+1.)-IP
Q=(1.-R)*X(ind(IP)) + R*X(ind(IP+1))
end if
Cuant=Q
end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Subroutine RFAST_H
!*********************************************************
subroutine rfast_h(X,Y,W,n,h,p,Xb,Pb,kbin,kernel,nh)
use iso_fortran_env
implicit none
INTEGER(kind=int32) n,i,j,kbin,p,nh,kernel
double precision x(n),y(n),W(n),Xb(kbin),Yb(kbin),Wb(kbin),&
  Pb(kbin,3),h,&
 rango,hmin,hmax,beta(10),xbb(kbin),pred(8)
double precision, allocatable::ls(:,:),li(:,:)

allocate(Li(kbin,3),Ls(kbin,3))

call Binning(X,Y,n,W,Xb,Yb,Wb,kbin)

rango=Xb(kbin)-Xb(1)
hmin= 0
hmax=1

if (h.eq.-1)  then
call Ventana1D(Xb,Yb,Wb,kbin,h,p,hmin,hmax,nh,rango,kernel)
end if

if(h.eq.0) then
call Reglineal (Xb,Yb,Wb,kbin,p,Beta)
do i=1,kbin
Pb(i,1)=beta(1)
Pb(i,2)=0
do j=1,p
pb(i,1)=pb(i,1)+beta(j+1)*Xb(i)**j
pb(i,2)=pb(i,2)+p*beta(j+1)*Xb(i)**(j-1)
end do
end do
goto 1
end if

if (h.eq.-2) then
Pb=0
goto 1
end if

xbb=xb

do i=1,kbin
call Reg1D(Xb,Yb,Wb,kbin,h,p,xbb(i),pred,rango,kernel,0)
pb(i,1)=pred(1)
pb(i,2)=pred(2)
pb(i,3)=pred(3)
end do


1 continue
end subroutine





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	BINNING LINEAL
!***************************************************
subroutine Binning(X,Y,n,W,Xb,Yb,Wb,kbin)
use iso_fortran_env
implicit none
INTEGER(kind=int32) n,i,j,kbin
double precision x(n),y(n),W(n),Xb(kbin),Yb(kbin),Wb(kbin),&
  Area(2),dis1,dis2

Wb=0
Yb=0
do i=1,n
if (W(i).gt.0) then
if (X(i).lt.Xb(1)) then
Wb(1)=wb(1)+W(i)
yb(1)=yb(1)+W(i)*Y(i)
elseif (X(i).gt.Xb(kbin)) then
Wb(kbin)=wb(kbin)+W(i)
yb(kbin)=yb(kbin)+W(i)*Y(i)
else
  do j=1,kbin-1
if (Xb(j).le.X(i).and.X(i).le.Xb(j+1)) then
dis1=X(i)-Xb(j)
dis2=Xb(j+1)-X(i)
Area(1)=dis2/(dis1+dis2)
Area(2)=dis1/(dis1+dis2)
Wb(j)=Wb(j)+W(i)*Area(1)
Yb(j)=Yb(j)+Y(i)*W(i)*Area(1)
Wb(j+1)=Wb(j+1)+W(i)*Area(2)
Yb(j+1)=Yb(j+1)+Y(i)*W(i)*Area(2)
end if
end do
end if
end if
end do
do i=1,kbin
if (Wb(i).gt.0) Yb(i)=Yb(i)/Wb(i)
end do
end subroutine






!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !	VENTANA1D
!**********************************************************
subroutine Ventana1D(X,Y,W,n,h,p,hmin,hmax,nh,rango,kernel)
use iso_fortran_env
implicit none
INTEGER(kind=int32) i,p,n,ih,nh,ih2,Err(nh),kernel
double precision x(n),x2(n),Y(n),h,hmin,hmax,&
  pred(8),W(n),hgrid(nh),ErrH(5000),sumW,rango,&
  VT,sumy,sumy2,maxr2,minr2
double precision,allocatable::ErrCV(:,:),Predh(:,:),WCV(:)
INTEGER(kind=int32),external::which_min
allocate(ErrCV(n,nh),Predh(n,nh),WCV(n))

do ih=1,nh
hgrid(ih)=hmin+(ih-1)*(hmax-hmin)/(nh-1)
end do

Err=0
do ih=nh,1,-1
do i=1,n
WCV=W
if(i.ne.1) WCV(i-1)=0
WCV(i)=0
if(i.ne.n) WCV(i+1)=0
x2(i)=x(i)
call Reg1D(X,Y,WCV,n,hgrid(ih),p,X2(i),pred,rango,kernel,1)
PredH(i,ih)=pred(1)
if (pred(1).eq.-1.0) then
do ih2=ih,1,-1
Err(ih2)=1
end do
goto 3333
end if
end do
end do

3333 continue


ErrH=9e9


do ih=1,nh
if (Err(ih).eq.0) then
sumw=0
ErrH(ih)=0
do i=1,n
sumw=sumw+W(i)
ErrH(ih)=ErrH(ih)+W(i)*(Y(i)-PredH(i,ih))**2
end do
ErrH(ih)=ErrH(ih)/sumw
end if
end do


ih=which_min(ErrH,nh)
h= hgrid(ih)

sumy=0
sumy2=0
sumw=0
do i=1,n
sumw=sumw+W(i)
sumy=sumy+W(i)*Y(i)
sumy2=sumy2+W(i)*Y(i)**2
end do

vt=(sumy2/sumw)-(sumy/sumw)**2


do i=1,nh
if (errh(i).ne.9e9) then
Errh(i)=(vT-ErrH(i))/Vt
else
  Errh(i)=0
end if
end do

minr2=9e9
maxr2=-minr2
do i=1,nh
if (Errh(i).gt.0) then
minr2=min(Errh(i),minr2)
maxr2=max(Errh(i),maxr2)
end if
end do

!11 continue

do i=ih,nh
if (Errh(i)+((maxr2-minr2)*0.5).ge.Errh(ih)) then
h=hgrid(i)
end if
end do


!33 continue

deallocate(ErrCV,Predh,WCV)
end subroutine






!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  Subroutine  REG1D
!****************************************************
subroutine Reg1D(X,Y,W,n,h,p,x0,pred,rango,kernel,ifcv)
use iso_fortran_env
implicit none
INTEGER(kind=int32) i,j,icont,p,iopt,ier,n,kernel,ifcv
double precision x(n),Y(n),h,waux,Beta(10),Sterr(20),se,r2,&
  pred(8),W(n),x0,rango,h2,u,pred2
double precision,allocatable::Vx(:),Vy(:),WW(:),XX(:,:)
allocate(Vx(n),Vy(n),WW(n))

pred=-1

h2=h
!345 continue

if (h2.gt.0) then
!3 continue
icont=0
do i=1,n
u=((X(i)-x0)/rango) /h2
if (W(i).gt.0) then

if(ifcv.eq.1.and.u.eq.0) then
waux=0
else
  if(kernel.eq.1.and.abs(u).le.1) then
waux=W(i) * ( (0.75* (1-(u**2))))
elseif(kernel.eq.2.and.abs(u).le.1) then
waux= W(i) * (1-(abs(u)))
elseif(kernel.eq.3) then
waux=W(i)*( (1/sqrt(2*3.1415927)) * dexp( -0.5*(u**2) ) )
else
  waux=W(i)*0.0
end if
end if



if (waux.gt.0) then
icont=icont+1
Vx(icont)=X(i)-x0
Vy(icont)=Y(i)
WW(icont)=waux
end if
end if
end do

if (icont.gt.6) then
allocate (XX(icont,4))
do i=1,icont
do j=1,p
XX(i,j)=Vx(i)**j
end do
end do
iopt=1
call WRegresion_Javier(XX,Vy,WW,icont,p,beta,sterr,se,r2,iopt,ier)
pred(1)=beta(1)
pred(2)=beta(2)
pred(3)=beta(3)
pred(4)=sterr(1)
pred(5)=sterr(2)
pred(6)=sterr(3)
pred(7)=r2
pred(8)=ier
deallocate(XX)
if (ier.ne.0)  then
pred=-1
goto 445
end if
else
  pred=-1
goto 445
end if
else
  continue
end if


445 continue
deallocate(Vx,Vy,WW)
pred2=pred(2)
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !		REGRESION LINEAL
!***************************************************
subroutine Reglineal (X,Y,W,n,p,Beta)
use iso_fortran_env
implicit none
INTEGER(kind=int32) i,n,j,p,iopt
double precision X(n),Y(n),W(n),beta(p+1),&
  sterr(p+1),se,r2,X2(n,p+1)

do i=1,n
do j=1,p
X2(i,j)=X(i)**j
end do
end do
iopt=0
call WRegresion(X2,Y,W,n,p,beta,sterr,se,r2,iopt)


end subroutine



!***************************************************
  !			WREGRESION_Javier  (ier)
!***************************************************
subroutine WRegresion_Javier(X,Y,W,n,nvar,beta,sterr,se,r2,iopt,ier)
USE lsq
use iso_fortran_env
IMPLICIT NONE

INTEGER(kind=int32)             :: i, ier, j, m, n,nvar,iopt

double precision          :: x(n,nvar), y(n),W(n), xrow(0:nvar+1),&

beta(0:nvar+1),var, covmat(231), sterr(0:nvar+1), &

  totalSS,se,r2

LOGICAL             :: fit_const = .TRUE., lindep(0:20)

! Least-squares calculations

m=nvar

CALL startup(m, fit_const)

DO i = 1, n

xrow(0) = 1.0_dp

DO j = 1, m

xrow(j) = x(i,j)

END DO

CALL includ(W(i), xrow, y(i))

END DO





if (iopt.gt.0) then

CALL sing(lindep, ier)

IF (ier /= 0) THEN

DO i = 0, m

!    IF (lindep(i)) WRITE(*, '(a, i3)') ' Singularity detected for power: ', i

!    IF (lindep(i)) WRITE(9, '(a, i3)') ' Singularity detected for power: ', i

END DO

END IF

end if



! Calculate progressive residual sums of squares

CALL ss()

var = rss(m+1) / (n - m - 1)



! Calculate least-squares regn. coeffs.

CALL regcf(beta, m+1, ier)



if (iopt.gt.0) then

! Calculate covariance matrix, and hence std. errors of coeffs.

CALL cov(m+1, var, covmat, 231, sterr, ier)



!WRITE(*, *) 'Least-squares coefficients & std. errors'

!WRITE(9, *) 'Least-squares coefficients & std. errors'

!WRITE(*, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'

!WRITE(9, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'

DO i = 0, m

!  WRITE(*, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &

  !        i, beta(i), sterr(i), rss(i+1)

!  WRITE(9, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &

  !        i, beta(i), sterr(i), rss(i+1)

END DO



!WRITE(*, *)



!WRITE(*, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)



se=SQRT(var)

totalSS = rss(1)

!WRITE(*, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS



r2=(totalSS - rss(m+1))/totalSS



end if

END







!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE qsortd(x,ind,n)

! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18  Time: 11:55:47
use iso_fortran_env
IMPLICIT NONE
INTEGER(kind=int32), PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER(kind=int32) n,ind(n)
double precision x(n)



!***************************************************************************

  !                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A double precision
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

  ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

  INTEGER(kind=int32)   :: iu(21), il(21)
INTEGER(kind=int32)   :: m, i, j, k, l, ij, it, itt, indx
double precision     :: r
double precision :: t

! LOCAL PARAMETERS -

  ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
r = r + .0390625
ELSE
r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = int(i + r*(j-i))
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
il(m) = i
iu(m) = l
i = k
m = m + 1
GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90
END SUBROUTINE qsortd


subroutine GRID(X,W,n,Xbb,nb)
use iso_fortran_env
implicit none
integer(kind=int32) i,nb,n
double precision X(n),W(n),xmin,xmax,Xbb(nb)
xmin=9e9
xmax=-xmin
do i=1,n
if (W(i).gt.0) then
xmin=min(X(i),xmin)
xmax=max(X(i),xmax)
end if
end do

do i=1,nb
Xbb(i)=xmin+(i-1)*(xmax-xmin)/(nb-1)
end do
end

subroutine Interpola (Xgrid,Pgrid,kbin,X0,P0,P1,n)
! Fit a quintic spline with user control of knot positions.
! If the knots are at tk1, tk2,..., then the fitted spline is
! b0 + b1.t + b2.t^2 + b3.t^3 + b4.t^4 + b5.t^5    for t <= tk1
! b0 + ... + b5.t^5 + b6.(t-tk1)^5                 for tk1 < t <= tk2
! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5  for tk2 < t <= tk3
! b0 + ... + b5.t^5 + b6.(t-tk1)^5 + b7.(t-tk2)^5 + b8.(t-tk3)^5
!                                                  for tk3 < t <= tk4, etc.

! In this version, the knots are evenly spaced.
! Also calculates first & 2nd derivatives of the spline.

! Uses the author's least-squares package in file lsq.f90
! Latest revision - 2 November 2003
! Alan Miller (amiller @ bigpond.net.au)

USE lsq
use iso_fortran_env
IMPLICIT NONE

INTEGER(kind=int32)                :: i, ier, j, n, nk,next_knot,kbin,icont
double precision               :: t, t1, y, dist,Xgrid(kbin),Pgrid(kbin),X0(n),P0(n),P1(n),P2(n)
double precision, PARAMETER    :: one = 1.0_dp,cero=0.0_dp
double precision, ALLOCATABLE  :: knot(:), xrow(:), b(:)



icont=0
do i=1,kbin
if (pgrid(i).ne.-1.0) icont=icont+1
end do

if (icont.gt.5) then
nk=icont/5


!numero de nodos
!if (nk>kbin/5) stop ! '** Too many knots requested - TRY AGAIN'

ALLOCATE ( knot(nk),xrow(0:5+nk), b(0:5+nk) )


! Calculate knot positions, evenly spaced.

dist = (Xgrid(kbin) - Xgrid(1)) / (nk + 1)
t1=Xgrid(1)
DO i = 1, nk
knot(i) = t1 + dist * i
END DO

! WRITE(9, '(a, i4)') 'Number of knots = ', nk


next_knot = 1

! Initialize the least-squares calculations
CALL startup(6+nk, .FALSE.)

DO i=1,kbin
t=Xgrid(i)
y=Pgrid(i)
xrow(0) = one
xrow(1) = (t - t1)
xrow(2) = (t - t1) * xrow(1)
xrow(3) = (t - t1) * xrow(2)
xrow(4) = (t - t1) * xrow(3)
xrow(5) = (t - t1) * xrow(4)
IF (t > knot(next_knot)) next_knot = MIN(nk, next_knot + 1)
DO j = 1, next_knot-1
xrow(5+j) = (t - knot(j))**5
END DO
xrow(5+next_knot:5+nk) = 0.0_dp
if (y.ne.-1.0_dp) CALL includ(one, xrow, y)

END DO

CALL regcf(b, 6+nk, ier)

!WRITE(*, *) ' Coefficient   Value'
!WRITE(*, '(a, g13.5)') ' Constant   ', b(0)
!WRITE(*, '(a, g13.5)') ' Linear     ', b(1)
!WRITE(*, '(a, g13.5)') ' Quadratic  ', b(2)
!WRITE(*, '(a, g13.5)') ' Cubic      ', b(3)
!WRITE(*, '(a, g13.5)') ' Quartic    ', b(4)
!WRITE(*, '(a, g13.5)') ' Quintic    ', b(5)
!WRITE(*, *) ' Knot position   Quintic Coefficient'
!DO j = 1, nk
!  WRITE(*, '(g13.5, t17, g13.5)') knot(j), b(5+j)
!END DO

! Calculate fitted values and derivatives


!call Ordena(X0,n,II)
next_knot = 1
DO i = 1, n
next_knot = 1
t=X0(i)
xrow(0) = one
xrow(1) = (t - t1)
xrow(2) = (t - t1) * xrow(1)
xrow(3) = (t - t1) * xrow(2)
xrow(4) = (t - t1) * xrow(3)
xrow(5) = (t - t1) * xrow(4)
if (i.eq.45) then
continue
end if
55 continue
IF (t > knot(next_knot)) THEN
next_knot = next_knot + 1
IF (next_knot <= nk) THEN
!      WRITE(9, '(a, g13.5)') 'New knot at t = ', knot(next_knot-1)
goto 55
ELSE
next_knot = nk + 1
goto 56
END IF
END IF

56 continue
DO j = 1, next_knot-1
xrow(5+j) = (t - knot(j))**5
END DO
p0(i) = DOT_PRODUCT( b(0:5+next_knot-1), xrow(0:5+next_knot-1) )
p2(i) = ((20*b(5)*(t-t1) + 12*b(4))*(t-t1) + 6*b(3))*(t-t1) +2*b(2)
p1(i) = (((5*b(5)*(t-t1) + 4*b(4))*(t-t1) + 3*b(3))*(t-t1) +2*b(2))*(t-t1) + b(1)
DO j = 1, next_knot-1
p1(i) = p1(i) + 5*b(j+5)*(t - knot(j))**4
p2(i) = p2(i) + 20*b(j+5)*(t - knot(j))**3
END DO
!  WRITE(9, '(f8.3, 4g13.4)') t, d1, d2, fitted, y

!  write (*,*) i
END DO
deallocate ( knot,xrow, b )
else
  p0=-1
p1=-1
p2=-1
end if


end subroutine
