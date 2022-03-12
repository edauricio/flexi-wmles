!=================================================================================================================================
! This file is NOT part of the original FLEXI project. Instead, it is an addition made by the LES group at
! Instituto Tecnológico de Aeronáutica, DCTA/IAE, Brazil, in order to implement wall modeling for LES computations.
!
! For more information on the original project, see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! ORIGINAL DISCLAIMER:
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================

#include "flexi.h"

!==================================================================================================================================
!> Module to implement wall modeling for LES computations, hence dealing with boundary condition operations
!==================================================================================================================================
MODULE MOD_WMLES_Utils
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
#if USE_MPI
USE MPI 
#endif

IMPLICIT NONE
PRIVATE

ABSTRACT INTERFACE
  SUBROUTINE RHS(t, y, res)
    REAL,INTENT(IN)  :: t
    REAL,INTENT(IN)  :: y(:)
    REAL,INTENT(OUT) :: res(:)
  END SUBROUTINE
END INTERFACE

REAL, PARAMETER              :: MAX_FACTOR = 10
REAL, PARAMETER              :: MIN_FACTOR = 0.2
REAL, PARAMETER              :: SAFETY = 0.9
!----------------------------------------------------------------------------------------------
!-------------------- Coefficients for Dormand-Prince RK method -------------------------------
!----------------------------------------------------------------------------------------------
REAL, DIMENSION(6,5), PARAMETER :: A_DP = TRANSPOSE(RESHAPE(&
(/&
0.,                 0.,             0.,                 0.,                 0.,&
1./5.,              0.,             0.,                 0.,                 0.,&
3./40.,             9./40.,         0.,                 0.,                 0.,&
44./45.,            -56./15.,       32./9.,             0.,                 0.,&
19372./6561.,       -25360./2187.,  64448./6561.,       -212./729.,         0.,&
9017./3168.,        -355./33.,      46732./5247.,       49./176.,           -5103./18656.&
/),&
(/5,6/)))

REAL, DIMENSION(6), PARAMETER :: B_DP = (/&
35./384.,     0.,       500./1113.,   125./192.,    -2187./6784.,     11./84.&
/)

REAL, DIMENSION(6), PARAMETER :: C_DP = (/&
0.,     1./5.,      3./10.,     4./5.,      8./9.,      1.&
/)

! Error estimate = y_n - y(hat)_n = b_i - b(hat)_i --- resulting in the following vector:
REAL, DIMENSION(7), PARAMETER :: E_DP = (/&
-71./57600.,  0., 71./16695., -71./1920., 17253./339200., -22./525., 1./40.&
/)

!PROCEDURE(RHS),POINTER :: RHS_Pointer     !< Point to the filter routine to be used

PUBLIC :: SOLVE_ODE, RHS ! Generic interfaces / subroutines
PUBLIC :: RHS_TEST, RHS16, RHS25, RHS27 ! Currently implemented right-hand side routines
PUBLIC :: GetParams
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine performs a single time-step for a self-starting, explicit Runge-Kutta method,
!> according to the coefficients passed in as inputs
!==================================================================================================================================
SUBROUTINE RK_STEP(fun, t, y, f, h, A, B, C, K, y_new)
! MODULES 

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
PROCEDURE(RHS),POINTER,INTENT(IN)    :: fun      !< Pointer to the function that evaluates the RHS of the system, given (t,y)
REAL,INTENT(IN)                      :: t        !< Current time (independent variable)
REAL,INTENT(IN)                      :: y(:)     !< Current value of the dependet variables (1:nSystem)
REAL,INTENT(IN)                      :: f(:)     !< Current value of the derived depent variables (1:nSystem)
REAL,INTENT(IN)                      :: h        !< Time-step size
REAL,INTENT(IN)                      :: A(:,:)   !< Weight coefficients of each RK stage (1:nStages,1:nStages-1)
REAL,INTENT(IN)                      :: B(:)     !< Final stage weights (1:nStages)
REAL,INTENT(IN)                      :: C(:)     !< Coefficients of the intermediate time levels (1:nStages)
REAL,INTENT(OUT)                     :: K(:,:)   !< RHS at each stage, plus the computed new value
                                                 !< of the RHS at t+h (1:nStages+1,1:nSystem)
REAL,INTENT(OUT),OPTIONAL            :: y_new(:) !< Solution vector at the end of time-step, t+h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                     :: nStages, nSystem, i, j, m, n
REAL                        :: t_new, sum
REAL,ALLOCATABLE            :: y_new_local(:), f_new(:), temp_f(:)
!==================================================================================================================================
nStages = SIZE(B,1)
nSystem = SIZE(y,1)


ALLOCATE(y_new_local(nSystem),f_new(nSystem),temp_f(nSystem))
y_new_local = 0.
f_new = 0.
temp_f = 0.
t_new = t

! Write the RHS for the current time (first stage, i.e., k=1)
DO n=1,nSystem
    K(1,n) = f(n)
END DO

! Subsequent stages
DO i=2,nStages
    t_new = t_new + C(i)*h
    DO n=1,nSystem
        sum = 0.
        DO m=1,i-1
            sum = sum + A(i,m)*K(m,n)
        END DO
        y_new_local(n) = y(n) + h*sum
    END DO

    CALL fun(t_new, y_new_local, temp_f)
    DO n=1,nSystem
        K(i,n) = temp_f(n)
    END DO
END DO

DO n=1,nSystem
    sum = 0
    DO j=1,nStages
        sum = sum + B(j)*K(j,n)
    END DO
    y_new_local(n) = y(n) + h*sum
END DO
IF (PRESENT(y_new)) y_new(:) = y_new_local(:)

CALL fun(t_new, y_new_local, temp_f)
DO n=1,nSystem
    K(nStages+1,n) = temp_f(n)
END DO

END SUBROUTINE RK_STEP

!==================================================================================================================================
!> This routine solves an ordinary differencial equation using a Runge-Kutta scheme, given the appropriate inputs
!==================================================================================================================================
SUBROUTINE SOLVE_ODE(fun, t0, y0, tmax, res, ts, max_h, rtol, atol, RKMethod, resFullInt, preStep)
! MODULES 
USE MOD_StringTools,        ONLY: STRICMP
USE MOD_Globals,            ONLY: Abort

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
PROCEDURE(RHS),POINTER,INTENT(IN)           :: fun                  !< Pointer to the function that evaluates the RHS of the system, given (t,y)
REAL,INTENT(IN)                             :: t0                   !< Initial integration time (independent variable)
REAL,INTENT(IN)                             :: y0(:)                !< Initial value of the solution vector (initial conditions)
REAL,INTENT(IN)                             :: tmax                 !< Final integration time (independent variable)
REAL,INTENT(OUT),ALLOCATABLE                :: res(:,:)             !< Results at the final integration time, tmax, or every timestep
REAL,INTENT(OUT),ALLOCATABLE                :: ts(:)                !< Vector containing the discrete integration time / points
REAL,INTENT(IN),OPTIONAL                    :: max_h                !< Maximum value of the time-step h
REAL,INTENT(IN),OPTIONAL                    :: atol                 !< Absolute tolerance
REAL,INTENT(IN),OPTIONAL                    :: rtol                 !< Relative tolerance
CHARACTER(LEN=255),INTENT(IN),OPTIONAL      :: RKMethod             !< Name of the Runge-Kutta method to be used
LOGICAL,INTENT(IN),OPTIONAL                 :: resFullInt           !< Indicates whether res should contain the value at every time step or only the last one
LOGICAL,INTENT(IN),OPTIONAL                 :: preStep              !< Indicates whether there is a pre-RK-step computation that should be done via GetParams(t)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)          :: RK
INTEGER                     :: RK_ID,nStages,nSystem,RK_ErrOrder,i,n,sizeRes,nSteps
REAL                        :: h_max, abs_tol, rel_tol, tstep, t_new, t, err_norm, err_exp, factor
REAL, ALLOCATABLE           :: A(:,:), B(:), C(:), E(:), K(:,:), f(:), y(:), y_new(:), err_est(:)
REAL, ALLOCATABLE           :: ts_buf(:), res_buf(:,:), scale(:)
LOGICAL                     :: fullResult, stepRejected, stepAccepted, preproc
INTEGER :: deb_step
!==================================================================================================================================
nSystem = SIZE(y0,1)
h_max = HUGE(1.)
abs_tol = 1E-6
rel_tol = 0.001
fullResult=.TRUE.
preproc=.FALSE.
IF (PRESENT(max_h)) h_max = max_h
IF (PRESENT(atol)) abs_tol = atol
IF (PRESENT(rtol)) rel_tol = rtol
IF (PRESENT(resFullInt)) fullResult=resFullInt
IF (PRESENT(preStep)) preproc=preStep

! Allocate an initial size for the vector
! Adjust on-the-fly if needed
IF (fullResult) THEN
    IF (ALLOCATED(ts)) DEALLOCATE(ts)
    IF (ALLOCATED(res)) DEALLOCATE(res)
    sizeRes = 100. ! initial estimate
    ALLOCATE(ts(sizeRes))
    ALLOCATE(res(nSystem,sizeRes))
ELSE
    sizeRes = 1
    ALLOCATE(ts(1))
    ALLOCATE(res(nSystem,1))
END IF

RK = "Dormand-Prince"
IF (PRESENT(RKMethod)) RK=RKMethod

! Check which RK method was chosen
RK_ID = 0
IF(STRICMP(RK,"Dormand-Prince")) THEN
    RK_ID=1
    RK_ErrOrder=4
    err_exp = -1. / (RK_ErrOrder + 1)
    nStages=6
ELSE IF (STRICMP(RK,"Fehlberg")) THEN
    ! NOTE: TO BE IMPLEMENTED
    RK_ID=2
    nStages=4 ! ??
END IF

ALLOCATE(A(nStages,nStages-1), B(nStages), C(nStages), E(nStages+1), K(nStages+1,nSystem))
K = 0.
SELECT CASE(RK_ID)
    CASE(1)
        A = A_DP
        B = B_DP
        C = C_DP
        E = E_DP

    CASE DEFAULT
        CALL Abort(__STAMP__, "Runge-Kutta Method Not Implemented")

END SELECT

ALLOCATE(scale(nSystem))
ALLOCATE(f(nSystem),y(nSystem),y_new(nSystem),err_est(nSystem))
f = 0.
y = 0.
y_new = 0.
CALL fun(t0, y0, f)
tstep = InitialStep(fun, t0, y0, f, RK_ErrOrder, rel_tol, abs_tol)
IF (tstep.GT.h_max) tstep = h_max

t = t0
y(:) = y0(:)

nSteps = 1
ts(1) = t
res(:,1) = y(:)
deb_step=0
DO WHILE (t.LT.tmax)
    stepRejected = .FALSE.
    stepAccepted = .FALSE.

    IF (tstep.GT.h_max) tstep = h_max
    IF (preproc) CALL GetParams(t)

    DO WHILE (.NOT.stepAccepted)
        t_new = t + tstep
        IF ((t_new-tmax).GT.0) t_new=tmax
        tstep = t_new - t
        CALL RK_STEP(fun, t, y, f, tstep, A, B, C, K, y_new)
        scale = abs_tol + MAX(y, y_new) * rel_tol
        err_est = 0.
        DO n=1,nSystem
            DO i=1,nStages+1
                err_est(n) = err_est(n) + K(i,n)*E(i)
            END DO
        END DO
        err_est = err_est * tstep / scale
        err_norm = NORM2(err_est) / SQRT(REAL(nSystem))

        IF (err_norm.LT.1) THEN
            IF (ABS(err_norm).LT.1E-10) THEN
                factor = MAX_FACTOR
            ELSE
                factor = MIN(MAX_FACTOR, SAFETY * (err_norm ** err_exp))
            END IF

            IF (stepRejected) factor = MIN(1.,factor)

            tstep = factor * tstep
            stepAccepted = .TRUE.

        ELSE
            tstep = MAX(MIN_FACTOR, SAFETY * (err_norm ** err_exp)) * tstep
            stepRejected = .TRUE.
        END IF
        deb_step = deb_step + 1
    END DO
    t = t_new
    y(:) = y_new(:)
    CALL fun(t, y, f)

    
    IF (fullResult) THEN
        nSteps = nSteps + 1
        IF (nSteps.GT.sizeRes) THEN
            ! No more space in our vector, so we use our buffer to temporarily hold 
            ! the results, so we can allocate more memory for the results
            ALLOCATE(ts_buf(sizeRes))
            ALLOCATE(res_buf(nSystem,sizeRes))
            ts_buf(:) = ts(:)
            res_buf(:,:) = res(:,:)
            DEALLOCATE(ts)
            DEALLOCATE(res)
            ! double the memory space
            ALLOCATE(ts(2*sizeRes))
            ALLOCATE(res(nSystem,2*sizeRes))
            ts(1:sizeRes) = ts_buf(:)
            res(:,1:sizeRes) = res_buf(:,:)
            DEALLOCATE(ts_buf)
            DEALLOCATE(res_buf)
            sizeRes = 2*sizeRes
        END IF
    END IF

    ts(nSteps) = t
    res(:,nSteps) = y
END DO

IF (fullResult) THEN
    ! There might be unused space left in the results vector
    ! Trim them down
    ALLOCATE(ts_buf(sizeRes))
    ALLOCATE(res_buf(nSystem,sizeRes))
    ts_buf(:) = ts(:)
    res_buf(:,:) = res(:,:)
    DEALLOCATE(ts)
    DEALLOCATE(res)
    ! double the memory space
    ALLOCATE(ts(nSteps))
    ALLOCATE(res(nSystem,nSteps))
    ts(:) = ts_buf(1:nSteps)
    res(:,:) = res_buf(:,1:nSteps)
    DEALLOCATE(ts_buf)
    DEALLOCATE(res_buf)
END IF


END SUBROUTINE SOLVE_ODE

!==================================================================================================================================
!> Function to calculate a good initial time step for RK integration
!==================================================================================================================================
FUNCTION InitialStep(fun, t0, y0, f0, order, rtol, atol)
! MODULES 

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
PROCEDURE(RHS),POINTER,INTENT(IN)           :: fun                  !< Pointer to the function that evaluates the RHS of the system, given (t,y)
REAL,INTENT(IN)                             :: t0                   !< Initial integration time (independent variable)
REAL,INTENT(IN)                             :: y0(:)                !< Initial value of the solution vector
REAL,INTENT(IN)                             :: f0(:)                !< Initial value of the derivatives
INTEGER,INTENT(IN)                          :: order                !< Order of the error approximation
REAL,INTENT(IN)                             :: atol                 !< Absolute tolerance
REAL,INTENT(IN)                             :: rtol                 !< Relative tolerance
REAL                                        :: InitialStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                        :: n
REAL                           :: d0,d1,d2,h0,h1
REAL, ALLOCATABLE              :: y1(:),f1(:),scale(:)
!==================================================================================================================================
n = SIZE(y0,1)
ALLOCATE(y1(n),f1(n),scale(n))

scale = atol + ABS(y0)*rtol
d0 = NORM2(y0 / scale) / SQRT(REAL(n))
d1 = NORM2(f0 / scale) / SQRT(REAL(n))

IF ((d0.LT.1E-5).OR.(d1.LT.1E-5)) THEN
    h0 = 1E-6
ELSE
    h0 = 0.01 * d0 / d1
END IF

y1 = y0 + h0 * f0
CALL fun(t0 + h0, y1, f1)
d2 = NORM2((f1-f0) / scale) / SQRT(REAL(n)) / h0

IF ((d1.LE.1E-15) .AND. (d2.LE.1E-15)) THEN
    h1 = MAX(1E-6, h0 * 1E-3)
ELSE
    h1 = (0.01 / MAX(d1,d2)) ** (1. / (order + 1))
END IF

InitialStep = MIN(100. * h0, h1)

END FUNCTION InitialStep


!==================================================================================================================================
!> An implementation of a system of equations
!> dy/dx = z
!> dz/dx = 6y - z
!==================================================================================================================================
SUBROUTINE RHS_TEST(t, y, res)
! MODULES 

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                             :: t                   !< Evaluation time (independent variable)
REAL,INTENT(IN)                             :: y(:)                !< Value of the solution vector at t
REAL,INTENT(OUT)                            :: res(:)              !< Evaluated right-hand side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
res(1) = y(2)
res(2) = 6.*y(1) - y(2)

END SUBROUTINE RHS_TEST

!==================================================================================================================================
!> Implementation of the system of ODEs of Eq. 16 from Zhang & Chen - An Iterative Method .... Falkner-Skan Equation
!==================================================================================================================================
SUBROUTINE RHS16(t, y, res)
! MODULES
USE MOD_WMLES_Vars,             ONLY: eta_inf, beta, beta_0

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                             :: t                   !< Evaluation time (independent variable)
REAL,INTENT(IN)                             :: y(:)                !< Value of the solution vector at t
REAL,INTENT(OUT)                            :: res(:)              !< Evaluated right-hand side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
res(1) = eta_inf*y(2)
res(2) = eta_inf*y(3)
res(3) = -eta_inf*(beta_0*y(1)*y(3) + beta*(1. - y(2)**2))

END SUBROUTINE RHS16

!==================================================================================================================================
!> Auxiliar subroutine to calculate parameters for input in subsequent RHS functions during Falkner-Skan solver
!==================================================================================================================================
SUBROUTINE GetParams(xi_i)
! MODULES
USE MOD_WMLES_Vars,             ONLY: xi_16, sol_16, fs_f, fs_u, fs_v

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                             :: xi_i                !< Time point needed for value of y
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                     :: i, n, sz
REAL, DIMENSION(3)          :: m
!==================================================================================================================================
! Naive, linear search will be used
! Obviously, this can get much better, although vector xi(:) tends to be small
sz = SIZE(xi_16,1)
DO i=1,sz
    IF (xi_16(i).GE.xi_i) EXIT
END DO

IF (i.EQ.1) THEN
    fs_f = sol_16(1,1)
    fs_u = sol_16(2,1)
    fs_v = sol_16(3,1)
ELSE
    m(1) = (xi_i - xi_16(i-1)) * (sol_16(1,i) - sol_16(1,i-1)) / (xi_16(i) - xi_16(i-1))
    m(2) = (xi_i - xi_16(i-1)) * (sol_16(2,i) - sol_16(2,i-1)) / (xi_16(i) - xi_16(i-1))
    m(3) = (xi_i - xi_16(i-1)) * (sol_16(3,i) - sol_16(3,i-1)) / (xi_16(i) - xi_16(i-1))

    fs_f = m(1) + sol_16(1,i-1)
    fs_u = m(2) + sol_16(2,i-1)
    fs_v = m(3) + sol_16(3,i-1)
END IF

END SUBROUTINE GetParams

!==================================================================================================================================
!> Implementation of the system of ODEs of Eq. 16 from Zhang & Chen - An Iterative Method .... Falkner-Skan Equation
!==================================================================================================================================
SUBROUTINE RHS25(t, y, res)
! MODULES
USE MOD_WMLES_Vars,             ONLY: eta_inf, beta, beta_0, fs_f, fs_u, fs_v

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                             :: t                   !< Evaluation time (independent variable)
REAL,INTENT(IN)                             :: y(:)                !< Value of the solution vector at t
REAL,INTENT(OUT)                            :: res(:)              !< Evaluated right-hand side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
res(1) = eta_inf*y(2)
res(2) = eta_inf*y(3)
res(3) = -eta_inf*(beta_0*(y(1)*fs_v + y(3)*fs_f) - 2.*beta*fs_u*y(2))

END SUBROUTINE RHS25

!==================================================================================================================================
!> Implementation of the system of ODEs of Eq. 16 from Zhang & Chen - An Iterative Method .... Falkner-Skan Equation
!==================================================================================================================================
SUBROUTINE RHS27(t, y, res)
! MODULES
USE MOD_WMLES_Vars,             ONLY: eta_inf, beta, beta_0, fs_f, fs_u, fs_v

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                             :: t                   !< Evaluation time (independent variable)
REAL,INTENT(IN)                             :: y(:)                !< Value of the solution vector at t
REAL,INTENT(OUT)                            :: res(:)              !< Evaluated right-hand side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
res(1) = fs_u + eta_inf*y(2)
res(2) = fs_v + eta_inf*y(3)
res(3) = -1. * (beta_0*fs_f*fs_v + beta*(1. - fs_u**2)) - eta_inf*(beta_0*(y(1)*fs_v + fs_f*y(3)) - 2.*beta*fs_u*y(2))

END SUBROUTINE RHS27

END MODULE MOD_WMLES_UTILS