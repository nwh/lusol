!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lusol_constants.f90
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lusol_constants
  use  lusol_precision

  implicit none

  public

  intrinsic ::  epsilon

  ! To indicate undefined parameters:
  character(8), parameter :: unDefC = '-1111111'
  integer(ip),  parameter :: unDefI =  -11111
  real(rp),     parameter :: unDefR =  -11111.0


  ! General parameters:
  integer(ip),  parameter :: COLD   = 0, BASIS  = 1, WARM  = 2, HOT = 3
  integer(ip),  parameter :: Rowtyp = 0, AStats = 1
  integer(ip),  parameter :: Init   = 0, Optml  = 1, Cycle = 2  ! s5dgen
  integer(ip),  parameter :: Intern = 0, Extern = 1
  integer(ip),  parameter :: SaveB  = 0, PrintS = 1, Wrap  = 1  ! s4savB
  integer(ip),  parameter :: FP     = 0, LP     = 1, QP    = 2, &
                             FPE    = 3, FPS    = 4, QPS   = 5
  integer(ip),  parameter :: DefltF = 0, OpenF = 1, StdIn = 2

  integer(ip), parameter  :: NoSubOpt = -1, YesSubOpt = 0


  ! For Hessian:
  integer(ip),  parameter :: INDEF = -1, SEMDEF = 0, POSDEF = 1, &
                             CHOL  =  0, MDCHOL = 1

  integer(ip),  parameter :: UnLim =  0, violLim = 1, UserLim = 2, UserRej = 3

  ! For quasi-Newton:
  integer(ip),  parameter :: LM     =  0, FM     = 1, LMCBD  = 3
  integer(ip),  parameter :: HNorml =  0, HDiag  = 1, HUnit  = 2, HUnset = -1  ! HDinfo
  integer(ip),  parameter :: noBFGS = -1, cnBFGS = 0, ssBFGS = 1               ! QNinfo
  integer(ip),  parameter :: modA   =  1, modAB  = 2                           ! MDinfo


  ! For LUSOL:
  integer(ip),  parameter :: WithL  = 0, WithB  = 1, WithBt = 2, &
                             B      = 0, BR     = 1, BS     = 2, &
                             BT     = 3, WithR  = 0, WithRt = 1, &
                             Normal = 0, Transp = 1
  integer(ip),  parameter :: No     =-1


  ! Miscellaneous numbers:
  integer(ip),  parameter :: i0     = 0,   i1    = 1,   i2   =  2
  real(rp),     parameter :: zero   = 0.0, half  = 0.5, one  = 1.0
  real(rp),     parameter :: two    = 2.0, three = 3.0


  ! For calls to fgwrap functions:
  integer(ip),  parameter :: getF   = 0, &
                             getG   = 1, &
                             getFG  = 2, &
                             getFH  = 3, &
                             getGH  = 4, &
                             getFGH = 5, &
                             getH   = 6

  ! For QP modes:
  integer(ip),  parameter :: genQP  = 10, blkQP = 11, varQP  = 12, &
                             regQP  = 20, pdrQP = 21, dualQP = 30, &
                             QPChol = 0,  CG    = 1,  QNN    = 2

  ! For linear solvers:
  integer(ip),  parameter :: iLUSOL  = 1, &
                             iMA57   = 2, &
                             iUMF    = 3, &
                             iPARD   = 4, &
                             iSUPER  = 5, &
                             iMA97   = 6, &
                             iMUMPS  = 7

  character(8), parameter :: LUnames(7) =    &
                               (/'lusol   ', &
                                 'hsl_ma57', &
                                 'umfpack ', &
                                 'pardiso ', &
                                 'superlu ', &
                                 'hsl_ma97', &
                                 'mumps   '/)

  ! Machine precision:
  real(rp),     parameter :: eps  = epsilon(eps),   &
                             eps0 = eps**(0.80_rp), &  ! IEEE DP  3.00e-13
                             eps1 = eps**(0.67_rp), &  ! IEEE DP  3.67e-11
                             eps2 = eps**(0.50_rp), &  ! IEEE DP  1.49e-08
                             eps3 = eps**(0.33_rp), &  ! IEEE DP  6.05e-06
                             eps4 = eps**(0.25_rp), &  ! IEEE DP  1.22e-04
                             eps5 = eps**(0.20_rp)

  ! Exit codes
  integer(ip), public, parameter :: exit_optimal       = 1, &
                                    exit_feasiblept    = 2, &
                                    exit_accuracy      = 3, &
                                    exit_weakminimizer = 4, &
                                    exit_elasticmin    = 5, &
                                    exit_elasticinf    = 6, &
                                    exit_deadpoint     = 7, &
                                    inf_linearcons     = 11, &
                                    inf_lineareq       = 12, &
                                    inf_nlninfmin      = 13, &
                                    inf_lninfmin       = 14, &
                                    inf_qplincon       = 15, &
                                    inf_nonelastics    = 16, &
                                    exit_unbounded     = 21, &
                                    exit_itnlimit      = 31, &
                                    exit_mjritnlimit   = 32, &
                                    exit_superbasic    = 33, &
                                    exit_timelimit     = 34, &
                                    num_noimprove      = 41, &
                                    num_singularb      = 42, &
                                    num_generalcon     = 43, &
                                    num_illcond        = 44, &
                                    err_dobjective     = 51, &
                                    err_dconstraint    = 52, &
                                    err_QPindef        = 53, &
                                    err_dsecond        = 54, &
                                    err_derivs         = 55, &
                                    usrf_firstfeas     = 61, &
                                    usrf_initpt        = 62, &
                                    usrf_undef         = 63, &
                                    ustp_feval         = 71, &
                                    ustp_ceval         = 72, &
                                    ustp_oeval         = 73, &
                                    err_input          = 91, &
!                                    err_basis          = 92, &
                                    err_indef          = 93, &
                                    specs_read         = 101, &
                                    jac_estimated      = 102, &
                                    mps_read           = 103, &
                                    spec_keywords      = 107, &
                                    err_variables      = 141, &
                                    err_basis          = 142, &
                                    lu_singular        = 151, &
                                    lu_numerical       = 152, &
                                    lu_input           = 153, &
                                    lu_memory          = 154, &
                                    lu_fatal           = 155, &
                                    lu_unsupport       = 156

end module lusol_constants
