!$Id: IASI_o3_mod.f,v 1.3 2011/02/23 00:08:48 daven Exp $
MODULE IASI_CO_OBS_MOD
  
  IMPLICIT NONE 
  
  !mkeller
#include "CMN_SIZE"
  !#include 'netcdf.inc'

  !=================================================================
  ! MODULE VARIABLES
  !=================================================================   

  PRIVATE
      
  PUBLIC READ_IASI_CO_OBS
  PUBLIC CALC_IASI_CO_FORCE
 
  ! Parameters
  INTEGER, PARAMETER           :: N_IASI_NLA = 19
  INTEGER, PARAMETER           :: N_IASI_NPR = 20
  INTEGER, PARAMETER           :: MAXIASI = 2000000
  INTEGER, PARAMETER           :: IASI_COL = 60!59

  ! Module variables

  ! IASI data
  REAL*8, ALLOCATABLE :: IASI_LON(:)
  REAL*8, ALLOCATABLE :: IASI_LAT(:)
  REAL*8, ALLOCATABLE :: IASI_TIME(:)
  REAL*8, ALLOCATABLE :: IASI_CO(:)
  REAL*8, ALLOCATABLE :: IASI_CO_STD(:)
  REAL*8, ALLOCATABLE :: IASI_ALT(:)
  REAL*8, ALLOCATABLE :: IASI_CO_APR(:,:)
  REAL*8, ALLOCATABLE :: IASI_CO_AVK(:,:)
  REAL*8, ALLOCATABLE :: IASI_QUAL_FLAG(:)
  REAL*8, ALLOCATABLE :: IASI_CLOUD_COVER(:)
  REAL*8, ALLOCATABLE :: IASI_SOLAR_ZENITH(:)
  REAL*8, ALLOCATABLE :: IASI_CO_LVL(:)
  !REAL*8, ALLOCATABLE :: TIME_FRAC(:)
  !REAL*8, ALLOCATABLE :: TEMPDATA(:,:)

  ! MLS grid specification
  INTEGER :: N_IASI_NOB, N_IASI_DATA


  ! mkeller: logical flag to check whether data is available for given day
  LOGICAL :: DATA_PRESENT
CONTAINS
  !------------------------------------------------------------------------------

  SUBROUTINE READ_IASI_CO_OBS( YYYYMMDD, N_IASI_NOB )
    !
!******************************************************************************
!  Subroutine READ_IASI_CO_OBS reads the file and passes back info contained
!  therein. (dkh, 02/19/09) 
! 
!  Based on READ_IASI_NH3 OBS (dkh, 04/26/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    INTEGER : Current year-month-day
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NIASI      (INTEGER) : Number of IASI retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) IASI    (IASI_CO_OBS) : IASI retrieval for current day 
!     
!  NOIASI:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!******************************************************************************
!
    ! Reference to f90 modules
    USE DIRECTORY_MOD,          ONLY : DATA_DIR
    USE TIME_MOD,               ONLY : EXPAND_DATE
    USE FILE_MOD,               ONLY : IOERROR

#include "CMN_SIZE" ! size parameters


    ! Arguments
    INTEGER,            INTENT(IN)  :: YYYYMMDD
    
    ! local variables 
    INTEGER                         :: FID
    INTEGER                         :: LIASI
    INTEGER                         :: N_IASI_NOB
    INTEGER                         :: N, J
    
    CHARACTER(LEN=5)                :: TMP
    CHARACTER(LEN=255)              :: READ_FILENAME

    REAL*8, PARAMETER               :: FILL = -999.0D0
    REAL*8, PARAMETER               :: TOL  = 1d-04
    REAL*8                          :: TMP1
    REAL*8                          :: DUMMY(IASI_COL)
    REAL*8                          :: TEMPDATA(MAXIASI, IASI_COL)
    INTEGER                         :: I, II, III, K
    INTEGER                         :: IU_FILE, IU_DATA, IOS

    
    !=================================================================
    ! READ_IASI_CO_OBS begins here!
    !=================================================================
    CALL CLEANUP_IASI
    ! filename root 
    READ_FILENAME = TRIM( 'iasi_CO_LATMOS_ULB_YYYYMMDD_v20140922.txt' )

    ! Expand date tokens in filename 
    CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 ) 
    
    ! Construct complete filename 
    !      READ_FILENAME = TRIM( DATA_DIR ) // TRIM( 'IASI_CO/' ) // 
    !     &                TRIM( READ_FILENAME )
    READ_FILENAME = '/users/jk/15/xzhang/IASI_CO/' // TRIM( READ_FILENAME )
    WRITE(6,*) '    - READ_IASI_CO_OBS: reading file: ', READ_FILENAME
    IU_FILE = 67
    ! mkeller: check to see if file exists
    INQUIRE(FILE=READ_FILENAME, EXIST = DATA_PRESENT)

    IF (.NOT. DATA_PRESENT) THEN
       PRINT *,"IASI file", TRIM(READ_FILENAME), " not found, "// "assuming that there is no data for this day."
       RETURN
    ELSE
       PRINT *,"IASI file found!"
    ENDIF
    
    OPEN(IU_FILE, FILE=READ_FILENAME, IOSTAT=IOS)
    N_IASI_DATA = 0
    N_IASI_NOB = 0 
    DO
       READ(IU_FILE, *, IOSTAT=IOS) DUMMY
       IF (IOS /= 0) EXIT
       N_IASI_DATA = N_IASI_DATA + 1
       I = I + 1
       TEMPDATA(I,:) = DUMMY(:)
    ENDDO
    CLOSE(IU_FILE)

    DO I = 1, N_IASI_DATA
       IF (TEMPDATA(I,3) .EQ. REAL(YYYYMMDD)) THEN
          N_IASI_NOB = N_IASI_NOB + 1
       ENDIF
    ENDDO
    
    ALLOCATE(IASI_LAT(N_IASI_NOB))
    ALLOCATE(IASI_LON(N_IASI_NOB))
    ALLOCATE(IASI_TIME(N_IASI_NOB))
    ALLOCATE(IASI_SOLAR_ZENITH(N_IASI_NOB))
    ALLOCATE(IASI_QUAL_FLAG(N_IASI_NOB))
    ALLOCATE(IASI_CLOUD_COVER(N_IASI_NOB))
    ALLOCATE(IASI_CO(N_IASI_NOB))
    ALLOCATE(IASI_CO_STD(N_IASI_NOB))
    ALLOCATE(IASI_CO_APR(N_IASI_NOB, N_IASI_NLA))
    ALLOCATE(IASI_CO_AVK(N_IASI_NOB, N_IASI_NLA))
    ALLOCATE(IASI_ALT(N_IASI_NPR))
    ALLOCATE(IASI_CO_LVL(N_IASI_NOB))

    DO I = 1, N_IASI_NOB
       IASI_LAT(I) = TEMPDATA(I,1)
       IASI_LON(I) = TEMPDATA(I,2)
       IASI_TIME(I) = TEMPDATA(I,4)
       IASI_SOLAR_ZENITH(I) = TEMPDATA(I,5)
       IASI_QUAL_FLAG(I) = TEMPDATA(I,16)!15)
       !PRINT *, "IASI_QUAL_FLAG", IASI_QUAL_FLAG(I)
       IASI_CLOUD_COVER(I) = TEMPDATA(I,17)!16)
       !READ(IU_FILE, *) (TEMPDATA(K), K=19,IASI_COL)
       IASI_CO(I) = TEMPDATA(I,21)!20)
       IASI_CO_STD(I) = TEMPDATA(I,21)*TEMPDATA(I,22)
       DO J = 1, N_IASI_NLA
          IASI_CO_APR(I,J) = TEMPDATA(I,22+J) !21
          IASI_CO_AVK(I,J) = TEMPDATA(I,41+J) !40
          !IASI_ALT(J) = REAL(J)-1d0
       ENDDO
    !IASI_ALT(N_IASI_NLA+1) = 60d0
    ENDDO
    DO J = 1, N_IASI_NLA
       IASI_ALT(J) = REAL(J) -1d0
    ENDDO
    
    IASI_ALT(N_IASI_NPR) = 60d0
    !PRINT *, "TEMPDATA", TEMPDATA(:)
    !PRINT *, "IASI_TIME", TEMPDATA(N_IASI_NOB-10000:N_IASI_NOB,3)
    !PRINT *, "IASI_CO_APR", IASI_CO_APR(1,:), "LEVEL2", IASI_CO_APR(2,:)
    !PRINT *, "IASI_CO_AVK", IASI_CO_AVK(1,:), "LEVEL2", IASI_CO_AVK(2,:)
    !PRINT *, "IASI_CO", IASI_CO(1:10)
      !-------------------------------- 
      ! Calculate S_OER_INV
      !-------------------------------- 

      ! loop over records
    ! Now determine how many of the levels in CO are
    ! 'good' and how many are just FILL.
    !ALLOCATE(IASI_CO_LVL(N_IASI_NOB))
    DO N = 1, N_IASI_NOB
       J = 1
       DO WHILE ( J .le. N_IASI_NLA )
          ! check if the value is good
          IF (IASI_CO_APR(N,J) > FILL ) THEN
             ! save the number of good levels as LTES               
             IASI_CO_LVL(N) = N_IASI_NLA - J + 1
             ! and now we can exit the while loop
             J = N_IASI_NLA + 1
             ! otherwise this level is just filler               
          ELSE
             ! so proceed to the next one up
             J = J + 1
          ENDIF
       ENDDO
    ENDDO

      ! Return to calling program
  END SUBROUTINE READ_IASI_CO_OBS

!------------------------------------------------------------------------------

  SUBROUTINE CHECK( STATUS, LOCATION )
!
!******************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was made   
!     
!     
!  NOIASI:
!
!******************************************************************************
!
    ! Reference to f90 modules 
    USE ERROR_MOD,    ONLY  : ERROR_STOP
    USE NETCDF 
 
    ! Arguments
    INTEGER, INTENT(IN)    :: STATUS 
    INTEGER, INTENT(IN)    :: LOCATION
    
    !=================================================================
    ! CHECK begins here!
    !=================================================================

    IF ( STATUS /= NF90_NOERR ) THEN 
       WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
       WRITE(6,*) 'At location = ', LOCATION 
       CALL ERROR_STOP('netCDF error', 'IASI_o3_mod')
     ENDIF
     
     ! Return to calling program
   END SUBROUTINE CHECK

!------------------------------------------------------------------------------

   SUBROUTINE CALC_IASI_CO_FORCE( COST_FUNC )
        !
!******************************************************************************
!  Subroutine CALC_IASI_CO_FORCE calculaIASI the adjoint forcing from the IASI
!  CO observations and updates the cost function. (dkh, 02/15/09)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost function                        [unitless]
!     
!     
!  NOIASI:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (2 ) Add more diagnostics.  Now read and write doubled CO (dkh, 11/08/09) 
!  (3 ) Now use CSPEC_AFTER_CHEM and CSPEC_AFTER_CHEM_ADJ (dkh, 02/09/11) 
!******************************************************************************
!
     ! Reference to f90 modules
     USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ, ADJ_FORCE
     USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
     USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
     !USE ADJ_ARRAYS_MOD,     ONLY : CO_PROF_SAV
     USE ADJ_ARRAYS_MOD,     ONLY : ID2C
     USE CHECKPT_MOD,        ONLY : CHK_STT
     USE COMODE_MOD,         ONLY : JLOP, VOLUME
     USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
     USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
     USE DAO_MOD,            ONLY : AD, AIRDEN, AIRVOL
     USE DAO_MOD,            ONLY : BXHEIGHT
     USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
     USE GRID_MOD,           ONLY : GET_IJ, GET_XMID, GET_YMID
     USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
     USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
     USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_HOUR
     USE TRACER_MOD,         ONLY : XNUMOLAIR, XNUMOL
     USE TRACER_MOD,         ONLY : TCVV
     USE TRACERID_MOD,       ONLY : IDTCO
     USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
     
     
#     include      "CMN_SIZE"      ! Size params
     
     ! Arguments
     REAL*8, INTENT(INOUT)       :: COST_FUNC
     
     ! Local variables 
     INTEGER                     :: NTSTART, NTSTOP, NT, I_IASI 
     INTEGER                     :: IIJJ(2), I,      J
     INTEGER                     :: L,       LL,     LIASI,  L0
     INTEGER                     :: JLOOP
     INTEGER                     :: GC_HOUR
     REAL*8                      :: GC_PRES(LLPAR)
     REAL*8                      :: GC_CO_NATIVE(LLPAR)
     REAL*8                      :: GC_CO(N_IASI_NLA), IASI_APR_WEB(N_IASI_NLA)
     REAL*8                      :: GC_PSURF, IASI_PSURF, GC_ALT(IIPAR, JJPAR, LLPAR+1)
     REAL*8                      :: CO_HAT(N_IASI_NLA)
     REAL*8                      :: SOBS_CO_AVK(IIPAR,JJPAR,N_IASI_NLA)
     REAL*8                      :: CO_PERT(N_IASI_NLA), SOBS_AVK_TOT(IIPAR,JJPAR,N_IASI_NLA)
     REAL*8                      :: FORCE, IASI_RATIO(N_IASI_NLA)
     REAL*8                      :: DIFF, CO_COL_OBS, CO_COL_HAT
     REAL*8                      :: IASI_PCENTER(N_IASI_NLA)
     REAL*8                      :: IASI_APR(N_IASI_NLA)
     REAL*8                      :: NEW_COST(IIPAR,JJPAR)
     REAL*8                      :: OLD_COST
     REAL*8                      :: XNUAIR, SOBS_RAND
     REAL*8, SAVE                :: TIME_FRAC(MAXIASI)
     REAL*8                      :: ALT_SURF(IIPAR,JJPAR)
     REAL*8                      :: SOBS_CO_STD, SOBS_CO_COL
     REAL*8                      :: SOBS_CO_NORM(IIPAR,JJPAR), SOBS_CO_MEAN(IIPAR,JJPAR)
     REAL*8                      :: SOBS_CO_TOT(IIPAR,JJPAR), SOBS_CO_AVG(IIPAR,JJPAR)
     REAL*8                      :: SOBS_CO(IIPAR,JJPAR), SOBS_VAR
     REAL*8                      :: SOBS_CO_ERR(IIPAR,JJPAR), SOBS_CO_VAR(IIPAR,JJPAR)
     REAL*8                      :: SOBS_ERR, SOBS_VAR_LIMIT, SOBS_CO_HAT(IIPAR,JJPAR)
     REAL*8                      :: SOBS_HAT(IIPAR,JJPAR)
     
     
     REAL*8                      :: GC_CO_NATIVE_ADJ(LLPAR)
     REAL*8                      :: CO_HAT_ADJ
     REAL*8                      :: CO_PERT_ADJ(N_IASI_NLA)
     REAL*8                      :: GC_CO_ADJ(N_IASI_NLA)
     
     REAL*8                      :: GC_ADJ_TEMP(IIPAR,JJPAR,LLPAR)
     REAL*8                      :: GC_ADJ_TEMP_COST(IIPAR,JJPAR)
     ! arrays needed for superobservations
     LOGICAL                     :: SUPER_OBS = .TRUE.  ! do super observations?
     REAL*8                      :: SOBS_COUNT(IIPAR,JJPAR)
     
     LOGICAL, SAVE               :: FIRST = .TRUE. 
     LOGICAL, SAVE                   :: SECOND = .TRUE.
     INTEGER                         :: IU_FILE, IU_DATA, IOS
     CHARACTER(LEN=255)              :: FILENAME_ALT, ALT_PATH, FILE_ALT, FILENAME
    
     
     !=================================================================
     ! CALC_IASI_CO_FORCE begins here!
     !=================================================================
     
     print*, '     - CALC_IASI_CO_FORCE '
     CALL RAND
     ! Reset
     IASI_APR(1:N_IASI_NLA) = (/10.161773,9.485584,9.1867937,8.9939589,8.8739,8.8418514,8.7545818, &
                            8.6420669,8.3986495,8.0001663,7.5277656,6.9161481,6.3326457,5.7343264, &
                            5.1686031,4.6027417,4.1054126,3.6674835,12.211041/)
     XNUAIR = 28.964d-3
     NEW_COST = 0D0 
     SOBS_COUNT = 0d0
     ADJ_FORCE = 0d0
     IASI_RATIO = 0d0
     GC_ADJ_TEMP_COST = 0d0
     SOBS_CO_NORM = 0d0
     SOBS_CO_MEAN = 0d0
     SOBS_CO_TOT = 0d0
     SOBS_CO_AVG = 0d0
     SOBS_CO = 0d0
     SOBS_CO_ERR = 0d0
     SOBS_CO_VAR = 0d0

     GC_HOUR = GET_HOUR()
     SOBS_VAR_LIMIT = 2.5e17
     ! Save a value of the cost function first
     OLD_COST = COST_FUNC
     ALT_PATH = '/users/jk/07/xzhang/met_field/'
     FILE_ALT = '20000101.cn.4x5.dat'
     
     FILENAME_ALT = TRIM(ALT_PATH) // TRIM(FILE_ALT)
     OPEN(UNIT = 13, FILE = "/users/jk/07/xzhang/met_field/20000101.cn.4x5.dat", STATUS="old",ACTION="read")
     READ(13,*) ALT_SURF
     CLOSE(13)
     PRINT *, "GET_NHMS", GET_NHMS()
     IF ( SECOND ) THEN
        FILENAME = 'lat_orb_iasico.NN.m'
        CALL EXPAND_NAME( FILENAME, N_CALC )
        FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
        OPEN( 901,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        FILENAME = 'lon_orb_iasico.NN.m'
        CALL EXPAND_NAME( FILENAME, N_CALC )
        FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
        OPEN( 902,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        FILENAME = 'co_chi_sq_iasico.NN.m'
        CALL EXPAND_NAME( FILENAME, N_CALC )
        FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
        OPEN( 904,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        FILENAME = 'diff_iasico.NN.m'
        CALL EXPAND_NAME( FILENAME, N_CALC )
        FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
        OPEN( 912,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    ENDIF

     ! Check if it is the last hour of a day
     IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN
        
        ! Read the IASI CO file for this day
        CALL READ_IASI_CO_OBS( GET_NYMD(), N_IASI_NOB )
        !PRINT *, "N_IASI", N_IASI_NOB
        !PRINT *, "TIME", IASI_TIME(N_IASI_NOB-1000:N_IASI_NOB)
        ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
        TIME_FRAC(1:N_IASI_NOB) = IASI_TIME(1:N_IASI_NOB)/240000d0
     ENDIF
     
     !IF(.NOT. DATA_PRESENT) THEN
     !PRINT *,"No IASI data present for this day, nothing to do here."
     !RETURN
     !ENDIF
     
     ! Get the range of TES retrievals for the current hour
     CALL GET_NT_RANGE( N_IASI_NOB, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP )

     IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN
        
        print*, ' No matching IASI CO obs for this hour'
        RETURN
     ENDIF
     
     !print*, ' for hour range: ', GET_NHMS(), IASI_TIME(NTSTART), IASI_TIME(NTSTOP)
     !print*, ' found record range: ', NTSTART, NTSTOP
     
     
     DO I_IASI = NTSTART, NTSTOP, -1
        IF ( (IASI_QUAL_FLAG(I_IASI) .EQ. 0d0 ) .AND. &
             (IASI_CLOUD_COVER(I_IASI) .EQ. 0d0) .AND. &
             ( IASI_LAT(I_IASI) > -60d0 ) .AND. &
             ( IASI_LAT(I_IASI) < 75d0 ) .AND. &
             ( IASI_CO(I_IASI) > 0d0 ) ) THEN
           ! Get model grid coordinate indices that correspond to the observation
                  ! For safety, initialize these up to LLTES
           GC_CO       = 0d0
           IASI_APR_WEB= 0d0
           CO_PERT     = 0d0
           CO_HAT      = 0d0
           IASI_RATIO  = 0d0

           LIASI = IASI_CO_LVL(I_IASI)
           IIJJ = GET_IJ(REAL(IASI_LON(I_IASI),4),REAL(IASI_LAT(I_IASI),4))
           I = IIJJ(1)
           J = IIJJ(2)
           L0 = N_IASI_NLA - LIASI

           ! Get GC surface pressure (mbar)
           GC_PSURF = GET_PEDGE(I,J,1)
           GC_ALT(I,J,1) = ALT_SURF(I,J)*1d-3
           
           DO L = 1, LLPAR
              JLOOP = JLOP(I,J,L)
              GC_PRES(L) = GET_PCENTER(I,J,L)
              GC_ALT(I,J,L+1) = (SUM(BXHEIGHT(I,J,1:L)) + ALT_SURF(I,J))*1d-3
              GC_CO_NATIVE(L) = CHK_STT(I,J,L,IDTCO) * TCVV(IDTCO)* XNUMOLAIR* BXHEIGHT(I,J,L) *100d0 /(AIRVOL(I,J,L)* 1d6)
           ENDDO
           !PRINT *, "GC_ALT", GC_ALT(I,J,:)
           !PRINT *, "IASI_ALT", IASI_ALT(:)
           !PRINT *, "GC_CO_NATIVE", SUM(GC_CO_NATIVE(:))           
           CALL BIN_DATA_IASI(GC_ALT(I,J,:),IASI_ALT(L0+1:L0+LIASI),GC_CO_NATIVE(:), GC_CO, IASI_APR(L0+1:L0+LIASI),LIASI, 1)
           !PRINT *, "GC_ALT", GC_ALT(I,J,:)
           !PRINT *, "IASI_ALT", IASI_ALT(L0+1:L0+LIASI)
           !PRINT *, " GC_CO_NATIVE", GC_CO_NATIVE(:)
           !PRINT *, " GC_CO", GC_CO(:)
           !--------------------------------------------------------------
           ! Apply IASI CO observation operator
           !
           !   x_hat = x_a + A_k ( x_m - x_a )
           !
           !  where
           !    x_hat = GC modeled column as seen by IASI [lnvmr]
           !    x_a   = IASI apriori column               [lnvmr]
           !    x_m   = GC modeled column                [lnvmr]
           !    A_k   = IASI averaging kernel
           !
           !   OR
           ! 
           !   x_smoothed = A_k*x_m + (I-A_k)*x_a
           !  where
           !    x_smoothed = GC modeled column smoothed by IASI [vmr]
           !    x_a   = IASI apriori partial column               [vmr]
           !    x_m   = GC modeled profile                [vmr]
           !    A_k   = IASI averaging kernel
           !--------------------------------------------------------------
           ! x_m - x_a
           ! x_a + A_k * ( x_m - x_a )
           !PRINT *, "GC_CO", SUM(GC_CO(:))
           !PRINT *, "IASI_CO_AVK", IASI_CO_AVK(I_IASI,:)
           DO L = 1, LIASI
              IF ( (IASI_CO_APR(I_IASI,L+L0) > 0) .AND. &
                   (GC_CO(L) > 0) ) THEN
                 IASI_APR_WEB(L) = IASI_APR(L+L0) * AIRDEN(L,I,J) * XNUMOLAIR * 1d-14 * 1d5
                 !IASI_RATIO(L) = IASI_CO_APR(I_IASI,L+L0)/IASI_APR_WEB(L)
                 !CO_PERT(L) = IASI_RATIO(L)*GC_CO(L) - IASI_CO_APR(I_IASI,L+L0)
                 CO_PERT(L) = GC_CO(L) - IASI_CO_APR(I_IASI,L+L0)
                 CO_HAT(L)  = IASI_CO_APR(I_IASI,L+L0) + IASI_CO_AVK(I_IASI,L0+L) * CO_PERT(L)
              ENDIF
              ! actual comparison
           ENDDO
           CO_COL_HAT = SUM(CO_HAT(1:LIASI))
           SOBS_CO_STD = (1d0 - IASI_CO_STD(I_IASI)/IASI_CO(I_IASI))**2
           SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
           SOBS_CO_NORM(I,J) = SOBS_CO_NORM(I,J) + SOBS_CO_STD
           SOBS_CO_MEAN(I,J) = SOBS_CO_MEAN(I,J) + IASI_CO(I_IASI) * SOBS_CO_STD
           SOBS_CO_TOT(I,J) = SOBS_CO_TOT(I,J) + IASI_CO(I_IASI)
           SOBS_CO_HAT(I,J) = SOBS_CO_HAT(I,J) + CO_COL_HAT
           SOBS_AVK_TOT(I,J,:) = SOBS_AVK_TOT(I,J,:) + IASI_CO_AVK(I_IASI,:) * SOBS_CO_STD
           !PRINT *, "SOBS_CO_STD", SOBS_CO_STD
        ENDIF
     ENDDO
     DO I = 1, IIPAR
        DO J = 1, JJPAR
           IF ( ( SOBS_CO_NORM(I,J) > 0d0 ) .AND. &
                ( SOBS_COUNT(I,J) > 0d0 ) ) THEN
              SOBS_CO(I,J) = SOBS_CO_MEAN(I,J)/SOBS_CO_NORM(I,J)
              SOBS_CO_AVG(I,J) = SOBS_CO_TOT(I,J)/SOBS_COUNT(I,J)
              SOBS_CO_AVK(I,J,:) = SOBS_AVK_TOT(I,J,:)/SOBS_CO_NORM(I,J)
              SOBS_HAT(I,J) = SOBS_CO_HAT(I,J)/SOBS_COUNT(I,J)
              SOBS_CO_ERR(I,J) = SOBS_CO(I,J)*(1d0 - SQRT(SOBS_CO_NORM(I,J)/SOBS_COUNT(I,J)))
           ENDIF
        ENDDO
     ENDDO
     DO I_IASI = NTSTART, NTSTOP, -1
        IF ( (IASI_QUAL_FLAG(I_IASI) .EQ. 0d0 ) .AND. &
             (IASI_CLOUD_COVER(I_IASI) .EQ. 0d0) .AND. &
             ( IASI_LAT(I_IASI) > -60d0 ) .AND. &
             ( IASI_LAT(I_IASI) < 75d0 ) .AND. &
             ( IASI_CO(I_IASI) > 0d0 ) ) THEN
           ! Get model grid coordinate indices that correspond to the observation
                  ! For safety, initialize these up to LLTES

           IIJJ = GET_IJ(REAL(IASI_LON(I_IASI),4),REAL(IASI_LAT(I_IASI),4))
           I = IIJJ(1)
           J = IIJJ(2)

           IF ( ( SOBS_CO_NORM(I,J) > 0d0 ) .AND. &
                ( SOBS_COUNT(I,J) > 0d0 ) ) THEN
              SOBS_CO_VAR(I,J) = SOBS_CO_VAR(I,J) + (IASI_CO(I_IASI) - SOBS_CO_AVG(I,J))**2
           ENDIF
        ENDIF
     ENDDO

     DO I=1,IIPAR
        DO J=1,JJPAR
           FORCE = 0d0
           DIFF = 0d0
           CO_PERT_ADJ = 0d0
           IF ( SOBS_COUNT(I,J) > 3d0 ) THEN
              SOBS_VAR = SOBS_CO_VAR(I,J)/SOBS_COUNT(I,J)
              SOBS_VAR = MIN(SOBS_VAR,SOBS_VAR_LIMIT)                 
              SOBS_ERR = SQRT(SOBS_CO_ERR(I,J)**2/SOBS_COUNT(I,J)+SOBS_VAR)
              IF ( SOBS_ERR/SOBS_CO(I,J) < 0.025 ) THEN
                 SOBS_ERR = 0.025*SOBS_CO(I,J)
              ENDIF
              CALL RANDOM_NUMBER(SOBS_RAND)
              SOBS_CO_COL = SOBS_CO(I,J) !- SQRT(SOBS_VAR)*SOBS_RAND
              DIFF = SOBS_HAT(I,J) - SOBS_CO_COL
              FORCE = DIFF/SOBS_ERR**2
              NEW_COST(I,J) = 0.5 * DIFF * FORCE
              !PRINT *, "FORCE", FORCE
              !PRINT *, "CO_COL_HAT", SOBS_HAT(I,J)
              !PRINT *, "SOBS_CO_COL", SOBS_CO_COL
              !PRINT *, "SOBS_ERR", SOBS_ERR
              WRITE(912,110) (DIFF/1d12)
              WRITE(902,110) (GET_XMID(I))
              WRITE(901,110) (GET_YMID(J))
              WRITE(904,110) (2*NEW_COST(I,J))

              !DO L = 1, LIASI
                 ! adjoint of IASI operator
                 !CO_PERT_ADJ(L) = SOBS_CO_AVK(I,J,L) * FORCE
                 !CO_PERT_ADJ(L) = CO_PERT_ADJ(L)*IASI_RATIO(L)
              !ENDDO
              !CALL BIN_DATA_IASI(GC_ALT(I,J,:),IASI_ALT(:),GC_CO_NATIVE_ADJ(:), CO_PERT_ADJ, IASI_APR(1+L0:L0+LIASI), LIASI, -1)
              !PRINT *, "CO_PERT_ADJ", CO_PERT_ADJ(:)
              !PRINT *, "GC_CO_NATIVE_ADJ(:)", GC_CO_NATIVE_ADJ(:)

              DO L = 1, LLPAR
                 ! Adjoint of unit conversion
                 ADJ_FORCE(I,J,L,IDTCO) = FORCE * TCVV(IDTCO) * BXHEIGHT(I,J,L) * 100d0 * XNUMOLAIR / (1d6 *AIRVOL(I,J,L))
                 !PRINT *, "ADJ_FORCE", ADJ_FORCE(I,J,L,IDTCO)
                 !PRINT *, "BXHEIGHT", BXHEIGHT(I,J,L)/AIRVOL(I,J,L)
                 IF ( ITS_IN_THE_TROP(I,J,L) ) THEN
                    STT_ADJ(I,J,L,IDTCO) = STT_ADJ(I,J,L,IDTCO) + ADJ_FORCE(I,J,L,IDTCO)
                 ENDIF
              ENDDO
              !PRINT *, "ADJ_FORCE", ADJ_FORCE(I,J,:,IDTCO)
              COST_FUNC = COST_FUNC + NEW_COST(I,J)
           ENDIF
        ENDDO
     ENDDO

110  FORMAT(F18.6,1X)
     !PRINT *, "GC_ADJ_TEMP", GC_ADJ_TEMP(:,:,5)
      !PRINT *, "STT_ADJ BEFORE", STT_ADJ(:,:,5,IDTCO)
     !PRINT *, "GC_ADJ_TEMP", GC_ADJ_TEMP(:,:,5)/SOBS_COUNT(:,:)
     !PRINT *, "STT_ADJ AFTER",STT_ADJ(:,:,5,IDTCO)
     
     
     IF ( FIRST ) FIRST = .FALSE.
     
     !     Update cost function
     !COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))
     
     print*, ' Updated value of COST_FUNC = ', COST_FUNC
     print*, ' IASI CO contribution           = ', COST_FUNC - OLD_COST
     
     ! Return to calling program
   END SUBROUTINE CALC_IASI_CO_FORCE
    
!----------------------------------------------------------------------------
    SUBROUTINE CLEANUP_IASI

      IF(ALLOCATED(IASI_LON)) DEALLOCATE(IASI_LON)
      IF(ALLOCATED(IASI_LAT)) DEALLOCATE(IASI_LAT)
      IF(ALLOCATED(IASI_TIME)) DEALLOCATE(IASI_TIME)
      IF(ALLOCATED(IASI_CO_AVK)) DEALLOCATE(IASI_CO_AVK)
      IF(ALLOCATED(IASI_CO_APR)) DEALLOCATE(IASI_CO_APR)
      IF(ALLOCATED(IASI_CO_STD)) DEALLOCATE(IASI_CO_STD)
      IF(ALLOCATED(IASI_CO)) DEALLOCATE(IASI_CO)
      IF(ALLOCATED(IASI_ALT)) DEALLOCATE(IASI_ALT)
      IF(ALLOCATED(IASI_SOLAR_ZENITH)) DEALLOCATE(IASI_SOLAR_ZENITH)
      IF(ALLOCATED(IASI_CLOUD_COVER)) DEALLOCATE(IASI_CLOUD_COVER)
      IF(ALLOCATED(IASI_QUAL_FLAG)) DEALLOCATE(IASI_QUAL_FLAG)
      IF(ALLOCATED(IASI_CO_LVL)) DEALLOCATE(IASI_CO_LVL)
      !IF(ALLOCATED(TEMPDATA)) DEALLOCATE(TEMPDATA)
      !IF(ALLOCATED(TIME_FRAC)) DEALLOCATE(TIME_FRAC)
    END SUBROUTINE CLEANUP_IASI
!------------------------------------------------------------------------------
      SUBROUTINE GET_NT_RANGE( N_IASI_NOB, HHMMSS, TIME_FRAC, NTSTART, NTSTOP)
!
!******************************************************************************
!  Subroutine GET_NT_RANGE retuns the range of retrieval records for the
!  current model hour
!
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTES   (INTEGER) : Number of TES retrievals in this day
!  (2 ) HHMMSS (INTEGER) : Current model time
!  (3 ) TIME_FRAC (REAL) : Vector of times (frac-of-day) for the TES retrievals
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP  (INTEGER) : TES record number at which to stop
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : YMD_EXTRACT

      ! Arguments
      INTEGER, INTENT(IN)   :: N_IASI_NOB
      INTEGER, INTENT(IN)   :: HHMMSS
      REAL*8,  INTENT(IN)   :: TIME_FRAC(N_IASI_NOB)
      INTEGER, INTENT(OUT)  :: NTSTART
      INTEGER, INTENT(OUT)  :: NTSTOP

      ! Local variables
      INTEGER, SAVE         :: NTSAVE
      LOGICAL               :: FOUND_ALL_RECORDS, FOUND_BAD_RECORDS
      INTEGER               :: NTEST
      INTEGER               :: HH, MM, SS
      REAL*8                :: GC_HH_FRAC
      REAL*8                :: H1_FRAC

      !=================================================================
      ! GET_NT_RANGE begins here!
      !=================================================================
      ! Initialize
      FOUND_ALL_RECORDS  = .FALSE.
      FOUND_BAD_RECORDS  = .TRUE.
      NTSTART            = 0
      NTSTOP             = 0

      ! set NTSAVE to NTES every time we start with a new file
      IF ( HHMMSS == 230000 ) THEN 
         NTSAVE = N_IASI_NOB
         IF ( NTSAVE > 1 ) THEN
            DO WHILE (FOUND_BAD_RECORDS == .TRUE.)                  
               IF (TIME_FRAC(NTSAVE) > 0.5d0) THEN
                  FOUND_BAD_RECORDS = .FALSE.
               ELSE
                  NTSAVE = NTSAVE - 1
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      DO WHILE (TIME_FRAC(NTSAVE) < 0)
         NTSAVE = NTSAVE -1
         IF (NTSAVE == 0) EXIT
      ENDDO
      !print*, ' GET_NT_RANGE for ', HHMMSS
      !print*, ' NTSAVE ', NTSAVE
      !print*, ' N_IASI_NOB   ', N_IASI_NOB

      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


      ! Convert HH from hour to fraction of day
      GC_HH_FRAC = REAL(HH,8) / 24d0
      ! one hour as a fraction of day
      H1_FRAC    = 0d0 / 24d0


      ! All records have been read already
      IF ( NTSAVE == 0 ) THEN

         !print*, 'All records have been read already '
         RETURN
      ! No records reached yet
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC < GC_HH_FRAC ) THEN
         !PRINT *, "TIME_FRAC", TIME_FRAC(NTSAVE)
         
         print*, 'No records reached yet'
         RETURN

      !
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC >=  GC_HH_FRAC ) THEN
         ! Starting record found
         NTSTART = NTSAVE

         !print*, ' Starting : TIME_FRAC(NTSTART) ', TIME_FRAC(NTSTART), NTSTART

         ! Now search forward to find stopping record
         NTEST = NTSTART

         DO WHILE ( FOUND_ALL_RECORDS == .FALSE. )

            ! Advance to the next record
            NTEST = NTEST - 1

            ! Stop if we reach the earliest available record
            IF ( NTEST == 0 ) THEN

               NTSTOP            = NTEST + 1
               FOUND_ALL_RECORDS = .TRUE.

               !print*, ' Records found '
               !print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP
               ! Reset NTSAVE
               NTSAVE = NTEST

            ! When the combined test date rounded up to the nearest
            ! half hour is smaller than the current model date, the
            ! stopping record has been passed.
            ELSEIF (  TIME_FRAC(NTEST) + H1_FRAC <  GC_HH_FRAC ) THEN

               !print*, ' Testing : TIME_FRAC ', TIME_FRAC(NTEST), NTEST

               NTSTOP            = NTEST + 1
               FOUND_ALL_RECORDS = .TRUE.

               !print*, ' Records found '
               !print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP

               ! Reset NTSAVE
               NTSAVE = NTEST
            !ELSE
               !print*, ' still looking ', NTEST

            ENDIF

         ENDDO

      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' )

      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE

!------------------------------------------------------------------------------------

      SUBROUTINE BIN_DATA_IASI( GC_EDGE, OBS_EDGE, DATA_MODEL, DATA_IASI, OBS_IASI, LIASI, FB )

!******************************************************************************
!Based on the code from Monika.  (zhe 1/19/11)
!FB = 1 for forward
!FB = -1 for adjoint
!******************************************************************************

      INTEGER :: L, LL, FB
      INTEGER :: LIASI, NB, LVL_CRT1, LVL_CRT2
      REAL*8  :: ALT_MODEL(LLPAR), GC_EDGE(LLPAR+1), HI, LOW
      REAL*8  :: DATA_MODEL(LLPAR), BIN_IASI(LLPAR,LIASI), DIFF_BIN
      REAL*8  :: OBS_EDGE(LIASI+1), DATA_IASI(LIASI), OBS_IASI(LIASI)
      BIN_IASI(:,:) = 0d0
      LVL_CRT1 = 0
      LVL_CRT2 = 0
      !=================================================================
      ! BIN_DATA_V4 begins here!
      !=================================================================
      DO L = 1, LLPAR
         ALT_MODEL(L) = 0.5d0*(GC_EDGE(L)+GC_EDGE(L+1))
      ENDDO
      IF (FB > 0) THEN

         !DO L = 1, LIASI
            !DO LL = 1, LLPAR
               !IF ( ALT_MODEL(LL) >= OBS_EDGE(L) ) THEN
                  !DATA_IASI(L) = DATA_MODEL(LL)
                  !EXIT
               !ENDIF
            !ENDDO
         !ENDDO
         
         DO L = 1, LIASI
            DO LL = 1, LLPAR
               LOW = GC_EDGE(LL)
               HI = GC_EDGE(LL+1)
               IF ( GC_EDGE(LL) >= OBS_EDGE(L)) THEN
                  IF ( GC_EDGE(LL+1) <= OBS_EDGE(L+1) ) THEN
                     BIN_IASI(LL,L) = 1d0 
                     !NB = NB + 1
                  ELSEIF (GC_EDGE(LL) <= OBS_EDGE(L+1)) THEN
                     DIFF_BIN = HI - LOW
                     BIN_IASI(LL,L) = ( OBS_EDGE(L+1) - LOW)/DIFF_BIN
                     BIN_IASI(LL,L+1) = ( HI - OBS_EDGE(L+1))/DIFF_BIN
                  ELSEIF (GC_EDGE(LL) > OBS_EDGE(LIASI+1)) THEN
                     BIN_IASI(LL,LIASI) = 1d0
                  ENDIF
               ELSEIF (GC_EDGE(LL) < OBS_EDGE(1) ) THEN
                  BIN_IASI(LL,1) = 1D0
               ENDIF
            ENDDO
            !IF (NB > 0) DATA_IASI(L) = DATA_TEM !/ NB
         ENDDO
         
         DO L = 1, LIASI
            DATA_IASI(L) = 0d0
            DO LL = 1, LLPAR
               DATA_IASI(L) = DATA_IASI(L) + BIN_IASI(LL,L) * DATA_MODEL(LL)
            ENDDO
         ENDDO
         DO L = 2, LIASI-1
            IF (DATA_IASI(L) == 0d0) THEN
               IF ( DATA_IASI(L-1) > 0d0) THEN
                  LVL_CRT1 = L-1
                  !PRINT *, "DATA_IASI1", DATA_IASI(L-1)
               ENDIF
               IF (DATA_IASI(L+1) > 0d0) THEN
                  LVL_CRT2 = L+1
                  !PRINT *, "DATA_IASI2", DATA_IASI(L+1)
               ENDIF
               IF (REAL(LVL_CRT1)*REAL(LVL_CRT2) > 0D0) THEN
                  DO LL = LVL_CRT1, LVL_CRT2
                     DATA_IASI(LL) = ((DATA_IASI(LVL_CRT1)+DATA_IASI(LVL_CRT2))/SUM(OBS_IASI(LVL_CRT1:LVL_CRT2)))*OBS_IASI(LL)
                  ENDDO
                  !PRINT *, "OBS_IASI", SUM(OBS_IASI(LVL_CRT1:LVL_CRT2))
                  LVL_CRT1 = 0
                  LVL_CRT2 = 0
               ENDIF
            ENDIF
         ENDDO

      ELSE

         DATA_MODEL(:) = 0.
         DO L = 1, LLPAR
            DO LL = 1, LIASI
               IF ( ( ALT_MODEL(L) >= OBS_EDGE(LL)) .and. ( ALT_MODEL(L) < OBS_EDGE(LL+1)) ) THEN
                  DATA_MODEL(L) = DATA_IASI(LL)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      ! Return to calling program
      END SUBROUTINE BIN_DATA_IASI




      FUNCTION GET_IJ_2x25( LON, LAT ) RESULT ( IIJJ )

!
!******************************************************************************
!  Subroutine GET_IJ_2x25 returns I and J index from the 2 x 2.5 grid for a 
!  LON, LAT coord. (dkh, 11/08/09) 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LON (REAL*8) : Longitude                          [degrees]
!  (2 ) LAT (REAL*8) : Latitude                           [degrees]
!     
!  Function result
!  ============================================================================
!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!     
!  NOIASI:
!
!******************************************************************************
!     
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP

      ! Arguments
      REAL*4    :: LAT, LON
      
      ! Return
      INTEGER :: I, J, IIJJ(2)
      
      ! Local variables 
      REAL*8              :: TLON, TLAT, DLON, DLAT
      REAL*8,  PARAMETER  :: DISIZE = 2.5d0
      REAL*8,  PARAMETER  :: DJSIZE = 2.0d0
      INTEGER, PARAMETER  :: IIMAX  = 144
      INTEGER, PARAMETER  :: JJMAX  = 91
      
      
      !=================================================================
      ! GET_IJ_2x25 begins here!
      !=================================================================

      TLON = 180d0 + LON + DISIZE
      TLAT =  90d0 + LAT + DJSIZE
      
      I = TLON / DISIZE
      J = TLAT / DJSIZE

      
      IF ( TLON / DISIZE - REAL(I)  >= 0.5d0 ) THEN
         I = I + 1
      ENDIF
      
      IF ( TLAT / DJSIZE - REAL(J)  >= 0.5d0 ) THEN
         J = J + 1
      ENDIF

      
      ! Longitude wraps around
      !IF ( I == 73 ) I = 1 
      IF ( I == ( IIMAX + 1 ) ) I = 1
      
      ! Check for impossible values 
      IF ( I > IIMAX .or. J > JJMAX .or. I < 1     .or. J < 1 ) THEN
         CALL ERROR_STOP('Error finding grid box', 'GET_IJ_2x25')
      ENDIF
      
      IIJJ(1) = I
      IIJJ(2) = J
      
      ! Return to calling program
      END FUNCTION GET_IJ_2x25
!-------------------------------------------------------------------
      
      END MODULE IASI_CO_OBS_MOD
