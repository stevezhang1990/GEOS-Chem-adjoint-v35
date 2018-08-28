! $Id: osiris_obs_mod.f,v 1.0 2012/02/23 06:47:07 twalker Exp $
MODULE OSIRIS_OBS_MOD

!******************************************************************************
! Module OSIRIS_OBS_MOD contains subroutines necessary to 
! 1. Read OSIRIS (ASCII) file with O3 observations, preprocessed onto model grid
! 2. Determine when OSIRIS O3 obs are available
! 3. Compute adjoint forcing in model space
!
!  Module Variables:
!  ============================================================================
!  (1 ) OSIRIS_DATA       : OSIRIS data pre-averaged onto 4x5 grid
!  (2 ) ADJ_FORCE_OSIRIS  : Adjoint forcing
!  (3 ) LOCATION_DATA     : Array of locations and times of observations
!  (4 ) FILEDATE          : Date associated with file currently read
!  (5 ) COUNT_TOTAL       : Number of observations
!  (6 ) OSIRIS_ERR        : OSIRIS error values pre-averaged onto 4x5 grid
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_OSIRIS_FILE           : Read OSIRIS ASCII file
!  (2 ) ITS_TIME_FOR_OSIRIS_OBS    : Checks model time
!  (3 ) CALC_OSIRIS_FORCE          : Calculates cost fnc and ADJ_STT increments
!  (4 ) INIT_OSIRIS                : Allocates memory of arrays
!  (5 ) CLEANUP_OSIRIS             : Deallocates memory of arrays
!  (6 ) IS_OSIRIS_NONZERO         : Determines if OSIRIS observed a grid
!
!  ============================================================================
!  NOTES: 
!  (1 ) Based on obs operators implemented in v8 adjoint. (tww, 20120223)
!
!******************************************************************************

  IMPLICIT NONE

#include "CMN_SIZE"   ! Size parameters

  ! Everything PRIVATE unless specified otherwise
  ! PRIVATE module variables
  ! PRIVATE module routines
  PRIVATE 

  PUBLIC :: READ_OSIRIS_FILE
  PUBLIC :: ITS_TIME_FOR_OSIRIS_OBS
  PUBLIC :: CALC_OSIRIS_FORCE
  PUBLIC :: COUNT_TOTAL
  PUBLIC :: IS_OSIRIS_NONZERO
  PUBLIC :: CALC_GC_O3

  REAL*8, ALLOCATABLE   :: OSIRIS_DATA(:,:)
  REAL*8, ALLOCATABLE   :: OSIRIS_ERR(:,:)
  REAL*8, ALLOCATABLE   :: ADJ_FORCE_OSIRIS(:,:,:)
  INTEGER, ALLOCATABLE  :: LOCATION_DATA(:,:)
  INTEGER, ALLOCATABLE  :: TIME_DATA(:)
  
  INTEGER               :: FILEDATE
  INTEGER               :: NOCC
  REAL*8                :: COUNT_TOTAL
  INTEGER, PARAMETER    :: MAX_OBS_PER_DAY=1000

CONTAINS


!----------------------------------------------------------------------

  SUBROUTINE READ_OSIRIS_FILE( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_OSIRIS_FILE reads an ASCII file containing OSIRIS data
!  for the given day.  Data are already averaged onto 4x5 grid, but are given
!  every km altitude from 0.5km to 64.5km. (tww, 20120223)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   : Hour-Min-Sec
!

    USE DAO_MOD,              ONLY : BXHEIGHT
    USE FILE_MOD,             ONLY : IOERROR
    USE TIME_MOD,             ONLY : EXPAND_DATE
    USE DIRECTORY_ADJ_MOD,    ONLY : DIAGADJ_DIR
    USE ADJ_ARRAYS_MOD,       ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,       ONLY : EXPAND_NAME

#include "CMN_SIZE" ! size parameters

    ! Arguments
    INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

    ! Local Variables
    CHARACTER(LEN=255)  :: DIR_OSIRIS
    CHARACTER(LEN=255)  :: DIR_MONTH_OSIRIS
    CHARACTER(LEN=255)  :: FILENAME_OSIRIS
    CHARACTER(LEN=255)  :: FILENAME

    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(I8,I5,I5,I5,F8.2,F8.2,F8.2,I9)"
    INTEGER,          PARAMETER  :: OSIRIS_LVL=65
  
    REAL*8              :: TEMPALT(OSIRIS_LVL)
    REAL*8              :: TEMPO3DATA(OSIRIS_LVL)
    REAL*8              :: TEMPO3ERR(OSIRIS_LVL)
    REAL*8              :: TEMPLAT, TEMPLON, TEMPH, Z1, Z2
    INTEGER             :: IU_FILE, IOS
    INTEGER             :: NOBS, I, J, L, K
    INTEGER             :: OCCID, OLMAX, LATIND, LONIND
    INTEGER             :: TEMPYMD
    INTEGER             :: NDAT
    LOGICAL, SAVE       :: FIRST = .TRUE.
    CHARACTER(LEN=255)  :: FILENAME_DIAG


    !=================================================================
    ! READ_OSIRIS_FILE begins here!
    !=================================================================

    !=============================
    ! FIRST CLEANUP IF NECESSARY:
    !=============================
    CALL CLEANUP_OSIRIS
    CALL INIT_OSIRIS

    IF ( FIRST ) THEN
       COUNT_TOTAL = 0

       FILENAME_DIAG = 'lat_orb_osi.NN.m'
       CALL EXPAND_NAME( FILENAME_DIAG, N_CALC )
       FILENAME_DIAG =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME_DIAG )
       OPEN( 405,      FILE=TRIM( FILENAME_DIAG  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME_DIAG = 'lon_orb_osi.NN.m'
       CALL EXPAND_NAME( FILENAME_DIAG, N_CALC )
       FILENAME_DIAG =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME_DIAG )
       OPEN( 406,      FILE=TRIM( FILENAME_DIAG  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FIRST = .FALSE.
    ENDIF
  
    !========================
    ! FILENAME
    !=========================

    DIR_OSIRIS = '/users/jk/15/xzhang/OSIRIS_O3/'
    DIR_MONTH_OSIRIS = '/YYYY/MM/'
    FILENAME_OSIRIS = TRIM( 'OSIRIS_507_4x5_YYYYMMDD_O3.data' )
    !      FILENAME_OSIRIS = TRIM( 'osiristestfile_YYYYMMDD.data' )
    IU_FILE = 15
    ! EXPAND_DATE replaces tokens like 'YYYY' with the year
    CALL EXPAND_DATE( DIR_MONTH_OSIRIS, YYYYMMDD, 9999 )
    CALL EXPAND_DATE( FILENAME_OSIRIS, YYYYMMDD, 9999 )
    FILENAME = TRIM( DIR_OSIRIS ) // TRIM(DIR_MONTH_OSIRIS) // FILENAME_OSIRIS

    WRITE(6,*) ' - READ_OSIRIS_FILE: reading: ', FILENAME

    NOCC = 0

    ! Open file
    OPEN( IU_FILE, FILE=FILENAME, IOSTAT=IOS )
    IF ( IOS == 0 ) RETURN
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'osiris:1' )
      
    !========================
    ! Read in data blocks
    !========================
    ! Needs to be true until end of file
    DO I = 1, MAX_OBS_PER_DAY
       
       IF ( IOS < 0 ) EXIT
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'osiris:2' )

       ! Read meta data line (OCCID, RETLVLS, LONIND, 
       !                      LATIND, LON, LAT, HH, YMD)
       READ( IU_FILE, FMT1, IOSTAT=IOS) OCCID, OLMAX, LONIND, LATIND, TEMPLON, TEMPLAT, TEMPH, TEMPYMD

       IF ( IOS < 0 ) EXIT
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'osiris:3' )

       !DEBUG, tww
       !print*, 'METADATA is: ', OCCID, OLMAX, LONIND, LATIND, TEMPLON, TEMPLAT, TEMPH, TEMPYMD
 
       LOCATION_DATA(I,1) = LONIND
       LOCATION_DATA(I,2) = LATIND
       TIME_DATA(I) = FLOOR(TEMPH)*10000d0
 
       ! Read altitude, ozone and error profiles
       READ( IU_FILE, * )(TEMPALT(K), K=1,OLMAX)
       READ( IU_FILE, * )(TEMPO3DATA(K), K=1,OLMAX)
       READ( IU_FILE, * )(TEMPO3ERR(K), K=1,OLMAX)
 
       ! OSIRIS data comes every km
       ! Compute average of data on each GEOS-Chem vertical level
       Z1 = 0d0
       Z2 = 0d0
       DO L = 1, LLPAR
          Z1 = Z2
          Z2 = Z1 + BXHEIGHT(LONIND,LATIND,L)/1000d0
          NDAT = 0
          DO K = 1, OSIRIS_LVL
             IF (TEMPALT(K) > Z1 .AND. TEMPALT(K) <=Z2) THEN
                IF (TEMPO3DATA(K)>0d0) THEN
                   OSIRIS_DATA( I, L ) = OSIRIS_DATA( I, L ) + TEMPO3DATA( K )
                   OSIRIS_ERR( I, L ) = OSIRIS_ERR( I, L ) + TEMPO3ERR( K ) ** 2
                   NDAT = NDAT + 1
                ENDIF
             ELSEIF (TEMPALT(K) > Z2) THEN
                EXIT
             ENDIF
          ENDDO
          IF (NDAT>0) THEN
             OSIRIS_DATA(I,L) = OSIRIS_DATA(I,L)/NDAT
             OSIRIS_ERR(I,L) = SQRT(OSIRIS_ERR(I,L))/NDAT
          ENDIF
       ENDDO

       NOCC = NOCC + 1
       !debug, tww
       !print*, 'o3data read: ', TEMPO3DATA
       !print*, 'o3err read: ', TEMPO3ERR
       WRITE(405,118) ( TEMPLAT             )
       WRITE(406,118) ( TEMPLON             )
118       FORMAT(F18.6,1X)
    ENDDO

    FILEDATE = TEMPYMD !- 10000
    
    ! Close file
    CLOSE( IU_FILE )

    ! DEBUG, tww
    !print*, 'FINISHED READING OSIRIS FILE'
    !print*, 'NOCC = ', NOCC
    !print*, 'FILEDATE = ', FILEDATE
    !print*, 'LOCATIONS = ', LOCATION_DATA
    !print*, 'TIMES = ', TIME_DATA
    !print*, 'OSIRIS_DATA = ', OSIRIS_DATA(1:3,:)
    !print*, 'OSIRIS_ERR = ', OSIRIS_ERR(1:3,:)

  END SUBROUTINE READ_OSIRIS_FILE


!----------------------------------------------------------------------

  FUNCTION ITS_TIME_FOR_OSIRIS_OBS( ) RESULT( FLAG )

!******************************************************************************
!  Function ITS_TIME_FOR_OSIRIS_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day). (tww, 20120223)
!

    USE TIME_MOD, ONLY : GET_NYMD, GET_NHMS

#include "CMN_SIZE"   ! Size parameters

    ! Function value
    LOGICAL :: FLAG
    
    INTEGER :: N

    !=================================================================
    ! ITS_TIME_FOR_OSIRIS_OBS begins here!
    !=================================================================

    ! Default to false
    FLAG = .FALSE.
  
    ! Observation times on this day are stored in TIME_DATA
    DO N = 1, NOCC
       IF( TIME_DATA( N ) == GET_NHMS() ) THEN
          ! DEBUG
          !WRITE(6,*) ' - ITS_TIME_FOR_OSIRIS_OBS found at: ',N
          !WRITE(6,*) ' - ITS_TIME_FOR_OSIRIS_OBS found: ', GET_NHMS()
          FLAG = .TRUE.
       ENDIF
    ENDDO
  
    ! If we have the wrong day
    IF( GET_NYMD() /= FILEDATE ) THEN
       WRITE(6,*) ' - ITS_TIME_FOR_OSIRIS_OBS wrong day: ', FILEDATE
       WRITE(6,*) ' - ITS_TIME_FOR_OSIRIS_OBS wrong day: ', GET_NYMD()
       FLAG = .FALSE.
    ENDIF

  END FUNCTION ITS_TIME_FOR_OSIRIS_OBS


!----------------------------------------------------------------------

  FUNCTION IS_OSIRIS_NONZERO( I,J,L ) RESULT ( FLAG )

!******************************************************************************
!  Function IS_OSIRIS_NONZERO returns TRUE if there are observations
!  available for particular location. (tww, 20120229)
!

    USE TIME_MOD, ONLY : GET_NYMD, GET_NHMS

#include "CMN_SIZE"   ! Size parameters

    ! Function value
    LOGICAL :: FLAG

    ! Arguments
    INTEGER :: I, J, L

    INTEGER :: N

    !=================================================================
    ! IS_OSIRIS_NONZERO begins here!
    !=================================================================

    ! Default to false
    FLAG = .FALSE.
      
    DO N = 1, NOCC
       IF( TIME_DATA( N ) == GET_NHMS() ) THEN
          IF ( ( LOCATION_DATA(N,1)==I ) .AND. &
               ( LOCATION_DATA(N,2)==J ) .AND. &
               ( OSIRIS_DATA(N,L) > 0 )) THEN
             WRITE(6,*) 'IS_OSIRIS_NONZERO - yes, at: ', I, J, L
             FLAG = .TRUE.
          ENDIF
       ENDIF
    ENDDO

    ! If we have the wrong day
    IF( GET_NYMD() /= FILEDATE ) THEN
       WRITE(6,*) ' - IS_OSIRIS_NONZERO wrong day: ', FILEDATE
       WRITE(6,*) ' - IS_OSIRIS_NONZERO wrong day: ', GET_NYMD()
       FLAG = .FALSE.
    ENDIF

  END FUNCTION IS_OSIRIS_NONZERO


!----------------------------------------------------------------------

  SUBROUTINE CALC_OSIRIS_FORCE( COST_FUNC )

!******************************************************************************
!  Subroutine CALC_OSIRIS_FORCE (tww, 20120223)
!
    USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ
    USE ADJ_ARRAYS_MOD,       ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,       ONLY : EXPAND_NAME
    USE CHECKPT_MOD,          ONLY : CHK_STT
    USE DAO_MOD,              ONLY : AD
    USE ERROR_MOD,            ONLY : IT_IS_NAN, ERROR_STOP
    USE TIME_MOD,             ONLY : GET_NHMS
    USE TRACERID_MOD,         ONLY : IDTOX
    USE PRESSURE_MOD,         ONLY : GET_PCENTER
    USE DIRECTORY_ADJ_MOD,    ONLY : DIAGADJ_DIR

#include "CMN_SIZE"    ! Size parameters

    ! Arguments
    REAL*8, INTENT(INOUT)       :: COST_FUNC
    
    ! Parameters
    REAL*8, PARAMETER           :: ADJ_TCVVOX = 28.97d0/48.d0
    
    REAL*8                      :: NEW_COST(IIPAR,JJPAR)
    REAL*8                      :: GC_PRES(LLPAR), DIFF(LLPAR)
    REAL*8                      :: OBS_ERRCOV(LLPAR)
    REAL*8                      :: GC_O3(LLPAR)
    REAL*8                      :: GC_ADJ_COUNT(IIPAR,JJPAR,LLPAR), COST_CONT(LLPAR)
    INTEGER                     :: I, J, L, N
    INTEGER                     :: THISYMD
    
    LOGICAL :: SUPER_OBS = .TRUE.
    REAL*8 :: SOBS_COUNT(IIPAR,JJPAR)
    REAL*8 :: SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)
    LOGICAL, SAVE               :: FIRST = .TRUE.
    INTEGER                     :: IOS
    CHARACTER(LEN=255)          :: FILENAME

    IF ( FIRST ) THEN
       FILENAME = 'gc_press_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 401,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'diff_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 402,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
         
       FILENAME = 'gc_o3_stt_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 403,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'gc_o3_stt_adj_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 404,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'obs_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 409,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'err_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 412,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'cfn_l_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 413,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )


    ENDIF
   
    !================================================================ 
    ! CALC_OSIRIS_FORCE begins here!
    !================================================================
      
! this comes from old CALC_ADJ_FORCE, same for v7 and v8 (tww, 20101027)

    ! Initialize to be safe
    NEW_COST(:,:) = 0d0 
    SOBS_COUNT(:,:) = 0d0
    SOBS_ADJ_FORCE(:,:,:) = 0d0
    GC_ADJ_COUNT(:,:,:) = 0d0
    
    


    DO N = 1, NOCC

       DIFF(:) = 0d0
       OBS_ERRCOV(:) = 0d0
       GC_PRES(:) = 0d0
       GC_O3(:) = 0d0
       COST_CONT = 0d0

       ! Only get obs at this time
       IF( TIME_DATA(N) == GET_NHMS() ) THEN

          I = LOCATION_DATA(N,1)
          J = LOCATION_DATA(N,2)
          DO L = 1, LLPAR
             GC_PRES(L) = GET_PCENTER(I,J,L)
             ! DEBUG, tww
             !print*, 'ADDING FORCING AT: ', I, J, L
             GC_O3(L) = CHK_STT(I,J,L,IDTOX) * ADJ_TCVVOX / AD(I,J,L) * 1d9
             ! Make sure data and error are not zero or fill
             IF( ( OSIRIS_DATA(N,L) > 0d0    ) .AND. &
                 ( OSIRIS_ERR(N,L) > 0d0   ) ) THEN
                ! Only force below level 40 (about 25km)
                ! Try without this condition (tww, 20120317)
                !IF( L < 40 ) THEN
                ! Condition data with very small errors
                ! errors smaller than 1% on a single obs are removed

                IF( OSIRIS_DATA(N,L) < 100d0 * OSIRIS_ERR(N,L) ) THEN 
                   ! CHK_STT is in units of [kg/box] here. Convert to ppb
                   DIFF(L) = ( GC_O3(L) - OSIRIS_DATA(N,L) * 1d9 )
                   
                   ! Get obs error covariance
                   OBS_ERRCOV(L) =  OSIRIS_ERR(N,L) * OSIRIS_ERR(N,L) * 1d18

                   COST_CONT(L) = 0.5d0 * (DIFF(L)**2)/ OBS_ERRCOV(L)
                   ! Calculate new additions to cost function
                   IF ( ( COST_CONT(L) > 0d0) .AND. &
                        ( GET_PCENTER(I,J,L) < 300d0) ) THEN
                      IF (SUPER_OBS) THEN
                         NEW_COST(I,J)  = NEW_COST(I,J) + COST_CONT(L)
                         SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                      ELSE
                         COST_FUNC = COST_FUNC + COST_CONT(L)
                      ENDIF
 
                      ! Force the adjoint variables x with dJ/dx
                      ! Change to get units right [kg/box]
                      ADJ_FORCE_OSIRIS(I,J,L) = DIFF(L) / OBS_ERRCOV(L) * ADJ_TCVVOX / AD(I,J,L) * 1d9
                      IF (SUPER_OBS) THEN
                         SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + ADJ_FORCE_OSIRIS(I,J,L)
                         GC_ADJ_COUNT(I,J,L) = GC_ADJ_COUNT(I,J,L) + 1
                      ELSE
                         STT_ADJ(I,J,L,IDTOX) = STT_ADJ(I,J,L,IDTOX) + ADJ_FORCE_OSIRIS(I,J,L)
                      ENDIF
                   ENDIF
                ELSE
                   !print*, 'removed outlier at ', N, L
                   !print*, 'outlier is ', OSIRIS_DATA(N,L)
                   !print*, 'outlier error is ', OSIRIS_ERR(N,L)
                ENDIF
                !ENDIF
             ENDIF
         
          ENDDO
          !SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
          !WRITE(6,102) (DIFF(L), L=LLPAR,1,-1)
102       FORMAT(1X,d14.6)
          !WRITE(401,110) ( GC_PRES(L),       L=LLPAR,1,-1 )
          !WRITE(402,110) ( DIFF(L),          L=LLPAR,1,-1 )
          !WRITE(403,110) ( GC_O3(L),         L=LLPAR,1,-1 )
          !WRITE(404,110) ( GC_STT_ADJ(L),    L=LLPAR,1,-1 )
          !WRITE(409,110) ( OSIRIS_DATA(N,L) * 1d9,    L=LLPAR,1,-1 )
          !WRITE(412,110) ( OSIRIS_ERR(N,L) * 1d9,     L=LLPAR,1,-1 )
          !WRITE(413,110) ( COST_CONT(L),     L=LLPAR,1,-1 )

110       FORMAT(F18.6,1X)
       ENDIF
    ENDDO
    IF (SUPER_OBS) THEN
       DO I=1,IIPAR
          DO J=1,JJPAR
             IF ( SOBS_COUNT(I,J) > 0d0 ) THEN
                DO L=1,LLPAR
                   IF ( ( GET_PCENTER(I,J,L) < 300d0 )  .AND. &
                        ( GC_ADJ_COUNT(I,J,L) > 0d0 )  ) THEN
                      STT_ADJ(I,J,L,IDTOX) = STT_ADJ(I,J,L,IDTOX) + SOBS_ADJ_FORCE(I,J,L)/GC_ADJ_COUNT(I,J,L)
                   ENDIF
                   
                ENDDO
                COST_FUNC = COST_FUNC + NEW_COST(I,J)/SOBS_COUNT(I,J)
             ENDIF
          ENDDO
       ENDDO
    ELSE
       WRITE(6,*) ' CALC_OSIRIS_FORCE: NEW_COST = ', SUM( NEW_COST(:,:))
       ! Update cost function
       COST_FUNC = COST_FUNC + SUM ( NEW_COST(:,:) )
       COUNT_TOTAL = COUNT_TOTAL + SUM (SOBS_COUNT(:,:)  )
       print*, ' Total observation number = ', COUNT_TOTAL
    ENDIF
    ! dkh debug
    !print*, ' CHK_STT  = ', CHK_STT(25,35,1:8,IDTOX)
    !print*, ' AD = ', AD(25,35,1:8)
    !print*, ' OSIRIS_DATA = ', OSIRIS_DATA(25,35,1:8)
    !print*, ' OSIRIS_ERR = ', OSIRIS_ERR(25,35,1:8)
    !print*, ' NEW_COST = ', NEW_COST(25,35,1:8)

    ! Error check
    IF ( IT_IS_NAN( COST_FUNC ) ) THEN
       CALL ERROR_STOP( 'COST_FUNC IS NaN', 'CALC_OSIRIS_FORCE')
    ENDIF

  END SUBROUTINE CALC_OSIRIS_FORCE


!-----------------------------------------------------------------------------

  SUBROUTINE INIT_OSIRIS
!
!*****************************************************************************
!  Subroutine INIT_OSIRIS allocates all module arrays. (dkh, 11/16/06)
!
!  NOTES:
!
!******************************************************************************
!
    USE ERROR_MOD,  ONLY : ALLOC_ERR

#include "CMN_SIZE"   ! IIPAR, JJPAR, LLPAR

    ! Local variables
    INTEGER :: AS

    !=================================================================
    ! INIT_OSIRIS begins here
    !=================================================================
    ALLOCATE( LOCATION_DATA( MAX_OBS_PER_DAY, 2 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOCATION_DATA' )
    LOCATION_DATA = 0

    ALLOCATE( TIME_DATA( MAX_OBS_PER_DAY ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TIME_DATA' )
    TIME_DATA = 0

    ALLOCATE( OSIRIS_DATA( MAX_OBS_PER_DAY, LLPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OSIRIS_DATA' )
    OSIRIS_DATA = 0d0

    ALLOCATE( OSIRIS_ERR( MAX_OBS_PER_DAY, LLPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OSIRIS_ERR' )
    OSIRIS_ERR = 0d0

    ALLOCATE( ADJ_FORCE_OSIRIS( IIPAR, JJPAR, LLPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ADJ_FORCE_OSIRIS' )
    ADJ_FORCE_OSIRIS = 0d0
  
    FILEDATE = 0
    NOCC = 0

  END SUBROUTINE INIT_OSIRIS


!------------------------------------------------------------------------------

  SUBROUTINE CLEANUP_OSIRIS
  
!******************************************************************************
!  Deallocate all memory (done before reading each monthly file)
!

    ! Deallocate
    IF ( ALLOCATED( LOCATION_DATA    ) ) DEALLOCATE( LOCATION_DATA    )
    IF ( ALLOCATED( TIME_DATA        ) ) DEALLOCATE( TIME_DATA        )
    IF ( ALLOCATED( OSIRIS_DATA      ) ) DEALLOCATE( OSIRIS_DATA      )
    IF ( ALLOCATED( OSIRIS_ERR       ) ) DEALLOCATE( OSIRIS_ERR       )
    IF ( ALLOCATED( ADJ_FORCE_OSIRIS ) ) DEALLOCATE( ADJ_FORCE_OSIRIS )

  END SUBROUTINE CLEANUP_OSIRIS

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  SUBROUTINE CALC_GC_O3

!!
!! Subroutine CALC_OMI_O3_FORCE computes the O3 adjoint forcing from OMI column data
!!

    USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
    USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
    USE CHECKPT_MOD,        ONLY : CHK_STT
    USE DAO_MOD,            ONLY : AD
    USE TIME_MOD,           ONLY : GET_HOUR
    USE TRACERID_MOD,       ONLY : IDO3, IDTOX
    USE TRACER_MOD,         ONLY : TCVV
    USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR

#include "CMN_SIZE" ! size parameters

    INTEGER :: I,J,L

    INTEGER :: GC_HOUR

    ! variables for observation operator and adjoint thereof

    REAL*8 :: GC_O3(LLPAR)
    LOGICAL, SAVE               :: FIRST = .TRUE.
    INTEGER                     :: IOS
    CHARACTER(LEN=255)          :: FILENAME

    IF ( FIRST ) THEN
       FILENAME = 'stt_o3_osi.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 407,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    ENDIF

    GC_HOUR = GET_HOUR()
    GC_O3 = 0d0

    DO I = 1, IIPAR
       DO J = 1, JJPAR
          GC_O3 = 0d0
          DO L = 1, LLPAR
             GC_O3(L) =  CHK_STT(I,J,L,IDTOX) * TCVV(IDTOX) / AD(I,J,L) * 1d9
          ENDDO
          !WRITE(407,117) ( GC_O3(L),      L=LLPAR,1,-1 )
117       FORMAT(F18.6,1X)
       ENDDO
    ENDDO

  END SUBROUTINE CALC_GC_O3

!---------------------------------------------------------------------------------------------------------------------------

END MODULE OSIRIS_OBS_MOD
