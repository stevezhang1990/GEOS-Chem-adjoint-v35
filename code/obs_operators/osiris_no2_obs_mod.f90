MODULE OSIRIS_NO2_OBS_MOD

  !
  ! Module OSIRIS_NO2_OBS contains all subroutines and variables needed to assimilate OSIRIS NO2 tropospheric column data
  !
  ! Module Routines:
  !
  ! (1) CALC_OSIRIS_NO2_FORCE      : calculates adjoint forcing and cost function contribution for OSIRIS tropospheric NO2 columns
  ! (2) TAI2UTC                 : converts TAI93 (seconds since 1.1.1993) to UTC
  ! (3) MAKE_OSIRIS_BIAS_FILE_HDF5 : writes OSIRIS satellite diagnostics in satellite diagnostic HDF5 file
  !

  IMPLICIT NONE

#include "CMN_SIZE"

  PRIVATE

  PUBLIC CALC_OSIRIS_NO2_FORCE

  ! Module variables

  ! Arrays for diagnostic output

  TYPE FLEX_REAL                                ! Type to store information for "flexible" arrays. Think of this as a cruddy
     INTEGER :: CURRENT_N, MAX_N                ! implementation of some of the features of the C++ std::vector<T> container
     REAL*8,ALLOCATABLE :: DATA(:)              ! This only works in Fortran 2003, would have to use a pointer in Fortran 95
  ENDTYPE FLEX_REAL

  ! arrays to store diagnostic information
  REAL*4:: OSIRIS_NO2_MEAN(IIPAR,JJPAR) = 0d0      ! Mean OSIRIS columns
  REAL*4:: OSIRIS_GEOS_NO2_MEAN(IIPAR,JJPAR) = 0d0 ! Mean GEOS-Chem columns
  REAL*4:: OSIRIS_NO2_ERR_MEAN(IIPAR,JJPAR) = 0d0  ! Mean OSIRIS observation errors
  REAL*4:: OSIRIS_BIAS(IIPAR,JJPAR)=0d0            ! Model biases
  REAL*4:: OSIRIS_VAR(IIPAR,JJPAR)=0d0             ! Model variances
  REAL*4:: OSIRIS_DELTA=0d0                        ! temporary storage variable
  REAL*4:: OSIRIS_BIAS_COUNT(IIPAR,JJPAR) = 0d0    ! counter for number of observations in grid box
  REAL*4:: OSIRIS_CHISQUARED(IIPAR,JJPAR) = 0d0   ! Chi-squared values
  LOGICAL :: FIRST = .TRUE.
  TYPE(FLEX_REAL) :: FLEX_LON, FLEX_LAT, FLEX_TIME, FLEX_OSIRIS_NO2, FLEX_GC_NO2 ! flex arrays to store satellite diagnostics sequentially

CONTAINS

  !-----------------------------------------------------------------------------!

  SUBROUTINE CALC_OSIRIS_NO2_FORCE


    !!
    !! Subroutine CALC_OSIRIS_NO2_FORCE computes the NO2 adjoint forcing and cost function contribution from OSIRIS column data
    !!
    !! References:
    !!
    !! Bucsela2013:
    !! "A new stratospheric and tropospheric NO2 retrieval algorithm for nadir-viewing satellite instruments: applications to OSIRIS"
    !! E.J. Bucsela et.al
    !! Atmos. Meas. Tech., 6, 2607-2626, 2013
    !! www.atmos-meas-tech.net/6/2607/2013/
    !! doi:10.5194/amt-6-2607-2013
    !!
    !! Chan83
    !! "Algorithms for Computing the Sample Variance: Analysis and Recommendations"
    !! Tony F. Chan, Gene H. Golub, Randall J. LeVeque
    !! The American Statistician
    !! Vol. 37, No. 3 (Aug. 1983), pp. 242-247
    !!

    USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
    USE CHECKPT_MOD,        ONLY : CHK_STT
    USE COMODE_MOD,         ONLY : JLOP
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
    USE GRID_MOD,           ONLY : GET_IJ
    USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
    USE TIME_MOD,           ONLY : GET_HOUR, GET_DAY, GET_YEAR, GET_MONTH
    USE DAO_MOD,            ONLY : BXHEIGHT, AD, AIRDEN
    USE FILE_MOD,           ONLY : IOERROR
    USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
    USE ADJ_ARRAYS_MOD,     ONLY : ID2C
    USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
    USE TRACER_MOD,         ONLY : TCVV, XNUMOLAIR
    USE TRACERID_MOD,       ONLY : IDTNOX, IDNO2
    USE ADJ_ARRAYS_MOD,     ONLY : COST_FUNC
    USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
    USE TRACER_MOD,         ONlY : XNUMOLAIR
    USE DAO_MOD,            ONLY : T, AIRDEN
    USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
    USE ERROR_MOD,          ONLY : IT_IS_NAN


    INTEGER :: I,J,L,K,KK,MM,DD,CC
    INTEGER :: I_OSIRIS, J_OSIRIS, K_OSIRIS, JLOOP
    INTEGER :: IIJJ(2)
    INTEGER :: DAY, YEAR, MONTH, MJD_CALC, MJD_INI, MJD_FIN
    INTEGER :: UNIT, MAX_OBS_PER_DAY
    CHARACTER(255) :: ORBIT_PATH, FILE_ORBIT, ALT_PATH, FILE_ALT 
    CHARACTER(2) :: I_CHAR, I_CHAR_1, I_CHAR_2
    INTEGER :: IO_ORBIT_STATUS
    INTEGER, PARAMETER :: OLMAX = 100
    CHARACTER(LEN=255) :: filename_orbit, dsetname

    REAL*8, ALLOCATABLE :: LON_ORBIT(:), LAT_ORBIT(:), LON_ORBIT2(:)
    REAL*8, ALLOCATABLE :: TIME_ORBIT(:)
    REAL*8, ALLOCATABLE :: NO2_STRAT(:,:), NO2_STRAT_STD(:), NO2_STRAT2(:)
    REAL*8, ALLOCATABLE :: ALT_OSIRIS_NO2(:)
    REAL*8              :: TIME_ORBIT_START, TIME_ORBIT_END

    ! variables for time unit conversion
    REAL*8 :: tai93
    INTEGER :: iy,im,id,ih,imin
    REAL*8 :: sec
    INTEGER :: GC_HOUR, MIN_HOUR, MAX_HOUR, MIN_DAY, MAX_DAY

    ! variables for observation operator and adjoint thereof

    REAL*8 :: NO2_STRAT_GC(LLPAR), NO2_STRAT_GC_STD(LLPAR)
    REAL*8 :: GC_NO2_NATIVE(LLPAR), NCP(LLPAR), GC_ALT_NO2(IIPAR, JJPAR, LLPAR+1)
    REAL*8 :: GC_NO2_COL
    REAL*8 :: GC_NO2(OLMAX), DIFF(OLMAX), DIFF_ADJ(OLMAX)
    REAL*8 :: OBS_ERROR(OLMAX)
    REAL*8  :: COUNT_GRID(IIPAR,JJPAR), ALT_SURF(IIPAR,JJPAR)
    REAL*8  :: COST_CONTRIB(IIPAR,JJPAR)
    REAL*8 :: ADJ_FORCING(LLPAR)
    REAL*8  :: OSIRIS_NO2_STD(27)

    ! arrays needed for superobservations 
    LOGICAL :: SUPER_OBS = .TRUE.                       ! do super observations?
    REAL*8 ::  SOBS_COUNT(IIPAR,JJPAR)                   ! super observation count
    REAL*8 ::  SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)         ! super observation adjoint forcing
    REAL*8 ::  SOBS_COST_CONTRIBUTION(IIPAR,JJPAR)       ! super observation cost function contribution
    REAL*8 ::  SOBS_GC(IIPAR,JJPAR)       
    REAL*8 ::  SOBS_OSIRIS(IIPAR,JJPAR)       
    REAL*8 ::  SOBS_BIAS(IIPAR,JJPAR)
    REAL*8 ::  SOBS_CHISQUARED(IIPAR,JJPAR)

    LOGICAL, SAVE               :: SECOND = .TRUE.
    INTEGER                     :: IU_FILE, IU_DATA, IOS
    CHARACTER(LEN=255)          :: FILENAME
    CHARACTER(LEN=255)          :: FILENAME_ALT

    IF ( SECOND ) THEN
       FILENAME = 'lat_orb_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 601,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'lon_orb_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 602,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'amf_gc_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 603,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'amf_obs_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 604,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'gc_press_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 605,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'diff_osirisno2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 612,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    ENDIF


    !=================================================================
    ! CALC_OSIRIS_NO2_FORCE begins here!
    !=================================================================

    IF(FIRST) THEN

       ! initialize flexible arrays
       CALL INIT_FLEX_REAL(FLEX_LON)
       CALL INIT_FLEX_REAL(FLEX_LAT)
       CALL INIT_FLEX_REAL(FLEX_TIME)
       CALL INIT_FLEX_REAL(FLEX_OSIRIS_NO2)
       CALL INIT_FLEX_REAL(FLEX_GC_NO2)

       FIRST = .FALSE.

    ENDIF

    ! initialize arrays
    COUNT_GRID = 0d0
    COST_CONTRIB = 0d0
    GC_NO2 = 0d0
    GC_NO2_COL = 0d0
    GC_ALT_NO2 = 0d0
    OSIRIS_NO2_STD(1:27) = (/6.0,5.14,4.9,3.91,3.23,3.48,3.64,4.15,4.12,3.85,3.75,3.92,4.67,4.81,4.95,5.18,5.54,6.75,4.08,4.08,3.19,3.19,3.22,2.4,1.58,1.22,10/)
    SOBS_COUNT = 0d0
    SOBS_ADJ_FORCE = 0d0
    SOBS_COST_CONTRIBUTION = 0d0
    SOBS_GC = 0d0
    KK = 0
    !DD = MAX_OBS_PER_DAY
    ! Loop through data to find observations
    !PRINT *, "ID2C(IDNO2)", ID2C(IDNO2)
    GC_HOUR = GET_HOUR()
    DAY = GET_DAY()
    MONTH = GET_MONTH()
    YEAR = GET_YEAR()
    !USE REAL OR FLOOR
    MJD_CALC = FLOOR(2-FLOOR(REAL(YEAR/100))+FLOOR(REAL(FLOOR(REAL(YEAR)/100))/4)+REAL(DAY)+365.25*REAL(YEAR+4716)+30.6001*REAL(MONTH+1)-1524.5-2400000.5)
    MJD_INI = MJD_CALC - 50337
    MJD_FIN = MJD_INI + 2

    ORBIT_PATH = '/users/jk/07/xzhang/OSIRIS_NO2/ORBIT_DATA/2009/05/'
    ALT_PATH = '/users/jk/07/xzhang/met_field/'
    FILE_ALT = '20000101.cn.4x5.dat'
    
    FILENAME_ALT = TRIM(ALT_PATH) // TRIM(FILE_ALT)
    
    WRITE(I_CHAR,'(I2.2)') DAY

    CALL SYSTEM("ls "//TRIM(ORBIT_PATH)//"OMPS_L2_200905"//I_CHAR//"* > osiris_file_list"//I_CHAR//".txt")

    OPEN(UNIT = 18, FILE = "/users/jk/07/xzhang/met_field/20000101.cn.4x5.dat", STATUS="old",ACTION="read")
    READ(18,*) ALT_SURF
    CLOSE(18)
    IU_DATA = 20
    IU_FILE = 13
    !PRINT *, "IU_FILE IS INITIALIZED"
    CLOSE(IU_FILE) ! ugly...

    OPEN(IU_FILE,FILE="osiris_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")

271 DO
       READ(IU_FILE,'(A)',IOSTAT=IO_ORBIT_STATUS) FILE_ORBIT
       IF(IO_ORBIT_STATUS < 0) EXIT
       WRITE(6,*) ' - READ_OSIRIS_NO2_FILE: reading: ', FILE_ORBIT
       IU_DATA = IU_DATA + KK
       KK = KK + 1
       OPEN(IU_DATA, FILE=FILE_ORBIT,IOSTAT=IOS)
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_DATA, 'osiris:1' )
       
       !========================
       ! Read in data blocks
       !========================
       ! Needs to be true until end of file
       
       IF ( IOS < 0 ) EXIT
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_DATA, 'osiris:2' )
       READ( IU_DATA, * ) MAX_OBS_PER_DAY
       !PRINT *, "MAX_OBS_PER_DAY", MAX_OBS_PER_DAY
       MM = 1
       DD = MAX_OBS_PER_DAY
       ! Read altitude, ozone and error profiles
       ALLOCATE(NO2_STRAT2(OLMAX*MAX_OBS_PER_DAY))
       READ( IU_DATA, * )(NO2_STRAT2(K), K=1,OLMAX*MAX_OBS_PER_DAY)
       ALLOCATE(NO2_STRAT(OLMAX,MAX_OBS_PER_DAY))
       NO2_STRAT = RESHAPE(NO2_STRAT2,(/OLMAX,MAX_OBS_PER_DAY/))
       !PRINT *, "NO2_STRAT", NO2_STRAT

       ALLOCATE(ALT_OSIRIS_NO2(OLMAX))
       READ( IU_DATA, * )(ALT_OSIRIS_NO2(K), K=1,OLMAX)
       !PRINT *, "ALT_OSI_NO2", ALT_OSIRIS_NO2

       ALLOCATE(LAT_ORBIT(MAX_OBS_PER_DAY))       
       READ( IU_DATA, * )(LAT_ORBIT(K), K=1,MAX_OBS_PER_DAY)

       ALLOCATE(LON_ORBIT2(MAX_OBS_PER_DAY))
       READ( IU_DATA, * )(LON_ORBIT2(K), K=1,MAX_OBS_PER_DAY)
       ALLOCATE(LON_ORBIT(MAX_OBS_PER_DAY))
       LON_ORBIT = MOD(REAL(LON_ORBIT2+180d0),360d0) -180d0

       ALLOCATE(TIME_ORBIT(MAX_OBS_PER_DAY))
       READ( IU_DATA, * )(TIME_ORBIT(K), K=1,MAX_OBS_PER_DAY)
       !PRINT *, "TIME_ORBIT", TIME_ORBIT
       CLOSE( IU_DATA )

289    IF ( IT_IS_NAN(TIME_ORBIT(MM))) THEN
          MM = MM + 1
          GO TO 289
       ELSE
          !PRINT *, "MM", MM
          TIME_ORBIT_START = TIME_ORBIT(MM)
          !PRINT *, "TIME_ORBIT_START", TIME_ORBIT_START
       ENDIF
       
296    IF ( IT_IS_NAN(TIME_ORBIT(DD)) ) THEN
          DD = DD - 1
          GO TO 296
       ELSE
          !PRINT *, "DD", DD
          !PRINT *, "TIME_ORBIT_END", TIME_ORBIT_END
          TIME_ORBIT_END = TIME_ORBIT(DD)
       ENDIF
       
       ! check if current hour is in dataset
       CALL MJD2UTC(TIME_ORBIT_START,IY,IM,ID,IH,IMIN,SEC)
       MIN_HOUR = IH
       MIN_DAY = ID
       !PRINT *, "TIME_START", IY, IM, ID, IH, IMIN, SEC
       CALL MJD2UTC(TIME_ORBIT_END,IY,IM,ID,IH,IMIN,SEC)
       MAX_HOUR = IH
       MAX_DAY = ID
       !PRINT *, "TIME_FINISH", IY, IM, ID, IH, IMIN, SEC
       ! go to next dataset if current hour is not contained in dataset
       IF ( (GC_HOUR<MIN_HOUR .OR. GC_HOUR>MAX_HOUR) .OR. &
            (DAY<MIN_DAY .OR. DAY>MAX_DAY) ) THEN
          DEALLOCATE(TIME_ORBIT)          
          DEALLOCATE(NO2_STRAT)
          DEALLOCATE(NO2_STRAT2)
          DEALLOCATE(ALT_OSIRIS_NO2)
          DEALLOCATE(LAT_ORBIT)
          DEALLOCATE(LON_ORBIT)
          DEALLOCATE(LON_ORBIT2)
          !ALLOCATE(NO2_STRAT_STD(OLMAX))
          !GO TO 263
          CYCLE
       ENDIF
       ALLOCATE(NO2_STRAT_STD(OLMAX))
       !PRINT *, "data_dims_orbit", (/DATA_MAX_OBS_PER_DAY,0/)
       !! close file
       DO J_OSIRIS = 1, OLMAX
          IF ( IT_IS_NAN(ALT_OSIRIS_NO2(J_OSIRIS)) ) THEN
             NO2_STRAT_STD(J_OSIRIS) = 0d0
          ELSEIF ( ALT_OSIRIS_NO2(J_OSIRIS) < 12 ) THEN
             NO2_STRAT_STD(J_OSIRIS) = OSIRIS_NO2_STD(1) * 1d8
          ELSEIF (ALT_OSIRIS_NO2(J_OSIRIS) > 36) THEN
             NO2_STRAT_STD(J_OSIRIS) = OSIRIS_NO2_STD(27) * 1d8
          ELSE
             NO2_STRAT_STD(J_OSIRIS) = OSIRIS_NO2_STD(J_OSIRIS-11) * 1d8
          ENDIF
       ENDDO
       
       !! loop over data
       
       DO I_OSIRIS= MM,DD
          
          IF(TIME_ORBIT(I_OSIRIS)>0d0) THEN ! very basic quality check, most likely not needed anymore
             ! Convert TAI93 to UTC
             CALL MJD2UTC(TIME_ORBIT(I_OSIRIS),IY,IM,ID,IH,IMIN,SEC)
             ! A number of conditions have to be met for OSIRIS NO2 data to actually be assimilated

             IF ( ( GC_HOUR .EQ. ih ) .AND. &
                  (REAL(LAT_ORBIT(I_OSIRIS),4) < 60d0) .AND. &
                  (REAL(LAT_ORBIT(I_OSIRIS),4) > -60d0) .AND. &
                  ( DAY .EQ. id ) ) THEN
                
                ! Get model grid coordinate indices that correspond to the observation
                IIJJ = GET_IJ(REAL(LON_ORBIT(I_OSIRIS),4), REAL(LAT_ORBIT(I_OSIRIS),4))
                
                I = IIJJ(1)
                J = IIJJ(2)
                ! initialize variables & arrays
                
                GC_NO2_NATIVE = 0d0
                GC_NO2_COL = 0d0
                ! Get GEOS-CHEM NO2 values [#/cm3]
                GC_ALT_NO2(I,J,1) = ALT_SURF(I,J)*1d-3
                DO L = 1, LLPAR
                   GC_ALT_NO2(I,J,L+1) = (SUM(BXHEIGHT(I,J,1:L)) + ALT_SURF(I,J))*1d-3
                   IF (ITS_IN_THE_TROP(I,J,L)) THEN
                      JLOOP=JLOP(I,J,L)
                      GC_NO2_NATIVE(L) = CSPEC_AFTER_CHEM(JLOOP,ID2C(IDNO2))
                   ENDIF
                ENDDO
                CALL BIN_OSIRIS_NO2(GC_ALT_NO2(I,J,:), ALT_OSIRIS_NO2, GC_NO2_NATIVE(:), GC_NO2, OLMAX, 1)
                DIFF_ADJ(:) = 0
                DO J_OSIRIS = 1, OLMAX
                   L = J_OSIRIS
                   OBS_ERROR(L) = 2*NO2_STRAT_STD(J_OSIRIS)
                   !PRINT *, "GC_NO2", GC_NO2(L)
                   IF( ( ALT_OSIRIS_NO2(J_OSIRIS) < 37d0 ) .AND. &
                       ( ALT_OSIRIS_NO2(J_OSIRIS) > 6d0 ) .AND. &
                       ( NO2_STRAT(J_OSIRIS,I_OSIRIS ) > 0d0 ) .AND. &
                       ( GC_NO2(L) > 0d0                   )  ) THEN
                      DIFF(L) = GC_NO2(L) - NO2_STRAT(J_OSIRIS,I_OSIRIS)
                      !PRINT *, "GC_NO2", GC_NO2(L)
                      !PRINT *, "NO2_STRAT", NO2_STRAT(J_OSIRIS,I_OSIRIS)
                   ELSE
                      DIFF(L) = 0d0
                   ENDIF
                   !IF (GC_NO2(L) < 50*NO2_STRAT(J_OSIRIS,I_OSIRIS)) THEN
                   IF (SUPER_OBS) THEN
                      SOBS_COST_CONTRIBUTION(I,J) = SOBS_COST_CONTRIBUTION(I,J) + 0.5 * (DIFF(L)/OBS_ERROR(L)) ** 2
                               !OBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + DIFF(L)/(OBS_ERROR(L)**2) * TCVV(IDTNOX) * 1d-6 * XNUMOLAIR * AIRDEN(L,I,J)/AD(I,J,L)*BXHEIGHT(I,J,L)*100d0
                         !SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + DIFF(L)/(OBS_ERROR(L)**2)
                      DIFF_ADJ(L) = DIFF(L)/OBS_ERROR(L)**2
                   ELSE
                         !CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) + DIFF(L)/(OBS_ERROR(L)**2)
                               !TT_ADJ(I,J,L,IDTNOX) = STT_ADJ(I,J,L,IDTNOX) + DIFF(L)/(OBS_ERROR(L)**2) * TCVV(IDTNOX) * 1d-6 * XNUMOLAIR * AIRDEN(L,I,J)/AD(I,J,L)*BXHEIGHT(I,J,L)*100d0
                      COST_FUNC = COST_FUNC + 0.5 * (DIFF(L)/OBS_ERROR(L))**2
                   ENDIF
                   !ENDIF
                ENDDO
                !PRINT *, "DIFF_ADJ", DIFF_ADJ(:)
                CALL BIN_OSIRIS_NO2(GC_ALT_NO2(I,J,:), ALT_OSIRIS_NO2, ADJ_FORCING, DIFF_ADJ, OLMAX, -1)
                !PRINT *, "ADJ_FORCING", ADJ_FORCING(:)
                DO L = 1, LLPAR
                   IF (SUPER_OBS) THEN
                      SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + ADJ_FORCING(L)
                   ELSE
                      CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) + ADJ_FORCING(L)
                   ENDIF
                ENDDO

                WRITE(601,110) (LAT_ORBIT(I_OSIRIS))
                WRITE(602,110) (LON_ORBIT(I_OSIRIS))
                WRITE(612,110) (DIFF(L), L = LLPAR-1,1,-1)
110             FORMAT(F18.6,1X)
                !PRINT *, "SUPER_OBS_FORCE", SOBS_ADJ_FORCE(I,J,:)
                ! update cost function
                IF (SUPER_OBS) THEN
                   
                   !SOBS_COST_CONTRIBUTION(I,J) = SOBS_COST_CONTRIBUTION(I,J) + 0.5 * (DIFF/OBS_ERROR)**2
                   SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       !PRINT *, "CHK_STT", CHK_STT(:,:,:,IDNO2)
       !PRINT *, "SUPER_OBS_FORCE", SOBS_ADJ_FORCE
       !PRINT *, "SUPER_COST_CONTRIBUTION", SOBS_COST_CONTRIBUTION
       
       IF (SUPER_OBS) THEN
             
          DO J = 1,JJPAR
             DO I = 1, IIPAR
                   
                IF(SOBS_COUNT(I,J) > 0d0) THEN
                   DO L = 1,LLPAR
                      IF( (GC_ALT_NO2(I,J,L) < 37d0 ) .AND. ( GC_ALT_NO2(I,J,L) > 6d0 ) ) THEN
                         !JLOOP = JLOP(I,J,L)
                         !ELSE
                         !PRINT *, "STT_ADJ BEF",  STT_ADJ(I,J,L,IDNO2)
                         !TT_ADJ(I,J,L,IDTNOX) = STT_ADJ(I,J,L,IDTNOX) + SOBS_ADJ_FORCE(I,J,L)/SOBS_COUNT(I,J)
                         IF(ITS_IN_THE_TROP(I,J,L)) THEN
                            JLOOP = JLOP(I,J,L)
                            CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) + SOBS_ADJ_FORCE(I,J,L)/SOBS_COUNT(I,J)
                            !RINT *, "STT_ADJ AFT", STT_ADJ(I,J,L,IDTNOX)
                            !IF (SOBS_ADJ_FORCE(I,J,L)/SOBS_COUNT(I,J) > 0) THEN
                            !PRINT *, "STT_ADJ_FORCING", SOBS_ADJ_FORCE(I,J,L)/SOBS_COUNT(I,J)
                            !ENDIF
                         ENDIF
                      ENDIF
                   ENDDO
                   
                   COST_FUNC = COST_FUNC + SOBS_COST_CONTRIBUTION(I,J)/SOBS_COUNT(I,J)
                   !WRITE(104,110) ( GC_STT_ADJ(L), L=LLPAR,1,-1 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       
       !PRINT *, "STT_ADJ AFTER MLS O3", SOBS_ADJ_FORCE
       PRINT *, "COST FUNCTION OF OSIRIS NO2", COST_FUNC
       
       ! deallocate OSIRIS arrays
          
483    IF(ALLOCATED(LON_ORBIT)) DEALLOCATE(LON_ORBIT)
       IF(ALLOCATED(LON_ORBIT2)) DEALLOCATE(LON_ORBIT2)
       IF(ALLOCATED(LAT_ORBIT)) DEALLOCATE(LAT_ORBIT)
       IF(ALLOCATED(TIME_ORBIT)) DEALLOCATE(TIME_ORBIT)
       IF(ALLOCATED(NO2_STRAT)) DEALLOCATE(NO2_STRAT)
       IF(ALLOCATED(NO2_STRAT2)) DEALLOCATE(NO2_STRAT2)
       IF(ALLOCATED(NO2_STRAT_STD)) DEALLOCATE(NO2_STRAT_STD)
       IF(ALLOCATED(ALT_OSIRIS_NO2)) DEALLOCATE(ALT_OSIRIS_NO2)
       !KK = KK + 1
    ENDDO ! loop over OSIRIS files
    
    CLOSE(IU_FILE)
    !IF ( CC < 2 ) THEN
   !    CC = CC  + 1
  !     GO TO 262
       !PRINT *, "CC is updated to", CC
 !   ENDIF
    !ENDDO
    
  END SUBROUTINE CALC_OSIRIS_NO2_FORCE

  !-----------------------------------------------------------------------------!

  SUBROUTINE MJD2UTC(mjd_input,iy,im,id,ih,imin,sec)

    !!
    !! SUBROUTINE MJD2UTC converts MJD time (seconds since 1.1.1993) to UTC
    !!
    !! adapted from
    !! http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    !!
    
    IMPLICIT NONE

    REAL*8,INTENT(IN) :: MJD_INPUT
    INTEGER,INTENT(OUT) :: IY,IM,ID,IH,IMIN
    REAL*8,INTENT(OUT) :: SEC
    REAL*8 :: JD,JDF,JDL,JDK,JDH,JDMIN,JDSEC
    INTEGER :: JDI,JDJ,JDN
    !PRINT *, "MJD_INPUT", MJD_INPUT
    JD = MJD_INPUT + 2400001
    JDF = JD - FLOOR(JD)
    JDL = FLOOR(JD) + 68569
    JDN = FLOOR(REAL(4*JDL)/146097)
    JDL = JDL - FLOOR((REAL(146097*JDN)+3)/4)
    JDI = FLOOR(4000*REAL(JDL+1)/1461001)
    JDL = JDL - FLOOR(REAL(1461*JDI)/4) + 31
    JDJ = FLOOR(REAL(80*JDL)/2447)
    JDK = JDL - FLOOR(REAL(2447*JDJ)/80)
    JDL = FLOOR(REAL(JDJ)/11)
    JDJ = JDJ +2 - REAL(12*JDL)
    JDI = 100*REAL(JDN-49)+JDI+JDL
    JDH = REAL(JDF)*24
    IY = JDI
    IM = JDJ
    ID = JDK
    IH = FLOOR(JDH)
    JDMIN = (JDH-FLOOR(JDH))*60
    IMIN = FLOOR(JDMIN)
    JDSEC = (JDMIN-FLOOR(JDMIN))*60
    SEC = FLOOR(JDSEC)

    RETURN

  END SUBROUTINE MJD2UTC
  !--------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  SUBROUTINE BIN_OSIRIS_NO2( GC_EDGE, OBS_MODEL, DATA_MODEL, DATA_OSI, LOSIRIS, FB )
    
    !******************************************************************************
    !Based on the code from Monika.  (zhe 1/19/11)
    !FB = 1 for forward
    !FB = -1 for adjoint
    !******************************************************************************
    
    INTEGER :: L, LL, FB
    INTEGER :: LOSIRIS, NB
    REAL*8  :: ALT_MODEL(LLPAR), GC_EDGE(LLPAR+1)
    REAL*8  :: DATA_MODEL(LLPAR), DATA_OSI(LOSIRIS), DATA_TEM
    REAL*8  :: OBS_MODEL(LOSIRIS)
    
    !=================================================================
    ! BIN_DATA_V4 begins here!
    !=================================================================
    IF (FB > 0) THEN
       
       DO L = 1, LOSIRIS
          DO LL = 1, LLPAR
             IF ( GC_EDGE(LL) >= OBS_MODEL(L) ) THEN
                DATA_OSI(L) = DATA_MODEL(LL)
                EXIT
             ENDIF
          ENDDO
       ENDDO
       !PRINT *, "DATA_MODEL", DATA_MODEL(:)
       !PRINT *, "GC_EDGE", GC_EDGE(:)
       !PRINT *, "OBS_MODEL", OBS_MODEL(:)
       DO L = 1, LOSIRIS-1
          NB = 0
          DATA_TEM = 0
          DO LL = 1, LLPAR
             IF ( ( GC_EDGE(LL) >= OBS_MODEL(L)) .and. ( GC_EDGE(LL) < OBS_MODEL(L+1)) ) THEN
                !PRINT *, "GC_EDGE", GC_EDGE(LL)
                !PRINT *, "OBS_MODEL", OBS_MODEL(L)
                DATA_TEM = DATA_TEM + DATA_MODEL(LL)
                NB = NB + 1
             ENDIF
          ENDDO
          IF (NB > 0) DATA_OSI(L) = DATA_TEM / NB
       ENDDO
       
    ELSE
       DATA_MODEL(:) = 0
       DO L = 1, LOSIRIS-1
          DO LL = 1, LLPAR
             IF ( ( GC_EDGE(LL) >= OBS_MODEL(L)) .and. ( GC_EDGE(LL) < OBS_MODEL(L+1)) ) THEN
                DATA_MODEL(LL) = DATA_OSI(L)
                !PRINT *, "DATA_MODEL", DATA_MODEL(LL)
             ENDIF
          ENDDO
       ENDDO
       
    ENDIF
    
    
    ! Return to calling program
  END SUBROUTINE BIN_OSIRIS_NO2
!-------------------------------------------------------------------------------------------------
  !mkeller: helper routines for managing flexible arrays
  !         reinventing the wheel here, but hey...

  SUBROUTINE INIT_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL):: INPUT
    INPUT%CURRENT_N = 0
    INPUT%MAX_N = 1000
    IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA) ! safety first
    ALLOCATE(INPUT%DATA(INPUT%MAX_N))

  END SUBROUTINE INIT_FLEX_REAL

  SUBROUTINE GROW_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL) :: INPUT
    REAL*8, ALLOCATABLE :: TEMP_ARRAY(:)
    ALLOCATE(TEMP_ARRAY(INPUT%MAX_N * 2))
    TEMP_ARRAY(1:INPUT%MAX_N) = INPUT%DATA
    DEALLOCATE(INPUT%DATA)
    ALLOCATE(INPUT%DATA(INPUT%MAX_N * 2))
    INPUT%DATA = TEMP_ARRAY
    DEALLOCATE(TEMP_ARRAY)
    INPUT%MAX_N = INPUT%MAX_N * 2

  END SUBROUTINE GROW_FLEX_REAL

  SUBROUTINE PUSH_FLEX_REAL(INPUT, NEW_VAL)

    TYPE(FLEX_REAL) :: INPUT
    REAL*8 :: NEW_VAL
    IF(INPUT%CURRENT_N == INPUT%MAX_N) THEN
       CALL GROW_FLEX_REAL(INPUT)
    ENDIF
    INPUT%CURRENT_N = INPUT%CURRENT_N + 1
    INPUT%DATA(INPUT%CURRENT_N) = NEW_VAL

  END SUBROUTINE PUSH_FLEX_REAL

  SUBROUTINE CLEAR_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL) :: INPUT
    IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA)

  END SUBROUTINE CLEAR_FLEX_REAL

END MODULE OSIRIS_NO2_OBS_MOD
