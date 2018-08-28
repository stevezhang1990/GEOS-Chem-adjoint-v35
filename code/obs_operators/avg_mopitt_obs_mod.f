      MODULE MOPITT_OBS_MOD
!*****************************************************************************
!  Module MOPITT_OBS_MOD contains all the subroutines for the using of MOPITT
!  observation (version 3 and version 4).(zhe 1/19/11)
!  Remove the support to MOPITT v3 and v4. Now support v5 and v6. (Zhe 1/20/14)
!  Module Routines:
!  ============================================================================
!  (1 ) READ_MOPITT_FILE       : Read MOPITT hdf file
!  (2 ) CALC_MOPITT_FORCE      : Calculates cost function and STT_ADJ increments
!  (3 ) CALC_AVGKER            : Construct the averging kernel matrix
!  (4 ) BIN_DATA               : Interpolation between different vertical resolutions
!  (5 ) INIT_DOMAIN            : Define the observation window
!  (6 ) CALC_OBS_HOUR          : Calculated hour of morning obs
!  (7 ) ITS_TIME_FOR_MOPITT_OBS: FUNCTION that checks time vs. OBS_HOUR array
!  (8 ) READ_MOP02             : Reads MOPITT data fields from the HDF-EOS file
!  (9)  INFO_MOP02             : Prints name, dims, type, etc. of MOPITT data fields
!  (10) CLEANUP_MOP02          : Deallocates all module arrays
!  =============================================================================

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "../adjoint/define_adj.h"

      PRIVATE

      PUBLIC OBS_HOUR_MOPITT
      PUBLIC COUNT_TOTAL
      PUBLIC ITS_TIME_FOR_MOPITT_OBS
      PUBLIC READ_MOPITT_FILE
      PUBLIC CALC_MOPITT_FORCE

      !=============================================================================
      ! MODULE VARIABLES
      !=============================================================================

      INTEGER   :: OBS_HOUR_MOPITT(IIPAR,JJPAR)
      INTEGER   :: DOMAIN_OBS(IIPAR,JJPAR)
      REAL*8    :: COUNT_TOTAL

      REAL*4    :: ERR_PERCENT(IIPAR,JJPAR)
      REAL*4,  ALLOCATABLE :: A(:,:)
      REAL*4,  ALLOCATABLE :: T(:)
      REAL*4,  ALLOCATABLE :: XA(:)
      REAL*8,  ALLOCATABLE :: AC(:)

      ! MOPITT dimension fields
      INTEGER   :: T_DIM, Z_DIM
      REAL*4, ALLOCATABLE :: LATITUDE(:)
      REAL*4, ALLOCATABLE :: LONGITUDE(:)
      REAL*4, ALLOCATABLE :: PRESSURE(:)
      REAL*4, ALLOCATABLE :: SECONDS_IN_DAY(:)
      REAL*4, ALLOCATABLE :: MOPITT_GMT(:)
      REAL*8, ALLOCATABLE :: TAU(:)

      ! MOPITT data quantities
      REAL*4, ALLOCATABLE :: BOTTOM_PRESSURE(:)
      REAL*4, ALLOCATABLE :: CO_MIXING_RATIO(:,:,:)
      REAL*4, ALLOCATABLE :: CO_RET_BOT_MIXING_RATIO(:,:)
      REAL*4, ALLOCATABLE :: CO_TOTAL_COLUMN(:,:)
      REAL*4, ALLOCATABLE :: AVGKER(:,:,:)
      REAL*4, ALLOCATABLE :: RET_ERR_COV(:,:,:)
      INTEGER, ALLOCATABLE :: CLOUD_DES(:)
      INTEGER, ALLOCATABLE :: SURFACE_INDEX(:)

      ! MOPITT a priori
      INTEGER :: NLEV_AP
      REAL*4, ALLOCATABLE :: PLEV_AP(:)
      REAL*4, ALLOCATABLE :: CO_MR_AP(:,:,:)
      REAL*4, ALLOCATABLE :: CO_MR_AP_BOTTOM(:,:)

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_MOPITT_FILE( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_MOPITT_FILE reads the MOPITT hdf file.
!  (mak, 7/12/07, zhe 1/19/11)
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!******************************************************************************

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE TIME_MOD,  ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, AS
      CHARACTER(LEN=255)  :: DIR_MOPITT
      CHARACTER(LEN=255)  :: DIR_MONTH
      CHARACTER(LEN=255)  :: FILENAMEM
      CHARACTER(LEN=255)  :: FILENAME2
      LOGICAL             :: IT_EXISTS
      LOGICAL, SAVE       :: FIRST = .TRUE.

      !=================================================================
      ! READ_MOPITT_FILE begins here!
      !=================================================================
#if defined( MOPITT_V5_CO_OBS )
      DIR_MOPITT = '/nobackupp8/zjiang2/mopitt/'
      DIR_MONTH = 'v5/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V10.1.3.beta.hdf'
#endif
#if defined( MOPITT_V6_CO_OBS )
      DIR_MOPITT = '/nobackupp8/zjiang2/mopitt/'
      DIR_MONTH = 'v6/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V16.2.3.he5'
#endif
#if defined( MOPITT_V7_CO_OBS )
      DIR_MOPITT = '/users/jk/15/xzhang/MOPITT/'
      DIR_MONTH = '/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V17.9.3.he5'
#endif

      IF ( FIRST ) THEN
         ERR_PERCENT(:,:) = 0.0
         COUNT_TOTAL = 0
         FIRST            = .FALSE.
      ENDIF

      OBS_HOUR_MOPITT(:,:) = -99

      CALL EXPAND_DATE( FILENAMEM, YYYYMMDD, 0 )
      CALL EXPAND_DATE( DIR_MONTH, YYYYMMDD, 0 )

      FILENAME2 = TRIM( DIR_MOPITT ) // TRIM( DIR_MONTH ) // FILENAMEM
      PRINT*, '=== Reading ===:', TRIM( FILENAME2 )
	  
	     INQUIRE( FILE = FILENAME2, EXIST = IT_EXISTS )
	     IF (IT_EXISTS) THEN

         !CALL INFO_MOP02(FILENAME2)

         CALL READ_MOP02( FILENAME2 )

         CALL INIT_DOMAIN

         ! Calculate hour of day when obs should be compared to model
         CALL CALC_OBS_HOUR
		
	     ENDIF

      !CALL READ_ERROR_VARIANCE
      !We assume 20% uniform observation error
      ERR_PERCENT(:,:) = 0.2/LOG(10d0)

      END SUBROUTINE READ_MOPITT_FILE
!-------------------------------------------------------------------------------------------------

      SUBROUTINE CALC_MOPITT_FORCE

!******************************************************************************
! CALC_MOPITT_FORCE calculate cost function and STT_ADJ increments
! (zhe 1/19/11)
!******************************************************************************

      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_AP, GET_BP
      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_YEAR
      USE TIME_MOD,     ONLY : GET_HOUR
      USE CHECKPT_MOD,  ONLY : CHK_STT
      USE TRACER_MOD,   ONLY : TCVV
      USE TRACERID_MOD,   ONLY : IDTCO
      USE DAO_MOD,      ONLY : AD, IS_LAND
      USE ADJ_ARRAYS_MOD,   ONLY : SET_FORCING, SET_MOP_MOD_DIFF,
     &                  SET_MODEL_BIAS, SET_MODEL, SET_OBS,
     &                  COST_ARRAY, DAY_OF_SIM, IFD, JFD, LFD, NFD,
     &                  COST_FUNC, ADJ_FORCE, STT_ADJ
      USE LOGICAL_ADJ_MOD,  ONLY : LPRINTFD, LDCOSAT
      USE ERROR_MOD,        ONLY : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,         ONLY : GET_IJ, GET_XMID, GET_YMID
      USE DIRECTORY_ADJ_MOD, ONLY: DIAGADJ_DIR
      USE ADJ_ARRAYS_MOD,     ONLY: N_CALC, EXPAND_NAME
      USE TROPOPAUSE_MOD, ONLY: ITS_IN_THE_TROP

      LOGICAL, SAVE       :: SECOND = .TRUE.
      CHARACTER(LEN=255)  :: FILENAME

      ! Local Variables
      INTEGER ::  W, I, J, Z, ZZ, L,LL,IOS
      INTEGER ::  LON15, IIJJ(2)
      INTEGER ::  NLEV_RET

      REAL*4  ::  RETLEV(Z_DIM+1)
      REAL*8  ::  P_EDGE(Z_DIM+2), MODEL_COL, MOPITT_COL
      REAL*8  ::  UTC, TAU0
      REAL*8  ::  MODEL_P(LLPAR), MODEL_CO_MR(LLPAR)
      REAL*8  ::  COUNT_GRID(IIPAR,JJPAR)
      REAL*8  ::  GC_ADJ_COUNT(IIPAR,JJPAR,LLPAR)
      REAL*8  ::  COUNT(IIPAR,JJPAR)
      REAL*8  ::  MOP_COL_GRID(IIPAR,JJPAR)
      REAL*8  ::  MODEL_COL_GRID(IIPAR,JJPAR)
      REAL*8  ::  NEW_COST(IIPAR,JJPAR)
      REAL*8  ::  
      REAL*8  ::  ADJ_F(LLPAR)
      REAL*8  ::  SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)
      REAL*8  ::  SY
      REAL*8  ::  COST_CONTRIB(LLPAR)
      REAL*8  ::  MODEL_P_EDGE(LLPAR+1)
      !REAL*8  ::  DIFF_COST
      REAL*4  ::  MOP_CO_BIAS(IIPAR,JJPAR,11)
      REAL*4  ::  MOP_BIAS_COUNT(IIPAR,JJPAR,11)
      REAL*4  ::  MOP_CO_CHI_SQ(IIPAR,JJPAR,11)
      REAL*4  ::  MOP_CO_BIAS_SOBS(IIPAR,JJPAR,11)
      REAL*4  ::  MOP_CO_CHI_SQ_SOBS(IIPAR,JJPAR,11)      

      REAL*8, ALLOCATABLE :: GEOS_RAW(:)
      REAL*8, ALLOCATABLE :: MOP_CO(:)
      REAL*8, ALLOCATABLE :: DIFF_ADJ(:)
      REAL*8, ALLOCATABLE :: GEOS_CO(:)
      REAL*8, ALLOCATABLE :: DIFF_COST(:)


      !=================================================================
      ! CALC_MOPITT_FORCE begins here!
      !=================================================================

      TAU0 = GET_TAU0( GET_MONTH(), GET_DAY(), GET_YEAR() )

      COUNT_GRID(:,:)     = 0d0
      COUNT(:,:)          = 0d0
      SOBS_COUNT(:,:)     = 0d0
      MOP_COL_GRID(:,:)   = -999.0
      MODEL_COL_GRID(:,:) = -999.0
      ADJ_FORCE(:,:,:,:)  = 0d0
      NEW_COST(:,:)       = 0d0
      MOP_CO_BIAS(:,:,:)  = 0d0
      MOP_BIAS_COUNT(:,:,:) = 0d0
      MOP_CO_CHI_SQ(:,:,:) = 0d0
      MOP_CO_BIAS_SOBS(:,:,:) = 0d0
      MOP_CO_CHI_SQ_SOBS(:,:,:) = 0d0
      GC_ADJ_COUNT(:,:,:) = 0d0
      SOBS_ADJ_FORCE(:,:,:) = 0d0
      
      IF ( SECOND ) THEN
         FILENAME = 'co_bias_mop.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 201,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &         IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
         FILENAME = 'co_chi_square_mop.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 202,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &         IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
         FILENAME = 'lat_orb_mop.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 203,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &         IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
         FILENAME = 'lon_orb_mop.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 204,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN',
     &         IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      ENDIF

      !=================================================================
      ! Loop over MOPITT data
      !=================================================================
      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         LON15 = LONGITUDE(W) / 15.
         UTC   = TAU(W) - TAU0 + LON15
         IF (UTC < 0. )  UTC = UTC + 24
         IF (UTC > 24.)  UTC = UTC - 24


         !Only consider day time MOPITT measurements
         ! am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)
#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box
         IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
         I = IIJJ(1)
         J = IIJJ(2)

         !=================================================================
         ! Data selection
         !=================================================================
         IF( GET_HOUR() == OBS_HOUR_MOPITT(I,J) .and.
     &      CLOUD_DES(W) == 2.0 .and.
     &      CO_TOTAL_COLUMN(1,W) > 5E17 .and.
     &      DOMAIN_OBS(I,J) == 1 ) THEN
     
!            IF ( (IS_LAND(I,J) .AND. 
!     &         LATITUDE(W) .GE. -52 .AND. LATITUDE(W) .LE. 52 ) .OR.   !52S-52N
!     &        (LATITUDE(W) .GE. -40 .AND. LATITUDE(W) .LE. 40) ) THEN  !40S-40N

            RETLEV(:) = -999.0
            MODEL_COL = 0D0
            MOPITT_COL = 0D0
            
            ! Create pressure profile
            RETLEV(1) = BOTTOM_PRESSURE(W)

            ZZ = 0
            ! Loop over Mopitt levels
            DO Z = 1, Z_DIM
            ! Always start from the bottom pressure,
            ! even if it means skipping a MOPITT pressure level
               IF ( PRESSURE(Z) >= RETLEV(1) ) THEN
                  ZZ = ZZ + 1
                  CYCLE
               ENDIF
               ! Save into profile
               RETLEV(Z+1-ZZ) = PRESSURE(Z)
            ENDDO
            NLEV_RET = Z_DIM+1 - ZZ

            DO L = 1, NLEV_RET
               P_EDGE(L) = RETLEV(L)
            ENDDO
            P_EDGE(NLEV_RET+1) = 36

            ALLOCATE( XA( NLEV_RET ) )
            ALLOCATE( T( NLEV_RET ) )
            ALLOCATE( A( NLEV_RET,NLEV_RET ) )
            ALLOCATE( AC( NLEV_RET ) )
            ALLOCATE( MOP_CO( NLEV_RET ) )
            ALLOCATE( GEOS_RAW( NLEV_RET ) )
            ALLOCATE( DIFF_ADJ( NLEV_RET ) )
            ALLOCATE( GEOS_CO( NLEV_RET ) )
            ALLOCATE( DIFF_COST( NLEV_RET ) )

            ! MOPITT CO vertical profile
            MOP_CO(1) = CO_RET_BOT_MIXING_RATIO(1,W)
            MOP_CO(2:NLEV_RET) = CO_MIXING_RATIO(1,11-NLEV_RET:9,W)
            MOP_CO = MOP_CO * 1E-9

            ! COMPUTE AVERAGING KERNEL
            CALL CALC_AVGKER(NLEV_RET, W, RETLEV, MOP_CO)

            !USE MOPITT SURFACE PRESSURE
            !DO L=1, LLPAR + 1
            !   MODEL_P_EDGE(L) = GET_AP(L) + GET_BP(L) * RETLEV(1)
            !ENDDO

            DO L = 1, LLPAR
               !MOPITT PRESSURE LEVEL
               !MODEL_P(L) = (MODEL_P_EDGE(L) + MODEL_P_EDGE(L+1)) / 2

               ! Get GC pressure levels (mbar)
               MODEL_P(L) = GET_PCENTER(I,J,L)

               ! Obtain archieved forward model results
               ! kg -> v/v
               MODEL_CO_MR(L) = CHK_STT(I,J,L,IDTCO) *
     &                            TCVV(IDTCO) / AD(I,J,L)
            ENDDO

            ! Interplote the model to MOPITT vertical grids
            CALL BIN_DATA(MODEL_P, P_EDGE, MODEL_CO_MR(:),
     &            GEOS_RAW, NLEV_RET, 1)

            !=================================================================
            ! Apply MOPITT observation operator
            !=================================================================

            ! Total Column: C = T * XA  + AC * ( Xm - XA )
            ! Stratosphere Levels are removed
            !DO L = 1, NLEV_RET
            DO L = 1, NLEV_RET - 2
               MODEL_COL = MODEL_COL
     &                   + T(L) * XA(L)
     &                   + AC(L) * (LOG10(GEOS_RAW(L))
     &                   - LOG10(XA(L)))
               !MOPITT_COL = MOPITT_COL + T(L) * MOP_CO(L)
            ENDDO
            
            MOPITT_COL = CO_TOTAL_COLUMN(1,W)

            GEOS_CO(:) = 0d0
            ! Smoothed Profile: X_hat = XA  + A * ( Xm - XA )
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  GEOS_CO(L) = GEOS_CO(L)
     &                       + A(L,LL)
     &                       * (LOG10( GEOS_RAW(LL) ) - LOG10( XA(LL) ))
               ENDDO
               GEOS_CO(L) = LOG10( XA(L) ) + GEOS_CO(L)
            ENDDO

            !=================================================================
            ! COST FUNCTION
            !=================================================================
            !SY = ( ERR_PERCENT(I,J) * MOPITT_COL )**2
            !DIFF_COST     = MODEL_COL - MOPITT_COL
            !NEW_COST(I,J) = NEW_COST(I,J) + 0.5 * (DIFF_COST ** 2) / SY
            !COUNT(I,J) = COUNT(I,J) +1
            DIFF_COST(:) = 0D0
            COST_CONTRIB(:) = 0D0
            SY = ERR_PERCENT(I,J) **2
            DO L = 1, NLEV_RET - 2
               DIFF_COST(L) = GEOS_CO(L) - LOG10( MOP_CO(L) )
               COST_CONTRIB(L) = 0.5d0 * ( DIFF_COST(L)**2 ) / SY
               COUNT(I,J) = COUNT(I,J) + 1
               MOP_CO_BIAS(I,J,L) = MOP_CO_BIAS(I,J,L) +
     &             10**(GEOS_CO(L)) - MOP_CO(L)
               MOP_CO_CHI_SQ(I,J,L) = MOP_CO_CHI_SQ(I,J,L) +
     &         (DIFF_COST(L))**2/ SY
               MOP_BIAS_COUNT(I,J,L) = MOP_BIAS_COUNT(I,J,L) + 1d0
            ENDDO
            !=================================================================
            ! adjoint operator
            !=================================================================
            DIFF_ADJ(:) = 0D0
            !DO L = 1, NLEV_RET
               !DIFF_ADJ(L) = DIFF_COST *  AC(L) / SY
               !DIFF_ADJ(L) = DIFF_ADJ(L) / (GEOS_RAW(L) * LOG(10.0))
            !ENDDO
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  DIFF_ADJ(L) = DIFF_ADJ(L)
     &                        + A(LL,L) * DIFF_COST(LL) / SY
               ENDDO
               ! fwd code: LOG(GEOS_RAW) - LOG(XA)
               ! mkeller: this is just plain wrong!
               !          the forward code is LOG10(GEOS_RAW) - LOG10(XA)
               !          a factor of 1/LOG(10) is missing
               !DIFF_ADJ(L) = DIFF_ADJ(L) / GEOS_RAW(L)
               DIFF_ADJ(L) = DIFF_ADJ(L) / (GEOS_RAW(L) * LOG(10d0))
            ENDDO


            CALL BIN_DATA( MODEL_P,  P_EDGE, ADJ_F,
     &                     DIFF_ADJ, NLEV_RET, -1   )

            ! adjoint FORCE
            DO L = 1, LLPAR
               IF (ITS_IN_THE_TROP(I,J,L)) THEN
               !v/v->kg
                  ADJ_FORCE(I,J,L,IDTCO) = ADJ_FORCE(I,J,L,IDTCO)
     &                 + ADJ_F(L) * TCVV(IDTCO)/ AD(I,J,L)
                  GC_ADJ_COUNT(I,J,L) = GC_ADJ_COUNT(I,J,L) + 1
                  NEW_COST(I,J) = NEW_COST(I,J) + COST_CONTRIB(L)
                  SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1
               ENDIF
            ENDDO

            COUNT_GRID(I,J)     = COUNT_GRID(I,J) + 1.d0
            MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J) + MOPITT_COL
            MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J) + MODEL_COL

            IF ( ALLOCATED( GEOS_RAW ) ) DEALLOCATE( GEOS_RAW )
            IF ( ALLOCATED( MOP_CO   ) ) DEALLOCATE( MOP_CO   )
            IF ( ALLOCATED( DIFF_ADJ ) ) DEALLOCATE( DIFF_ADJ )
            IF ( ALLOCATED( A        ) ) DEALLOCATE( A        )
            IF ( ALLOCATED( AC       ) ) DEALLOCATE( AC       )
            IF ( ALLOCATED( T        ) ) DEALLOCATE( T        )
            IF ( ALLOCATED( XA       ) ) DEALLOCATE( XA       )
            IF ( ALLOCATED( GEOS_CO  ) ) DEALLOCATE( GEOS_CO  )
            IF ( ALLOCATED( DIFF_COST) ) DEALLOCATE( DIFF_COST)

         !ENDIF !IS_LAND

         ENDIF !OBS_HOUR

         ENDIF !local time

      ENDDO  !loop over MOPITT data

      !=================================================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=================================================================
      DO I = 1, IIPAR
         DO J = 1, JJPAR

            IF ( COUNT_GRID(I,J) > 0d0 ) THEN

               !The mean value in the grid
               MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               DO L = 1, LLPAR
                  IF (ITS_IN_THE_TROP(I,J,L)) THEN
                     COST_FUNC = COST_FUNC + 
     &                    NEW_COST(I,J,L)/GC_ADJ_COUNT(I,J,L)
                     STT_ADJ(I,J,L,IDTCO) = STT_ADJ(I,J,L,IDTCO) +
     &                    ADJ_FORCE(I,J,L,IDTCO)/GC_ADJ_COUNT(I,J,L)
                  ENDIF
               ENDDO
               ! Diagnostic stuff: FORCING, MOP_MOD_DIFF, MODEL_BIAS
               IF( LDCOSAT )THEN

                  CALL SET_FORCING( I, J, DAY_OF_SIM,
     &                              ADJ_FORCE(I,J,1,IDTCO) )
                  CALL SET_MOP_MOD_DIFF( I, J, DAY_OF_SIM,
     &               MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) )

                  CALL SET_MODEL_BIAS( I, J, DAY_OF_SIM, 1,
     &              ( MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) ) /
     &                                 MOP_COL_GRID(I,J)           )
                  CALL SET_MODEL     ( I, J, DAY_OF_SIM, 1,
     &                                 MODEL_COL_GRID(I,J)         )
                  CALL SET_OBS       ( I, J, DAY_OF_SIM, 1,
     &                                 MOP_COL_GRID(I,J)           )

                  COST_ARRAY(I,J,DAY_OF_SIM) =
     &               COST_ARRAY(I,J,DAY_OF_SIM) + NEW_COST(I,J)

               ENDIF

               IF ( IT_IS_NAN( NEW_COST(I,J) ) ) THEN
                  PRINT*, 'I=', I, 'J=', J
                  CALL ERROR_STOP( 'NEW_COST is NaN',
     &                             'CALC_MOPITT_FORCE')
               ENDIF

            ENDIF !COUNT_GRID
            !DO L=1,NLEV_RET
            IF (MOP_BIAS_COUNT(I,J,6) > 0d0) THEN
               MOP_CO_BIAS_SOBS(I,J,6) = 
     &              MOP_CO_BIAS(I,J,6)/MOP_BIAS_COUNT(I,J,6)
               !PRINT *, "MOP_CO_BIAS", MOP_CO_BIAS_SOBS(I,J,6)
               MOP_CO_CHI_SQ_SOBS(I,J,6) = 
     &              MOP_CO_CHI_SQ(I,J,6)/MOP_BIAS_COUNT(I,J,6)
               WRITE(201,110) (1e12*MOP_CO_BIAS_SOBS(I,J,6))
               WRITE(202,110) (MOP_CO_CHI_SQ_SOBS(I,J,6))
               WRITE(203,110) (GET_XMID(I))
               WRITE(204,110) (GET_YMID(J))
            ENDIF
            !ENDDO
         ENDDO
      ENDDO
 110  FORMAT(F18.6,1X)
      IF (LPRINTFD)  THEN
         PRINT*, 'IFD, JFD= ', IFD, JFD
         PRINT*, 'MODEL_STT:', MODEL_COL_GRID(IFD,JFD)
         PRINT*, 'OBS_STT:', MOP_COL_GRID(IFD,JFD)
         PRINT*, 'NEW_COST', NEW_COST(IFD,JFD)
         PRINT*, 'ADJ_FORCE:', ADJ_FORCE(IFD,JFD,:,IDTCO)
         PRINT*, 'STT_ADJ:', STT_ADJ(IFD,JFD,:,IDTCO)
      ENDIF

      ! Update cost function
      !PRINT*, 'TOTAL NEW_COST = ', SUM(NEW_COST)
      !PRINT*, 'COST_FUNC BEFORE ADDING NEW_COST=', COST_FUNC
      !COST_FUNC   = COST_FUNC   + SUM ( NEW_COST )
      !COUNT_TOTAL = COUNT_TOTAL + SUM ( COUNT    )
      !PRINT*, 'Total observation number:', COUNT_TOTAL

      ! Return to calling program
      END SUBROUTINE CALC_MOPITT_FORCE
!--------------------------------------------------------------------------------------------

      SUBROUTINE CALC_AVGKER( NLEV_RET, W, RETLEV, MOP_CO )

!******************************************************************************
! SUBROUTINE CALC_AVGKER construct the averging kernel matrix
! (zhe 1/19/11)
!******************************************************************************

      INTEGER :: ILEV, JLEV, ILEV2, JLEV2, Z, W
      INTEGER :: NLEV_RET
      REAL*4  :: DELP(NLEV_RET)
      REAL*4  :: RETLEV(NLEV_RET)
      REAL*8  :: MOP_CO(NLEV_RET)
      REAL*8, PARAMETER :: log10e = LOG10(2.71828183)

      !=================================================================
      ! CALC_AVGKER begins here!
      !=================================================================

      A(:,:) = 0d0
      AC(:)  = 0d0

      XA(1) = CO_MR_AP_BOTTOM(1, W)
      XA(2:NLEV_RET) = CO_MR_AP(1,11-NLEV_RET:9,W)
      XA = XA * 1E-9

      !Remove bad levels from averging kernel matrix
      IF ( NLEV_RET < 10 ) THEN
         DO ILEV = 1, NLEV_RET
            ILEV2 = ILEV + ( 10 - NLEV_RET )
            DO JLEV =1, NLEV_RET
               JLEV2 = JLEV + ( 10 - NLEV_RET)
               A(ILEV,JLEV) =
     &              AVGKER(ILEV2,JLEV2,W)
            ENDDO
         ENDDO
      ELSE
         A(:,:) = AVGKER(:,:,W)
      ENDIF

      DELP(1) = RETLEV(1) - RETLEV(2)
      DELP(2:NLEV_RET-1) = 100D0
      DELP(NLEV_RET) = 74D0

      ! transfer function [v/v -> molec/cm2]
      T = 2.12E+22 * DELP

      ! Convert to column averaging kernel
      DO JLEV = 1, NLEV_RET
         DO ILEV = 1, NLEV_RET
            AC(JLEV) = AC(JLEV) + DELP(ILEV) * MOP_CO(ILEV)
     &               * A(ILEV,JLEV)
         ENDDO
         AC(JLEV) = (2.12E+22 / log10e ) * AC(JLEV)
      ENDDO


      END SUBROUTINE CALC_AVGKER
!------------------------------------------------------------------------------------

      SUBROUTINE BIN_DATA( P_MODEL,  P_EDGE, DATA_MODEL, DATA_MOP,
     &                        NLEV_RET, FB                            )

!******************************************************************************
!Based on the code from Monika.  (zhe 1/19/11)
!FB = 1 for forward
!FB = -1 for adjoint
!******************************************************************************

      INTEGER :: L, LL, FB
      INTEGER :: NLEV_RET, NB
      REAL*8  :: P_MODEL(LLPAR)
      REAL*8  :: DATA_MODEL(LLPAR), DATA_MOP(NLEV_RET), DATA_TEM
      REAL*8  :: P_EDGE(NLEV_RET+1)

      !=================================================================
      ! BIN_DATA begins here!
      !=================================================================

      IF (FB > 0) THEN
         
         DO L = 1, NLEV_RET
            DO LL = 1, LLPAR
               IF ( P_MODEL(LL) <= P_EDGE(L) ) THEN
                  DATA_MOP(L) = DATA_MODEL(LL)
                  EXIT
               ENDIF
            ENDDO
         ENDDO
                  
         DO L = 1, NLEV_RET
            NB = 0
            DATA_TEM = 0
            DO LL = 1, LLPAR
               IF ( ( P_MODEL(LL) <= P_EDGE(L)) .and.
     &              ( P_MODEL(LL) > P_EDGE(L+1)) ) THEN
                  DATA_TEM = DATA_TEM + DATA_MODEL(LL)
                  NB = NB + 1
               ENDIF
            ENDDO
            IF (NB > 0) DATA_MOP(L) = DATA_TEM / NB
         ENDDO

      ELSE

         DATA_MODEL(:) = 0.
         DO L = 1, LLPAR
            DO LL = 1, NLEV_RET
               IF ( ( P_MODEL(L) <= P_EDGE(LL)) .and.
     &              ( P_MODEL(L) > P_EDGE(LL+1)) ) THEN
                  DATA_MODEL(L) = DATA_MOP(LL)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      ! Return to calling program
      END SUBROUTINE BIN_DATA
!-----------------------------------------------------------------------------------

      SUBROUTINE INIT_DOMAIN

!******************************************************************************
!Define the observatio region
!******************************************************************************
#     include "CMN_SIZE"   ! Size parameters

      !local variables
      INTEGER :: I, J

      !=================================================================
      ! INIT_DOMAIN begins here!
      !=================================================================

      DOMAIN_OBS(:,:) = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR

#if   defined( GRID05x0666 )
!     The surrounding region is used as cushion
!     (zhe 11/28/10)
         IF ( J >= 8 .and. J <= JJPAR-7 .and.
     &        I >= 7 .and. I <= IIPAR-6
#elif defined( GRID2x25 )
         IF ( J >= 16 .and. J <= 76   !60S-60N
#elif defined( GRID4x5 )
         IF ( J >= 9 .and. J <= 39    !60S-60N
#endif
     &       ) DOMAIN_OBS(I,J) = 1d0

      ENDDO
      ENDDO

      PRINT*, sum(DOMAIN_obs), 'MAX observations today'

      END SUBROUTINE INIT_DOMAIN

!-----------------------------------------------------------------------------

      SUBROUTINE CALC_OBS_HOUR

!***************************************************************************
! Subroutine CALC_OBS_HOUR computes an array of hours for each day of obs.
! If there is an obs in a particular gridbox on that day, it assigns the
! hour (0..23). If there isn't, OBS_HOUR stays initialized to -1.
! (mak, 12/14/05)
!***************************************************************************

      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY,
     &                         GET_YEAR, GET_HOUR
      USE GRID_MOD,     ONLY : GET_IJ

#     include "CMN_SIZE"

      REAL*4   :: OBS_HOUR(IIPAR,JJPAR)
      REAL*8   :: TAU0, UTC
      INTEGER  :: W, I, J
      INTEGER  :: LON15, IIJJ(2)
      INTEGER  :: COUNT_GRID(IIPAR,JJPAR)

      !=================================================================
      ! CALC_OBS_HOUR begins here!
      !=================================================================

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = GET_TAU0(GET_MONTH(), GET_DAY(), GET_YEAR())

      OBS_HOUR_MOPITT(:,:) = -1
      OBS_HOUR(:,:)        = 0
      COUNT_GRID(:,:)      = 0

      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         !============================================================
         LON15 = LONGITUDE(W) / 15d0
         UTC   = TAU(W) - TAU0 + LON15
         IF ( UTC < 0d0  )  UTC = UTC + 24
         IF ( UTC > 24d0 )  UTC = UTC - 24

         !Only consider day time MOPITT measurements
         !am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)

#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box of current record
            IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
            I = IIJJ(1)
            J = IIJJ(2)

            ! If there's an obs, calculate the time
            IF ( CO_TOTAL_COLUMN(1,W) > 0d0 ) THEN

               COUNT_GRID(I,J) = COUNT_GRID(I,J) + 1d0
               !Add the time of obs, to be averaged and floored later
               OBS_HOUR(I,J) = OBS_HOUR(I,J) + MOPITT_GMT(W)

            ENDIF
         ENDIF
      ENDDO

      ! average obs_hour on the grid
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID(I,J) > 0d0 ) THEN

            OBS_HOUR_MOPITT(I,J) =
     &         FLOOR( OBS_HOUR(I,J) / COUNT_GRID(I,J) )

         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE CALC_OBS_HOUR

!----------------------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_MOPITT_OBS( ) RESULT( FLAG )

!******************************************************************************
!  Function ITS_TIME_FOR_MOPITT_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day) based on
!  the OBS_HOUR_MOPITT array which holds the hour of obs in each gridbox
!  (computed when file read in mop02_mod.f) (mak, 7/12/07)
!******************************************************************************

      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE

#     include "CMN_SIZE"  ! Size params

      ! Function value
      LOGICAL :: FLAG

      INTEGER :: I,J

      !=================================================================
      ! ITS_TIME_FOR_MOPITT_OBS begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.

      DO J = 1,JJPAR
      DO I = 1,IIPAR
         IF( GET_HOUR()   == OBS_HOUR_MOPITT(I,J) .and.
     &       GET_MINUTE() == 0                          ) THEN

               PRINT*, 'obs_hour was', get_hour(), 'in box', I, J
               FLAG = .TRUE.

               !GOTO 11
               RETURN

         ENDIF
      ENDDO
      ENDDO

      END FUNCTION ITS_TIME_FOR_MOPITT_OBS

!----------------------------------------------------------------------------

      SUBROUTINE READ_MOP02( FILENAME )

!******************************************************************************
!  Subroutine READ_MOP02 allocates all module arrays and reads data into
!  them from the HDF file. (bmy, 7/2/03, zhe 1/19/11)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
#if defined( MOPITT_V5_CO_OBS )
      USE HdfSdModule
      USE HdfVdModule
#endif
      USE BPCH2_MOD,  ONLY : GET_TAU0
      USE ERROR_MOD,  ONLY : ALLOC_ERR

      ! Local variables
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      INTEGER  :: as, i, year, month, day
      REAL*8   :: TAU0
      
#if defined( MOPITT_V5_CO_OBS )
      INTEGER  :: sId, vId, vSize, nDims, dims(4)
#endif
#if defined( MOPITT_V6_CO_OBS ) .or. defined( MOPITT_V7_CO_OBS )

      INTEGER  :: he5_swopen, he5_swattach, he5_swfldinfo, 
     &            he5_swrdfld, he5_swdetach, he5_swclose
     
      INTEGER    :: N, fId, swathid, rank
      INTEGER    :: ntype(4)
      INTEGER*8  :: dims8(4)
      INTEGER*8  :: START1(1), STRIDE1(1), EDGE1(1)
      INTEGER*8  :: START2(2), STRIDE2(2), EDGE2(2)
      INTEGER*8  :: START3(3), STRIDE3(3), EDGE3(3)
      INTEGER, PARAMETER  :: HE5F_ACC_RDONLY=101
      character*72 dimlist, maxdimlist      
      
#endif
                        
      !=================================================================
      ! Mop02Read begins here!
      !=================================================================

      ! Deallocate arrays
      CALL CLEANUP_MOP02

      ! Get date from filename (next to the '-' character)
      i    = INDEX( FILENAME, '-' )
      READ( FILENAME(i+1:i+4), '(i4)' ) year
      READ( FILENAME(i+5:i+6), '(i2)' ) month
      READ( FILENAME(i+7:i+8), '(i2)' ) day

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = GET_TAU0( month, day, year )
      
#if defined( MOPITT_V6_CO_OBS ) .or. defined( MOPITT_V7_CO_OBS )

      ! Opening an HDF-EOS5 swath file
      fId = he5_swopen(FILENAME, HE5F_ACC_RDONLY)
      
      ! Attaching to a swath object
      swathid = he5_swattach(fId, 'MOP02' )
      
      !=================================================================
      ! Seconds in day (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SecondsinDay", rank, dims8,
     &                   ntype, dimlist, maxdimlist)
     
      ! Allocate arrays
      ALLOCATE( SECONDS_IN_DAY( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SECONDS_IN_DAY' )

      ALLOCATE( TAU( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'TAU' )

      ALLOCATE( MOPITT_GMT( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'MOPITT_GMT' )
      
      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SecondsinDay', 
     &     START1, STRIDE1, EDGE1, SECONDS_IN_DAY)
      
      ! Compute GMT of MOPITT observations
      MOPITT_GMT = ( DBLE( SECONDS_IN_DAY ) / 3600d0 )

      ! Compute TAU values for GAMAP from SECONDS_IN_DAY
      TAU        = MOPITT_GMT + TAU0

      ! Save time dimension in T_DIM
      T_DIM      = dims8(1)
      
      !=================================================================
      ! LONGITUDE (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Longitude", rank, dims8,
     &                   ntype, dimlist, maxdimlist)
     
      ! Allocate array
      ALLOCATE( LONGITUDE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LONGITUDE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Longitude', 
     &     START1, STRIDE1, EDGE1, LONGITUDE)
      
      !=================================================================
      ! LATITUDE (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Latitude", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( LATITUDE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LATITUDE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Latitude', 
     &     START1, STRIDE1, EDGE1, LATITUDE)
      
      !=================================================================
      ! PRESSURE (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Pressure", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(PRESSURE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'PRESSURE' )
      
      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Pressure', 
     &     START1, STRIDE1, EDGE1, PRESSURE)

      ! Save PRESSURE dimension in Z_DIM
      Z_DIM = dims8(1)
      
      !=================================================================
      ! Cloud Description (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "CloudDescription", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CLOUD_DES( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CLOUD_DES' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'CloudDescription', 
     &     START1, STRIDE1, EDGE1, CLOUD_DES)
      
      !=================================================================
      ! Surface Index (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SurfaceIndex", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(SURFACE_INDEX( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SURFACE_INDEX' )
      
      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SurfaceIndex', 
     &     START1, STRIDE1, EDGE1, SURFACE_INDEX)
      
      !=================================================================
      ! Retrieval Bottom Pressure (1-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SurfacePressure", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(BOTTOM_PRESSURE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE' )
      
      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SurfacePressure', 
     &     START1, STRIDE1, EDGE1, BOTTOM_PRESSURE)
      
      !=================================================================
      ! CO Mixing Ratio (3-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOMixingRatioProfile", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_MIXING_RATIO( dims8(1), dims8(2), dims8(3) ), 
     &          stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MIXING_RATIO' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)
        
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOMixingRatioProfile', 
     &     START3, STRIDE3, EDGE3, CO_MIXING_RATIO)
      
      !=================================================================
      ! SDATA field: CO Retrieval Bottom Mixing Ratio (2-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOSurfaceMixingRatio", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_RET_BOT_MIXING_RATIO( dims8(1), dims8(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_RET_BOT_MIXING_RATIO' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOSurfaceMixingRatio', 
     &     START2, STRIDE2, EDGE2, CO_RET_BOT_MIXING_RATIO)
      
      !=================================================================
      ! CO Total Column (2-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOTotalColumn", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_TOTAL_COLUMN( dims8(1), dims8(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_TOTAL_COLUMN' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOTotalColumn', 
     &     START2, STRIDE2, EDGE2, CO_TOTAL_COLUMN)
      
      !=================================================================
      ! Retrieval Averaging Kernel Matrix (3-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievalAveragingKernelMatrix", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( AVGKER( dims8(1), dims8(2), dims8(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'AVGKER' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievalAveragingKernelMatrix', 
     &     START3, STRIDE3, EDGE3, AVGKER)
      
      !=================================================================
      ! A Priori CO Mixing Ratio Profile (3-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "APrioriCOMixingRatioProfile", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_MR_AP( dims8(1), dims8(2), dims8(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'APrioriCOMixingRatioProfile', 
     &     START3, STRIDE3, EDGE3, CO_MR_AP)

      !=================================================================
      ! A Priori CO Surface Mixing Ratio (2-D)
      !=================================================================
      
      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "APrioriCOSurfaceMixingRatio", 
     &                   rank, dims8, ntype, dimlist, maxdimlist)
     
      ! Allocate array
      ALLOCATE( CO_MR_AP_BOTTOM( dims8(1), dims8(2)), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_BOTTOM' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)
      
      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'APrioriCOSurfaceMixingRatio', 
     &     START2, STRIDE2, EDGE2, CO_MR_AP_BOTTOM)
     
      ! Detaching from the swath object
      as = he5_swdetach(swathid)
      
      ! Closing the file
      as = he5_swclose(fId)


#endif  !MOPITT v6

#if defined( MOPITT_V5_CO_OBS )

      ! Open file for HDF-VDATA interface
      CALL vdOpen( FILENAME )

      !=================================================================
      ! VDATA field: Time (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Seconds in Day', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate arrays
      ALLOCATE( SECONDS_IN_DAY( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SECONDS_IN_DAY' )

      ALLOCATE( TAU( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'TAU' )

      ALLOCATE( MOPITT_GMT( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'MOPITT_GMT' )

      ! Read data
      CALL vdGetData( vId, vSize, SECONDS_IN_DAY )

      ! Close field
      CALL vdCloseField( vId )

      ! Compute GMT of MOPITT observations
      MOPITT_GMT = ( DBLE( SECONDS_IN_DAY ) / 3600d0 )

      ! Compute TAU values for GAMAP from SECONDS_IN_DAY
      TAU       = MOPITT_GMT + TAU0

      ! Save time dimension in T_DIM
      T_DIM      = vSize

      !=================================================================
      ! VDATA field: LONGITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Longitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LONGITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LONGITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LONGITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: LATITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Latitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LATITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LATITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LATITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Cloud Description (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Cloud Description', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( CLOUD_DES( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CLOUD_DES' )

      ! Read data
      CALL vdGetData( vId, vSize, CLOUD_DES )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Surface Index (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Surface Index', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( SURFACE_INDEX( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SURFACE_INDEX' )

      ! Read data
      CALL vdGetData( vId, vSize, SURFACE_INDEX )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: PRESSURE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Pressure Grid', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( PRESSURE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'PRESSURE' )

      ! Read data
      CALL vdGetData( vId, vSize, PRESSURE )

      ! Close field
      CALL vdCloseField( vId )

      ! Save PRESSURE dimension in Z_DIM
      Z_DIM = vSize

      !=================================================================
      ! VDATA field: Retrieval Bottom Pressure (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Surface Pressure', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( BOTTOM_PRESSURE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE' )

      ! Read data
      CALL vdGetData( vId, vSize, BOTTOM_PRESSURE )

      ! Close field
      CALL vdCloseField( vId )

      ! Close HDF-VDATA interface
      CALL vdClose( FILENAME )
      


      ! Open file for HDF-SDATA interface
      CALL sdOpen( FILENAME )

      !=================================================================
      ! SDATA field: CO Mixing Ratio (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieved CO Mixing Ratio Profile', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MIXING_RATIO( dims(1), dims(2), dims(3) ), 
     &          stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), 
     &                CO_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Retrieval Bottom Mixing Ratio (2-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieved CO Surface Mixing Ratio', sId )
      
      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_RET_BOT_MIXING_RATIO( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_RET_BOT_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_RET_BOT_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Total Column (2-D)
      !=================================================================

      ! Open field

      CALL sdOpenFieldByName( 'Retrieved CO Total Column', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_TOTAL_COLUMN( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_TOTAL_COLUMN' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_TOTAL_COLUMN )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: Retrieval Averaging Kernel Matrix (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieval Averaging Kernel Matrix', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( AVGKER( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'AVGKER' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), AVGKER )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Mixing Ratio Profile (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Mixing Ratio Profile', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), CO_MR_AP )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Surface Mixing Ratio (2-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Surface Mixing Ratio', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP_BOTTOM( dims(1), dims(2)), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_BOTTOM' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_MR_AP_BOTTOM )

      ! Close field
      CALL sdCloseField( sId )

      ! Close file and quit
      CALL sdClose( FILENAME )
      
#endif  !MOPITT v5

      ! Return to calling program
      END SUBROUTINE READ_MOP02

!------------------------------------------------------------------------------------

      SUBROUTINE READ_ERROR_VARIANCE
!
!******************************************************************************
!  Subroutine READ_ERROR_VARIANCE reads observation error from binary punch files
!  (zhe 4/20/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,   ONLY : GET_TAUb

      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_ERROR_VARIANCE begins here!
      !=================================================================

      ! Filename
        FILENAME = TRIM( 'OBS_ERR_' ) // GET_RES_EXT()

      ! Echo some information to the standard output
        WRITE( 6, 110 ) TRIM( FILENAME )
 110    FORMAT( '     - READ_ERROR_VARIANCE: Reading ERR_PERCENT
     &                from: ', a )

      ! Read data from the binary punch file
        CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,
     &           GET_TAUb(),    IGLOB,     JGLOB,
     &           1,  ERR_PERCENT,  QUIET=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_ERROR_VARIANCE
  
!------------------------------------------------------------------------------

      SUBROUTINE INFO_MOP02( FILENAME )
!
!******************************************************************************
!  Subroutine INFO_MOP02 Info prints info about all VDATA and SDATA fields
!  contained within the MOPITT HDF file. (bmy, 7/3/03, 4/27/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE HdfSdModule
      USE HdfVdModule

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      !=================================================================
      ! INFO_MOP02 begins here!
      !=================================================================

      ! Print HDF-VDATA variables
      CALL vdOpen( FILENAME )
      CALL vdPrintInfo
      CALL vdClose( FILENAME )

      ! Print HDF-SDATA variables
      CALL sdOpen( FILENAME )
      CALL sdPrintInfo
      CALL sdClose( FILENAME )

      ! Return to calling program
      END SUBROUTINE INFO_MOP02

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MOP02
!
!******************************************************************************
!  Subroutine CLEANUP_MOP02 deallocates all module arrays (bmy, 4/27/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_MOP02 begins here!
      !=================================================================
      IF ( ALLOCATED( LATITUDE        ) ) DEALLOCATE( LATITUDE        )
      IF ( ALLOCATED( LONGITUDE       ) ) DEALLOCATE( LONGITUDE       )
      IF ( ALLOCATED( PRESSURE        ) ) DEALLOCATE( PRESSURE        )
      IF ( ALLOCATED( CLOUD_DES       ) ) DEALLOCATE( CLOUD_DES       )
      IF ( ALLOCATED( SURFACE_INDEX   ) ) DEALLOCATE( SURFACE_INDEX   )
      IF ( ALLOCATED( TAU             ) ) DEALLOCATE( TAU             )
      IF ( ALLOCATED( SECONDS_IN_DAY  ) ) DEALLOCATE( SECONDS_IN_DAY  )
      IF ( ALLOCATED( MOPITT_GMT      ) ) DEALLOCATE( MOPITT_GMT      )
      IF ( ALLOCATED( BOTTOM_PRESSURE ) ) DEALLOCATE( BOTTOM_PRESSURE )
      IF ( ALLOCATED( CO_MIXING_RATIO ) ) DEALLOCATE( CO_MIXING_RATIO )
      
      IF ( ALLOCATED( CO_RET_BOT_MIXING_RATIO)) THEN
         DEALLOCATE( CO_RET_BOT_MIXING_RATIO )
      ENDIF
      
      IF ( ALLOCATED( CO_TOTAL_COLUMN ) ) DEALLOCATE( CO_TOTAL_COLUMN )
      IF ( ALLOCATED( AVGKER          ) ) DEALLOCATE( AVGKER          )
      IF ( ALLOCATED( PLEV_AP         ) ) DEALLOCATE( PLEV_AP         )
      IF ( ALLOCATED( CO_MR_AP        ) ) DEALLOCATE( CO_MR_AP        )
      IF ( ALLOCATED( CO_MR_AP_BOTTOM ) ) DEALLOCATE( CO_MR_AP_BOTTOM )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MOP02

!---------------------------------------------------------------------------------------------------


      END MODULE MOPITT_OBS_MOD
