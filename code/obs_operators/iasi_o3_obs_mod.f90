!$Id: IASI_o3_mod.f,v 1.3 2011/02/23 00:08:48 daven Exp $
MODULE IASI_O3_OBS_MOD
  
  IMPLICIT NONE 
  
  !mkeller
#include "CMN_SIZE"
  !#include 'netcdf.inc'

  !=================================================================
  ! MODULE VARIABLES
  !=================================================================   

  PRIVATE
      
  PUBLIC READ_IASI_O3_OBS
  PUBLIC CALC_IASI_O3_FORCE
 
  ! Parameters
  INTEGER, PARAMETER           :: MAXLEV = 41
  INTEGER, PARAMETER           :: MAXIASI = 500000

  ! Module variables

  ! IASI data
  REAL*8, ALLOCATABLE :: IASI_LON(:)
  REAL*8, ALLOCATABLE :: IASI_LAT(:)
  REAL*8, ALLOCATABLE :: IASI_TIME(:)
  REAL*8, ALLOCATABLE :: IASI_O3(:,:)
  REAL*8, ALLOCATABLE :: IASI_O3_STD(:,:)
  REAL*8, ALLOCATABLE :: IASI_O3_APR(:,:)
  REAL*8, ALLOCATABLE :: IASI_AIR_DEN(:,:)
  REAL*8, ALLOCATABLE :: IASI_O3_AVK(:,:,:)
  REAL*8, ALLOCATABLE :: IASI_QUAL_FLAG(:)
  REAL*8, ALLOCATABLE :: IASI_CLOUD_COVER(:)
  REAL*8, ALLOCATABLE :: IASI_SATELLITE_ZENITH(:)
  REAL*8, ALLOCATABLE :: IASI_SOLAR_ZENITH(:)
  REAL*8, ALLOCATABLE :: IASI_PRESSURE(:,:)
  REAL*8, ALLOCATABLE :: IASI_DOFS(:)
  REAL*8, ALLOCATABLE :: IASI_O3_LVL(:)

  ! MLS grid specification
  INTEGER :: N_IASI_NLA
  INTEGER :: N_IASI_NPR
  INTEGER :: N_IASI_NOB


  ! mkeller: logical flag to check whether data is available for given day
  LOGICAL :: DATA_PRESENT
CONTAINS
  !------------------------------------------------------------------------------

  SUBROUTINE READ_IASI_O3_OBS( YYYYMMDD, N_IASI_NOB )
    !
!******************************************************************************
!  Subroutine READ_IASI_O3_OBS reads the file and passes back info contained
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
!  (1 ) IASI    (IASI_O3_OBS) : IASI retrieval for current day 
!     
!  NOIASI:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!******************************************************************************
!
    ! Reference to f90 modules
    USE DIRECTORY_MOD,          ONLY : DATA_DIR
    USE NETCDF 
    USE TIME_MOD,               ONLY : EXPAND_DATE

#include "CMN_SIZE" ! size parameters


    ! Arguments
    INTEGER,            INTENT(IN)  :: YYYYMMDD
    
    ! local variables 
    INTEGER                         :: FID
    INTEGER                         :: LIASI
    INTEGER                         :: N_IASI_NOB
    INTEGER                         :: N, J
    INTEGER                         :: NLA_ID, NPR_ID, NOB_ID
    INTEGER                         :: O3_ID, AVK_ID, PRE_ID, O3E_ID, APR_ID
    INTEGER                         :: LAT_ID, LON_ID, TIM_ID, DOF_ID
    INTEGER                         :: SAZ_ID, SOZ_ID, CLR_ID, QUA_ID, APC_ID
    
    CHARACTER(LEN=5)                :: TMP
    CHARACTER(LEN=255)              :: READ_FILENAME
    CHARACTER(LEN=255)              :: FILENAME_IASIO3
    CHARACTER(LEN=255)              :: DIR_IASIO3
    CHARACTER(LEN=255)              :: DIR_MONTH_IASIO3

    REAL*8, PARAMETER               :: FILL = -999.0D0
    REAL*8, PARAMETER               :: TOL  = 1d-04
    REAL*8                          :: U(MAXLEV,MAXLEV)
    REAL*8                          :: VT(MAXLEV,MAXLEV)
    REAL*8                          :: S(MAXLEV)
    REAL*8                          :: TMP1
    REAL*8                          :: TEST(MAXLEV,MAXLEV)
    INTEGER                         :: I, II, III
    
    !=================================================================
    ! READ_IASI_O3_OBS begins here!
    !=================================================================
    CALL CLEANUP_IASI
    ! filename root 
    DIR_IASIO3 = '/users/jk/15/xzhang/IASI_O3/'
    DIR_MONTH_IASIO3 = '/YYYY/MM/'
    FILENAME_IASIO3 = 'IASI_FORLI_O3_metopa_YYYYMMDD_v20151001.nc'

    ! Expand date tokens in filename 
    CALL EXPAND_DATE( DIR_MONTH_IASIO3, YYYYMMDD, 9999)
    CALL EXPAND_DATE( FILENAME_IASIO3, YYYYMMDD, 9999 ) 
    
    ! Construct complete filename 
    !      READ_FILENAME = TRIM( DATA_DIR ) // TRIM( 'IASI_O3/' ) // 
    !     &                TRIM( READ_FILENAME )
    READ_FILENAME = TRIM(DIR_IASIO3) // TRIM(DIR_MONTH_IASIO3) // TRIM( FILENAME_IASIO3 )
    
    WRITE(6,*) '    - READ_IASI_O3_OBS: reading file: ', READ_FILENAME
    
    ! mkeller: check to see if file exists
    INQUIRE(FILE=READ_FILENAME, EXIST = DATA_PRESENT)
    
    IF (.NOT. DATA_PRESENT) THEN
       PRINT *,"IASI file '", TRIM(READ_FILENAME), " not found, "// "assuming that there is no data for this day."
       RETURN
    ELSE
       PRINT *,"IASI file found!"
    ENDIF
    
    ! Open file and assign file id (FID)
    CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 800 )
    
    !-------------------------------- 
    ! Get data record IDs
    !-------------------------------- 
    CALL CHECK( NF90_INQ_DIMID( FID, "nlayers",         NLA_ID),  811 )
    CALL CHECK( NF90_INQ_DIMID( FID, "npressures",         NPR_ID),  812 )
    CALL CHECK( NF90_INQ_DIMID( FID, "nobservations",         NOB_ID),  813 )

    CALL CHECK( NF90_INQ_VARID( FID, "ozone_partial_column_profile",         O3_ID ), 821 )
    CALL CHECK( NF90_INQ_VARID( FID, "averaging_kernels_matrix", AVK_ID ), 822 )
    CALL CHECK( NF90_INQ_VARID( FID, "atmosphere_pressure_grid",        PRE_ID ), 823 )
    CALL CHECK( NF90_INQ_VARID( FID, "ozone_partial_column_error", O3E_ID ), 824 )  
    CALL CHECK( NF90_INQ_VARID( FID, "ozone_apriori_partial_column_profile",APR_ID ), 825 )
    CALL CHECK( NF90_INQ_VARID( FID, "latitude",        LAT_ID ), 826 )
    CALL CHECK( NF90_INQ_VARID( FID, "longitude",       LON_ID ), 827 )
    CALL CHECK( NF90_INQ_VARID( FID, "time",        TIM_ID ), 828 )
    CALL CHECK( NF90_INQ_VARID( FID, "sun_zen_angle", SOZ_ID ), 829 )
    CALL CHECK( NF90_INQ_VARID( FID, "satellite_zen_angle", SAZ_ID ), 830 )
    CALL CHECK( NF90_INQ_VARID( FID, "cloud_cover", CLR_ID ), 831 )
    CALL CHECK( NF90_INQ_VARID( FID, "retrieval_quality_flag", QUA_ID ), 832 )
    CALL CHECK( NF90_INQ_VARID( FID, "air_partial_column_profile", APC_ID ), 833 )
    CALL CHECK( NF90_INQ_VARID( FID, "dofs", DOF_ID ), 834 )
    
      ! READ number of retrievals, NIASI 
    CALL CHECK( NF90_INQUIRE_DIMENSION( FID, NLA_ID, TMP, N_IASI_NLA),  841 )
    CALL CHECK( NF90_INQUIRE_DIMENSION( FID, NPR_ID, TMP, N_IASI_NPR),  842 )
    CALL CHECK( NF90_INQUIRE_DIMENSION( FID, NOB_ID, TMP, N_IASI_NOB),  843 )

    
    print*, 'IASI dimensions on layers, pressures and observations are ', N_IASI_NLA, N_IASI_NPR, N_IASI_NOB

    ALLOCATE(IASI_O3(N_IASI_NLA, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, O3_ID, IASI_O3), 851)
    ALLOCATE(IASI_O3_STD(N_IASI_NLA, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, O3E_ID, IASI_O3_STD), 852)
    ALLOCATE(IASI_O3_APR(N_IASI_NLA, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, APR_ID, IASI_O3_APR), 858)
    ALLOCATE(IASI_AIR_DEN(N_IASI_NLA, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, APC_ID, IASI_AIR_DEN), 859)
    ALLOCATE(IASI_SOLAR_ZENITH(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, SOZ_ID, IASI_SOLAR_ZENITH), 853)
    ALLOCATE(IASI_SATELLITE_ZENITH(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, SAZ_ID, IASI_SATELLITE_ZENITH), 861)
    ALLOCATE(IASI_PRESSURE(N_IASI_NPR, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, PRE_ID, IASI_PRESSURE), 854)
    ALLOCATE(IASI_LAT(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, LAT_ID, IASI_LAT), 855)
    ALLOCATE(IASI_LON(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, LON_ID, IASI_LON), 856)
    ALLOCATE(IASI_TIME(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, TIM_ID, IASI_TIME), 857)
    !PRINT *, "IASI_TIME", IASI_TIME(1:N_IASI_NOB)
    ALLOCATE(IASI_O3_AVK(N_IASI_NLA, N_IASI_NLA, N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, AVK_ID, IASI_O3_AVK), 860)
    ALLOCATE(IASI_QUAL_FLAG(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, QUA_ID, IASI_QUAL_FLAG), 862)
    ALLOCATE(IASI_CLOUD_COVER(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, CLR_ID, IASI_CLOUD_COVER), 863)
    ALLOCATE(IASI_DOFS(N_IASI_NOB))
    CALL CHECK( NF90_GET_VAR(FID, DOF_ID, IASI_DOFS), 864)
    ALLOCATE(IASI_O3_LVL(N_IASI_NOB))
    CALL CHECK( NF90_CLOSE( FID ), 899)
      !-------------------------------- 
      ! Calculate S_OER_INV
      !-------------------------------- 

      ! loop over records
    ! Now determine how many of the levels in O3 are
    ! 'good' and how many are just FILL.
    !ALLOCATE(IASI_O3_LVL(N_IASI_NOB))
    DO N = 1, N_IASI_NOB
       J = 1
       DO WHILE ( J .le. MAXLEV )
          ! check if the value is good
          IF (IASI_O3(J,N) > FILL ) THEN
             ! save the number of good levels as LTES               
             IASI_O3_LVL(N) = MAXLEV - J + 1
             ! and now we can exit the while loop
             J = MAXLEV + 1
             ! otherwise this level is just filler               
          ELSE
             ! so proceed to the next one up
             J = J + 1
          ENDIF
       ENDDO
    ENDDO

         !print*,  ' IASI TEST ', IASI(N)%O3 
         !print*,  ' IASI good ', IASI(N)%LIASI
         !print*,  ' IASI pres ', IASI(N)%PRES(1:J)

         ! Add a bit to the diagonal to regularize the inversion
         ! (ks, ml, dkh, 11/18/10) 
         ! mkeller: this makes no sense to me.
         !DO II=1,J
         !   IASI(N)%S_OER(II,II) = IASI(N)%S_OER(II,II)+ 0.001D0
         !ENDDO

         !CALL SVD( IASI(N)%S_OER(1:J,1:J), J, 
     !&                        U(1:J,1:J), S(1:J), 
     !&                       VT(1:J,1:J)          )

         ! U = S^-1 * U^T  
         !TEST = 0d0
         !DO I = 1, J

            ! mkeller: regularize matrix inverse by ignoring all singular values below a certain cutoff.
            !          This is horrendously inefficient, but should work for now. In the 
            !          future, Thikonov regularization should be implemented instead.            
            ! xzhang: svd TEST critical value changes from 1e-2 to 5e-2
            !IF ( S(I)/S(1) < 1e-2 ) THEN
               !S(I) = 1e-2 * S(1)
            !ENDIF
            !DO II = 1, J
               !TEST(I,II) = U(II,I) / S(I)
            !ENDDO              
         !ENDDO

         !TEST = 0d0
         !U    = TEST
         !TEST = 0d0

   
         ! S_OER_INV = V * S^-1 * U^T
         !DO I = 1, J
            !DO II = 1, J
               !TMP1 = 0d0
               !DO III = 1, J
                  !TMP1 = TMP1 + VT(III,I) * U(III,II)
               !ENDDO
               !IASI(N)%S_OER_INV(I,II) = TMP1
            !ENDDO
         !ENDDO
         
         ! TEST: calculate 2-norm of I - S_OER_INV * S_OER
         ! mkeller: comment this out for now; pointless given the regularization
         !          performed above. 
         !          Need to come up with an alternative TEST in the future.
         !DO I = 1, J
         !   DO II = 1, J
         !     TMP1 = 0d0
         !      DO III = 1, J
         !         TMP1 = TMP1 
         !&                 + IASI(N)%S_OER_INV(III,I) * IASI(N)%S_OER(III,II)
         !ENDDO
         !TEST(I,II) = - TMP1
         !ENDDO
         !TEST(I,I) = ( TEST(I,I) + 1 ) ** 2
         !ENDDO

         !IF ( SUM(TEST(1:J,1:J)) > TOL ) THEN  
         !   print*, ' WARNING: inversion error for retv N = ', 
         !&                SUM(TEST(1:J,1:J)), N 
         !       print*, '   in IASI obs ', READ_FILENAME 
         !    ENDIF 

   !ENDDO  ! N

      ! Return to calling program
  END SUBROUTINE READ_IASI_O3_OBS

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

   SUBROUTINE CALC_IASI_O3_FORCE( COST_FUNC )
        !
!******************************************************************************
!  Subroutine CALC_IASI_O3_FORCE calculaIASI the adjoint forcing from the IASI
!  O3 observations and updaIASI the cost function. (dkh, 02/15/09)
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost function                        [unitless]
!     
!     
!  NOIASI:
!  (1 ) Updated to GCv8 (dkh, 10/07/09) 
!  (2 ) Add more diagnostics.  Now read and write doubled O3 (dkh, 11/08/09) 
!  (3 ) Now use CSPEC_AFTER_CHEM and CSPEC_AFTER_CHEM_ADJ (dkh, 02/09/11) 
!******************************************************************************
!
     ! Reference to f90 modules
     USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
     USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
     USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
     USE ADJ_ARRAYS_MOD,     ONLY : O3_PROF_SAV
     USE ADJ_ARRAYS_MOD,     ONLY : ID2C
     USE CHECKPT_MOD,        ONLY : CHK_STT
     USE COMODE_MOD,         ONLY : JLOP
     USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
     USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
     USE DAO_MOD,            ONLY : AD, AIRDEN, AIRVOL
     USE DAO_MOD,            ONLY : BXHEIGHT
     USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
     USE GRID_MOD,           ONLY : GET_IJ, GET_XMID, GET_YMID
     USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
     USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
     USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_HOUR
     USE TRACER_MOD,         ONLY : XNUMOLAIR, STT
     USE TRACER_MOD,         ONLY : TCVV
     USE TRACERID_MOD,       ONLY : IDO3, IDTOX
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
      INTEGER                     :: LVL_CRT1, LVL_CRT2
      REAL*8                      :: GC_PRES(LLPAR+1)
      REAL*8                      :: GC_O3_NATIVE(LLPAR)
      REAL*8                      :: GC_O3(MAXLEV)
      REAL*8                      :: GC_PSURF, IASI_PSURF
      REAL*8                      :: MAP(LLPAR,MAXLEV), GC_O3_CONTRIB(LLPAR,MAXLEV)
      REAL*8                      :: O3_HAT(MAXLEV)
      REAL*8                      :: O3_PERT(MAXLEV)
      REAL*8                      :: FORCE(MAXLEV)
      REAL*8                      :: DIFF(MAXLEV), COST_CONTRIB(MAXLEV), COST_CONTRIB_COL
      REAL*8                      :: DIFF_COL, FORCE_COL, IASI_O3_STD_COL
      REAL*8                      :: IASI_PCENTER(MAXLEV)
      REAL*8                      :: NEW_COST(IIPAR,JJPAR)
      REAL*8                      :: OLD_COST
      REAL*8                      :: XNUAIR
      REAL*8, SAVE                :: TIME_FRAC(MAXIASI)

      REAL*8                      :: GC_O3_NATIVE_ADJ(LLPAR)

      REAL*8                      :: GC_ADJ_TEMP(IIPAR,JJPAR,LLPAR)
      REAL*8                      :: GC_ADJ_COUNT(IIPAR,JJPAR,LLPAR)
      REAL*8                      :: GC_ADJ_TEMP_COST(IIPAR,JJPAR)
      REAL*8                      :: IASI_O3_BIAS(IIPAR,JJPAR), IASI_O3_BIAS_FILT(IIPAR,JJPAR)
      REAL*8                      :: IASI_O3_CHI_SQ(IIPAR,JJPAR), IASI_O3_CHI_SQ_FILT(IIPAR,JJPAR)
      REAL*8                      :: IASI_O3_BIAS_SOBS(IIPAR,JJPAR), IASI_O3_BIAS_FILT_SOBS(IIPAR,JJPAR)
      REAL*8                      :: IASI_O3_CHI_SQ_SOBS(IIPAR,JJPAR), IASI_O3_CHI_SQ_FILT_SOBS(IIPAR,JJPAR)

      ! arrays needed for superobservations
      LOGICAL                     :: SUPER_OBS = .TRUE.  ! do super observations?
      REAL*8                      :: SOBS_COUNT(IIPAR,JJPAR)
   
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      
      !=================================================================
      ! CALC_IASI_O3_FORCE begins here!
      !=================================================================

      print*, '     - CALC_IASI_O3_FORCE '
      
      ! Reset 
      XNUAIR = 28.9644d-3
      NEW_COST = 0d0
      SOBS_COUNT = 0d0
      GC_ADJ_TEMP = 0d0
      GC_ADJ_TEMP_COST = 0d0
      GC_ADJ_COUNT = 0d0
      IASI_O3_BIAS = 0d0
      IASI_O3_BIAS_FILT = 0d0
      IASI_O3_CHI_SQ = 0d0
      IASI_O3_CHI_SQ_FILT = 0d0
      IASI_O3_BIAS_SOBS = 0d0
      IASI_O3_BIAS_FILT_SOBS = 0d0
      IASI_O3_CHI_SQ_SOBS = 0d0
      IASI_O3_CHI_SQ_FILT_SOBS = 0d0
      GC_HOUR = GET_HOUR()
      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      PRINT *, "GET_NHMS", GET_NHMS()

      ! Check if it is the last hour of a day
      IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN

         ! Read the IASI O3 file for this day
         CALL READ_IASI_O3_OBS( GET_NYMD(), N_IASI_NOB )

         ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
         TIME_FRAC(1:N_IASI_NOB) = IASI_TIME(1:N_IASI_NOB)/240000d0
      ENDIF
      IF ( FIRST ) THEN
         FILENAME = 'chi_sq2_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 801, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'sobs_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 802, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'chi_sq2_filt_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 803, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'diff2_filt_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 804, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'sobs_count_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 805, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'sobs_count_filt_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 806, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

         FILENAME = 'lat_orb_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 813, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )
         
         FILENAME = 'lon_orb_iasio3.NN.m'
         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
         OPEN( 814, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED', ACCESS='SEQUENTIAL' )

      ENDIF

      ! Get the range of TES retrievals for the current hour
      CALL GET_NT_RANGE( N_IASI_NOB, GET_NHMS(), TIME_FRAC, NTSTART, NTSTOP )

      IF ( NTSTART == 0 .and. NTSTOP == 0 ) THEN

         print*, ' No matching IASI O3 obs for this hour'
        RETURN
      ENDIF

      print*, ' for hour range: ', GET_NHMS(), IASI_TIME(NTSTART), IASI_TIME(NTSTOP)
      print*, ' found record range: ', NTSTART, NTSTOP


      DO I_IASI = NTSTART, NTSTOP, -1
         !PRINT *, "IASI_TIME", IASI_TIME(I_IASI)
         IF ( ( IASI_LAT(I_IASI) > -60d0 ) .AND. &
              ( IASI_LAT(I_IASI) < 75d0 ) .AND. &
              ( IASI_QUAL_FLAG(I_IASI) == 0d0 ) .AND. &
              ( IASI_SOLAR_ZENITH(I_IASI) < 80d0 ) .AND. &
              ( IASI_DOFS(I_IASI) > 2.75 ) .AND. &
              ( IASI_CLOUD_COVER(I_IASI) < 13d0 ) ) THEN 

                 !( IASI_SATELLITE_ZENITH(I_IASI) > 0d0 ) ) THEN
               ! Get model grid coordinate indices that correspond to the observation
                  ! For safety, initialize these up to LLTES
            GC_O3(:)       = 0d0
            MAP(:,:)       = 0d0
            GC_O3_CONTRIB(:,:) = 0d0
            FORCE(:)       = 0d0
            DIFF(:)        = 0d0
            COST_CONTRIB(:) = 0d0
            COST_CONTRIB_COL = 0d0
            DIFF_COL       = 0d0
            IASI_O3_STD_COL = 0d0
            FORCE_COL = 0d0
            LVL_CRT1 = 0
            LVL_CRT2 = 0
            LIASI = IASI_O3_LVL(I_IASI)
            IIJJ = GET_IJ(REAL(IASI_LON(I_IASI),4),REAL(IASI_LAT(I_IASI),4))
            IASI_PCENTER = 0d0
            I = IIJJ(1)
            J = IIJJ(2)

            L0 = MAXLEV - IASI_O3_LVL(I_IASI)

            ! Get GC pressure levels (mbar)
            DO L = 1, LLPAR+1
               GC_PRES(L) = GET_PEDGE(I,J,L)
            ENDDO
            DO L = 1, LIASI
               IASI_PCENTER(L) = 0.5d0*(IASI_PRESSURE(L+L0,I_IASI)+IASI_PRESSURE(L+L0+1,I_IASI))
            ENDDO
            ! Get GC surface pressure (mbar)
            GC_PSURF = GET_PEDGE(I,J,1)

            IASI_PSURF = IASI_PRESSURE(1+L0,I_IASI)
            ! Calculate the interpolation weight matrix
            MAP(1:LLPAR,1:LIASI) = GET_INTMAP( LLPAR, GC_PRES, GC_PSURF, &
                 LIASI,  IASI_PRESSURE(1+L0:LIASI+L0,I_IASI), IASI_PSURF )
            DO L = 1, LLPAR
               IF ( GC_PRES(L) > 300d0 ) THEN
                  GC_O3_NATIVE(L) = CHK_STT(I,J,L,IDTOX) * TCVV(IDTOX) * BXHEIGHT(I,J,L)/(AIRVOL(I,J,L)*XNUAIR)
               ELSE
                  GC_O3_NATIVE(L) = STT(I,J,L,IDTOX) * TCVV(IDTOX) * BXHEIGHT(I,J,L)/(AIRVOL(I,J,L)*XNUAIR)
               ENDIF
            ENDDO
            !PRINT *, "IASI_PCENTER", IASI_PCENTER
            !PRINT *, "GC_PRES", GC_PRES
            !PRINT *, " GC_O3_NATIVE", GC_O3_NATIVE
            ! Interpolate GC O3 column to TES grid
            DO L = 1, LIASI
               GC_O3(L) = 0d0
               DO LL = 1, LLPAR
                  GC_O3_CONTRIB(LL,L) = MAP(LL,L) * GC_O3_NATIVE(LL) 
                  GC_O3(L) = GC_O3(L) + MAP(LL,L) * GC_O3_NATIVE(LL)
               ENDDO
            ENDDO
            DO L = 2, LIASI-1
               IF (GC_O3(L) == 0d0) THEN
                  IF (GC_O3(L-1) > 0d0) THEN
                     LVL_CRT1 = L-1
                  ENDIF
                  IF (GC_O3(L+1) > 0d0) THEN
                     LVL_CRT2 = L+1
                  ENDIF
                  IF (REAL(LVL_CRT1)*REAL(LVL_CRT2)>0d0) THEN
                     DO LL = LVL_CRT1,LVL_CRT2
                        GC_O3(LL) = ((GC_O3(LVL_CRT1)+GC_O3(LVL_CRT2))/SUM(IASI_O3(LVL_CRT1+L0:LVL_CRT2+L0,I_IASI)))* &
                             IASI_O3(LL+L0,I_IASI)
                     ENDDO
                     LVL_CRT1 = 0
                     LVL_CRT2 = 0
                  ENDIF
               ENDIF
            ENDDO
            !PRINT *, " GC_O3", GC_O3
            !--------------------------------------------------------------
            ! Apply IASI O3 observation operator
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
            DO L = 1, LIASI
               IF (IASI_O3_APR(L+L0,I_IASI) > 0) THEN
                  GC_O3(L)   = MAX(GC_O3(L), 1d-4)
                  O3_PERT(L) = GC_O3(L) - IASI_O3_APR(L+L0,I_IASI)
               ENDIF
            ENDDO
            
            ! x_a + A_k * ( x_m - x_a )
            DO L = 1, LIASI
               O3_HAT(L) = 0d0
               FORCE(L) = 0d0
               IF (IASI_O3_APR(L+L0,I_IASI) > 0) THEN
                  DO LL = 1, LIASI
                     O3_HAT(L) = O3_HAT(L) + IASI_O3_AVK(LL+L0,L0+L,I_IASI) * O3_PERT(LL)
                  ENDDO
                  O3_HAT(L)    = O3_HAT(L) + IASI_O3_APR(L+L0,I_IASI)
               ENDIF
               ! actual comparison with bias correction information
               IF ( ( IASI_O3(L+L0,I_IASI) > 0d0 ) .AND. &
                    ( O3_HAT(L) > 0d0            ) .AND. &
                    ( IASI_O3(L+L0,I_IASI) < 100d0 * IASI_O3_STD(L+L0,I_IASI) ) ) THEN
                  IF ( REAL(IASI_LAT(I_IASI)) >= 60d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN                        
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.049
                     ELSEIF (IASI_PCENTER(L) >= 150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.805
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.975
                     ENDIF
                  ELSEIF ( REAL(IASI_LAT(I_IASI)) >= 30d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.132
                     ELSEIF (IASI_PCENTER(L) >=150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.018
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.921
                     ENDIF
                  ELSEIF ( REAL(IASI_LAT(I_IASI)) >= 0d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN                        
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.142
                     ELSEIF (IASI_PCENTER(L) >= 150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.934
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.876
                     ENDIF
                  ELSEIF ( REAL(IASI_LAT(I_IASI)) >= -30d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.122
                     ELSEIF (IASI_PCENTER(L) >= 150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.963
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.886
                     ENDIF
                  ELSEIF ( REAL(IASI_LAT(I_IASI)) >= -60d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*1.12
                     ELSEIF (IASI_PCENTER(L) >= 150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.994
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.94
                     ENDIF
                  ELSEIF ( REAL(IASI_LAT(I_IASI)) >= -90d0 ) THEN
                     IF (IASI_PCENTER(L) >= 300d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.98
                     ELSEIF (IASI_PCENTER(L) >= 150d0) THEN
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.765
                     ELSE
                        DIFF(L) = O3_HAT(L) - IASI_O3(L+L0,I_IASI)*0.71
                     ENDIF
                  ENDIF
                  !PRINT *, "DIFF", DIFF(L)
               ELSE
                  DIFF(L) = 0d0
               ENDIF
               !PRINT *, "O3_HAT", O3_HAT(L)
               !PRINT *, "IASI_O3", IASI_O3(L+L0,I_IASI)
               !PRINT *, "IASI_O3_STD", IASI_O3_STD(L+L0, I_IASI)
            ! Forcing is computed in profile and in columns
               FORCE(L) = DIFF(L)/(0.5*IASI_O3_STD(L+L0,I_IASI))**2
               COST_CONTRIB(L) = 0.5d0 * DIFF(L) * FORCE(L)
               IF ( ( IASI_PCENTER(L) >= 300d0 ) ) THEN
                  DIFF_COL = DIFF_COL + DIFF(L)
                  IASI_O3_STD_COL = IASI_O3_STD_COL + IASI_O3_STD(L+L0,I_IASI)**2
               ENDIF
            ENDDO

            FORCE_COL = DIFF_COL/IASI_O3_STD_COL
            ! Diagnostics of chi square and model-obs biases
            IF ( ( IASI_PCENTER(8) >= 300d0 ).AND. &
                 ( SQRT(COST_CONTRIB(8)) < 10d0 ) .AND. &
                 ( SQRT(COST_CONTRIB(8)) > 0d0  ) ) THEN
               IASI_O3_BIAS_FILT(I,J) = IASI_O3_BIAS_FILT(I,J)+ DIFF(8)
               IASI_O3_CHI_SQ_FILT(I,J) = IASI_O3_CHI_SQ_FILT(I,J) + 2*COST_CONTRIB(8)
            ENDIF
            COST_CONTRIB_COL = 0.5d0 * DIFF_COL * FORCE_COL
            IF ( ( COST_CONTRIB_COL > 0d0)  .AND. &
                 ( COST_CONTRIB_COL < 200d0) ) THEN
               ! adjoint of interpolation
               DO L = 1, LLPAR
                  ! Adjoint of unit conversion
                  GC_O3_NATIVE_ADJ(L) = FORCE_COL  * TCVV(IDTOX) * BXHEIGHT(I,J,L)/ (XNUAIR*AIRVOL(I,J,L))
                  GC_ADJ_TEMP(I,J,L) = GC_ADJ_TEMP(I,J,L) + GC_O3_NATIVE_ADJ(L)
                  GC_ADJ_COUNT(I,J,L) = GC_ADJ_COUNT(I,J,L) + 1d0
               ENDDO
            
               NEW_COST(I,J) = NEW_COST(I,J) + COST_CONTRIB_COL!0.5 * DIFF_COL * FORCE_COL
               IF (SUPER_OBS) THEN
                  !IF (COST_CONTRIB_COL > 0d0) THEN
                  SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                  IASI_O3_BIAS(I,J) = IASI_O3_BIAS(I,J) + DIFF(8)
                  IASI_O3_CHI_SQ(I,J) = IASI_O3_CHI_SQ(I,J) + 2*COST_CONTRIB(8)
                  !ENDIF
               ENDIF
            ENDIF
         
 110        FORMAT(F18.6,1X)
640      ENDIF
      ENDDO  ! NT

      ! Compute adjoint forcing and cost function in superobservation
      IF (SUPER_OBS) THEN
         DO I=1,IIPAR
            DO J=1,JJPAR
               IF ( SOBS_COUNT(I,J) > 0d0 ) THEN
                  DO L=1,LLPAR
                     IF ( ( GET_PCENTER(I,J,L) >= 300d0 )  .AND. &
                          ( GC_ADJ_COUNT(I,J,L) > 0d0 )  ) THEN
                        STT_ADJ(I,J,L,IDTOX) = STT_ADJ(I,J,L,IDTOX) + GC_ADJ_TEMP(I,J,L)/GC_ADJ_COUNT(I,J,L)
                     ENDIF
                  
                  ENDDO
                  COST_FUNC = COST_FUNC + NEW_COST(I,J)/SOBS_COUNT(I,J)
                  IASI_O3_BIAS_SOBS(I,J) = IASI_O3_BIAS(I,J)/SOBS_COUNT(I,J)
                  IASI_O3_CHI_SQ_SOBS(I,J) = IASI_O3_CHI_SQ(I,J)/SOBS_COUNT(I,J)
               ENDIF
               IF (GC_ADJ_COUNT(I,J,8) > 0d0) THEN
                  IASI_O3_BIAS_FILT_SOBS(I,J) = IASI_O3_BIAS_FILT(I,J)/GC_ADJ_COUNT(I,J,8)
                  IASI_O3_CHI_SQ_FILT_SOBS(I,J) = IASI_O3_CHI_SQ_FILT(I,J)/GC_ADJ_COUNT(I,J,8)
               ENDIF
               WRITE(801,110)(IASI_O3_CHI_SQ_SOBS(I,J))
               WRITE(802,110)(1e6*IASI_O3_BIAS_SOBS(I,J))
               WRITE(803,110)(IASI_O3_CHI_SQ_FILT_SOBS(I,J))
               WRITE(804,110)(1e6*IASI_O3_BIAS_FILT_SOBS(I,J))
               WRITE(805,110)(SOBS_COUNT(I,J))
               WRITE(806,110)(GC_ADJ_COUNT(I,J,8))
               WRITE(813,110)(GET_XMID(I))
               WRITE(814,110)(GET_YMID(J))
            ENDDO
         ENDDO
      ENDIF
   


      IF ( FIRST ) FIRST = .FALSE.
      
      print*, ' Updated value of COST_FUNC = ', COST_FUNC
      print*, ' IASI contribution           = ', COST_FUNC - OLD_COST

      ! Return to calling program
    END SUBROUTINE CALC_IASI_O3_FORCE

!----------------------------------------------------------------------------
    SUBROUTINE CLEANUP_IASI

      IF(ALLOCATED(IASI_LON)) DEALLOCATE(IASI_LON)
      IF(ALLOCATED(IASI_LAT)) DEALLOCATE(IASI_LAT)
      IF(ALLOCATED(IASI_TIME)) DEALLOCATE(IASI_TIME)
      IF(ALLOCATED(IASI_O3_AVK)) DEALLOCATE(IASI_O3_AVK)
      IF(ALLOCATED(IASI_O3_APR)) DEALLOCATE(IASI_O3_APR)
      IF(ALLOCATED(IASI_O3_STD)) DEALLOCATE(IASI_O3_STD)
      IF(ALLOCATED(IASI_O3)) DEALLOCATE(IASI_O3)
      IF(ALLOCATED(IASI_SATELLITE_ZENITH)) DEALLOCATE(IASI_SATELLITE_ZENITH)
      IF(ALLOCATED(IASI_SOLAR_ZENITH)) DEALLOCATE(IASI_SOLAR_ZENITH)
      IF(ALLOCATED(IASI_CLOUD_COVER)) DEALLOCATE(IASI_CLOUD_COVER)
      IF(ALLOCATED(IASI_QUAL_FLAG)) DEALLOCATE(IASI_QUAL_FLAG)
      IF(ALLOCATED(IASI_AIR_DEN)) DEALLOCATE(IASI_AIR_DEN)
      IF(ALLOCATED(IASI_PRESSURE)) DEALLOCATE(IASI_PRESSURE)
      IF(ALLOCATED(IASI_O3_LVL)) DEALLOCATE(IASI_O3_LVL)
      IF(ALLOCATED(IASI_DOFS)) DEALLOCATE(IASI_DOFS)
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
      LOGICAL               :: FOUND_ALL_RECORDS
      INTEGER               :: NTEST
      INTEGER               :: HH, MM, SS
      REAL*8                :: GC_HH_FRAC
      REAL*8                :: H1_FRAC

      !=================================================================
      ! GET_NT_RANGE begins here!
      !=================================================================
      ! Initialize
      FOUND_ALL_RECORDS  = .FALSE.
      NTSTART            = 0
      NTSTOP             = 0

      ! set NTSAVE to NTES every time we start with a new file
      IF ( HHMMSS == 230000 ) NTSAVE = N_IASI_NOB

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

         print*, 'All records have been read already '
         RETURN

      ! No records reached yet
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC < GC_HH_FRAC ) THEN


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

      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP, LTM_TOP, TM_PRESC, TM_SURFP ) RESULT( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!     
!  NOIASI:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      REAL*8             :: GC_PRESC(LGC_TOP+1)
      REAL*8             :: TM_PRESC(LTM_TOP+1) 
      REAL*8             :: GC_SURFP
      REAL*8             :: TM_SURFP
 
      ! Return value 
      REAL*8             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables 
      INTEGER  :: LGC, LTM
      REAL*8   :: DIFF, DELTA_SURFP
      REAL*8   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0 
  
!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN 
!            CALL ERROR_STOP( 'highly unlikey', 
!     &                       'read_sciano2_mod.f')
!         ENDIF 
!
!      ENDDO 
      

      ! Loop over each pressure level of TM grid
      DO LGC = 1, LGC_TOP
 
         ! Find the levels from GC that bracket level LTM
         DO LTM = 1, LTM_TOP

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP

            ! Linearly interpolate value on the LTM grid 
            IF ( GC_PRESC(LGC) <= TM_PRESC(LTM ) ) THEN
               IF (GC_PRESC(LGC+1)  >= TM_PRESC(LTM+1)) THEN
                  HINTERPZ(LGC,LTM) = 1D0
               ELSEIF (GC_PRESC(LGC) >= TM_PRESC(LTM+1)) THEN
                  DIFF                = HI - LOW  
                  HINTERPZ(LGC,LTM) = ( HI - TM_PRESC(LTM+1)  ) / DIFF
                  HINTERPZ(LGC,LTM+1) = ( TM_PRESC(LTM+1) - LOW ) / DIFF
               ENDIF
            ELSEIF (GC_PRESC(LGC) > TM_PRESC(1) ) THEN
               HINTERPZ(LGC,1) = 1D0
            ENDIF

 
            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO

       ! Correct for case where IASI pressure is higher than the
       ! highest GC pressure.  In this case, just 1:1 map. 
       !DO LTM = 1, LTM_TOP
          !IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             !HINTERPZ(:,LTM)   = 0D0
             !HINTERPZ(LTM,LTM) = 1D0
          !ENDIF
       !ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP



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
      
      END MODULE IASI_O3_OBS_MOD
