MODULE OMI_CH2O_OBS_MOD

!
!
! Module OMI_CH2O_OBS contains all subroutines and variables needed for OMH CH2O column data
!
!
! Module Routines:
!
! (1) READ_OMI_CH2O_FILE       : Read OMI hdf file


  IMPLICIT NONE
  
#include "CMN_SIZE"
  
  PRIVATE

  PUBLIC READ_OMI_CH2O_FILE
  PUBLIC CALC_OMI_CH2O_FORCE
  
  !Arrays for diagnostic outpit
  ! Module variables
      
  ! Diagnostic arrays for output in netCDF file
  REAL*8, ALLOCATABLE :: OMH_LON(:,:), OMH_LAT(:,:)
  REAL*8, ALLOCATABLE :: OMH_TIME(:), OMH_AMF_TROP(:,:)
  REAL*8, ALLOCATABLE :: OMH_CH2O_TROP(:,:), OMH_CH2O_TROP_STD(:,:)
  REAL*8, ALLOCATABLE :: OMH_VIEW_ZEN(:,:), OMH_SOLAR_ZEN(:,:)
  REAL*8, ALLOCATABLE :: OMH_CFR(:,:), OMH_SCW_PRE(:,:,:)
  REAL*8, ALLOCATABLE :: OMH_SCW(:,:,:)
  REAL*8, ALLOCATABLE :: OMH_X_Q_FLAG(:,:), OMH_M_Q_FLAG(:,:)
  REAL*8, ALLOCATABLE :: LON_ORB(:,:), LAT_ORB(:,:)
  REAL*8, ALLOCATABLE :: TIME_ORB(:), AMF_TROP_ORB(:,:)
  REAL*8, ALLOCATABLE :: CH2O_TROP(:,:), CH2O_TROP_STD(:,:)
  REAL*8, ALLOCATABLE :: VIEW_ZEN(:,:), SOLAR_ZEN(:,:)
  REAL*8, ALLOCATABLE :: CFR(:,:), SCW_PRE(:,:,:)
  REAL*8, ALLOCATABLE :: SCW(:,:,:)
  REAL*8, ALLOCATABLE :: X_Q_FLAG(:,:), M_Q_FLAG(:,:)

  INTEGER :: N_OMH_ORB
  INTEGER, PARAMETER :: MAX_ORB = 50000
  INTEGER, PARAMETER :: N_OMH_SWATHS = 60
  INTEGER, PARAMETER :: N_OMH_LEVELS = 47

  REAL*4:: OMH_CH2O_MEAN(IIPAR,JJPAR) = 0d0 ! Mean OMI column
  REAL*4:: OMH_GEOS_CH2O_MEAN(IIPAR,JJPAR) = 0d0 ! Mean GEOS-Chem column
  REAL*4:: OMH_CH2O_ERR_MEAN(IIPAR,JJPAR) = 0d0 ! Mean OMI observation error
     
  REAL*4:: OMH_BIAS(IIPAR,JJPAR)=0d0 ! Model bias 
  REAL*4:: OMH_VAR(IIPAR,JJPAR)=0d0 ! Model variance
  REAL*4:: OMH_DELTA=0d0 ! temporary storage variable
  REAL*4:: OMH_BIAS_COUNT(IIPAR,JJPAR) = 0d0 ! counter for number of observations
  REAL*4:: OMH_CHI_SQUARED(IIPAR,JJPAR) = 0d0 ! Chi-squared values
  
CONTAINS
      
!-----------------------------------------------------------------------------!

  SUBROUTINE READ_OMI_CH2O_FILE ( YYYYMMDD, N_OMH_ORB )

    USE ERROR_MOD, ONLY : ALLOC_ERR
    USE TIME_MOD,  ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR
    USE HDF5
    USE TIME_MOD,  ONLY : GET_HOUR, GET_DAY
    ! Arguments
    INTEGER, INTENT(IN) :: YYYYMMDD!, HHMMSS

    CHARACTER(LEN=255) :: DIR_OMH
    CHARACTER(LEN=255) :: FILENAME_OMH
    CHARACTER(LEN=255) :: FILENAME_FULL
    CHARACTER(LEN=255) :: DSET_NAME
    CHARACTER(255) :: ORB_PATH,FILE_ORB,FILE_OMH
    CHARACTER(2) :: I_CHAR
    INTEGER :: IO_ORB_STATUS, IO_STATUS
    INTEGER :: DAY
    INTEGER(HID_T) :: file_id, dset_id, dspace_id
    !INTEGER(HSIZE_T) :: data_dims
    INTEGER :: error, DIMS_COUNTER, DIMS_INDEX
    INTEGER :: N_OMH_ORB
    INTEGER(HID_T) :: file_id_orb
    INTEGER(HID_T) :: dset_id_orb
    INTEGER(HID_T) :: dspace_id_orb
    INTEGER :: error_orb
    CHARACTER(LEN=255) :: filename_orb, dsetname

    INTEGER :: rank_orb, rank_omh
    INTEGER(HSIZE_T) :: dims_orb(3), maxdims_orb(3)
    INTEGER(HSIZE_T) :: dims_omh(3), maxdims_omh(3)
    INTEGER(HSIZE_T) :: DATA_DIMS_ORB(3)
    INTEGER :: GC_HOUR
    LOGICAL          :: DATA_VALID

    CALL CLEANUP_OMH

    GC_HOUR = GET_HOUR()
    DAY = GET_DAY()
    DIR_OMH = '/users/jk/15/xzhang/OMI_HCHO/'
    ORB_PATH = '/YYYY/MM/'
    FILENAME_OMH = 'OMI-Aura_L2-OMHCHO_YYYYmMMDD'

    CALL EXPAND_DATE( ORB_PATH, YYYYMMDD, 9999 )
    CALL EXPAND_DATE( FILENAME_OMH, YYYYMMDD, 9999 )
    WRITE(I_CHAR,'(I2.2)') DAY

    !CALL SYSTEM("ls "//TRIM(ORB_PATH)//"OMH-Aura_L2-OMNO2_2016m08"//I_CHAR//"* > OMH_file_list"//I_CHAR//".txt")
    CALL SYSTEM("ls "//TRIM(DIR_OMH)//TRIM(ORB_PATH)//TRIM(FILENAME_OMH)//"* > omi_hcho_file_list"//I_CHAR//".txt")
    CLOSE(66) ! ugly...

    OPEN(66,FILE="omi_hcho_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")

    N_OMH_ORB = 0
    DIMS_COUNTER = 1
    DIMS_INDEX = 0
    DO ! loop over all available OMH NO2 files for the current day

       READ(66,'(A)',IOSTAT=IO_STATUS) FILE_OMH

       IF(IO_STATUS < 0) EXIT
       !FILE_ORB = TRIM(ORB_PATH) // FILE_ORB

       PRINT *,"Reading OMI HCHO file "//TRIM(FILE_OMH)
      !! open OMH ORB file

       CALL H5OPEN_F(error)

       CALL H5FOPEN_f (FILE_OMH, H5F_ACC_RDWR_F,file_id,error)

       !PRINT *,"OMH file status: ",error_orb

       ! Open an existing dataset.

       DSETNAME = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ScatteringWeights'

       CALL H5DOPEN_F(FILE_ID, DSETNAME, DSET_ID, ERROR)

       ! open dataspace

       CALL H5DGET_SPACE_F(DSET_ID, DSPACE_ID, ERROR)

       CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(DSPACE_ID, RANK_OMH,ERROR)

       CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DSPACE_ID, DIMS_OMH, MAXDIMS_OMH, ERROR)

       CALL H5DCLOSE_F(DSET_ID,ERROR)
       IF (ERROR < 0) DIMS_OMH(2) = 0

       N_OMH_ORB = N_OMH_ORB + DIMS_OMH(2)

       CALL H5FCLOSE_F(FILE_ID,ERROR)
       CALL H5CLOSE_F(ERROR)
       !PRINT *, "DIMS3", DIMS_ORB(3)
    ENDDO

    CLOSE(66)
    ALLOCATE(TIME_ORB(N_OMH_ORB))
    ALLOCATE(LON_ORB(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(LAT_ORB(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(AMF_TROP_ORB(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(CH2O_TROP(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(CH2O_TROP_STD(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(VIEW_ZEN(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(SOLAR_ZEN(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(CFR(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(X_Q_FLAG(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(M_Q_FLAG(N_OMH_SWATHS,N_OMH_ORB))
    ALLOCATE(SCW_PRE(N_OMH_SWATHS,N_OMH_ORB,N_OMH_LEVELS))
    ALLOCATE(SCW(N_OMH_SWATHS,N_OMH_ORB,N_OMH_LEVELS))
    CLOSE(76)
    OPEN(76,FILE="omi_hcho_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")
    DO
       READ(76,'(A)',IOSTAT=IO_ORB_STATUS) FILE_ORB

       IF(IO_ORB_STATUS < 0) EXIT
       DATA_VALID = .TRUE.
       !FILE_ORB = TRIM(ORB_PATH) // FILE_ORB

       !PRINT *,"Reading OMH file "//TRIM(FILE_ORB)

       !! open OMH ORB file

       CALL H5OPEN_F(error_orb)

       CALL H5FOPEN_f (FILE_ORB, H5F_ACC_RDWR_F,file_id_orb,error_orb)

       ! Open an existing dataset.

       DSETNAME = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ScatteringWeights'

       CALL H5DOPEN_F(FILE_ID_ORB, DSETNAME, DSET_ID_ORB, ERROR_ORB)
      ! open dataspace

       CALL h5dget_space_f(dset_id_orb, dspace_id_orb, error_orb)

       CALL h5sget_simple_extent_ndims_f(dspace_id_orb, rank_orb, error_orb)

       CALL h5sget_simple_extent_dims_f(dspace_id_orb, dims_orb, maxdims_orb, error_orb)

       CALL h5dclose_f(dset_id_orb,error_orb)

       IF ( REAL(DIMS_ORB(2)) == 0d0 ) THEN
          DATA_VALID = .FALSE.
       ELSEIF (REAL(DIMS_ORB(3)) .NE. REAL(N_OMH_LEVELS)) THEN
          DATA_VALID = .FALSE.
       ELSEIF (REAL(DIMS_ORB(1)) .NE. REAL(N_OMH_SWATHS)) THEN
          DATA_VALID = .FALSE.
       ELSEIF (ERROR_ORB < 0) THEN
          DATA_VALID = .FALSE.
       ENDIF

       IF (DATA_VALID == .FALSE.) THEN
          CALL H5FCLOSE_F(FILE_ID_ORB,ERROR_ORB)
          CALL H5CLOSE_F(ERROR_ORB)
          CYCLE
       ENDIF

       ALLOCATE(OMH_TIME(dims_orb(2)))
       ALLOCATE(OMH_LON(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_LAT(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_AMF_TROP(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_CH2O_TROP(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_CH2O_TROP_STD(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_VIEW_ZEN(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_SOLAR_ZEN(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_CFR(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_SCW_PRE(dims_orb(1),dims_orb(2),dims_orb(3)))
       ALLOCATE(OMH_SCW(dims_orb(1),dims_orb(2), dims_orb(3)) )
       ALLOCATE(OMH_M_Q_FLAG(dims_orb(1),dims_orb(2)))
       ALLOCATE(OMH_X_Q_FLAG(dims_orb(1),dims_orb(2)))

       DIMS_INDEX = DIMS_COUNTER+DIMS_ORB(2)-1
       !! read times
       !! open OMH ORB file

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Time'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_TIME,(/data_dims_orb(2),0/), error_orb)
       TIME_ORB(DIMS_COUNTER:DIMS_INDEX) = OMH_TIME
       CALL h5dclose_f(dset_id_orb,error_orb)

       !PRINT *,"Found matching OMH file for hour ", DAY, ",",GC_HOUR, ":", TRIM(FILE_ORB)
       !! read tropospheric air mass factors

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/AirMassFactor'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_AMF_TROP, data_dims_orb, error_orb)
       AMF_TROP_ORB(:,DIMS_COUNTER:DIMS_INDEX) = OMH_AMF_TROP
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read tropospheric CH2O column

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ReferenceSectorCorrectedVerticalColumn'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_CH2O_TROP, data_dims_orb(1:2), error_orb)
       CH2O_TROP(:,DIMS_COUNTER:DIMS_INDEX) = OMH_CH2O_TROP
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read tropospheric CH2O column

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ColumnUncertainty'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_CH2O_TROP_STD, data_dims_orb(1:2), error_orb)
       CH2O_TROP_STD(:,DIMS_COUNTER:DIMS_INDEX) = OMH_CH2O_TROP_STD
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read quality flag array

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/MainDataQualityFlag'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_INTEGER, OMH_M_Q_FLAG, data_dims_orb(2:3), error_orb)
       M_Q_FLAG(:,DIMS_COUNTER:DIMS_INDEX) = OMH_M_Q_FLAG
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read longitudes

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Longitude'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_LON, data_dims_orb(1:2), error_orb)
       LON_ORB(:,DIMS_COUNTER:DIMS_INDEX) = OMH_LON
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read latitudes

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/Latitude'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_LAT, data_dims_orb(1:2), error_orb)
       LAT_ORB(:,DIMS_COUNTER:DIMS_INDEX) = OMH_LAT
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read viewing zenith angles

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/ViewingZenithAngle'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_VIEW_ZEN, data_dims_orb(1:2), error_orb)
       VIEW_ZEN(:,DIMS_COUNTER:DIMS_INDEX) = OMH_VIEW_ZEN
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read solar zenith angles

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/SolarZenithAngle'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_SOLAR_ZEN, data_dims_orb(1:2), error_orb)
       SOLAR_ZEN(:,DIMS_COUNTER:DIMS_INDEX) = OMH_SOLAR_ZEN
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read cloud fraction

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/AMFCloudFraction'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_CFR, data_dims_orb(1:2), error_orb)
       CFR(:,DIMS_COUNTER:DIMS_INDEX) = OMH_CFR
       CALL h5dclose_f(dset_id_orb,error_orb)

       DSETNAME = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/XtrackQualityFlags'

       CALL H5DOPEN_F(FILE_ID_ORB, DSETNAME, DSET_ID_ORB, ERROR_ORB)

       CALL H5DREAD_F(DSET_ID_ORB, H5T_NATIVE_DOUBLE, OMH_X_Q_FLAG, DATA_DIMS_ORB(1:2), ERROR_ORB)
       X_Q_FLAG(:,DIMS_COUNTER:DIMS_INDEX) = OMH_X_Q_FLAG
       CALL H5DCLOSE_F(DSET_ID_ORB,ERROR_ORB)

       !! read scattering weight pressures

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ClimatologyLevels'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_SCW_PRE, data_dims_orb, error_orb)
       SCW_PRE(:,DIMS_COUNTER:DIMS_INDEX,:) = OMH_SCW_PRE
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! read scattering weights

       dsetname = '/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/ScatteringWeights'

       CALL h5dopen_f(file_id_orb, dsetname, dset_id_orb, error_orb)

       CALL h5dread_f(dset_id_orb, H5T_NATIVE_DOUBLE, OMH_SCW, data_dims_orb, error_orb)
       SCW(:,DIMS_COUNTER:DIMS_INDEX,:) = OMH_SCW
       CALL h5dclose_f(dset_id_orb,error_orb)

       !! close file

       CALL H5FCLOSE_F(file_id_orb,error_orb)

       DIMS_COUNTER = DIMS_COUNTER+DIMS_ORB(2)
       ! deallocate OMH arrays
       CALL H5CLOSE_F(ERROR_ORB)

       IF(ALLOCATED(OMH_LON)) DEALLOCATE(OMH_LON)
       IF(ALLOCATED(OMH_LAT)) DEALLOCATE(OMH_LAT)
       IF(ALLOCATED(OMH_TIME)) DEALLOCATE(OMH_TIME)
       IF(ALLOCATED(OMH_AMF_TROP)) DEALLOCATE(OMH_AMF_TROP)
       IF(ALLOCATED(OMH_CH2O_TROP)) DEALLOCATE(OMH_CH2O_TROP)
       IF(ALLOCATED(OMH_CH2O_TROP_STD)) DEALLOCATE(OMH_CH2O_TROP_STD)
       IF(ALLOCATED(OMH_VIEW_ZEN)) DEALLOCATE(OMH_VIEW_ZEN)
       IF(ALLOCATED(OMH_SOLAR_ZEN)) DEALLOCATE(OMH_SOLAR_ZEN)
       IF(ALLOCATED(OMH_CFR)) DEALLOCATE(OMH_CFR)
       IF(ALLOCATED(OMH_SCW_PRE)) DEALLOCATE(OMH_SCW_PRE)
       IF(ALLOCATED(OMH_X_Q_FLAG)) DEALLOCATE(OMH_X_Q_FLAG)
       IF(ALLOCATED(OMH_M_Q_FLAG)) DEALLOCATE(OMH_M_Q_FLAG) 
       IF(ALLOCATED(OMH_SCW)) DEALLOCATE(OMH_SCW)

    ENDDO
    CLOSE(76)

  END SUBROUTINE READ_OMI_CH2O_FILE

!=================================================================================================================
  SUBROUTINE CALC_OMI_CH2O_FORCE
    
    USE HDF5

    !!
    !! Subroutine CALC_OMI_CH2O_FORCE computes the O3 adjoint forcing from OMH column data
    !!

    USE COMODE_MOD,         ONLY : JLOP
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
    USE GRID_MOD,           ONLY : GET_IJ
    USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
    USE TIME_MOD,           ONLY : GET_HOUR, GET_DAY
    USE TIME_MOD,           ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR
    USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS, GET_TS_CHEM
    USE DAO_MOD,            ONLY : BXHEIGHT
    USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
    USE ADJ_ARRAYS_MOD,     ONLY : ID2C
    USE TRACERID_MOD,       ONLY : IDCH2O
    USE ADJ_ARRAYS_MOD,     ONLY : COST_FUNC
    
    USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
    USE TRACER_MOD,         ONlY : XNUMOLAIR 
    USE DAO_MOD,            ONLY : T, AIRDEN
    
#include "CMN_SIZE" ! size parameters
    
    INTEGER :: I,J,L
    INTEGER :: I_OMH, J_OMH, K_OMH, JLOOP
    INTEGER :: IIJJ(2)
    INTEGER :: NTSTART_OMH, NTSTOP_OMH
    
    !Arguments
    !INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS
    
    INTEGER :: DAY
    
    ! variables for time unit conversion
    REAL*8 :: tai93
    INTEGER :: iy,im,id,ih,imin
    REAL*8 :: sec
    INTEGER :: GC_HOUR, MIN_HOUR, MAX_HOUR
    
    ! variables for observation operator and adjoint thereof
    
    REAL*8 :: OMI_CH2O_GC(IIPAR,JJPAR)
    REAL*8 :: SCW_GC(LLPAR), DP(LLPAR)
    REAL*8 :: AMF_GC
    REAL*8 :: GC_CH2O(LLPAR)
    REAL*8 :: GC_CH2O_COL
    REAL*8 :: DIFF, FORCE_COL, COST_CONTRIB_COL
    REAL*8 :: OBS_ERROR
    REAL*8 :: MEAN_DIFF(IIPAR,JJPAR)
    REAL*8, SAVE :: OMH_HOUR(MAX_ORB)
    REAL*8 :: OLD_COST_OMH

    ! arrays needed for superobservations
    LOGICAL :: SUPER_OBS = .TRUE.                       ! do super observations?
    REAL*8 ::  SOBS_COUNT(IIPAR,JJPAR)                   ! super observation count
    REAL*8 ::  SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)         ! super observation adjoint forcing
    REAL*8 ::  GC_ADJ_COUNT(IIPAR,JJPAR,LLPAR)
    REAL*8 ::  NEW_COST(IIPAR,JJPAR)       ! super observation cost function contribution
    REAL*8 ::  SOBS_GC(IIPAR,JJPAR)
    REAL*8 ::  SOBS_OMH(IIPAR,JJPAR)
    REAL*8 ::  SOBS_BIAS(IIPAR,JJPAR)
    REAL*8 ::  SOBS_CHISQUARED(IIPAR,JJPAR)
    
    !=============================================================
    ! CALC_OMI_HCHO_FORCE begins here!
    !=============================================================
    
    PRINT *,"ID2C:",ID2C(IDCH2O)
    
    GC_HOUR = GET_HOUR()
    
    ! initialize arrays
    
    GC_CH2O = 0d0
    GC_CH2O_COL = 0d0
    MEAN_DIFF = 0d0
    OLD_COST_OMH = COST_FUNC
    OMI_CH2O_GC = 0d0
    GC_ADJ_COUNT = 0d0
    SOBS_COUNT = 0d0
    SOBS_ADJ_FORCE = 0d0
    NEW_COST = 0d0
    SOBS_GC = 0d0
    SOBS_OMH = 0d0
    SOBS_BIAS = 0d0
    SOBS_CHISQUARED = 0d0
    
    DAY = GET_DAY()
    GC_HOUR = GET_HOUR()

    IF ( GET_NHMS() == 236000 - GET_TS_CHEM()* 100 ) THEN
       CALL READ_OMI_CH2O_FILE(GET_NYMD(), N_OMH_ORB)!GET_NHMS())
       DO I_OMH = 1, N_OMH_ORB
          IF (TIME_ORB(I_OMH) > 0) THEN
             CALL TAI2UTC(TIME_ORB(I_OMH),IY,IM,ID,IH,IMIN,SEC)
             IF (ID == DAY) THEN
                OMH_HOUR(I_OMH) = IH
             ELSE
                OMH_HOUR(I_OMH) = -999
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    !! loop over data
    CALL GET_NT_RANGE_OMH(N_OMH_ORB, GET_NHMS(), OMH_HOUR(1:N_OMH_ORB), NTSTART_OMH, NTSTOP_OMH)
    IF ( NTSTART_OMH == 0 .and. NTSTOP_OMH == 0 ) THEN

       print*, ' No matching OMI HCHO obs for this hour'
       RETURN
    ENDIF
    PRINT *, 'found record range:', NTSTART_OMH, NTSTOP_OMH
    !PRINT *, "TIME_FRAC", TIME_FRAC(1:N_OMH_ORB)
    DO I_OMH=NTSTART_OMH,NTSTOP_OMH,-1
       DO J_OMH=1,N_OMH_SWATHS

                ! A number of conditions have to be met for OMH CH2O data to actually be assimilated
          IF ( ( TIME_ORB(I_OMH) > 0                  ) .AND. &
               ( REAL(M_Q_FLAG(J_OMH,I_OMH)) < 1d0      ) .AND. &
               ( REAL(X_Q_FLAG(J_OMH,I_OMH)) < 1d0      ) .AND. &
               ( ABS(LAT_ORB(J_OMH,I_OMH)) < 60d0     ) .AND. &
               ( CH2O_TROP(J_OMH,I_OMH) > 0d0           ) .AND. &
               ( ABS(SOLAR_ZEN(J_OMH, I_OMH)) < 75d0 ) .AND. &
               ( ABS(VIEW_ZEN(J_OMH,I_OMH)) < 65d0   ) .AND. &
               ( AMF_TROP_ORB(J_OMH,I_OMH) > 0d0      ) .AND. &
               ( CFR(J_OMH,I_OMH) < 0.4   ) ) THEN                
             ! Get model grid coordinate indices that correspond to the observation                     
             
             IIJJ = GET_IJ(REAL(LON_ORB(J_OMH,I_OMH),4), REAL(LAT_ORB(J_OMH,I_OMH),4))
             
             I = IIJJ(1)
             J = IIJJ(2)
             
             ! Get GEOS-CHEM CH2O values [#/cm3]
             
             GC_CH2O = 0d0
             GC_CH2O_COL = 0d0
             SCW_GC = 0d0
             DP = 0d0
             COST_CONTRIB_COL = 0d0
             FORCE_COL = 0d0

             DO L = 1, LLPAR
                
                IF( ITS_IN_THE_TROP(I,J,L) ) THEN
                   
                   JLOOP = JLOP(I,J,L) 
                   
                   GC_CH2O(L) = CSPEC_AFTER_CHEM(JLOOP,ID2C(IDCH2O))   
                ENDIF
                
             ENDDO
             
             ! Compute tropospheric CH2O column value [#/cm2]
             
             GC_CH2O_COL = SUM(GC_CH2O(:) * BXHEIGHT(I,J,:) * 100d0)
             
             ! interpolate scattering weights to GEOS-Chem grid
             
             DO L=1,LLPAR
                DO K_OMH = 2,N_OMH_LEVELS
                   
                   IF( GET_PCENTER(I,J,L) < SCW_PRE(J_OMH,I_OMH,K_OMH-1) .AND. GET_PCENTER(I,J,L) > SCW_PRE(J_OMH,I_OMH,K_OMH) ) THEN
                      
                      ! linearly interpolate scattering weights to GEOS-Chem center pressures
                      
                      SCW_GC(L) = SCW(J_OMH,I_OMH,K_OMH) + &
                           ( SCW(J_OMH,I_OMH,K_OMH-1) - SCW(J_OMH,I_OMH,K_OMH) ) * &
                           ( GET_PCENTER(I,J,L) - SCW_PRE(J_OMH,I_OMH,K_OMH) )/( SCW_PRE(J_OMH,I_OMH,K_OMH-1) - SCW_PRE(J_OMH,I_OMH,K_OMH) )
                      
                      ! save pressure differences
                      
                      DP(L) = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                      
                            !! convert CH2O concentrations from number density to mixing ratio
                      
                      GC_CH2O(L) = GC_CH2O(L) *1d6 / ( AIRDEN(L,I,J)   * XNUMOLAIR )
                      
                      !EXIT
                      
                   ENDIF
                   
                ENDDO
             ENDDO
             
             ! Use tropospheric air mass factors to convert vertical column to slant column
             !PRINT *, "SUM1", SUM(GC_CH2O * DP)
             !PRINT *, "SUM2", SUM(DP)
             
             AMF_GC = SUM(GC_CH2O * DP * SCW_GC)/SUM(GC_CH2O * DP)
             !PRINT *, "AMF_GC", AMF_GC
             GC_CH2O_COL = AMF_GC*GC_CH2O_COL
             
             ! The computation above is a little awkward, since the slant column can be computed directly from equation (2) in Bucsela2013 without 
             ! computing the airmass factors and CH2O column first.
             ! I chose to compute the slant column from the computed air mass factors (which already included the computation of the slant column)
             ! since the air mass factors might be of diagnostic interest (i.e. some reviewer might want to see them) and should be computed and saved
             ! alongside other observation operator diagnostics. Furthermore, this formulation makes the adjoint of the observation operator somewhat simpler to handle.
             
             DIFF = GC_CH2O_COL - CH2O_TROP(J_OMH,I_OMH) * AMF_TROP_ORB(J_OMH, I_OMH)
             
             !MEAN_DIFF(I,J) = MEAN_DIFF(I,J) + DIFF
             !PRINT *, "CHECK"
             OBS_ERROR = CH2O_TROP_STD(J_OMH,I_OMH) * AMF_TROP_ORB(J_OMH,I_OMH)
             IF (OBS_ERROR > 0d0) THEN
                FORCE_COL = DIFF/(OBS_ERROR**2)
                COST_CONTRIB_COL = 0.5d0 * DIFF * FORCE_COL
             ELSE
                FORCE_COL = 0d0
                COST_CONTRIB_COL = 0d0
             ENDIF
             IF ( ( COST_CONTRIB_COL > 0d0 ) .AND. &
                  ( COST_CONTRIB_COL <= 200d0 ) ) THEN   
                DO L=1,LLPAR
                   
                   IF(ITS_IN_THE_TROP(I,J,L)) THEN
                      IF (SUPER_OBS) THEN
                         SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + FORCE_COL * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                         GC_ADJ_COUNT(I,J,L) = GC_ADJ_COUNT(I,J,L) + 1
                      ELSE
                         JLOOP = JLOP(I,J,L)
                         CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDCH2O)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDCH2O)) +  &
                              FORCE_COL * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                      ENDIF
                   ENDIF
                   
                ENDDO
                      
                ! update cost function
                IF(SUPER_OBS) THEN
                   NEW_COST(I,J) = NEW_COST(I,J) + COST_CONTRIB_COL
                   SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                ELSE
                   COST_FUNC = COST_FUNC + COST_CONTRIB_COL
                   !PRINT *, "OBS_COST", COST_FUNC
                ENDIF
             ENDIF
             ! update dignostic arrays
             
             IF(SUPER_OBS) THEN
                SOBS_GC(I,J) = SOBS_GC(I,J) + GC_CH2O_COL
                SOBS_OMH(I,J) = SOBS_OMH(I,J) + CH2O_TROP(J_OMH,I_OMH)
                SOBS_BIAS(I,J) = SOBS_BIAS(I,J) + DIFF
                SOBS_CHISQUARED(I,J) = SOBS_CHISQUARED(I,J) + 0.5 * (DIFF/OBS_ERROR)**2
             ELSE
                ! calculate OMH bias and variance using Knuth's online algorithm
                OMH_BIAS_COUNT(I,J) = OMH_BIAS_COUNT(I,J) + 1d0
                OMH_CH2O_MEAN(I,J) = OMH_CH2O_MEAN(I,J) + (AMF_TROP_ORB(J_OMH,I_OMH)*CH2O_TROP(J_OMH,I_OMH))
                OMH_GEOS_CH2O_MEAN(I,J) = OMH_GEOS_CH2O_MEAN(I,J) + GC_CH2O_COL
                OMH_CH2O_ERR_MEAN(I,J) = OMH_CH2O_ERR_MEAN(I,J) + OBS_ERROR
                OMH_DELTA = DIFF - OMH_BIAS(I,J)
                OMH_BIAS(I,J) = OMH_BIAS(I,J) + OMH_DELTA/OMH_BIAS_COUNT(I,J) 
                OMH_VAR(I,J) = OMH_VAR(I,J) + OMH_DELTA*(DIFF-OMH_BIAS(I,J))
                OMH_CHI_SQUARED(I,J) = OMH_CHI_SQUARED(I,J) + ( DIFF/OBS_ERROR )**2
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !PRINT *, "SOBS_ADJ_FORCE", SOBS_ADJ_FORCE
    !PRINT *, "CSPEC_AFTER_CHEM", CSPEC_AFTER_CHEM_ADJ
    !PRINT *, "SOBS_COST_CONTRIBUTION", SOBS_COST_CONTRIBUTION
       
    IF(SUPER_OBS) THEN
       
       DO J=1,JJPAR
          DO I=1,IIPAR
             
             IF(SOBS_COUNT(I,J) > 0d0) THEN
                
                DO L=1,LLPAR
                   JLOOP = JLOP(I,J,L)
                   IF( ( ITS_IN_THE_TROP(I,J,L) ) .AND. &
                        ( GC_ADJ_COUNT(I,J,L) > 0 ) .AND. &
                        (JLOOP > 0) ) THEN
                      CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDCH2O)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDCH2O)) +  SOBS_ADJ_FORCE(I,J,L)/GC_ADJ_COUNT(I,J,L)
                      
                   ENDIF
                   
                ENDDO
                
                COST_FUNC = COST_FUNC + NEW_COST(I,J)/SOBS_COUNT(I,J)
                
                OMH_BIAS_COUNT(I,J) = OMH_BIAS_COUNT(I,J) + 1d0
                
                OMH_CH2O_MEAN(I,J) = OMH_CH2O_MEAN(I,J) + SOBS_OMH(I,J)/SOBS_COUNT(I,J)
                
                OMH_GEOS_CH2O_MEAN(I,J) = OMH_GEOS_CH2O_MEAN(I,J) + SOBS_GC(I,J)/SOBS_COUNT(I,J)
                
                OMH_CH2O_ERR_MEAN(I,J) = OMH_CH2O_ERR_MEAN(I,J) + OBS_ERROR  
                
                ! calculate bias and variance of GC-OMH bias using numerically stable one-pass algorithm (Chan83)
                
                OMH_DELTA = SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMH_BIAS(I,J)
                      
                OMH_BIAS(I,J) = OMH_BIAS(I,J) + OMH_DELTA/OMH_BIAS_COUNT(I,J)
                
                OMH_VAR(I,J) = OMH_VAR(I,J) + OMH_DELTA*(SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMH_BIAS(I,J))
                
                OMH_CHI_SQUARED(I,J) = OMH_CHI_SQUARED(I,J)  + SOBS_CHISQUARED(I,J)/SOBS_COUNT(I,J)
                
             ENDIF
             
          ENDDO
       ENDDO
       
    ENDIF
    
    PRINT *, "COST FUNCTION OF OMI HCHO", COST_FUNC-OLD_COST_OMH
      
  END SUBROUTINE CALC_OMI_CH2O_FORCE

!-----------------------------------------------------------------------------!
  SUBROUTINE CLEANUP_OMH
    ! deallocate OMH arrays

    IF(ALLOCATED(LON_ORB)) DEALLOCATE(LON_ORB)
    IF(ALLOCATED(LAT_ORB)) DEALLOCATE(LAT_ORB)
    IF(ALLOCATED(TIME_ORB)) DEALLOCATE(TIME_ORB)
    IF(ALLOCATED(AMF_TROP_ORB)) DEALLOCATE(AMF_TROP_ORB)
    IF(ALLOCATED(CH2O_TROP)) DEALLOCATE(CH2O_TROP)
    IF(ALLOCATED(CH2O_TROP_STD)) DEALLOCATE(CH2O_TROP_STD)
    IF(ALLOCATED(VIEW_ZEN)) DEALLOCATE(VIEW_ZEN)
    IF(ALLOCATED(SOLAR_ZEN)) DEALLOCATE(SOLAR_ZEN)
    IF(ALLOCATED(CFR)) DEALLOCATE(CFR)
    IF(ALLOCATED(SCW_PRE)) DEALLOCATE(SCW_PRE)
    IF(ALLOCATED(X_Q_FLAG)) DEALLOCATE(X_Q_FLAG)
    IF(ALLOCATED(M_Q_FLAG)) DEALLOCATE(M_Q_FLAG)
    IF(ALLOCATED(SCW)) DEALLOCATE(SCW)

  END SUBROUTINE CLEANUP_OMH
!----------------------------------------------------------------------------------
  SUBROUTINE GET_NT_RANGE_OMH( N_OMH_ORB, HHMMSS, OMH_HOUR, NTSTART_OMH, NTSTOP_OMH)
!
!******************************************************************************
!  Subroutine GET_NT_RANGE_OMH retuns the range of retrieval records for the
!  current model hour
!
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTES   (INTEGER) : Number of TES retrievals in this day
!  (2 ) HHMMSS (INTEGER) : Current model time
!  (3 ) TIME_FRAC (REAL) : Vector of times (frac-of-day) for the TES retrievals
!=====================================================================================
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART_OMH (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP_OMH  (INTEGER) : TES record number at which to stop
!
!  NOTES:
!
!******************************************************************************
!
    ! Reference to f90 modules
    USE ERROR_MOD,    ONLY : ERROR_STOP
    USE TIME_MOD,     ONLY : YMD_EXTRACT

    ! Arguments
    INTEGER, INTENT(IN)   :: N_OMH_ORB
    INTEGER, INTENT(IN)   :: HHMMSS
    REAL*8,  INTENT(IN)   :: OMH_HOUR(N_OMH_ORB)
    INTEGER, INTENT(OUT)  :: NTSTART_OMH
    INTEGER, INTENT(OUT)  :: NTSTOP_OMH

    ! Local variables
    INTEGER, SAVE         :: NTSAVE_OMH
    LOGICAL               :: FOUND_ALL_RECORDS
    INTEGER               :: NTEST_OMH
    INTEGER               :: HH, MM, SS
    REAL*8                :: GC_HH_FRAC
    REAL*8                :: H1_FRAC

    !=================================================================
    ! GET_NT_RANGE_OMH begins here!
    !=================================================================
    ! Initialize
    FOUND_ALL_RECORDS  = .FALSE.
    NTSTART_OMH            = 0
    NTSTOP_OMH             = 0

    ! set NTSAVE_OMH to NTES every time we start with a new file
    IF ( HHMMSS == 230000 ) NTSAVE_OMH = N_OMH_ORB
    !print*, ' GET_NT_RANGE_OMH for ', HHMMSS
    !print*, ' NTSAVE_OMH ', NTSAVE_OMH
    !print*, ' N_IASI_NOB   ', N_IASI_NOB
    DO WHILE (OMH_HOUR(NTSAVE_OMH) < 0 )
       NTSAVE_OMH = NTSAVE_OMH - 1
       IF (NTSAVE_OMH == 0) EXIT
    ENDDO

    !PRINT *, "TIME_FRAC", TIME_FRAC(NTSAVE_OMH-1000:NTSAVE_OMH)
    CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


    ! Convert HH from hour to fraction of day
    GC_HH_FRAC = REAL(HH,8)
    ! one hour as a fraction of day
    H1_FRAC    = 0d0


    ! All records have been read already
    IF ( NTSAVE_OMH == 0 ) THEN

       print*, 'All records have been read already '
       RETURN

       ! No records reached yet
    ELSEIF ( OMH_HOUR(NTSAVE_OMH) + H1_FRAC < GC_HH_FRAC ) THEN


       print*, 'No records reached yet'
       RETURN

       !
    ELSEIF ( OMH_HOUR(NTSAVE_OMH) + H1_FRAC >=  GC_HH_FRAC ) THEN
       ! Starting record found
       NTSTART_OMH = NTSAVE_OMH

       !print*, ' Starting : TIME_FRAC(NTSTART_OMH) ', TIME_FRAC(NTSTART_OMH), NTSTART_OMH

       ! Now search forward to find stopping record
       NTEST_OMH = NTSTART_OMH

       DO WHILE ( FOUND_ALL_RECORDS == .FALSE. )

          ! Advance to the next record
          NTEST_OMH = NTEST_OMH - 1

          ! Stop if we reach the earliest available record
          IF ( NTEST_OMH == 0 ) THEN

             NTSTOP_OMH            = NTEST_OMH + 1
             FOUND_ALL_RECORDS = .TRUE.

             !print*, ' Records found '
             !print*, ' NTSTART_OMH, NTSTOP_OMH = ', NTSTART_OMH, NTSTOP_OMH
             ! Reset NTSAVE_OMH
             NTSAVE_OMH = NTEST_OMH
             ! When the combined test date rounded up to the nearest
             ! half hour is smaller than the current model date, the
             ! stopping record has been passed.
        ELSEIF (  OMH_HOUR(NTEST_OMH) + H1_FRAC <  GC_HH_FRAC ) THEN

             !print*, ' Testing : TIME_FRAC ', TIME_FRAC(NTEST_OMH), NTEST_OMH

             NTSTOP_OMH            = NTEST_OMH + 1
             FOUND_ALL_RECORDS = .TRUE.

             !print*, ' Records found '
             !print*, ' NTSTART_OMH, NTSTOP_OMH = ', NTSTART_OMH, NTSTOP_OMH

             ! Reset NTSAVE_OMH
             NTSAVE_OMH = NTEST_OMH
             !ELSE
             !print*, ' still looking ', NTEST_OMH

          ENDIF

       ENDDO

    ELSE

       CALL ERROR_STOP('problem', 'GET_NT_RANGE_OMH' )

    ENDIF

    ! Return to calling program
  END SUBROUTINE GET_NT_RANGE_OMH
!================================================================================================================================

  SUBROUTINE TAI2UTC(tai93,iy,im,id,ih,imin,sec)
      
    !!
    !! SUBROUTINE TAI2UTC converts TAI93 time (seconds since 1.1.1993) to UTC 
    !!
    !! adapted from 
    !! code.google.com/p/miyoshi/source/browse/trunk/common/common.f90

    IMPLICIT NONE
      
    INTEGER,PARAMETER :: n=10  ! number of leap seconds after Jan. 1, 1993
    INTEGER,PARAMETER :: leapsec(n) = (/ 15638399, 47174400, 94608001, 141868802, 189302403, 410227204, 504921605, 615254406, 709862407, 757382408/)
    REAL*8,INTENT(IN) :: tai93
    INTEGER,INTENT(OUT) :: iy,im,id,ih,imin
    REAL*8,INTENT(OUT) :: sec
    REAL*8,PARAMETER :: mins = 60.0d0
    REAL*8,PARAMETER :: hour = 60.0d0*mins
    REAL*8,PARAMETER :: day = 24.0d0*hour
    REAL*8,PARAMETER :: year = 365.0d0*day
    INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    REAL*8 :: wk,tai
    INTEGER :: days,i,leap
       
    tai = tai93
    sec = 0.0d0
       
    DO i=1,n
       
       IF(FLOOR(tai93) == leapsec(i)+1) THEN
          sec = 60.0d0 + tai93-FLOOR(tai93)
       ENDIF

       IF(FLOOR(tai93) > leapsec(i)) tai = tai -1.0d0
       
    END DO

    iy = 1993 + FLOOR(tai /year)
    wk = tai - REAL(iy-1993)*year - FLOOR(REAL(iy-1993)/4.0)*day
    
    IF(wk < 0.0d0) THEN
       iy = iy -1
       wk = tai - REAL(iy-1993)*year - FLOOR(REAL(iy-1993)/4.0)*day
    END IF

    days = FLOOR(wk/day)
    wk = wk - REAL(days)*day
    im = 1
      
    DO i=1,12
       
       leap = 0
       IF(im == 2 .AND. MOD(iy,4)==0) leap=1
       IF(im == i .AND. days >= mdays(i)+leap) THEN
          im = im + 1
          days = days - mdays(i)-leap
       END IF
       
    END DO
    
    id = days +1
    
    ih = FLOOR(wk/hour)
    wk = wk - REAL(ih)*hour
    imin = FLOOR(wk/mins)

    IF(sec < 60.0d0) sec = wk - REAL(imin)*mins
    
    RETURN
       
  END SUBROUTINE TAI2UTC
  !------------------------------------------------------------------------------    
  SUBROUTINE OMH_HANDLE_ERR( RETVAL )

    INTEGER, INTENT(IN) :: RETVAL

    PRINT *,"AN ERROR OCCURED WHILE WRITING OUT OMH CH2O DIAGNOSTICS!"
    PRINT *,"NETCDF ERROR MESSAGE:"
    
  END SUBROUTINE OMH_HANDLE_ERR
      
  !--------------------------------------------------------------------------------
  
END MODULE OMI_CH2O_OBS_MOD
