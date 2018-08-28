MODULE OMI_NO2_OBS_MOD

  !!
  ! Module OMI_NO2_OBS contains all subroutines and variables needed to assimilate OMI NO2 tropospheric column data
  !
  ! Module Routines:
  !
  ! (1) CALC_OMI_NO2_FORCE      : calculates adjoint forcing and cost function contribution for OMI tropospheric NO2 columns
  ! (2) TAI2UTC                 : converts TAI93 (seconds since 1.1.1993) to UTC
  ! (3) MAKE_OMI_BIAS_FILE_HDF5 : writes OMI satellite diagnostics in satellite diagnostic HDF5 file
  !

  IMPLICIT NONE

#include "CMN_SIZE"

  PRIVATE

  PUBLIC READ_OMI_NO2_FILE
  PUBLIC CALC_OMI_NO2_FORCE

  ! Module variables

  ! Arrays for diagnostic output
  ! OMI data
  ! Parameters
  REAL*8, ALLOCATABLE :: OMI_LON(:,:), OMI_LAT(:,:)
  REAL*8, ALLOCATABLE :: OMI_TIME(:), OMI_AMF_TROP(:,:)
  REAL*8, ALLOCATABLE :: OMI_NO2_TROP(:,:), OMI_NO2_TROP_STD(:,:)
  REAL*8, ALLOCATABLE :: OMI_VIEW_ZENITH(:,:), OMI_SOLAR_ZENITH(:,:)
  REAL*8, ALLOCATABLE :: OMI_CLOUDFR(:,:), OMI_SCW_P(:)
  REAL*8, ALLOCATABLE :: OMI_SCATTERING_WEIGHTS(:,:,:)
  REAL*8, ALLOCATABLE :: OMI_X_QUAL_FLAG(:,:), OMI_M_QUAL_FLAG(:)
  REAL*8, ALLOCATABLE :: OMI_TROPO_PRESSURE(:,:)
  REAL*8, ALLOCATABLE :: LON_ORBIT(:,:), LAT_ORBIT(:,:)
  REAL*8, ALLOCATABLE :: TIME_ORBIT(:), AMF_TROP_ORBIT(:,:)
  REAL*8, ALLOCATABLE :: NO2_TROP(:,:), NO2_TROP_STD(:,:)
  REAL*8, ALLOCATABLE :: VIEW_ZENITH(:,:), SOLAR_ZENITH(:,:)
  REAL*8, ALLOCATABLE :: CLOUDFR(:,:), SCW_P(:)
  REAL*8, ALLOCATABLE :: SCATTERING_WEIGHTS(:,:,:)
  REAL*8, ALLOCATABLE :: X_QUAL_FLAG(:,:), M_QUAL_FLAG(:)
  REAL*8, ALLOCATABLE :: TROPO_PRESSURE(:,:)
  
  INTEGER :: N_OMI_ORBITS
  INTEGER, PARAMETER :: MAX_ORBITS = 50000
  INTEGER, PARAMETER :: N_OMI_SWATHS = 60
  INTEGER, PARAMETER :: N_OMI_LEVELS = 35

  ! arrays to store diagnostic information
  REAL*4:: OMI_NO2_MEAN(IIPAR,JJPAR) = 0d0      ! Mean OMI columns
  REAL*4:: OMI_GEOS_NO2_MEAN(IIPAR,JJPAR) = 0d0 ! Mean GEOS-Chem columns
  REAL*4:: OMI_NO2_ERR_MEAN(IIPAR,JJPAR) = 0d0  ! Mean OMI observation errors
  REAL*4:: OMI_BIAS(IIPAR,JJPAR)=0d0            ! Model biases
  REAL*4:: OMI_VAR(IIPAR,JJPAR)=0d0             ! Model variances
  REAL*4:: OMI_DELTA=0d0                        ! temporary storage variable
  REAL*4:: OMI_BIAS_COUNT(IIPAR,JJPAR) = 0d0    ! counter for number of observations in grid box
  REAL*4:: OMI_CHISQUARED(IIPAR,JJPAR) = 0d0   ! Chi-squared values
  LOGICAL :: FIRST = .TRUE.

CONTAINS

  !-----------------------------------------------------------------------------!

  SUBROUTINE READ_OMI_NO2_FILE ( YYYYMMDD, N_OMI_ORBITS )
    
    USE ERROR_MOD, ONLY : ALLOC_ERR
    USE TIME_MOD,  ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR
    USE HDF5
    USE TIME_MOD,  ONLY : GET_HOUR, GET_DAY
    ! Arguments
    INTEGER, INTENT(IN) :: YYYYMMDD!, HHMMSS
    
    CHARACTER(LEN=255) :: DIR_OMI
    CHARACTER(LEN=255) :: FILENAME_OMI
    CHARACTER(LEN=255) :: FILENAME_FULL
    CHARACTER(LEN=255) :: DSET_NAME
    CHARACTER(255) :: ORBIT_PATH,FILE_ORBIT, FILE_ORBIT2
    CHARACTER(2) :: I_CHAR
    INTEGER :: IO_ORBIT_STATUS, IO_ORBIT_STATUS2
    INTEGER :: DAY
    INTEGER(HID_T) :: file_id, dset_id, file_id2, dset_id2
    
    !INTEGER(HSIZE_T) :: data_dims
    INTEGER :: error, DIMS_COUNTER, DIMS_INDEX
    INTEGER :: N_OMI_ORBITS
    INTEGER(HID_T) :: file_id_orbit, file_id_orbit2
    INTEGER(HID_T) :: dset_id_orbit, dset_id_orbit2
    INTEGER(HID_T) :: dspace_id_orbit, dspace_id_orbit2
    INTEGER :: error_orbit, error_orbit2
    CHARACTER(LEN=255) :: filename_orbit, dsetname, filename_orbit2, dsetname2

    INTEGER :: rank_orbit, rank_orbit2
    INTEGER(HSIZE_T) :: dims_orbit(3), maxdims_orbit(3), dims_count(3), maxdims_orbit2(3)
    INTEGER(HSIZE_T) :: DATA_DIMS_ORBIT(3)
    INTEGER :: GC_HOUR
    LOGICAL          :: DATA_VALID

    CALL CLEANUP_OMI

    !PRINT *, "ID2C(IDNO2)", ID2C(IDNO2)
    GC_HOUR = GET_HOUR()
    DAY = GET_DAY()
    DIR_OMI = '/users/jk/15/xzhang/OMI_NO2/'
    ORBIT_PATH = '/YYYY/MM/'
    FILENAME_OMI = 'OMI-Aura_L2-OMNO2_YYYYmMMDD'

    CALL EXPAND_DATE( ORBIT_PATH, YYYYMMDD, 9999 )
    CALL EXPAND_DATE( FILENAME_OMI, YYYYMMDD, 9999 )
    WRITE(I_CHAR,'(I2.2)') DAY
    !WRITE(I_CHAR2, '(I4.4)') YEAR
    !WRITE(I_CHAR3, '(I2.2)') 
    !PRINT *, "PATH TO FILE", TRIM(ORBIT_PATH)//TRIM(FILENAME_OMI)
    !CALL SYSTEM("ls "//TRIM(ORBIT_PATH)//"OMI-Aura_L2-OMNO2_2016m08"//I_CHAR//"* > omi_file_list"//I_CHAR//".txt")
    CALL SYSTEM("ls "//TRIM(DIR_OMI)//TRIM(ORBIT_PATH)//TRIM(FILENAME_OMI)//"* > omi_file_list"//I_CHAR//".txt")
    CLOSE(65) ! ugly...

    OPEN(65,FILE="omi_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")

    N_OMI_ORBITS = 0
    DIMS_COUNTER = 1
    DIMS_INDEX = 0
    DO ! loop over all available OMI NO2 files for the current day

       READ(65,'(A)',IOSTAT=IO_ORBIT_STATUS2) FILE_ORBIT2

       IF(IO_ORBIT_STATUS2 < 0) EXIT
       !FILE_ORBIT = TRIM(ORBIT_PATH) // FILE_ORBIT

       PRINT *,"Reading OMI file "//TRIM(FILE_ORBIT2)

       !! open OMI orbit file

       CALL H5OPEN_F(error_orbit2)

       CALL H5FOPEN_f (FILE_ORBIT2, H5F_ACC_RDWR_F,file_id_orbit2,error_orbit2)

       !PRINT *,"OMI file status: ",error_orbit

       ! Open an existing dataset.

       DSETNAME2 = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight'

       CALL H5DOPEN_F(FILE_ID_ORBIT2, DSETNAME2, DSET_ID_ORBIT2, ERROR_ORBIT2)

       ! open dataspace

       CALL H5DGET_SPACE_F(DSET_ID_ORBIT2, DSPACE_ID_ORBIT2, ERROR_ORBIT2)

       CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(DSPACE_ID_ORBIT2, RANK_ORBIT2,ERROR_ORBIT2)

       CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DSPACE_ID_ORBIT2, DIMS_COUNT, MAXDIMS_ORBIT2, ERROR_ORBIT2)

       CALL H5DCLOSE_F(DSET_ID_ORBIT2,ERROR_ORBIT2)
       IF (ERROR_ORBIT2 < 0) THEN
          DIMS_COUNT(3) = 0
       ENDIF
       N_OMI_ORBITS = N_OMI_ORBITS + DIMS_COUNT(3)
       CALL H5FCLOSE_F(FILE_ID_ORBIT2,ERROR_ORBIT2)
       CALL H5CLOSE_F(ERROR_ORBIT2)
       !PRINT *, "DIMS3", DIMS_ORBIT(3)
    ENDDO
    !PRINT *, "DATA_DIMS", N_OMI_ORBITS
    CLOSE(65)
    !PRINT *, "N_OMI_ORBITS TOTAL", N_OMI_ORBITS
    ALLOCATE(TIME_ORBIT(N_OMI_ORBITS))
    ALLOCATE(LON_ORBIT(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(LAT_ORBIT(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(AMF_TROP_ORBIT(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(NO2_TROP(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(NO2_TROP_STD(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(VIEW_ZENITH(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(SOLAR_ZENITH(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(CLOUDFR(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(X_QUAL_FLAG(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(M_QUAL_FLAG(N_OMI_ORBITS))
    ALLOCATE(SCW_P(N_OMI_LEVELS))
    ALLOCATE(TROPO_PRESSURE(N_OMI_SWATHS,N_OMI_ORBITS))
    ALLOCATE(SCATTERING_WEIGHTS(N_OMI_LEVELS,N_OMI_SWATHS,N_OMI_ORBITS))
    CLOSE(75)
    OPEN(75,FILE="omi_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")
    DO
       READ(75,'(A)',IOSTAT=IO_ORBIT_STATUS) FILE_ORBIT

       IF(IO_ORBIT_STATUS < 0) EXIT
       DATA_VALID = .TRUE.
       !FILE_ORBIT = TRIM(ORBIT_PATH) // FILE_ORBIT

       !PRINT *,"Reading OMI file "//TRIM(FILE_ORBIT)

       !! open OMI orbit file

       CALL H5OPEN_F(error_orbit)

       CALL H5FOPEN_f (FILE_ORBIT, H5F_ACC_RDWR_F,file_id_orbit,error_orbit)

       ! Open an existing dataset.

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       ! open dataspace

       CALL H5DGET_SPACE_F(DSET_ID_ORBIT, DSPACE_ID_ORBIT, ERROR_ORBIT)

       CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(DSPACE_ID_ORBIT, RANK_ORBIT,ERROR_ORBIT)

       CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DSPACE_ID_ORBIT, DIMS_ORBIT, MAXDIMS_ORBIT, ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       IF ( REAL(DIMS_ORBIT(3)) == 0d0 ) THEN
          DATA_VALID = .FALSE.
       ELSEIF ( REAL(DIMS_ORBIT(1)) .NE. REAL(N_OMI_LEVELS)) THEN
          DATA_VALID = .FALSE.
       ELSEIF ( REAL(DIMS_ORBIT(2)) .NE. REAL(N_OMI_SWATHS)) THEN
          DATA_VALID = .FALSE.
       ELSEIF ( ERROR_ORBIT < 0 ) THEN
          DATA_VALID = .FALSE.
       ENDIF

       IF (DATA_VALID == .FALSE.) THEN
          CALL H5FCLOSE_F(FILE_ID_ORBIT,ERROR_ORBIT)
          CALL H5CLOSE_F(ERROR_ORBIT)
          CYCLE
       ENDIF
       !PRINT *,"Found matching OMI file for hour ", DAY, ",",GC_HOUR, ":", TRIM(FILE_ORBIT)
       !PRINT *,"Found matching OMI file for hour ", DAY, ",",GC_HOUR, ":", TRIM(FILE_ORBIT)

       ALLOCATE(OMI_TIME(DIMS_ORBIT(3)))
       ALLOCATE(OMI_LON(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_LAT(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_AMF_TROP(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_NO2_TROP(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_NO2_TROP_STD(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_VIEW_ZENITH(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_SOLAR_ZENITH(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_CLOUDFR(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_X_QUAL_FLAG(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_M_QUAL_FLAG(DIMS_ORBIT(3)))
       ALLOCATE(OMI_SCW_P(DIMS_ORBIT(1)))
       ALLOCATE(OMI_TROPO_PRESSURE(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(OMI_SCATTERING_WEIGHTS(DIMS_ORBIT(1),DIMS_ORBIT(2),DIMS_ORBIT(3)) )
       
       DIMS_INDEX = DIMS_COUNTER+DIMS_ORBIT(3)-1

       !! read in OMI data
       !! read time
       dsetname = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time'
       
       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)
       !PRINT *, "DATA_DIMS_ORBIT", DATA_DIMS_ORBIT(3)
       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_TIME, (/DATA_DIMS_ORBIT(3)/), ERROR_ORBIT)
       !PRINT *, "DIMS_COUNTER", DIMS_COUNTER
       !PRINT *, "DIMS_COUNTER2", DIMS_INDEX
       !PRINT *, "SHAPE1", SHAPE(TIME_ORBIT(DIMS_COUNTER:DIMS_INDEX))
       !PRINT *, "SHAPE2", SHAPE(OMI_TIME)
       TIME_ORBIT(DIMS_COUNTER:DIMS_INDEX) = OMI_TIME
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read tropospheric air mass factors

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTrop'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_AMF_TROP, DATA_DIMS_ORBIT, ERROR_ORBIT)
       AMF_TROP_ORBIT(:,DIMS_COUNTER:DIMS_INDEX) = OMI_AMF_TROP
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "AMF_TROP_ORBIT", SHAPE(AMF_TROP_ORBIT)

       !! read tropospheric NO2 column

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_NO2_TROP, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       NO2_TROP(:,DIMS_COUNTER:DIMS_INDEX) = OMI_NO2_TROP
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "NO2_TROP", NO2_TROP(:,DIMS_ORBIT(3))
       !! read tropospheric NO2 column
       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TropopausePressure'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_TROPO_PRESSURE, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       TROPO_PRESSURE(:,DIMS_COUNTER:DIMS_INDEX) = OMI_TROPO_PRESSURE
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !! read longitudes

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropStd'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_NO2_TROP_STD, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       NO2_TROP_STD(:,DIMS_COUNTER:DIMS_INDEX) = OMI_NO2_TROP_STD
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read longitudes

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_LON, DATA_DIMS_ORBIT, ERROR_ORBIT)
       LON_ORBIT(:,DIMS_COUNTER:DIMS_INDEX) = OMI_LON
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "LONG", LON_ORBIT(:,DIMS_ORBIT(3))


       !! read latitudes

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_LAT, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       LAT_ORBIT(:,DIMS_COUNTER:DIMS_INDEX) = OMI_LAT
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "LON", LAT_ORBIT(:,DIMS_ORBIT(3))

       !! read viewing zenith angles

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/ViewingZenithAngle'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_VIEW_ZENITH, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       VIEW_ZENITH(:,DIMS_COUNTER:DIMS_INDEX) = OMI_VIEW_ZENITH
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read solar zenith angles

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarZenithAngle'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_SOLAR_ZENITH, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       SOLAR_ZENITH(:,DIMS_COUNTER:DIMS_INDEX) = OMI_SOLAR_ZENITH
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

      !! read cloud fraction

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_CLOUDFR, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       CLOUDFR(:,DIMS_COUNTER:DIMS_INDEX) = OMI_CLOUDFR
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read scattering weight pressures

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWtPressure'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_SCW_P, (/DATA_DIMS_ORBIT(1),0/), ERROR_ORBIT)
       SCW_P = OMI_SCW_P
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "data_dims_orbit", (/DATA_DIMS_ORBIT(1),0/)
       !! read scattering weights

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_SCATTERING_WEIGHTS, DATA_DIMS_ORBIT, ERROR_ORBIT)
       SCATTERING_WEIGHTS(:,:,DIMS_COUNTER:DIMS_INDEX) = OMI_SCATTERING_WEIGHTS
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "SCATTERING WEIGHTS", SHAPE(SCATTERING_WEIGHTS)
       !! close file
       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/XTrackQualityFlags'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_X_QUAL_FLAG, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)
       X_QUAL_FLAG(:,DIMS_COUNTER:DIMS_INDEX) = OMI_X_QUAL_FLAG
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "X QUAL FLAG", SHAPE(OMI_X_QUAL_FLAG)
       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/MeasurementQualityFlags'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, OMI_M_QUAL_FLAG, (/DATA_DIMS_ORBIT(3)/), ERROR_ORBIT)
       M_QUAL_FLAG(DIMS_COUNTER:DIMS_INDEX) = OMI_M_QUAL_FLAG
       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)
       !PRINT *, "M_QUAL_FLAG", SHAPE(OMI_M_QUAL_FLAG)

       CALL H5FCLOSE_F(FILE_ID_ORBIT,ERROR_ORBIT)
       DIMS_COUNTER = DIMS_COUNTER+DIMS_ORBIT(3)
       CALL H5CLOSE_F(ERROR_ORBIT)
       ! deallocate OMI arrays
       IF(ALLOCATED(OMI_LON)) DEALLOCATE(OMI_LON)
       IF(ALLOCATED(OMI_LAT)) DEALLOCATE(OMI_LAT)
       IF(ALLOCATED(OMI_TIME)) DEALLOCATE(OMI_TIME)
       IF(ALLOCATED(OMI_AMF_TROP)) DEALLOCATE(OMI_AMF_TROP)
       IF(ALLOCATED(OMI_NO2_TROP)) DEALLOCATE(OMI_NO2_TROP)
       IF(ALLOCATED(OMI_NO2_TROP_STD)) DEALLOCATE(OMI_NO2_TROP_STD)
       IF(ALLOCATED(OMI_VIEW_ZENITH)) DEALLOCATE(OMI_VIEW_ZENITH)
       IF(ALLOCATED(OMI_SOLAR_ZENITH)) DEALLOCATE(OMI_SOLAR_ZENITH)
       IF(ALLOCATED(OMI_CLOUDFR)) DEALLOCATE(OMI_CLOUDFR)
       IF(ALLOCATED(OMI_SCW_P)) DEALLOCATE(OMI_SCW_P)
       IF(ALLOCATED(OMI_X_QUAL_FLAG)) DEALLOCATE(OMI_X_QUAL_FLAG)
       IF(ALLOCATED(OMI_M_QUAL_FLAG)) DEALLOCATE(OMI_M_QUAL_FLAG)
       IF(ALLOCATED(OMI_TROPO_PRESSURE)) DEALLOCATE(OMI_TROPO_PRESSURE)
       IF(ALLOCATED(OMI_SCATTERING_WEIGHTS)) DEALLOCATE(OMI_SCATTERING_WEIGHTS)
    ENDDO
    CLOSE(75)
    !N_OMI_ORBITS = DATA_DIMS
  END SUBROUTINE READ_OMI_NO2_FILE
!================================================================================================================================
    SUBROUTINE CALC_OMI_NO2_FORCE

    USE HDF5

    !!
    !! Subroutine CALC_OMI_NO2_FORCE computes the NO2 adjoint forcing and cost function contribution from OMI column data
    !!
    !! References:
    !!
    !! Bucsela2013:
    !! "A new stratospheric and tropospheric NO2 retrieval algorithm for nadir-viewing satellite instruments: applications to OMI"
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

    USE COMODE_MOD,         ONLY : JLOP
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM, CSPEC
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
    USE GRID_MOD,           ONLY : GET_IJ
    USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
    USE TIME_MOD,           ONLY : GET_HOUR, GET_DAY
    USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS, GET_TS_CHEM
    USE DAO_MOD,            ONLY : BXHEIGHT
    USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
    USE ADJ_ARRAYS_MOD,     ONLY : ID2C
    USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
    USE TRACERID_MOD,       ONLY : IDNO2
    USE ADJ_ARRAYS_MOD,     ONLY : COST_FUNC
    USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
    USE TRACER_MOD,         ONlY : XNUMOLAIR
    USE DAO_MOD,            ONLY : T, AIRDEN
    USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR


    INTEGER :: I,J,L
    INTEGER :: I_OMI, J_OMI, K_OMI, JLOOP
    INTEGER :: IIJJ(2)
    INTEGER :: DAY
    INTEGER :: NTSTART_OMI, NTSTOP_OMI

    ! variables for time unit conversion
    REAL*8 :: tai93
    INTEGER :: iy,im,id,ih,imin
    REAL*8 :: sec
    INTEGER :: GC_HOUR, MIN_HOUR, MAX_HOUR

    ! variables for observation operator and adjoint thereof

    REAL*8 :: OMI_NO2_GC(IIPAR,JJPAR)
    REAL*8 :: SCW_GC(LLPAR), DP(LLPAR)
    REAL*8 :: AMF_GC
    REAL*8 :: GC_NO2(LLPAR)
    REAL*8 :: GC_NO2_COL
    REAL*8 :: DIFF, FORCE_COL, COST_CONTRIB_COL
    REAL*8 :: OBS_ERROR
    REAL*8, SAVE :: OMI_HOUR(MAX_ORBITS)
    REAL*8 :: OLD_COST_OMI

    ! arrays needed for superobservations 
    LOGICAL :: SUPER_OBS = .TRUE.                       ! do super observations?
    REAL*8 ::  SOBS_COUNT(IIPAR,JJPAR)                   ! super observation count
    REAL*8 ::  SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)         ! super observation adjoint forcing
    REAL*8  :: GC_ADJ_COUNT(IIPAR,JJPAR,LLPAR)
    REAL*8 ::  NEW_COST(IIPAR,JJPAR)       ! super observation cost function contribution
    REAL*8 ::  SOBS_GC(IIPAR,JJPAR)       
    REAL*8 ::  SOBS_OMI(IIPAR,JJPAR)       
    REAL*8 ::  SOBS_BIAS(IIPAR,JJPAR)
    REAL*8 ::  SOBS_CHISQUARED(IIPAR,JJPAR)

    LOGICAL, SAVE               :: SECOND = .TRUE.
    INTEGER                     :: IOS
    CHARACTER(LEN=255)          :: FILENAME

    IF ( SECOND ) THEN
       FILENAME = 'lat_orb_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 301,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'lon_orb_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 302,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'amf_gc_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 303,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'amf_obs_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 304,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'sobs_count_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 305,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       FILENAME = 'diff_omino2.NN.m'
       CALL EXPAND_NAME( FILENAME, N_CALC )
       FILENAME =  TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
       OPEN( 312,      FILE=TRIM( FILENAME  ), STATUS='UNKNOWN', IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    ENDIF


    !=================================================================
    ! CALC_OMI_NO2_FORCE begins here!
    !=================================================================
    ! initialize arrays

    GC_NO2 = 0d0
    GC_NO2_COL = 0d0
    OMI_NO2_GC = 0d0
    OLD_COST_OMI = COST_FUNC
    SOBS_COUNT = 0d0
    SOBS_ADJ_FORCE = 0d0
    GC_ADJ_COUNT = 0d0
    NEW_COST = 0d0
    SOBS_GC = 0d0
    SOBS_OMI = 0d0
    SOBS_BIAS = 0d0
    SOBS_CHISQUARED = 0d0
    !DATA_DIMS_ORBIT(:) = 4565449537985793824
    ! Loop through data to find observations
    !PRINT *, "ID2C(IDNO2)", ID2C(IDNO2)
    GC_HOUR = GET_HOUR()
    DAY = GET_DAY()
    IF ( GET_NHMS() == 236000 - GET_TS_CHEM()* 100 ) THEN
       CALL READ_OMI_NO2_FILE(GET_NYMD(), N_OMI_ORBITS)!GET_NHMS())
       DO I_OMI = 1, N_OMI_ORBITS
          IF (TIME_ORBIT(I_OMI) > 0) THEN
             CALL TAI2UTC(TIME_ORBIT(I_OMI),IY,IM,ID,IH,IMIN,SEC)
             IF (ID == DAY) THEN
                OMI_HOUR(I_OMI) = IH
             ELSE
                OMI_HOUR(I_OMI) = -999
             ENDIF
          ENDIF
       ENDDO
    ENDIF
    !! loop over data
    CALL GET_NT_RANGE_OMI(N_OMI_ORBITS, GET_NHMS(), OMI_HOUR(1:N_OMI_ORBITS), NTSTART_OMI, NTSTOP_OMI)
    IF ( NTSTART_OMI == 0 .and. NTSTOP_OMI == 0 ) THEN

       print*, ' No matching OMI NO2 obs for this hour'
       RETURN
    ENDIF
    PRINT *, 'found record range:', NTSTART_OMI, NTSTOP_OMI
    !PRINT *, 'X_QUAL_FLAG',X_QUAL_FLAG(:,NTSTOP_OMI-1:NTSTART_OMI)
    !PRINT *, 'M_QUAL_FLAG',M_QUAL_FLAG(NTSTOP_OMI-1:NTSTART_OMI)
    !PRINT *, "TIME_FRAC", TIME_FRAC(1:N_OMI_ORBITS)
    DO I_OMI=NTSTART_OMI,NTSTOP_OMI,-1
       DO J_OMI=1,N_OMI_SWATHS
          !PRINT *, "QFLAG", M_QUAL_FLAG(I_OMI), X_QUAL_FLAG(J_OMI,I_OMI), NO2_TROP(J_OMI,I_OMI)
          ! A number of conditions have to be met for OMI NO2 data to actually be assimilated
          IF ( ( TIME_ORBIT(I_OMI) > 0 ) .AND. &
#if defined(NESTED_NA) || defined(NESTED_CH)
               ( LON_ORBIT(J_OMI,I_OMI) >= GET_XMID(1)    ) .AND. &
               ( LON_ORBIT(J_OMI,I_OMI) <= GET_XMID(IIPAR)) .AND. &
               ( LAT_ORBIT(J_OMI,I_OMI) >= GET_YMID(1)    ) .AND. &
               ( LAT_ORBIT(J_OMI,I_OMI) <= GET_YMID(JJPAR)) .AND. &
#endif
               ( ABS(LAT_ORBIT(J_OMI,I_OMI)) < 60d0       ) .AND. &
               ( NO2_TROP(J_OMI,I_OMI) > 0d0              ) .AND. &
               ( NO2_TROP_STD(J_OMI,I_OMI) > 0d0          ) .AND. &
               ( ABS(SOLAR_ZENITH(J_OMI,I_OMI)) < 75d0    ) .AND. &
               ( ABS(VIEW_ZENITH(J_OMI,I_OMI)) < 65d0     ) .AND. &
               ( AMF_TROP_ORBIT(J_OMI,I_OMI) > 0d0        ) .AND. &
               ( REAL(X_QUAL_FLAG(J_OMI,I_OMI)) < 1d0     ) .AND. &
               ( M_QUAL_FLAG(I_OMI) == 0d0          ) .AND. &
               ( CLOUDFR(J_OMI,I_OMI) >= 0d0              ) .AND. &
               ( CLOUDFR(J_OMI,I_OMI) < 200   ) ) THEN
             ! Get model grid coordinate indices that correspond to the observation
             
             IIJJ = GET_IJ(REAL(LON_ORBIT(J_OMI,I_OMI),4), REAL(LAT_ORBIT(J_OMI,I_OMI),4))
             
             I = IIJJ(1)
             J = IIJJ(2)
             
             ! initialize variables & arrays
             
             GC_NO2 = 0d0
             GC_NO2_COL = 0d0
             SCW_GC = 0d0
             DP = 0d0
             COST_CONTRIB_COL = 0d0
             FORCE_COL = 0d0
             
             ! Get GEOS-CHEM NO2 values [#/cm3]
             
             DO L = 1, LLPAR
                
                !IF( ITS_IN_THE_TROP(I,J,L) ) THEN
                IF ( GET_PEDGE(I,J,L) >= TROPO_PRESSURE(J_OMI,I_OMI) ) THEN
                   JLOOP = JLOP(I,J,L)
                   IF (GET_PEDGE(I,J,L) >= 400d0) THEN
                      GC_NO2(L) = CSPEC_AFTER_CHEM(JLOOP,ID2C(IDNO2))
                   ELSE
                      GC_NO2(L) = CSPEC(JLOOP,IDNO2)
                   ENDIF
                   !PRINT *, "GC_NO2", GC_NO2(2)
                ENDIF
                
             ENDDO
             
             ! Compute tropospheric NO2 vertical column [#/cm2]
             
             GC_NO2_COL = SUM(GC_NO2(:) * BXHEIGHT(I,J,:)*100d0)                   
             
             ! interpolate scattering weights to GEOS-Chem grid to compute GEOS-Chem air mass factors
             ! question: how do differences in surface pressures used in the retrieval and GEOS-Chem affect the computation below?
             
             DO L=1,LLPAR
                DO K_OMI = 2,N_OMI_LEVELS
                   
                   IF( GET_PCENTER(I,J,L) < SCW_P(K_OMI-1) .AND. GET_PCENTER(I,J,L) > SCW_P(K_OMI) ) THEN
                      
                      ! linearly interpolate scattering weights to GEOS-Chem center pressures
                      
                      SCW_GC(L) = SCATTERING_WEIGHTS(K_OMI,J_OMI,I_OMI) + &
                           ( SCATTERING_WEIGHTS(K_OMI-1,J_OMI,I_OMI) - SCATTERING_WEIGHTS(K_OMI,J_OMI,I_OMI) ) * &
                           ( GET_PCENTER(I,J,L) - SCW_P(K_OMI) ) / ( SCW_P(K_OMI-1) - SCW_P(K_OMI) )
                      
                      ! save pressure difference of edge pressures
                      
                      DP(L) = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                      
                      ! apply temperature correction, as in Bucsela2013, eq. (4)
                      
                      SCW_GC(L) = SCW_GC(L) * ( 1 - 0.003 * ( T(I,J,L) - 220 ) )
                      
                      ! convert NO2 concentrations from number density to mixing ratio, as required for the calculation of the air mass factors from scattering weights
                      
                      GC_NO2(L) = GC_NO2(L) *1d6 / ( AIRDEN(L,I,J)   * XNUMOLAIR )
                      
                      EXIT ! exit loop over K_OMI to go to next cycle in loop over L
                      
                   ENDIF
                   
                ENDDO
             ENDDO
             
             ! Use GEOS-Chem tropospheric air mass factors to convert vertical column to slant column
             AMF_GC = SUM(GC_NO2 * DP * SCW_GC)/SUM(GC_NO2 * DP)
             !PRINT *, "AMF_GC", AMF_GC
             GC_NO2_COL = AMF_GC*GC_NO2_COL
             
             ! The computation above is a little awkward, since the slant column can be computed directly from equation (2) in Bucsela2013 without
             ! computing the airmass factors and NO2 column first.
             ! I chose to compute the slant column from the computed air mass factors (which already included the computation of the slant column)
             ! since the air mass factors might be of diagnostic interest and can be computed and saved
             ! alongside other observation operator diagnostics. Furthermore, this formulation makes the adjoint of the observation operator somewhat simpler to handle.
             
             ! compute slant column difference
             
             DIFF = GC_NO2_COL - (NO2_TROP(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI))
             !PRINT *, "GC_NO2_COL", GC_NO2_COL
             !PRINT *, "NO2_TROP", GC_NO2_COL-DIFF
             !PRINT *, "AMF_GC", AMF_GC
             !PRINT *, "AMF_OBS", AMF_TROP_ORBIT(J_OMI,I_OMI)
             ! compute slant column standard deviation
             !PRINT *, "CHECK"
             OBS_ERROR =  0.05*NO2_TROP_STD(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI)
             IF (OBS_ERROR>0d0) THEN
                FORCE_COL = DIFF/(OBS_ERROR**2)
                COST_CONTRIB_COL = 0.5d0 * DIFF * FORCE_COL
             ELSE
                FORCE_COL = 0d0
                COST_CONTRIB_COL = 0d0
             ENDIF
             ! update adjoint NO2 concentration
             IF ( ( COST_CONTRIB_COL > 0d0 ) .AND. &
                  ( COST_CONTRIB_COL <= 80000d0) ) THEN
                DO L = 1, LLPAR
                   IF (ITS_IN_THE_TROP(I,J,L)) THEN
                      ! question: how do errors in retrieved surface pressure impact the NO2 column values?
                      ! question: how do errors in simulated surface pressures impact the NO2 column values?
                      
                      IF (SUPER_OBS) THEN                            
                         SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + FORCE_COL * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                         GC_ADJ_COUNT(I,J,L) = GC_ADJ_COUNT(I,J,L) + 1
                      ELSE
                         JLOOP = JLOP(I,J,L)
                         CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) &
                              + FORCE_COL * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                      ENDIF
                   ENDIF
                   
                ENDDO
                ! update cost function
                
                IF(SUPER_OBS) THEN
                   NEW_COST(I,J) = NEW_COST(I,J) + COST_CONTRIB_COL
                   SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                ELSE
                   COST_FUNC = COST_FUNC + COST_CONTRIB_COL
                ENDIF
                
             ENDIF
             
             WRITE(301,110) ( LAT_ORBIT(J_OMI,I_OMI)      )
             WRITE(302,110) ( LON_ORBIT(J_OMI,I_OMI)      )
             !WRITE(303,110) ( AMF_GC                      )
             !WRITE(304,110) ( AMF_TROP_ORBIT(J_OMI,I_OMI) )
             WRITE(312,110) ( DIFF/1e10                   )
110          FORMAT(F18.6,1X)
             
             ! update diagnostic arrays
             
             IF( SUPER_OBS) THEN                      
                SOBS_GC(I,J) = SOBS_GC(I,J) + GC_NO2_COL
                SOBS_OMI(I,J) = SOBS_OMI(I,J) + NO2_TROP(J_OMI,I_OMI)
                SOBS_BIAS(I,J) = SOBS_BIAS(I,J) + DIFF
                SOBS_CHISQUARED(I,J) = SOBS_CHISQUARED(I,J) + 0.5 * (DIFF/OBS_ERROR)**2
             ELSE
                
                OMI_BIAS_COUNT(I,J) = OMI_BIAS_COUNT(I,J) + 1d0
                
                OMI_NO2_MEAN(I,J) = OMI_NO2_MEAN(I,J) + NO2_TROP(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI)
                
                OMI_GEOS_NO2_MEAN(I,J) = OMI_GEOS_NO2_MEAN(I,J) + GC_NO2_COL
                
                OMI_NO2_ERR_MEAN(I,J) = OMI_NO2_ERR_MEAN(I,J) + OBS_ERROR
                
                OMI_DELTA = DIFF - OMI_BIAS(I,J)
                
                OMI_BIAS(I,J) = OMI_BIAS(I,J) + OMI_DELTA/OMI_BIAS_COUNT(I,J)
                
                OMI_VAR(I,J) = OMI_VAR(I,J) + OMI_DELTA*(DIFF-OMI_BIAS(I,J))
                
                OMI_CHISQUARED(I,J) = OMI_CHISQUARED(I,J) + ( DIFF/OBS_ERROR )**2
                
             ENDIF
          ENDIF ! data selection IF statement
       
       ENDDO ! J
    ENDDO ! I
    
    IF(SUPER_OBS) THEN
       
       DO J=1,JJPAR
          DO I=1,IIPAR
             
             IF(SOBS_COUNT(I,J) > 0d0) THEN
                
                DO L=1,LLPAR
                   JLOOP = JLOP(I,J,L)
                   IF( ( ITS_IN_THE_TROP(I,J,L) ).AND. &
                        ( GC_ADJ_COUNT(I,J,L) > 0 ) .AND. &
                        ( JLOOP > 0 ) ) THEN
                         
                      CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) &
                           + SOBS_ADJ_FORCE(I,J,L)/GC_ADJ_COUNT(I,J,L)
                      
                   ENDIF
                   
                ENDDO
                WRITE(305,110) (SOBS_COUNT(I,J))
                COST_FUNC = COST_FUNC + NEW_COST(I,J)/SOBS_COUNT(I,J)
                
                OMI_BIAS_COUNT(I,J) = OMI_BIAS_COUNT(I,J) + 1d0
                
                OMI_NO2_MEAN(I,J) = OMI_NO2_MEAN(I,J) + SOBS_OMI(I,J)/SOBS_COUNT(I,J)
                
                OMI_GEOS_NO2_MEAN(I,J) = OMI_GEOS_NO2_MEAN(I,J) + SOBS_GC(I,J)/SOBS_COUNT(I,J)
                
                OMI_NO2_ERR_MEAN(I,J) = OMI_NO2_ERR_MEAN(I,J) + OBS_ERROR  ! mkeller: need to change this to reflect super observation error, but how?
                
                ! calculate bias and variance of GC-OMI bias using numerically stable one-pass algorithm (Chan83)
                
                OMI_DELTA = SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMI_BIAS(I,J)
                
                OMI_BIAS(I,J) = OMI_BIAS(I,J) + OMI_DELTA/OMI_BIAS_COUNT(I,J)
                
                OMI_VAR(I,J) = OMI_VAR(I,J) + OMI_DELTA*(SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMI_BIAS(I,J))
                
                OMI_CHISQUARED(I,J) = OMI_CHISQUARED(I,J) + SOBS_CHISQUARED(I,J)/SOBS_COUNT(I,J)
                
             ENDIF
             
          ENDDO
       ENDDO
          
    ENDIF
    PRINT *, "OMI NO2 COST FUNCTION", COST_FUNC-OLD_COST_OMI
  END SUBROUTINE CALC_OMI_NO2_FORCE

  !-----------------------------------------------------------------------------!
  SUBROUTINE CLEANUP_OMI
    ! deallocate OMI arrays

    IF(ALLOCATED(LON_ORBIT)) DEALLOCATE(LON_ORBIT)
    IF(ALLOCATED(LAT_ORBIT)) DEALLOCATE(LAT_ORBIT)
    IF(ALLOCATED(TIME_ORBIT)) DEALLOCATE(TIME_ORBIT)
    IF(ALLOCATED(AMF_TROP_ORBIT)) DEALLOCATE(AMF_TROP_ORBIT)
    IF(ALLOCATED(NO2_TROP)) DEALLOCATE(NO2_TROP)
    IF(ALLOCATED(NO2_TROP_STD)) DEALLOCATE(NO2_TROP_STD)
    IF(ALLOCATED(VIEW_ZENITH)) DEALLOCATE(VIEW_ZENITH)
    IF(ALLOCATED(SOLAR_ZENITH)) DEALLOCATE(SOLAR_ZENITH)
    IF(ALLOCATED(CLOUDFR)) DEALLOCATE(CLOUDFR)
    IF(ALLOCATED(SCW_P)) DEALLOCATE(SCW_P)
    IF(ALLOCATED(X_QUAL_FLAG)) DEALLOCATE(X_QUAL_FLAG)
    IF(ALLOCATED(M_QUAL_FLAG)) DEALLOCATE(M_QUAL_FLAG)
    IF(ALLOCATED(TROPO_PRESSURE)) DEALLOCATE(TROPO_PRESSURE)
    IF(ALLOCATED(SCATTERING_WEIGHTS)) DEALLOCATE(SCATTERING_WEIGHTS)
    
  END SUBROUTINE CLEANUP_OMI
!----------------------------------------------------------------------------------
  SUBROUTINE GET_NT_RANGE_OMI( N_OMI_ORBITS, HHMMSS, OMI_HOUR, NTSTART_OMI, NTSTOP_OMI)
!
!******************************************************************************
!  Subroutine GET_NT_RANGE_OMI retuns the range of retrieval records for the
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
!  (1 ) NTSTART_OMI (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP_OMI  (INTEGER) : TES record number at which to stop
!
!  NOTES:
!
!******************************************************************************
!
    ! Reference to f90 modules
    USE ERROR_MOD,    ONLY : ERROR_STOP
    USE TIME_MOD,     ONLY : YMD_EXTRACT
    
    ! Arguments
    INTEGER, INTENT(IN)   :: N_OMI_ORBITS
    INTEGER, INTENT(IN)   :: HHMMSS
    REAL*8,  INTENT(IN)   :: OMI_HOUR(N_OMI_ORBITS)
    INTEGER, INTENT(OUT)  :: NTSTART_OMI
    INTEGER, INTENT(OUT)  :: NTSTOP_OMI
    
    ! Local variables
    INTEGER, SAVE         :: NTSAVE_OMI
    LOGICAL               :: FOUND_ALL_RECORDS
    INTEGER               :: NTEST_OMI
    INTEGER               :: HH, MM, SS
    REAL*8                :: GC_HH_FRAC
    REAL*8                :: H1_FRAC
    
    !=================================================================
    ! GET_NT_RANGE_OMI begins here!
    !=================================================================
    ! Initialize
    FOUND_ALL_RECORDS  = .FALSE.
    NTSTART_OMI            = 0
    NTSTOP_OMI             = 0
    
    ! set NTSAVE_OMI to NTES every time we start with a new file
    IF ( HHMMSS == 230000 ) NTSAVE_OMI = N_OMI_ORBITS       
    !print*, ' GET_NT_RANGE_OMI for ', HHMMSS
    !print*, ' NTSAVE_OMI ', NTSAVE_OMI
    !print*, ' N_IASI_NOB   ', N_IASI_NOB
    DO WHILE (OMI_HOUR(NTSAVE_OMI) < 0 )
       NTSAVE_OMI = NTSAVE_OMI - 1
       IF (NTSAVE_OMI == 0) EXIT       
    ENDDO

    !PRINT *, "TIME_FRAC", TIME_FRAC(NTSAVE_OMI-1000:NTSAVE_OMI)
    CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )
    
    
    ! Convert HH from hour to fraction of day
    GC_HH_FRAC = REAL(HH,8)
    ! one hour as a fraction of day
    H1_FRAC    = 0d0
    
    
    ! All records have been read already
    IF ( NTSAVE_OMI == 0 ) THEN
       
       print*, 'All records have been read already '
       RETURN
       
       ! No records reached yet
    ELSEIF ( OMI_HOUR(NTSAVE_OMI) + H1_FRAC < GC_HH_FRAC ) THEN
       
       
       print*, 'No records reached yet'
       RETURN
       
       !
    ELSEIF ( OMI_HOUR(NTSAVE_OMI) + H1_FRAC >=  GC_HH_FRAC ) THEN
       ! Starting record found
       NTSTART_OMI = NTSAVE_OMI
       
       !print*, ' Starting : TIME_FRAC(NTSTART_OMI) ', TIME_FRAC(NTSTART_OMI), NTSTART_OMI
       
       ! Now search forward to find stopping record
       NTEST_OMI = NTSTART_OMI
       
       DO WHILE ( FOUND_ALL_RECORDS == .FALSE. )
          
          ! Advance to the next record
          NTEST_OMI = NTEST_OMI - 1
          
          ! Stop if we reach the earliest available record
          IF ( NTEST_OMI == 0 ) THEN
             
             NTSTOP_OMI            = NTEST_OMI + 1
             FOUND_ALL_RECORDS = .TRUE.
             
             !print*, ' Records found '
             !print*, ' NTSTART_OMI, NTSTOP_OMI = ', NTSTART_OMI, NTSTOP_OMI
             ! Reset NTSAVE_OMI
             NTSAVE_OMI = NTEST_OMI
             ! When the combined test date rounded up to the nearest
             ! half hour is smaller than the current model date, the
             ! stopping record has been passed.
          ELSEIF (  OMI_HOUR(NTEST_OMI) + H1_FRAC <  GC_HH_FRAC ) THEN
             
             !print*, ' Testing : TIME_FRAC ', TIME_FRAC(NTEST_OMI), NTEST_OMI
             
             NTSTOP_OMI            = NTEST_OMI + 1
             FOUND_ALL_RECORDS = .TRUE.
             
             !print*, ' Records found '
             !print*, ' NTSTART_OMI, NTSTOP_OMI = ', NTSTART_OMI, NTSTOP_OMI
             
             ! Reset NTSAVE_OMI
             NTSAVE_OMI = NTEST_OMI
             !ELSE
             !print*, ' still looking ', NTEST_OMI
             
          ENDIF
          
       ENDDO
       
    ELSE
       
       CALL ERROR_STOP('problem', 'GET_NT_RANGE_OMI' )
       
    ENDIF
    
    ! Return to calling program
  END SUBROUTINE GET_NT_RANGE_OMI
!================================================================================================================================

  SUBROUTINE TAI2UTC(tai93,iy,im,id,ih,imin,sec)

    !!
    !! SUBROUTINE TAI2UTC converts TAI93 time (seconds since 1.1.1993) to UTC
    !!
    !! adapted from
    !! code.google.com/p/miyoshi/source/browse/trunk/common/common.f90
    !!
    
    IMPLICIT NONE

    INTEGER,PARAMETER :: N=10  ! number of leap seconds after Jan. 1, 1993
    !-----------------------------------93/06/30---94/06/30--95/12/31---97/06/30---98/12/31---05/12/31---08/12/31--12/06/30---15/06/30--16/12/31
    INTEGER,PARAMETER :: LEAPSEC(N) = (/ 15638399, 47174400, 94608001, 141868802, 189302403, 410227204, 504921605, 615254406, 709862407, 757382408/)
    REAL*8,INTENT(IN) :: TAI93
    INTEGER,INTENT(OUT) :: IY,IM,ID,IH,IMIN
    REAL*8,INTENT(OUT) :: SEC
    REAL*8,PARAMETER :: MINS = 60.0D0
    REAL*8,PARAMETER :: HOUR = 60.0D0*MINS
    REAL*8,PARAMETER :: DAY = 24.0D0*HOUR
    REAL*8,PARAMETER :: YEAR = 365.0D0*DAY
    INTEGER,PARAMETER :: MDAYS(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    REAL*8 :: WK,TAI
    INTEGER :: DAYS,I,LEAP

    TAI = TAI93
    SEC = 0.0D0

    DO I=1,N

       IF(FLOOR(TAI93) == LEAPSEC(I)+1) THEN
          SEC = 60.0D0 + TAI93-FLOOR(TAI93)
       ENDIF

       IF(FLOOR(TAI93) > LEAPSEC(I)) TAI = TAI -1.0D0

    END DO

    IY = 1993 + FLOOR(TAI /YEAR)
    WK = TAI - REAL(IY-1993)*YEAR - FLOOR(REAL(IY-1993)/4.0)*DAY

    IF(WK < 0.0D0) THEN
       IY = IY -1
       WK = TAI - REAL(IY-1993)*YEAR - FLOOR(REAL(IY-1993)/4.0)*DAY
    END IF

    DAYS = FLOOR(WK/DAY)
    WK = WK - REAL(DAYS)*DAY
    IM = 1

    DO I=1,12

       LEAP = 0
       IF(IM == 2 .AND. MOD(IY,4)==0) LEAP=1
       IF(IM == I .AND. DAYS >= MDAYS(I)+LEAP) THEN
          IM = IM + 1
          DAYS = DAYS - MDAYS(I)-LEAP
       END IF

    END DO

    ID = DAYS +1

    IH = FLOOR(WK/HOUR)
    WK = WK - REAL(IH)*HOUR
    IMIN = FLOOR(WK/MINS)

    IF(SEC < 60.0D0) SEC = WK - REAL(IMIN)*MINS

    RETURN

  END SUBROUTINE TAI2UTC

  !----------------------------------------------------------
END MODULE OMI_NO2_OBS_MOD