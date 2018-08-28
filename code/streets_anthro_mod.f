! $Id: streets_anthro_mod.f,v 1.3 2012/05/09 22:31:56 nicolas Exp $
      MODULE STREETS_ANTHRO_MOD
!
!******************************************************************************
!  Module STREETS_ANTHRO_MOD contains variables and routines to read the
!  David Streets et al Asian anthropogenic emissions for NOx and CO.
!  (yxw, bmy, 8/16/06, 3/11/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) A_CM2          (REAL*8 ) : Array for grid box surface area [cm2]
!  (2 ) MASK_CHINA_1x1 (INTEGER) : Mask for the China region at 1x1
!  (2 ) MASK_CHINA     (REAL*8)  : Mask for the China region (for 2001 CO)
!  (3 ) MASK_SE_ASIA   (REAL*8)  : Mask for the SE Asia region (for 2000 emiss)
!  (4 ) NOx            (REAL*8)  : Streets anthro NOx emissions [kg/yr]
!  (5 ) CO             (REAL*8)  : Streets anthro CO  emissions [kg/yr]
!  (6 ) SO2            (REAL*8)  : Streets anthro SO2 emissions [kg/yr]
!  (7 ) NH3            (REAL*8)  : Streets anthro NH3 emissions [kg/yr]
!  (8 ) CO2            (REAL*8)  : Streets anthro CO2 emissions [kg/yr]
!  (9 ) CH4            (REAL*8)  : Streets anthro CH4 emissions [kg/yr]
!  (10) following VOC in [atoms C/yr] or [molec/yr]: ACET, ALD2, ALK4, C2H6,
!       C3H8, CH2O, ISOP, MEK, PRPE
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_CHINA_MASK         : Gets the China mask value at (I,J)
!  (2 ) GET_SE_ASIA_MASK       : Gets the SE Asia mask value at (I,J)
!  (3 ) GET_STREETS_ANTHRO     : Gets emissions at (I,J) for emissions species
!  (4 ) EMISS_STREETS_ANTHRO   : Reads Streets' emissions from disk
!  (5 ) STREETS_SCALE_FUTURE   : Applies IPCC future scale factors to emissions
!  (6 ) READ_STREETS_MASKS     : Reads mask info from disk
!  (7 ) INIT_STREETS_ANTHRO    : Allocates and zeroes module arrays
!  (8 ) CLEANUP_STREETS_ANTHRO : Dealocates module arrays
!
!  GEOS-Chem modules referenced by "streets_anthro_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f        : Module w/ GEOS-Chem data & met field dirs
!  (3 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (4 ) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (5 ) grid_mod.f             : Module w/ horizontal grid information
!  (6 ) logical_mod.f          : Module w/ GEOS-Chem logical switches
!  (7 ) regrid_1x1_mod.f       : Module w/ routines to regrid 1x1 data
!  (8 ) time_mod.f             : Module w/ routines for computing time & date
!  (9 ) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!
!  References:
!  ============================================================================
!  (1 ) Streets, D.G, Q. Zhang, L. Wang, K. He, J. Hao, Y. Wu, Y. Tang,
!        and G.C. Carmichael, "Revisiting China's CO emissions after the
!        Transport and Chemical Evolution over the Pacific (TRACE-P) mission:
!        Synthesis of inventories, atmospheric modeling, and observations",
!        J. Geophys. Res, 111, D14306, doi:10.1029/2006JD007118, 2006.
!  (2 ) Streets, D.G., T.C. Bond, G.R. Carmichael, S.D. Fernandes, Q. Fu,
!        Z. Klimont, S.M. Nelson, N.Y. Tsai, M.Q. Wang, J-H. Woo, and
!        K.F. Yarber, "An inventory of gaseous and primary aerosol emissions
!        in Asia in the year 2000", J. Geophys. Res, 108, D21,
!        doi:10.1029/2002JD003093, 2003.
!  (3) Zhang, Q., Streets, D. G., Carmichael, G., He, K., Huo, H.,
!        Kannari, A., Klimont, Z., Park, I., Reddy, S., Chen, D., Duan, L.,
!        Lei, Y., Wang, L. and Yao, Z.: Asian emissions in 2006 for the
!        NASA INTEX-B mission, manuscript submitted to Atmospheric
!        Chemistry & Physics Discussions, 2009

!
!  NOTES:
!  (1 ) Modification: Now use 2001 CO over China, and 2000 CO over countries
!        other than China in the larger SE Asia region. (yxw, bmy, 9/5/06)
!  (2 ) Modifications for 0.5 x 0.667 nested grids (yxw, dan, bmy, 11/6/08)
!  (3 ) 2006 and 2020 inventories are now available. But species emitted
!        differ (phs, 3/7/08):
!                              2000/2001 = NOx, CO, SO2, NH3, CO2, CH4
!                              2006/2020 = NOx, CO, SO2, all VOC
!  (4 ) Now scale emissions using int'annual scale factors (amv, 08/24/07)
!  (5 ) Implemented monthly variations (phs, 4/12/08)
!  (6 ) Bug fix: call READ_STREETS_05x0666 in routine
!        EMISS_STREETS_ANTHRO_05x0666 (ccc, 3/11/09)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "streets_anthro_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_STREETS_ANTHRO
      PUBLIC :: EMISS_STREETS_ANTHRO
      PUBLIC :: EMISS_STREETS_ANTHRO_05x0666
      PUBLIC :: EMISS_STREETS_ANTHRO_025x03125    ! (lzh,02/01/2015)
      PUBLIC :: GET_CHINA_MASK
      PUBLIC :: GET_SE_ASIA_MASK
      PUBLIC :: GET_STREETS_ANTHRO

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Arrays
      INTEGER, ALLOCATABLE :: MASK_CHINA_1x1(:,:)
      INTEGER, ALLOCATABLE :: MASK_CHINA_05x0666(:,:)
      REAL*8,  ALLOCATABLE :: A_CM2(:)
      REAL*8,  ALLOCATABLE :: MASK_CHINA(:,:)
      REAL*8,  ALLOCATABLE :: MASK_SE_ASIA(:,:)
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)
      REAL*8,  ALLOCATABLE :: NH3(:,:)
      REAL*8,  ALLOCATABLE :: CO2(:,:)
      REAL*8,  ALLOCATABLE :: CH4(:,:)

      ! added VOC for 2006 inventory (phs, 3/7/08)
      ! Note ISOP is not used in GEOS-Chem but it is available for 2006.
      REAL*8,  ALLOCATABLE :: ALK4(:,:)
      REAL*8,  ALLOCATABLE :: ACET(:,:)
      REAL*8,  ALLOCATABLE :: MEK(:,:)
      REAL*8,  ALLOCATABLE :: PRPE(:,:)
      REAL*8,  ALLOCATABLE :: C2H6(:,:)
      REAL*8,  ALLOCATABLE :: C3H8(:,:)
      REAL*8,  ALLOCATABLE :: CH2O(:,:)
      REAL*8,  ALLOCATABLE :: ALD2(:,:)

      ! flag to denote if emission base year is 2006
      LOGICAL IS_2006

      ! month
      INTEGER MONTH


      ! Parameters
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_CHINA_MASK( I, J ) RESULT( THISMASK )
!
!******************************************************************************
!  Function GET_STREETS_MASK returns the value of the China mask for the David
!  Streets et al emissions at grid box (I,J).  MASK=1 if (I,J) is China, or
!  MASK=0 otherwise. (bmy, 8/16/06)
!
!  NOTE: The China Mask is used with the 2001 CO emissions.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-Chem longitude index
!  (2 ) J (INTEGER) : GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      REAL*8              :: THISMASK

      !=================================================================
      ! GET_CHINA_MASK begins here!
      !=================================================================
      THISMASK = MASK_CHINA(I,J)

      ! Return to calling program
      END FUNCTION GET_CHINA_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_SE_ASIA_MASK( I, J ) RESULT( THISMASK )
!
!******************************************************************************
!  Function GET_SE_ASIA_MASK returns the value of the China mask for the David
!  Streets et al emissions at grid box (I,J).  MASK=1 if (I,J) is China, or
!  MASK=0 otherwise. (bmy, 8/16/06)
!
!  NOTE: The SE Asia Mask is used with the 2000 emissions for
!        NOx, CO, CO2, SO2, NH3, and CH4.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-Chem longitude index
!  (2 ) J (INTEGER) : GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      REAL*8              :: THISMASK

      !=================================================================
      ! GET_SE_ASIA_MASK begins here!
      !=================================================================
      THISMASK = MASK_SE_ASIA(I,J)

      ! Return to calling program
      END FUNCTION GET_SE_ASIA_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_STREETS_ANTHRO( I,    J,     N,
     &                             MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_STREETS_ANTHRO returns the David Streets et al emission for
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].  (bmy, 8/16/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude index
!  (3 ) N           (INTEGER) : GEOS-Chem tracer number
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL -- return emissions in [molec/cm2/s]
!  (5 ) KG_S        (LOGICAL) : OPTIONAL -- return emissions in [kg/s]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD,           ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,           ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNOx, IDTCO, IDTSO2, IDTNH3
      USE TRACERID_MOD,         ONLY : IDTACET, IDTALK4, IDTC2H6
      USE TRACERID_MOD,         ONLY : IDTCH2O, IDTMEK,  IDTALD2
      USE TRACERID_MOD,         ONLY : IDTPRPE, IDTC3H8

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J, N
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S

      ! Local variables
      LOGICAL                       :: DO_KGS, DO_MCS, IS_NMVOC
      REAL*8                        :: VALUE

      !=================================================================
      ! GET_STREETS_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS   = .FALSE.
      DO_MCS   = .FALSE.
      IS_NMVOC = .TRUE.

      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Test for simulation type
      IF ( ITS_A_CH4_SIM() ) THEN

         !-------------------
         ! CH4 simulation
         !-------------------
         VALUE = CH4(I,J)

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !-------------------
         ! CO2 simulation
         !-------------------
         VALUE = CO2(I,J)

      ELSE

         !-------------------
         ! Other simulations
         !-------------------
         IF ( N == IDTNOx ) THEN

            ! NOx [kg/yr]
            VALUE = NOx(I,J)

            IS_NMVOC =.FALSE.       !PHS

         ELSE IF ( N == IDTCO ) THEN

            ! CO [kg/yr]
            VALUE = CO(I,J)

            IS_NMVOC =.FALSE.       !PHS

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2(I,J)

            IS_NMVOC =.FALSE.   !PHS

         ELSE IF ( N == IDTNH3 ) THEN

            ! NH3 [kg/yr]
            VALUE = NH3(I,J)

            IS_NMVOC =.FALSE.   !PHS (bug fix, 3/2/09)

         !========= start VOC modifications (phs, 3/7/08)
         ELSE IF ( N == IDTALK4 ) THEN

            ! SO2 [kg/yr]
            VALUE = ALK4(I,J)

         ELSE IF ( N == IDTALD2 ) THEN

            ! SO2 [kg/yr]
            VALUE = ALD2(I,J)

         ELSE IF ( N == IDTPRPE ) THEN

            ! SO2 [kg/yr]
            VALUE = PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN

            ! SO2 [kg/yr]
            VALUE = C3H8(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN

            ! SO2 [kg/yr]
            VALUE = C2H6(I,J)

         ELSE IF ( N == IDTMEK ) THEN

            ! SO2 [kg/yr]
            VALUE = MEK(I,J)

         ELSE IF ( N == IDTACET ) THEN

            ! SO2 [kg/yr]
            VALUE = ACET(I,J)

         ELSE IF ( N == IDTCH2O ) THEN

            ! SO2 [kg/yr]
            VALUE = CH2O(I,J)

         !========= end VOC modifications ================
         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no STREETS emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ENDIF

      ! Check if some species are missing
      IF ( VALUE .LT. 0D0 ) RETURN


      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN

         IF ( IS_NMVOC ) THEN
            ! Convert from [atom C/yr] to [kg/s] or from [molec/yr]
            ! to [kg/s]
            VALUE = VALUE / ( XNUMOL(N) * SEC_IN_YEAR )
         ELSE
            ! Convert from [kg/yr] to [kg/s]
            VALUE = VALUE / SEC_IN_YEAR
         ENDIF

      ELSE IF ( DO_MCS ) THEN

         IF ( IS_NMVOC ) THEN
            ! Convert from [atom C/yr] to [atom C/cm2/s] or
            ! from [molec/yr] to [molec/cm2/s]
            VALUE = VALUE / ( A_CM2(J) * SEC_IN_YEAR )
         ELSE
            ! Convert from [kg/yr] to [molec/cm2/s]
            VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_YEAR )
         ENDIF

      ENDIF

      ! Return to calling program
      END FUNCTION GET_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine EMISS_STREETS_ANTHRO reads the David Streets et al emission
!  fields at 1x1 resolution and regrids them to the current model resolution.
!  (bmy, 8/16/06, 9/5/06)
!
!  NOTES:
!  (1 ) Overwrite 2000 SE Asia CO with 2001 CO over China (bmy, 9/5/06)
!  (2 ) Now can use 2000(2001 for CO over CHINA), or 2006, or 2020 inventory
!       (phs,3/07/08)
!  (3 ) Added int'annual scale factors (amv, 08/24/07)
!  (4 ) Now accounts for FSCALYR and monthly variation (phs, 3/17/08)
!  (5 ) Now NH3 2000 is used for all simulation years (phs, 2/27/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,        ONLY : GET_TAU0,          READ_BPCH2
      USE DIRECTORY_MOD,    ONLY : DATA_DIR_1x1
      USE LOGICAL_MOD,      ONLY : LFUTURE
      USE TRACER_MOD,       ONLY : ITS_A_CO2_SIM,     ITS_A_CH4_SIM
      USE TIME_MOD,         ONLY : GET_YEAR,          GET_MONTH
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR

#     include "CMN_SIZE"         ! Size parameters
#     include "CMN_O3"           ! FSCALYR

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: BASE_YEAR,       SIM_YEAR
      CHARACTER(LEN=255)         :: FILENAME,        STREETS_DIR,
     $                              STREETS_DIR_2000

      ! to loop over the sources
      INTEGER, PARAMETER         :: NSRCE = 10
      INTEGER                    :: NSOURCE, NTSOURCE1, NTSOURCE2
      CHARACTER(LEN=3)           :: SOURCES(NSRCE)   !!
      REAL*8                     :: ONOFF(NSRCE)     !! To switch off/on each sources
      REAL*8                     :: SCALE2020(NSRCE) !! To scale 2006 to 2020

      ! to hold data and scale factors
      REAL*4                     :: SCALFAC( IIPAR, JJPAR )
      REAL*8                     :: TEMP( IIPAR, JJPAR )

      ! TAUs
      REAL*8                     :: TAU2000,       TAU2004,      TAU2006
      REAL*8                     :: TAUMONTH_2001, TAUMONTH_2004, TAU

      !=================================================================
      ! EMISS_STREETS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_STREETS_ANTHRO
         FIRST = .FALSE.
      ELSE
         NOx  = 0D0
         CO   = 0D0
         SO2  = 0D0
         ALK4 = 0D0
         ACET = 0D0
         MEK  = 0D0
         PRPE = 0D0
         C2H6 = 0D0
         C3H8 = 0D0
         CH2O = 0D0
         ALD2 = 0D0
      ENDIF

      ! TAU0 values for 2000, 2004 and 2006
      TAU2000 = GET_TAU0( 1, 1, 2000 )
      TAU2004 = GET_TAU0( 1, 1, 2004 )
      TAU2006 = GET_TAU0( 1, 1, 2006 )

      MONTH         = GET_MONTH()
      TAUMONTH_2001 = GET_TAU0( MONTH, 1, 2001 )
      TAUMONTH_2004 = GET_TAU0( MONTH, 1, 2004 )

      !-------------------------------------------------------------------------
      !     Base Year & Yearly Scale Factors used
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% To simulate 2020 :                                                 %
      ! %%%         you set BASE_YEAR = 2020 (hardwired, see below)            %
      ! %%%                                                                    %
      ! %%% To simulate 2006 and after, the program sets:                      %
      ! %%%      BASE_YEAR = 2006     for all species, except NH3              %
      ! %%%      BASE_YEAR = 2000     for NH3                                  %
      ! %%%                                                                    %
      ! %%% To simulate 2005 and before, it sets:                              %
      ! %%%     BASE_YEAR = 2004     for NOx                                   %
      ! %%%     BASE_YEAR = 2001     for CO in China                           %
      ! %%%     BASE_YEAR = 2000     for CO outside China, NH3, SO2, CH4 & CO2 %
      ! %%%     & VOC are not emitted -                                        %
      ! %%%                                                                    %
      ! %%%     & YEARLY SCALE FACTOR are applied to get 1985-2005 estimates   %
      ! %%%           of NOx, CO, SO2 if BASE_YEAR in 2000-4                   %
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !-------------------------------------------------------------------------

      ! select emissions year
      IF ( FSCALYR < 0 ) THEN
         SIM_YEAR = GET_YEAR()
      ELSE
         SIM_YEAR = FSCALYR
      ENDIF

      ! Pickup BASE_YEAR according to SIMulation YEAR
      IF ( SIM_YEAR >= 2006 ) THEN
         BASE_YEAR = 2006
      ELSE
         BASE_YEAR = 2000
      ENDIF

      ! set module flag
      IS_2006 = ( BASE_YEAR == 2006 )

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%%%    To simulate 2020 estimate, uncomment following two lines    %%%%
      !BASE_YEAR = 2020
      !IS_2006 = .TRUE.
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! define data directory and number of sources
      STREETS_DIR_2000 = TRIM( DATA_DIR_1x1 ) // 'Streets_200607/'

      IF ( IS_2006 ) THEN
         NTSOURCE1   = 4
         NTSOURCE2   = 10
         STREETS_DIR = TRIM( DATA_DIR_1x1 ) // 'Streets_200812/'
      ELSE
         NTSOURCE1   = 1
         NTSOURCE2   = 2
         STREETS_DIR = TRIM( DATA_DIR_1x1 ) // 'Streets_200607/'
      ENDIF


      !-----------------------------------------------------------------
      ! SOURCES = String array to identify emissions sources for 2006
      !           inventory. They correspond to :
      !
      !   Industry, Power, Residential, Transport for NOx/CO/SO2
      !
      !   Domestic Biofuel, Domestic Fossil Fuel, Domestic Non-Combustion,
      !   Industry, Power Plants, and Transportation for NMVOC
      !
      ! ONOFF = their orresponding switch
      !         ### MODIFY ONLY ONOFF  FOR SENSITIVITY STUDIES ##
      !-----------------------------------------------------------------
      SOURCES = (/ 'ind', 'pow', 'res', 'tra',                   ! for NOx/CO/SO2
     &             'dob', 'dof', 'dop', 'ind', 'pow', 'tra' /)   ! for VOC

      ONOFF = (/ 1D0, 1D0, 1D0, 1D0,
     &           1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)

      ! Corresponding scaling to get 2020 from 2006
      ! Note : first line only for NOx (no change in CO/SO2)
      IF ( BASE_YEAR == 2020 ) THEN
         SCALE2020 = (/ 2.36D0, 1.33D0, 1.02D0, 2.5D0,
     &                  1.02D0, 1.02D0, 1.02D0, 2.36D0, 1.33D0, 2.5D0 /)
      ELSE
         SCALE2020 = (/ 1D0, 1D0, 1D0, 1D0,
     &                  1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)
      ENDIF


      !-----------------------------------------------------------------
      !                      Test for simulation type
      !-----------------------------------------------------------------
      IF ( ITS_A_CH4_SIM() ) THEN

         !--------------------------
         ! Read CH4 and regrid
         ! (CH4 simulations only)
         !--------------------------

         ! File name
         FILENAME  = TRIM( STREETS_DIR ) //
     &              'Streets_CH4_FF_2000.generic.1x1'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS( FILENAME, 'CH4-EMIS', 1, TAU2000, CH4, 
     &                      IS_MASS=1 )

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------------
         ! Read CO2 and regrid
         ! (CO2 simulations only)
         !--------------------------

         ! File name
         FILENAME  = TRIM( STREETS_DIR ) //
     &              'Streets_CO2_FF_2000.generic.1x1'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS( FILENAME, 'CO2-SRCE', 1, TAU2000, CO2,
     $                      IS_MASS=1 )

      ELSE

         !--------------------------------------------------------------
         !                   Other simulations
         !--------------------------------------------------------------

         !--------------------------
         ! Read NOx and regrid
         !--------------------------
         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( IS_2006 ) THEN

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_NOx_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               TAU = TAU2006
            ELSE

!--- prior to 7/1/09 (has only Chinese data)
!
!               FILENAME  = TRIM( STREETS_DIR ) //
!     &              'Streets_NOx_FF_2004_monthly.generic.1x1'
!
!               TAU = TAUMONTH_2004
!

               FILENAME  = TRIM( STREETS_DIR ) //
     &              'Streets_NOx_FF_2000.generic.1x1'

               TAU = TAU2000

            ENDIF

            ! Read data
            CALL READ_STREETS( FILENAME, 'ANTHSRCE', 1, TAU, TEMP,
     &                         IS_MASS=1 )

            NOX = NOX + TEMP * ONOFF( NSOURCE ) * SCALE2020( NSOURCE )

         ENDDO


         !--------------------------
         ! Scale NOx
         !--------------------------

!----- prior to 7/1/09 (phs)
!         IF ( IS_2006 ) THEN
!
!            ! Monthly Variability for 2006 Base Year. Variability has
!            ! been obtained from 2004 data (phs, 12/2/08)
!            FILENAME  = TRIM( STREETS_DIR ) //
!     &                  'Streets_2004_NOx_MonthFctr_total.generic.1x1'
!
!            CALL READ_STREETS( FILENAME,     'RATIO-2D', 71,
!     $                         TAUMONTH_2004, TEMP,     'unitless' )
!
!            NOX = NOX * TEMP
!
!
!         ELSE
!
!            ! Annual scalar factor for NOx 2004 (amv, phs, 3/10/08)
!            CALL GET_ANNUAL_SCALAR( 71, 2004, SIM_YEAR, SCALFAC )
!
!            NOX = NOX * SCALFAC
!
!         ENDIF

         ! Annual scalar factor (phs, 3/10/08)
         !--------------------------
         IF ( BASE_YEAR == 2000 ) THEN

            CALL GET_ANNUAL_SCALAR( 71, 2000, SIM_YEAR, SCALFAC )

            NOX = NOX * SCALFAC
         ENDIF

         ! Seasonal Variation for NOx
         !--------------------------
         ! Monthly Variability for any year. Variability has
         ! been obtained from 2004 1x1 data, for two cases:
         ! FF seasonality for 2000, and TOTAL (BF+FF) for 2006
         ! inventories (phs, 12/2/08)

         IF ( IS_2006 ) THEN
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_2004_NOx_MonthFctr_total.generic.1x1'

         ELSE

            ! redefine the entire path here
            FILENAME = TRIM( DATA_DIR_1x1 ) // 'Streets_200812/'
     &                 // 'Streets_2004_NOx_MonthFctr_FF.generic.1x1'

         ENDIF

         CALL READ_STREETS( FILENAME,     'RATIO-2D', 71,
     $                      TAUMONTH_2004, TEMP, IS_MASS=0 )

         NOX = NOX * TEMP

         !--------------------------
         ! Read CO and scale CO
         !--------------------------

         ! Base year = 2006
         IF ( IS_2006 ) THEN

            DO NSOURCE = 1, NTSOURCE1

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CO_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 4,
     $                            TAU2006,  TEMP, IS_MASS=1 )

               ! No scaling for 2006-2020
               CO = CO + TEMP * ONOFF( NSOURCE )

            ENDDO

            ! Monthly Variability for 2006 Base Year. Variability has
            ! been obtained from 2001 data, and thus affects only China
            ! (phs, 12/2/08)
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_2001_CO_MonthFctr_total.generic.1x1'

            CALL READ_STREETS( FILENAME,     'RATIO-2D', 72,
     $                         TAUMONTH_2001, TEMP, IS_MASS=0 )

            CO = CO * TEMP

         ! Base year = 2000 (2001 for China)
         ELSE

            !-- PART 1 -- File name for 2000 CO over SE Asia
            FILENAME  = TRIM( STREETS_DIR ) //
     &                 'Streets_CO_FF_2000.generic.1x1'

            ! Read data
            CALL READ_STREETS( FILENAME, 'ANTHSRCE', 4, TAU2000, CO,
     &                         IS_MASS=1 )

            ! Annual scalar factor (amv, phs, 3/10/08)
            CALL GET_ANNUAL_SCALAR( 72, 2000, SIM_YEAR, SCALFAC )
            CO = CO * SCALFAC



            !-- PART 2 -- File name for 2001 CO over China only
            FILENAME  = TRIM( STREETS_DIR ) //
     &                 'Streets_CO_FF_2001_monthly.generic.1x1'

            ! Read data
            CALL READ_STREETS( FILENAME,      'ANTHSRCE', 4,
     $                         TAUMONTH_2001, TEMP, IS_MASS=1)

            ! Annual scalar factor (amv, phs, 3/10/08)
            CALL GET_ANNUAL_SCALAR( 72, 2001, SIM_YEAR, SCALFAC )
            TEMP = TEMP * SCALFAC


            !-- PART 3 -- Replace SE Asia CO for 2000 with China CO for 2001
            WHERE ( MASK_CHINA > 0 ) CO = TEMP

            ! switch and scale
            CO = CO * ONOFF( 1 )

         ENDIF


         !--------------------------
         ! Read SO2 and regrid
         !--------------------------
         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( BASE_YEAR == 2000 ) THEN
               FILENAME  = TRIM( STREETS_DIR ) //
     &                    'Streets_SO2_FF_2000.generic.1x1'

               TAU = TAU2000

            ELSE
               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_SO2_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               TAU = TAU2006

            ENDIF

            ! Read data
            CALL READ_STREETS( FILENAME, 'ANTHSRCE', 26, TAU, TEMP,
     &                         IS_MASS=1 )

            SO2 = SO2 + TEMP * ONOFF( NSOURCE )

         ENDDO

         ! Annual scalar factor (amv, phs, 3/10/08)
         IF ( .NOT. IS_2006  ) THEN
            CALL GET_ANNUAL_SCALAR( 73, 2000, SIM_YEAR, SCALFAC )
            SO2 = SO2 * SCALFAC
         ENDIF


         !---------------------------------------------
         ! Read NH3 only available for base year 2000
         !---------------------------------------------
!         IF ( IS_2006 ) THEN
!
!            NH3 = -1D0
!
!         ELSE

            ! File name
            FILENAME  = TRIM( STREETS_DIR_2000 ) //
     &                 'Streets_NH3_FF_2000.generic.1x1'

            ! Old file has NH3 as tracer #30
            CALL READ_STREETS( FILENAME, 'ANTHSRCE', 30, TAU2000, NH3,
     &                         IS_MASS=1 )

            ! switch and scale
            NH3 = NH3 * ONOFF( 1 )

!         ENDIF


         !---------------------------------------------
         ! Read VOC only if base year is 2006
         !---------------------------------------------
         IF ( IS_2006 ) THEN

            TAU = TAU2006

            !--------------------------
            ! Read ACET and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ACET_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 9, TAU, TEMP,
     &                            IS_MASS=1 )

               ACET = ACET + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read C2H6 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C2H6_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 21, TAU, TEMP,
     &                            IS_MASS=1 )

               C2H6 = C2H6 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read CH2O and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CH2O_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 20, TAU, TEMP,
     &                            IS_MASS=0 )

               CH2O = CH2O + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read C3H8 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C3H8_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 19, TAU, TEMP,
     &                            IS_MASS=1 )

               C3H8 = C3H8 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read PRPE and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_PRPE_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 18, TAU, TEMP,
     &                            IS_MASS=1 )

               PRPE = PRPE + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read ALD2 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALD2_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 11, TAU, TEMP,
     &                            IS_MASS=1 )

               ALD2 = ALD2 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO



            !--------------------------
            ! Read MEK and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_MEK_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 10, TAU, TEMP,
     &                            IS_MASS=1 )

               MEK = MEK + TEMP * ONOFF( NSOURCE )
     &                          * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read ALK4 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALK4_'//
     &                     SOURCES( NSOURCE ) // '_2006.generic.1x1'

               ! Read data [atom C/yr]
               CALL READ_STREETS( FILENAME, 'ANTHSRCE', 5, TAU, TEMP,
     &                            IS_MASS=1 )

               ALK4 = ALK4 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


         ELSE

            ! Set VOC to -1
            ALK4 = -1D0
            ACET = -1D0
            MEK  = -1D0
            PRPE = -1D0
            C2H6 = -1D0
            C3H8 = -1D0
            CH2O = -1D0
            ALD2 = -1D0


         ENDIF ! end VOCs

      ENDIF    ! end other simulations


      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN
         CALL STREETS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( SIM_YEAR, BASE_YEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS( FILENAME, CATEGORY, TRACERN, TAU, ARR,
     $                         IS_MASS )
!
!******************************************************************************
!     Subroutine READ_STREETS reads data from one STREETS data file
!     from disk, at GENERIC 1x1 resolution and regrids them to the
!     current model resolution.  (phs, 3/7/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of anthro or biomass file to read
!  (2 ) CATEGORY (CHARACTER) : Category name
!  (3 ) TRACERN  (INTEGER  ) : Tracer number
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2 ) ARR      (REAL*8   ) : Array to hold emissions
!
!  NOTES:
!  (1) UNIT argument in DO_REGRID_... is 'kg/yr' for VOCs
!        because 'atom C/yr' and 'molec/yr' are not recognized. The result is
!        still correct.
!  (2) Now inlcude seasonal scaling of NH3 emissions (jaf, 3/2/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0,      READ_BPCH2
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1
      USE TIME_MOD,       ONLY : GET_MONTH
      USE TRACERID_MOD,   ONLY : IDTNH3

      ! (lzh,02/01/2015) update regridding
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

#     include "CMN_SIZE"        ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME, CATEGORY
      INTEGER,          INTENT(IN)    :: TRACERN
      REAL*8,           INTENT(INOUT) :: ARR(IIPAR,JJPAR)
      REAL*8,           INTENT(IN)    :: TAU
      INTEGER,          INTENT(IN)    :: IS_MASS           ! For MAP_A2A regrid

      ! Local variables
      REAL*4            :: ARRAY(I1x1,J1x1-1,1)
      !REAL*8            :: GEN_1x1(I1x1,J1x1-1)
      REAL*8, TARGET    :: GEN_1x1(I1x1,J1x1-1)    !(lzh)
      REAL*8            :: GEOS_1x1(I1x1,J1x1,1)

      ! (lzh, 02/01/2015)
      CHARACTER(LEN=255)              :: LLFILENAME
      REAL*8, POINTER                 :: INGRID(:,:) => NULL()

      ! Variables for seasonal scaling of NH3 (jaf, 3/2/11)
      REAL*4            :: SCALAR_1x1(I1x1,J1x1-1,1)
      REAL*8            :: TAU1995
      INTEGER           :: RATIOID
      CHARACTER(LEN=250):: FILENAME_S

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )

      ARR   = 0d0
      ARRAY = 0.

      ! Read data
      CALL READ_BPCH2( FILENAME,  CATEGORY,  TRACERN,
     &                 TAU,       I1x1,      J1x1-1,
     &                 1,         ARRAY,     QUIET=.TRUE. )
      !=================================================================
      ! Apply seasonal variation to NH3 based on seasonality from
      ! Lex Bouwman.  Follow methodology in emep_mod.f (jaf, 3/2/11)
      !=================================================================

      ! Get TAU value for 1995, since the data is timestamped w/ this
      TAU1995 = GET_TAU0( GET_MONTH(), 1, 1995 )

! (lzh, 06/23/2014) disable NH3 seasonality
c$$$      ! For NH3 only ...
c$$$      IF ( TRACERN == 30 ) THEN
c$$$
c$$$         ! File name containing scaling factors
c$$$         FILENAME_S = TRIM( DATA_DIR_1x1 ) //
c$$$     &        'Streets_200607/NH3-Streets-SeasonalScalar.generic.1x1'
c$$$
c$$$         ! Tracer number for scale factor data
c$$$         RATIOID = 74
c$$$
c$$$         ! Echo info
c$$$         WRITE( 6, 101 ) TRIM( FILENAME_S )
c$$$101      FORMAT( '     - READ_STREETS: Reading ', a )
c$$$
c$$$         ! Read scaling factors
c$$$         CALL READ_BPCH2( FILENAME_S, 'RATIO-2D',  RATIOID,
c$$$     &                    TAU1995,     I1x1,       J1x1-1,
c$$$     &                    1,           SCALAR_1x1, QUIET=.TRUE. )
c$$$
c$$$         ! Apply seasonal scalar to NH3 emissions
c$$$         ARRAY(:,:,1) = ARRAY(:,:,1) * SCALAR_1x1(:,:,1)
c$$$
c$$$      ENDIF

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 --> GEOS 1x1
!      CALL DO_REGRID_G2G_1x1( THISUNIT, GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 --> current model resolution
!      CALL DO_REGRID_1x1( THISUNIT, GEOS_1x1, ARR )

      ! (lzh,02/01/2015)
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid from GENERIC 1x1 --> current model resolution
      INGRID => GEN_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1, J1x1-1,
     &                    INGRID,     ARR,  IS_MASS,
     &                    netCDF=.TRUE.              )

      ! Free pointer
      NULLIFY( INGRID )

      END SUBROUTINE READ_STREETS

!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS_05x0666( FILENAME, CATEGORY, TRACERN,
     $                                 TAU,      ARR )
!
!******************************************************************************
!     Subroutine READ_STREETS_05x0666 reads data from one STREETS data file
!     from disk, at 05x0666 resolution and cut them to the CHINA nested
!     window.  (phs, 12/2/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of anthro or biomass file to read
!  (2 ) CATEGORY (CHARACTER) : Category name
!  (3 ) TRACERN  (INTEGER  ) : Tracer number
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2 ) ARR      (REAL*8   ) : Array to hold emissions
!
!  NOTES:
!  (1) UNIT argument in DO_REGRID_... is 'kg/yr' for VOCs
!        because 'atom C/yr' and 'molec/yr' are not recognized. The result is
!        still correct.
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0,         READ_BPCH2
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_05X0666

#     include "CMN_SIZE"        ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME, CATEGORY
      INTEGER,          INTENT(IN)    :: TRACERN
      REAL*8,           INTENT(INOUT) :: ARR(IIPAR,JJPAR)
      REAL*8,           INTENT(IN)    :: TAU

      ! Local variables
      REAL*4            :: ARRAY(I05x0666,J05x0666,1)
      REAL*8            :: GEOS_05x0666(I05x0666,J05x0666,1)


      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME,  CATEGORY,  TRACERN,
     &                 TAU,       I05x0666,  J05x0666,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Cut to China nested simulation window
      CALL DO_REGRID_05x0666( 1,'kg/yr', GEOS_05x0666, ARR )


      END SUBROUTINE READ_STREETS_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_STREETS_ANTHRO_05x0666
!
!******************************************************************************
!  Subroutine EMISS_STREETS_ANTHRO_05x0666 reads the David Streets et al
!  emission fields at 0.5 x 0.666 resolution and regrids them to the current
!  nested-grid model resolution. (yxw, dan, bmy, 11/6/08)
!
!  NOTES:
!     (1 ) For now, disable the monthly CO emissions and just read the
!     same emissions as we do for the global simulations.  Update
!     emissions in a future release. (bmy, 11/6/08)
!     (2) Now read 2006 inventory (including VOCs) if needed. Apply monthly
!     variations for NOx
!     (3) Bug fixe : we call only read_streets_05x0666. (ccc, 3/11/09)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,        ONLY : GET_TAU0,         READ_BPCH2
      USE LOGICAL_MOD,      ONLY : LFUTURE
      USE TRACER_MOD,       ONLY : ITS_A_CO2_SIM,    ITS_A_CH4_SIM
      USE DIRECTORY_MOD,    ONLY : DATA_DIR
      USE TIME_MOD,         ONLY : GET_MONTH, GET_YEAR
      USE TRACER_MOD,       ONLY : XNUMOL
!      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR_05x0666

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_O3"         ! FSCALYR

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: BASE_YEAR,       SIM_YEAR
      CHARACTER(LEN=255)         :: FILENAME,        STREETS_DIR

      ! to loop over the sources
      INTEGER, PARAMETER         :: NSRCE = 10
      INTEGER                    :: NSOURCE, NTSOURCE1, NTSOURCE2
      CHARACTER(LEN=3)           :: SOURCES(NSRCE)   !!
      REAL*8                     :: ONOFF(NSRCE)     !! To switch off/on each sources
      REAL*8                     :: SCALE2020(NSRCE) !! To scale 2006 to 2020

      ! to hold temporary data and scale factors
      REAL*4                     :: SCALFAC( IIPAR, JJPAR )
      REAL*8                     :: TEMP( IIPAR, JJPAR )

      ! TAUs
      REAL*8                     :: TAU2000,       TAU2004,      TAU2006
      REAL*8                     :: TAUMONTH_2001, TAUMONTH_2004, TAU


      !=================================================================
      ! EMISS_STREETS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_STREETS_ANTHRO
         FIRST = .FALSE.
      ELSE
         NOx  = 0D0
         CO   = 0D0
         SO2  = 0D0
         ALK4 = 0D0
         ACET = 0D0
         MEK  = 0D0
         PRPE = 0D0
         C2H6 = 0D0
         C3H8 = 0D0
         CH2O = 0D0
         ALD2 = 0D0
      ENDIF

      ! TAU0 values for 2000, 2004 and 2006
      TAU2000 = GET_TAU0( 1, 1, 2000 )
      TAU2004 = GET_TAU0( 1, 1, 2004 )
      TAU2006 = GET_TAU0( 1, 1, 2006 )

      MONTH         = GET_MONTH()
      TAUMONTH_2001 = GET_TAU0( MONTH, 1, 2001 )
      TAUMONTH_2004 = GET_TAU0( MONTH, 1, 2004 )

      !-------------------------------------------------------------------------
      !     Base Year & Yearly Scale Factors
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% To simulate 2020 :                                                 %
      ! %%%         you must set BASE_YEAR = 2020 (hardwired, see below)       %
      ! %%%                                                                    %
      ! %%% To simulate 2006 and after, we use:                                %
      ! %%%      BASE_YEAR = 2006     for all species. NH3 is not emitted      %
      ! %%%                                                                    %
      ! %%% To simulate 2005 and before, we use:                               %
      ! %%%      BASE_YEAR = 2001   for CO in China                            %
      ! %%%      BASE_YEAR = 2000   for CO outside China, NOx, SO2, CH4, & CO2 %
      ! %%%     & VOC are not emitted -                                        %
      ! %%%                                                                    %
      ! %%% YEARLY SCALE FACTOR  (**** NOT AVAILABLE YET ****)                 %
      ! %%%      to be applied to get 1985-2005 estimates                      %
      ! %%%              of NOx, CO and SO2 if BASE_YEAR is 2000/1             %
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !-------------------------------------------------------------------------

      ! select emissions year
      IF ( FSCALYR < 0 ) THEN
         SIM_YEAR = GET_YEAR()
      ELSE
         SIM_YEAR = FSCALYR
      ENDIF

      ! Pickup BASE_YEAR according to SIMulation YEAR
      IF ( SIM_YEAR >= 2006 ) THEN
         IS_2006 = .TRUE.
         BASE_YEAR = 2006
      ELSE
         BASE_YEAR = 2000
      ENDIF

      ! %%%%  To simulate 2020 estimate, uncomment following line %%%%
      !BASE_YEAR = 2020
      !IS_2006 = .TRUE.

      ! define data directory and number of sources
      IF ( IS_2006 ) THEN
         NTSOURCE1   = 4
         NTSOURCE2   = 10
         STREETS_DIR = TRIM( DATA_DIR ) // 'Streets_200812/'
      ELSE
         NTSOURCE1   = 1
         NTSOURCE2   = 2
         STREETS_DIR = TRIM( DATA_DIR ) // 'Streets_200607/'
      ENDIF


      !-----------------------------------------------------------------
      ! SOURCES = String array to identify emissions sources for 2006
      !           inventory. They correspond to :
      !
      !   Industry, Power, Residential, Transport for NOx/CO/SO2
      !
      !   Domestic Biofuel, Domestic Fossil Fuel, Domestic Non-Combustion,
      !   Industry, Power Plants, and Transportation for NMVOC
      !
      ! ONOFF = their orresponding switch
      !         ### MODIFY ONLY ONOFF  FOR SENSITIVITY STUDIES ##
      !-----------------------------------------------------------------
      SOURCES = (/ 'ind', 'pow', 'res', 'tra',                   ! for NOx/CO/SO2
     &             'dob', 'dof', 'dop', 'ind', 'pow', 'tra' /)   ! for VOC

      ONOFF = (/ 1D0, 1D0, 1D0, 1D0,
     &           1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)

      ! Corresponding scaling to get 2020 from 2006
      ! Note : first line only for NOx (no change in CO/SO2)
      IF ( BASE_YEAR == 2020 ) THEN
         SCALE2020 = (/ 2.36D0, 1.33D0, 1.02D0, 2.5D0,
     &                  1.02D0, 1.02D0, 1.02D0, 2.36D0, 1.33D0, 2.5D0 /)
      ELSE
         SCALE2020 = (/ 1D0, 1D0, 1D0, 1D0,
     &                  1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)
      ENDIF


      !-----------------------------------------------------------------
      !                      Test for simulation type
      !-----------------------------------------------------------------
      IF ( ITS_A_CH4_SIM() ) THEN

         !--------------------------
         ! Read CH4 and regrid
         ! (CH4 simulations only)
         !--------------------------

         ! File name for 2000 CO over SE Asia
         FILENAME  = TRIM( STREETS_DIR ) //
     &             'Streets_CH4_FF_2000.geos5.05x0666'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS_05x0666( FILENAME, 'CH4-EMIS', 1,
     $                              TAU2000,   CH4 )


      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------------
         ! Read CO2 and regrid
         ! (CH2 simulations only)
         !--------------------------

         ! File name for 2000 CO over SE Asia
         FILENAME  = TRIM( STREETS_DIR ) //
     &             'Streets_CO2_FF_2000.geos5.05x0666'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS_05x0666( FILENAME, 'CO2-SRCE', 1,
     $                              TAU2000,   CO2 )

      ELSE

         !--------------------------------------------------------------
         !                   Other simulations
         !--------------------------------------------------------------

         !--------------------------
         ! Read NOx
         !--------------------------
         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( BASE_YEAR == 2000 ) THEN

               ! File name for 2000 CO over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) //
     &                     'Streets_NOx_FF_2000.geos5.05x0666'

               TAU = TAU2000

            ELSE

               ! File name for 2000 NOx over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_NOx_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               TAU = TAU2006

            ENDIF

            ! Read data
            CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 1,
     $                                 TAU,       TEMP )

            NOX = NOX + TEMP * ONOFF( NSOURCE ) * SCALE2020( NSOURCE )

         ENDDO

!------------------------------------------------------------------------
! Not available yet
!         ! Annual scalar factor (phs, 3/10/08)
!         IF ( BASE_YEAR == 2000 ) THEN
!            CALL GET_ANNUAL_SCALAR_05x0666( 71,       2000,
!     &                                      SIM_YEAR, SCALFAC )
!            NOX = NOX * SCALFAC
!         ENDIF
!------------------------------------------------------------------------
         !--------------------------
         ! Seasonal Variation for NOx
         !--------------------------

         ! Monthly Variability for any year. Variability has
         ! been obtained from 2004 1x1 data, for two cases:
         ! FF seasonality for 2000, and TOTAL (BF+FF) for 2006
         ! inventories (phs, 12/2/08)

         IF ( IS_2006 ) THEN
            FILENAME  = TRIM( STREETS_DIR ) //
     &           'Streets_2004_NOx_MonthFctr_total.geos5.05x0666'

         ELSE

            ! we need to redefine the entire path here
            !-----------------------------------------
            FILENAME = TRIM( DATA_DIR ) // 'Streets_200812/'
     &           // 'Streets_2004_NOx_MonthFctr_FF.geos5.05x0666'

         ENDIF

         CALL READ_STREETS_05x0666( FILENAME, 'RATIO-2D',
     $                              71,        TAUMONTH_2004, TEMP )

         NOX = NOX * TEMP


         !--------------------------
         ! Read CO 2006 (SE Asia)
         !--------------------------
         IF ( IS_2006 ) THEN

            DO NSOURCE = 1, NTSOURCE1

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CO_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 4,
     $                                    TAU2006,   TEMP )

               ! No scaling for 2006-2020
               CO = CO + TEMP * ONOFF( NSOURCE )

            ENDDO

            ! Monthly Variability for 2006 BF+FF. Variability has
            ! been obtained from 2001 05x0666 data, and like those
            ! those data, only China features variability (phs, 12/2/08)
            FILENAME  = TRIM( STREETS_DIR ) //
     &                 'Streets_2001_CO_MonthFctr_total.geos5.05x0666'

            CALL READ_STREETS_05x0666( FILENAME,      'RATIO-2D',  72,
     $                                 TAUMONTH_2001,  TEMP     )

            CO = CO * TEMP

         !------------------------------
         ! Read CO 2000 (2001 for China)
         !-----------------------------
         ELSE

            !-- PART 1 -- File name for 2000 CO over SE Asia
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_CO_FF_2000.geos5.05x0666'

            ! Read data
            CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 4,
     $                                 TAU2000,   CO )

!------------------------------------------------------------------------
! Not available yet
!            CALL GET_ANNUAL_SCALAR_05x0666( 72,       2000,
!     &                                      SIM_YEAR, SCALFAC )
!            CO = CO * SCALFAC
!------------------------------------------------------------------------


            !-- PART 2 -- File name for 2001 CO over China only
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_2001CO_monthly_ff.geos5.05x0666'

            ! Read data
            CALL READ_STREETS_05x0666( FILENAME,      'ANTHSRCE', 4,
     $                                 TAUMONTH_2001, TEMP )

!------------------------------------------------------------------------
! Not available yet
!            CALL GET_ANNUAL_SCALAR_05x0666( 72,       2001,
!     &                                      SIM_YEAR, SCALFAC )
!            TEMP = TEMP * SCALFAC
!------------------------------------------------------------------------


            !-- PART 3 -- Replace SE Asia CO for 2000 with China CO for 2001

            WHERE ( MASK_CHINA > 0 ) CO = TEMP


         ENDIF



         !--------------------------
         ! Read SO2 and regrid
         !--------------------------

         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( BASE_YEAR == 2000 ) THEN

               ! File name for 2000 SO2 over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) //
     &                    'Streets_SO2_FF_2000.geos5.05x0666'

               TAU = TAU2000

            ELSE

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_SO2_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               TAU = TAU2006

            ENDIF

            ! Read data
            CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 26,
     $                                 TAU,       TEMP )

            SO2 = SO2 + TEMP * ONOFF( NSOURCE )

         ENDDO

!------------------------------------------------------------------------
! Not available yet
!         ! Annual scalar factor (amv, phs, 3/10/08)
!         IF ( .NOT. IS_2006  ) THEN
!            CALL GET_ANNUAL_SCALAR_05x0666( 73,       2000,
!     $                                      SIM_YEAR, SCALFAC )
!            SO2 = SO2 * SCALFAC
!         ENDIF
!------------------------------------------------------------------------



         !---------------------------------------------
         ! Read NH3 only if base year is 2000
         !---------------------------------------------
         IF ( IS_2006 ) THEN

            NH3 = -1D0

         ELSE

            ! File name
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_NH3_FF_2000.geos5.05x0666'

            ! Old file has NH3 as tracer #30
            CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 30,
     $                                 TAU2000,   NH3 )

            ! switch and scale
            NH3 = NH3 * ONOFF( 1 )

         ENDIF



         !---------------------------------------------
         ! Read VOC only if base year is 2006
         !---------------------------------------------
         IF ( IS_2006 ) THEN

            TAU = TAU2006

            !--------------------------
            ! Read ACET and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ACET_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 9,
     &                                    TAU, TEMP )

               ACET = ACET + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read C2H6 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C2H6_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 21,
     &                                    TAU, TEMP )

               C2H6 = C2H6 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read CH2O and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CH2O_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 20,
     &                                    TAU, TEMP )

               CH2O = CH2O + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read C3H8 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C3H8_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 19,
     &                                    TAU, TEMP )

               C3H8 = C3H8 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read PRPE and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_PRPE_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 18,
     &                                    TAU, TEMP )

               PRPE = PRPE + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read ALD2 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALD2_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 11,
     &                                    TAU, TEMP )

               ALD2 = ALD2 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO



            !--------------------------
            ! Read MEK and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_MEK_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 10,
     &                                    TAU, TEMP )

               MEK = MEK + TEMP * ONOFF( NSOURCE )
     &                          * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read ALK4 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALK4_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.05x0666'

               ! Read data [atom C/yr]
               CALL READ_STREETS_05x0666( FILENAME, 'ANTHSRCE', 5,
     &                                    TAU, TEMP )

               ALK4 = ALK4 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO


         ELSE

            ! Set VOC to -1
            ALK4 = -1D0
            ACET = -1D0
            MEK  = -1D0
            PRPE = -1D0
            C2H6 = -1D0
            C3H8 = -1D0
            CH2O = -1D0
            ALD2 = -1D0


         ENDIF ! end VOCs

      ENDIF    ! end other simulations


!------------------------------------------------------------------------
! Not available yet
!      !--------------------------
!      ! Compute future emissions
!      !--------------------------
!      IF ( LFUTURE ) THEN
!         CALL STREETS_SCALE_FUTURE
!      ENDIF
!------------------------------------------------------------------------

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( SIM_YEAR, BASE_YEAR )

      END SUBROUTINE EMISS_STREETS_ANTHRO_05x0666

!-----------------------------------------------------------------------------
!===== (lzh, 04/25/2014) =====
!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS_025x03125( FILENAME, CATEGORY, TRACERN,
     $                                 TAU,      ARR )
!
!******************************************************************************
!     Subroutine READ_STREETS_025x03125 reads data from one STREETS data file
!     from disk, at 025x03125 resolution and cut them to the CHINA nested
!     window.  (phs, 12/2/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of anthro or biomass file to read
!  (2 ) CATEGORY (CHARACTER) : Category name
!  (3 ) TRACERN  (INTEGER  ) : Tracer number
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2 ) ARR      (REAL*8   ) : Array to hold emissions
!
!  NOTES:
!  (1) UNIT argument in DO_REGRID_... is 'kg/yr' for VOCs
!        because 'atom C/yr' and 'molec/yr' are not recognized. The result is
!        still correct.
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0,         READ_BPCH2
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_025X03125

#     include "CMN_SIZE"        ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME, CATEGORY
      INTEGER,          INTENT(IN)    :: TRACERN
      REAL*8,           INTENT(INOUT) :: ARR(IIPAR,JJPAR)
      REAL*8,           INTENT(IN)    :: TAU

      ! Local variables
      REAL*4            :: ARRAY(I025x031,J025x031,1)
      REAL*8            :: GEOS_025x03125(I025x031,J025x031,1)

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME,  CATEGORY,  TRACERN,
     &                 TAU,       I025x031,  J025x031,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_025x03125(:,:,1) = ARRAY(:,:,1)

      ! Cut to China nested simulation window
      CALL DO_REGRID_025x03125( 1,'kg/yr', GEOS_025x03125, ARR )

      END SUBROUTINE READ_STREETS_025x03125

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_STREETS_ANTHRO_025x03125
!
!******************************************************************************
!  Subroutine EMISS_STREETS_ANTHRO_05x0666 reads the David Streets et al
!  emission fields at 0.5 x 0.666 resolution and regrids them to the current
!  nested-grid model resolution. (yxw, dan, bmy, 11/6/08)
!
!  NOTES:
!     (1 ) For now, disable the monthly CO emissions and just read the
!     same emissions as we do for the global simulations.  Update
!     emissions in a future release. (bmy, 11/6/08)
!     (2) Now read 2006 inventory (including VOCs) if needed. Apply monthly
!     variations for NOx
!     (3) Bug fixe : we call only read_streets_05x0666. (ccc, 3/11/09)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,        ONLY : GET_TAU0,         READ_BPCH2
      USE LOGICAL_MOD,      ONLY : LFUTURE
      USE TRACER_MOD,       ONLY : ITS_A_CO2_SIM,    ITS_A_CH4_SIM
      USE DIRECTORY_MOD,    ONLY : DATA_DIR
      USE TIME_MOD,         ONLY : GET_MONTH, GET_YEAR
      USE TRACER_MOD,       ONLY : XNUMOL

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_O3"         ! FSCALYR

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: BASE_YEAR,       SIM_YEAR
      CHARACTER(LEN=255)         :: FILENAME,        STREETS_DIR

      ! to loop over the sources
      INTEGER, PARAMETER         :: NSRCE = 10
      INTEGER                    :: NSOURCE, NTSOURCE1, NTSOURCE2
      CHARACTER(LEN=3)           :: SOURCES(NSRCE)   !!
      REAL*8                     :: ONOFF(NSRCE)     !! To switch off/on each sources
      REAL*8                     :: SCALE2020(NSRCE) !! To scale 2006 to 2020

      ! to hold temporary data and scale factors
      REAL*4                     :: SCALFAC( IIPAR, JJPAR )
      REAL*8                     :: TEMP( IIPAR, JJPAR )

      ! TAUs
      REAL*8                     :: TAU2000,       TAU2004,      TAU2006
      REAL*8                     :: TAUMONTH_2001, TAUMONTH_2004, TAU

      !=================================================================
      ! EMISS_STREETS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_STREETS_ANTHRO
         FIRST = .FALSE.
      ELSE
         NOx  = 0D0
         CO   = 0D0
         SO2  = 0D0
         ALK4 = 0D0
         ACET = 0D0
         MEK  = 0D0
         PRPE = 0D0
         C2H6 = 0D0
         C3H8 = 0D0
         CH2O = 0D0
         ALD2 = 0D0

         ! (lzh, 09/15/2014)
         NH3  = 0D0
      ENDIF

      ! TAU0 values for 2000, 2004 and 2006
      TAU2000 = GET_TAU0( 1, 1, 2000 )
      TAU2004 = GET_TAU0( 1, 1, 2004 )
      TAU2006 = GET_TAU0( 1, 1, 2006 )

      MONTH         = GET_MONTH()
      TAUMONTH_2001 = GET_TAU0( MONTH, 1, 2001 )
      TAUMONTH_2004 = GET_TAU0( MONTH, 1, 2004 )

      !-------------------------------------------------------------------------
      !     Base Year & Yearly Scale Factors
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% To simulate 2020 :                                                 %
      ! %%%         you must set BASE_YEAR = 2020 (hardwired, see below)       %
      ! %%%                                                                    %
      ! %%% To simulate 2006 and after, we use:                                %
      ! %%%      BASE_YEAR = 2006     for all species. NH3 is not emitted      %
      ! %%%                                                                    %
      ! %%% To simulate 2005 and before, we use:                               %
      ! %%%      BASE_YEAR = 2001   for CO in China                            %
      ! %%%      BASE_YEAR = 2000   for CO outside China, NOx, SO2, CH4, & CO2 %
      ! %%%     & VOC are not emitted -                                        %
      ! %%%                                                                    %
      ! %%% YEARLY SCALE FACTOR  (**** NOT AVAILABLE YET ****)                 %
      ! %%%      to be applied to get 1985-2005 estimates                      %
      ! %%%              of NOx, CO and SO2 if BASE_YEAR is 2000/1             %
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !-------------------------------------------------------------------------

      ! select emissions year
      IF ( FSCALYR < 0 ) THEN
         SIM_YEAR = GET_YEAR()
      ELSE
         SIM_YEAR = FSCALYR
      ENDIF

      ! Pickup BASE_YEAR according to SIMulation YEAR
      IF ( SIM_YEAR >= 2006 ) THEN
         IS_2006 = .TRUE.
         BASE_YEAR = 2006
      ELSE
         BASE_YEAR = 2000
      ENDIF

      ! %%%%  To simulate 2020 estimate, uncomment following line %%%%
      !BASE_YEAR = 2020
      !IS_2006 = .TRUE.

      ! define data directory and number of sources
      IF ( IS_2006 ) THEN
         NTSOURCE1   = 4
         NTSOURCE2   = 10
         STREETS_DIR = TRIM( DATA_DIR ) // 'Streets_200812/'
      ELSE
         NTSOURCE1   = 1
         NTSOURCE2   = 2
         STREETS_DIR = TRIM( DATA_DIR ) // 'Streets_200607/'
      ENDIF

      !-----------------------------------------------------------------
      ! SOURCES = String array to identify emissions sources for 2006
      !           inventory. They correspond to :
      !
      !   Industry, Power, Residential, Transport for NOx/CO/SO2
      !
      !   Domestic Biofuel, Domestic Fossil Fuel, Domestic Non-Combustion,
      !   Industry, Power Plants, and Transportation for NMVOC
      !
      ! ONOFF = their orresponding switch
      !         ### MODIFY ONLY ONOFF  FOR SENSITIVITY STUDIES ##
      !-----------------------------------------------------------------
      SOURCES = (/ 'ind', 'pow', 'res', 'tra',                   ! for NOx/CO/SO2
     &             'dob', 'dof', 'dop', 'ind', 'pow', 'tra' /)   ! for VOC

      ONOFF = (/ 1D0, 1D0, 1D0, 1D0,
     &           1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)

      ! Corresponding scaling to get 2020 from 2006
      ! Note : first line only for NOx (no change in CO/SO2)
      IF ( BASE_YEAR == 2020 ) THEN
         SCALE2020 = (/ 2.36D0, 1.33D0, 1.02D0, 2.5D0,
     &                  1.02D0, 1.02D0, 1.02D0, 2.36D0, 1.33D0, 2.5D0 /)
      ELSE
         SCALE2020 = (/ 1D0, 1D0, 1D0, 1D0,
     &                  1D0, 1D0, 1D0, 1D0, 1D0, 1D0 /)
      ENDIF

      !-----------------------------------------------------------------
      !                      Test for simulation type
      !-----------------------------------------------------------------
      IF ( ITS_A_CH4_SIM() ) THEN

         !--------------------------
         ! Read CH4 and regrid
         ! (CH4 simulations only)
         !--------------------------

         ! File name for 2000 CO over SE Asia
         FILENAME  = TRIM( STREETS_DIR ) //
     &             'Streets_CH4_FF_2000.geos5.025x03125'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS_025x03125( FILENAME, 'CH4-EMIS', 1,
     $                              TAU2000,   CH4 )


      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------------
         ! Read CO2 and regrid
         ! (CH2 simulations only)
         !--------------------------

         ! File name for 2000 CO over SE Asia
         FILENAME  = TRIM( STREETS_DIR ) //
     &             'Streets_CO2_FF_2000.geos5.025x03125'

         BASE_YEAR = 2000

         ! Read data
         CALL READ_STREETS_025x03125( FILENAME, 'CO2-SRCE', 1,
     $                              TAU2000,   CO2 )

      ELSE

         !--------------------------------------------------------------
         !                   Other simulations
         !--------------------------------------------------------------

         !--------------------------
         ! Read NOx
         !--------------------------
         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( BASE_YEAR == 2000 ) THEN

               ! File name for 2000 CO over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) //
     &                     'Streets_NOx_FF_2000.geos5.025x03125'

               TAU = TAU2000

            ELSE

               ! File name for 2000 NOx over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_NOx_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               TAU = TAU2006

            ENDIF

            ! Read data
            CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 1,
     $                                 TAU,       TEMP )

            NOX = NOX + TEMP * ONOFF( NSOURCE ) * SCALE2020( NSOURCE )

         ENDDO

!------------------------------------------------------------------------
! Not available yet
!         ! Annual scalar factor (phs, 3/10/08)
!         IF ( BASE_YEAR == 2000 ) THEN
!            CALL GET_ANNUAL_SCALAR_05x0666( 71,       2000,
!     &                                      SIM_YEAR, SCALFAC )
!            NOX = NOX * SCALFAC
!         ENDIF
!------------------------------------------------------------------------
c$$$         !--------------------------
c$$$         ! Seasonal Variation for NOx
c$$$         !--------------------------
c$$$
c$$$         ! Monthly Variability for any year. Variability has
c$$$         ! been obtained from 2004 1x1 data, for two cases:
c$$$         ! FF seasonality for 2000, and TOTAL (BF+FF) for 2006
c$$$         ! inventories (phs, 12/2/08)
c$$$
c$$$         IF ( IS_2006 ) THEN
c$$$            FILENAME  = TRIM( STREETS_DIR ) //
c$$$     &           'Streets_2004_NOx_MonthFctr_total.geos5.05x0666'
c$$$
c$$$         ELSE
c$$$
c$$$            ! we need to redefine the entire path here
c$$$            !-----------------------------------------
c$$$            FILENAME = TRIM( DATA_DIR ) // 'Streets_200812/'
c$$$     &           // 'Streets_2004_NOx_MonthFctr_FF.geos5.05x0666'
c$$$
c$$$         ENDIF
c$$$
c$$$         CALL READ_STREETS_05x0666( FILENAME, 'RATIO-2D',
c$$$     $                              71,        TAUMONTH_2004, TEMP )
c$$$
c$$$         NOX = NOX * TEMP

         !--------------------------
         ! Read CO 2006 (SE Asia)
         !--------------------------
         IF ( IS_2006 ) THEN

            DO NSOURCE = 1, NTSOURCE1

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CO_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 4,
     $                                    TAU2006,   TEMP )

               ! No scaling for 2006-2020
               CO = CO + TEMP * ONOFF( NSOURCE )

            ENDDO

            ! Monthly Variability for 2006 BF+FF. Variability has
            ! been obtained from 2001 05x0666 data, and like those
            ! those data, only China features variability (phs, 12/2/08)
            FILENAME  = TRIM( STREETS_DIR ) //
     &                 'Streets_2001_CO_MonthFctr_total.geos5.025x03125'

            CALL READ_STREETS_025x03125( FILENAME,      'RATIO-2D',  72,
     $                                 TAUMONTH_2001,  TEMP     )

            CO = CO * TEMP

         !------------------------------
         ! Read CO 2000 (2001 for China)
         !-----------------------------
         ELSE

            !-- PART 1 -- File name for 2000 CO over SE Asia
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_CO_FF_2000.geos5.025x03125'

            ! Read data
            CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 4,
     $                                 TAU2000,   CO )

!------------------------------------------------------------------------
! Not available yet
!            CALL GET_ANNUAL_SCALAR_05x0666( 72,       2000,
!     &                                      SIM_YEAR, SCALFAC )
!            CO = CO * SCALFAC
!------------------------------------------------------------------------


            !-- PART 2 -- File name for 2001 CO over China only
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_2001CO_monthly_ff.geos5.025x03125'

            ! Read data
            CALL READ_STREETS_025x03125( FILENAME,      'ANTHSRCE', 4,
     $                                 TAUMONTH_2001, TEMP )

!------------------------------------------------------------------------
! Not available yet
!            CALL GET_ANNUAL_SCALAR_05x0666( 72,       2001,
!     &                                      SIM_YEAR, SCALFAC )
!            TEMP = TEMP * SCALFAC
!------------------------------------------------------------------------


            !-- PART 3 -- Replace SE Asia CO for 2000 with China CO for 2001

            WHERE ( MASK_CHINA > 0 ) CO = TEMP


         ENDIF

         !--------------------------
         ! Read SO2 and regrid
         !--------------------------

         DO NSOURCE = 1, NTSOURCE1

            ! File name
            IF ( BASE_YEAR == 2000 ) THEN

               ! File name for 2000 SO2 over SE Asia
               FILENAME  = TRIM( STREETS_DIR ) //
     &                    'Streets_SO2_FF_2000.geos5.025x03125'

               TAU = TAU2000

            ELSE

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_SO2_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               TAU = TAU2006

            ENDIF

            ! Read data
            CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 26,
     $                                 TAU,       TEMP )

            SO2 = SO2 + TEMP * ONOFF( NSOURCE )

         ENDDO

!------------------------------------------------------------------------
! Not available yet
!         ! Annual scalar factor (amv, phs, 3/10/08)
!         IF ( .NOT. IS_2006  ) THEN
!            CALL GET_ANNUAL_SCALAR_05x0666( 73,       2000,
!     $                                      SIM_YEAR, SCALFAC )
!            SO2 = SO2 * SCALFAC
!         ENDIF
!------------------------------------------------------------------------

         !---------------------------------------------
         ! Read NH3 only if base year is 2000
         !---------------------------------------------
         IF ( IS_2006 ) THEN

            NH3 = -1D0

         ELSE

            ! File name
            FILENAME  = TRIM( STREETS_DIR ) //
     &                  'Streets_NH3_FF_2000.geos5.025x03125'

            ! Old file has NH3 as tracer #30
            CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 30,
     $                                 TAU2000,   NH3 )

            ! switch and scale
            NH3 = NH3 * ONOFF( 1 )

         ENDIF

         !---------------------------------------------
         ! Read VOC only if base year is 2006
         !---------------------------------------------
         IF ( IS_2006 ) THEN

            TAU = TAU2006

            !--------------------------
            ! Read ACET and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ACET_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 9,
     &                                    TAU, TEMP )

               ACET = ACET + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read C2H6 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C2H6_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 21,
     &                                    TAU, TEMP )

               C2H6 = C2H6 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read CH2O and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_CH2O_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 20,
     &                                    TAU, TEMP )

               CH2O = CH2O + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read C3H8 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_C3H8_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 19,
     &                                    TAU, TEMP )

               C3H8 = C3H8 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read PRPE and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_PRPE_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 18,
     &                                    TAU, TEMP )

               PRPE = PRPE + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read ALD2 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALD2_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 11,
     &                                    TAU, TEMP )

               ALD2 = ALD2 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

            !--------------------------
            ! Read MEK and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_MEK_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 10,
     &                                    TAU, TEMP )

               MEK = MEK + TEMP * ONOFF( NSOURCE )
     &                          * SCALE2020( NSOURCE )
            ENDDO


            !--------------------------
            ! Read ALK4 and regrid
            !--------------------------
            DO NSOURCE = NTSOURCE1 + 1, NTSOURCE2

               FILENAME  = TRIM( STREETS_DIR ) // 'Streets_ALK4_'//
     &                     SOURCES( NSOURCE ) // '_2006.geos5.025x03125'

               ! Read data [atom C/yr]
               CALL READ_STREETS_025x03125( FILENAME, 'ANTHSRCE', 5,
     &                                    TAU, TEMP )

               ALK4 = ALK4 + TEMP * ONOFF( NSOURCE )
     &                            * SCALE2020( NSOURCE )
            ENDDO

         ELSE

            ! Set VOC to -1
            ALK4 = -1D0
            ACET = -1D0
            MEK  = -1D0
            PRPE = -1D0
            C2H6 = -1D0
            C3H8 = -1D0
            CH2O = -1D0
            ALD2 = -1D0


         ENDIF ! end VOCs

      ENDIF    ! end other simulations


!------------------------------------------------------------------------
! Not available yet
!      !--------------------------
!      ! Compute future emissions
!      !--------------------------
!      IF ( LFUTURE ) THEN
!         CALL STREETS_SCALE_FUTURE
!      ENDIF
!------------------------------------------------------------------------

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( SIM_YEAR, BASE_YEAR )

      END SUBROUTINE EMISS_STREETS_ANTHRO_025x03125

!-----------------------------------------------------------------------------

      SUBROUTINE STREETS_SCALE_FUTURE
!
!******************************************************************************
!  Subroutine STREETS_SCALE_FUTURE applies the IPCC future scale factors to
!  the David Streets' anthropogenic emissions. (swu, bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      INTEGER                       :: I, J

      !=================================================================
      ! STREETS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /yr]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/yr]
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future SO2 [kg SO2/yr]
         NH3(I,J)  = NH3(I,J) * GET_FUTURE_SCALE_NH3an( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE STREETS_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG ( YEAR, BASE_YEAR )
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the totals for the anthropogenic
!  emissions of NOx and CO. (bmy, 8/16/06)
!
!  NOTES:
!     (1 ) Now both simulation and base years are input. Output totals in
!          Tg/month instead of Tg/yr, except for CO2 and CH4 offline
!          simulations (phs, 12/9/08)
!     (2 ) Updated information output. Account for NH3 2000 used for all
!          simulation years (phs, 2/27/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : ITS_A_CH4_SIM, ITS_A_CO2_SIM

#     include "CMN_SIZE"   ! Size parameters

      ! argument
      INTEGER, INTENT(IN) :: YEAR, BASE_YEAR

      ! Local variables
      INTEGER             :: I,     J
      REAL*8              :: T_NOX, T_CO,  T_SO2
      REAL*8              :: T_NH3, T_CH4, T_CO2
      REAL*8              :: T_ACET, T_ALD2, T_ALK4, T_C2H6
      REAL*8              :: T_C3H8, T_CH2O, T_MEK,  T_PRPE
      REAL*8              :: AFACTOR

      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  ) BASE_YEAR
 100  FORMAT( 'M O N T H L Y   S T R E E T S   A S I A N',
     $     '   E M I S S I O N S', /, 'Scaled from base year : ', i4)


      ! Test for simulation type
      IF ( ITS_A_CH4_SIM() ) THEN

         !-----------------------
         ! CH4 simulation
         !-----------------------

         ! Total CH4 [Tg CH4]
         T_CH4 = SUM( CH4 ) * 1d-9

         ! Print totals
         WRITE( 6, 120 ) 'CH4 ', 2000, T_NOx,  ' CH4'

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !-----------------------
         ! CO2 simulation
         !-----------------------

         ! Total CO2 [Tg CO2]
         T_CH4 = SUM( CO2 ) * 1d-9

         ! Print totals
         WRITE( 6, 120 ) 'CO2 ', 2000, T_NOx,  ' CO2'

      ELSE

         !-----------------------
         ! Other simulations
         !-----------------------
#if   !defined( GRID05x0666 )
         IF ( .NOT. IS_2006 ) THEN
            WRITE( 6, * ) 'NOTES: '
            WRITE( 6, * ) '(1) Base year for NOx : 2004'
            WRITE( 6, * ) '(2) Annual scale factors applied to' //
     $                       ' NOx, CO & SO2'
            WRITE( 6, * ) '(3) Monthly variations applied to NOx & CO'
         ELSE
            WRITE( 6, * ) 'NOTES: '
            WRITE( 6, * ) '(1) Include ANTH and BIOFUEL'
            WRITE( 6, * ) '(2) Base year for NH3 : 2000'
            WRITE( 6, * ) '(3) Monthly variations applied to NOx & CO'
         ENDIF
#endif

         ! Total NOx [Tg N]
         T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / ( 12d0 * 46d0 ) )

         ! Total CO  [Tg CO]
         T_CO  = SUM( CO  ) * 1d-9 / 12d0

         ! Total SO2 [Tg S]
         T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / ( 12d0 * 64d0 ) )

         ! Total NH3 [Tg NH3]
         T_NH3 = SUM( NH3 ) * 1d-9 / 12d0

         IF ( IS_2006 ) THEN

            AFACTOR = 12d-12 / 6.0225d23 ! for C atom, units=Tg/yr

            AFACTOR = AFACTOR / 12d0 ! convert from Tg/yr to Tg/month

            T_ACET = SUM(ACET)  * AFACTOR
            T_ALD2 = SUM(ALD2)  * AFACTOR
            T_ALK4 = SUM(ALK4)  * AFACTOR
            T_C2H6 = SUM(C2H6)  * AFACTOR
            T_C3H8 = SUM(C3H8)  * AFACTOR
            T_CH2O = SUM(CH2O)  * 30d-12 / (6.0225d23 * 12d0)
            T_MEK  = SUM(MEK)   * AFACTOR
            T_PRPE = SUM(PRPE)  * AFACTOR

!-- prior 2/27/09
!         ELSE
!
!            ! Total NH3 [Tg NH3]
!            T_NH3 = SUM( NH3 ) * 1d-9 / 12d0

         ENDIF


         ! Print totals in [kg/month]
         WRITE( 6, 110 ) 'NOx  ', YEAR, MONTH, T_NOx,  '[Tg N   ]'
         WRITE( 6, 110 ) 'CO   ', YEAR, MONTH, T_CO,   '[Tg CO  ]'
         WRITE( 6, 110 ) 'SO2  ', YEAR, MONTH, T_SO2,  '[Tg S   ]'
         WRITE( 6, 110 ) 'NH3  ', YEAR, MONTH, T_NH3,  '[Tg NH3 ]'

         IF ( IS_2006 ) THEN

            WRITE( 6, 110 ) 'ALK4 ', YEAR, MONTH, T_ALK4, '[Tg C   ]'
            WRITE( 6, 110 ) 'ACET ', YEAR, MONTH, T_ACET, '[Tg C   ]'
            WRITE( 6, 110 ) 'MEK  ', YEAR, MONTH, T_MEK,  '[Tg C   ]'
            WRITE( 6, 110 ) 'PRPE ', YEAR, MONTH, T_PRPE, '[Tg C   ]'
            WRITE( 6, 110 ) 'C3H8 ', YEAR, MONTH, T_C3H8, '[Tg C   ]'
            WRITE( 6, 110 ) 'CH2O ', YEAR, MONTH, T_CH2O, '[Tg Ch2O]'
            WRITE( 6, 110 ) 'C2H6 ', YEAR, MONTH, T_C2H6, '[Tg C   ]'
            WRITE( 6, 110 ) 'ALD2 ', YEAR, MONTH, T_ALD2, '[Tg C   ]'

!-- prior 2/27/09
!         ELSE
!
!            WRITE( 6, 110 ) 'NH3  ', YEAR, MONTH, T_NH3,  '[Tg NH3 ]'

         ENDIF

      ENDIF


      ! Format statement
 110  FORMAT( 'David Streets anthro ', a5, 'for year ', i4,
     $     ' and month ', i2.2 ,': ', f11.4, 1x, a9 )

 120  FORMAT( 'David Streets anthro ', a5, 'for year ', i4,
     $     ': ', f11.4, 1x, a9 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg

!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS_MASKS
!
!******************************************************************************
!  Subroutine READ_STREETS_MASKS reads and regrids the China and SE Asia masks
!  that define the David Streets' emission regions (bmy, 8/16/06, 9/5/06)
!
!  NOTES:
!  (1 ) Now also save 1x1 CHINA MASK for use in other routines. (bmy, 9/5/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1

      ! (lzh,02/01/2015) update regridding
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME  !(lzh,02/01/2015)

      !=================================================================
      ! READ_STREETS_MASKS begins here!
      !=================================================================

      !------------------------------------
      ! China Mask (for 2001 CO emisisons)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) //
     &            'Streets_200607/China_mask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_STREETS_MASKS: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 0d0,       I1x1,     J1x1-1,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Save the 1x1 China mask for future use
      MASK_CHINA_1x1(:,:) = GEN_1x1(:,:)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
!      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
!      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK_CHINA )

      ! (lzh,02/01/2015)
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid from GENERIC 1x1 to current model resolution
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,       J1x1-1,
     &                    GEN_1x1,    MASK_CHINA, IS_MASS=0,
     &                    netCDF=.TRUE.                      )
      WHERE( MASK_CHINA > 0d0 ) MASK_CHINA = 1d0

      !------------------------------------
      ! SE Asia Mask (for 2000 emissions)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) //
     &            'Streets_200607/SE_Asia_mask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 0d0,       I1x1,     J1x1-1,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
!      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
!      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK_SE_ASIA )

      ! (lzh,02/01/2015)
      ! Regrid from GENERIC 1x1 GRID to current model resolution
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,         J1x1-1,
     &                    GEN_1x1,    MASK_SE_ASIA, IS_MASS=0,
     &                    netCDF=.TRUE.                        )
      WHERE( MASK_SE_ASIA > 0d0 ) MASK_SE_ASIA = 1d0

      ! Return to calling program
      END SUBROUTINE READ_STREETS_MASKS

!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS_MASKS_05x0666
!
!******************************************************************************
!  Subroutine READ_STREETS_MASKS reads and regrids the China and SE Asia
!  masks that define the David Streets' emission regions.  Specially modified
!  for the GEOS-5 0.5 x 0.666 nested grid simulations.
!  (yxw, dan, bmy, 11/6/08)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_05x0666

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4                  :: ARRAY(I05x0666,J05x0666,1)
      REAL*8                  :: GEOS_05x0666(I05x0666,J05x0666,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_STREETS_MASKS begins here!
      !=================================================================

      ! Zero arrays
      ARRAY        = 0d0
      GEOS_05x0666 = 0d0

      !------------------------------------
      ! China Mask (for 2001 CO emisisons)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR ) //
     &            'Streets_200607/China_mask.geos5.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_STREETS_MASKS: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 0d0,       I05x0666, J05x0666,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding


      GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Save the 1x1 China mask for future use
      MASK_CHINA_05x0666(:,:) = GEOS_05x0666(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_05x0666( 1, 'unitless', GEOS_05x0666,  MASK_CHINA )

      !------------------------------------
      ! SE Asia Mask (for 2000 emissions)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR ) //
     &            'Streets_200607/SE_Asia_mask.geos5.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 0d0,       I05x0666, J05x0666,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_05x0666( 1, 'unitless', GEOS_05x0666, MASK_SE_ASIA)

      ! Return to calling program
      END SUBROUTINE READ_STREETS_MASKS_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE INIT_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine INIT_STREETS_ANTHRO allocates and zeroes all module arrays.
!  (bmy, 8/16/06, 11/6/08)
!
!  NOTES:
!  (1 ) Now allocate MASK_CHINA_1x1 (bmy, 9/5/06)
!  (2 ) Now calls READ_STREETS_MASKS_05x0666 for the GEOS-5 0.5 x 0.666
!        nested-grid simulations (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LSTREETS

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_STREETS begins here!
      !=================================================================

      ! Return if LSTREETS is false
      IF ( .not. LSTREETS ) RETURN

      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0

      ALLOCATE( CO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO2' )
      CO2 = 0d0

      ALLOCATE( CH4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4' )
      CH4 = 0d0

      ! Now allocate VOCs (phs, 3/7/08)
      ALLOCATE( ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACET' )
      ACET = 0d0


      ALLOCATE( ALD2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALD2' )
      ALD2 = 0d0


      ALLOCATE( C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C2H6' )
      C2H6 = 0d0


      ALLOCATE( C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C3H8' )
      C3H8 = 0d0


      ALLOCATE( PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRPE' )
      PRPE = 0d0


      ALLOCATE( ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK4' )
      ALK4 = 0d0


      ALLOCATE( CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH2O' )
      CH2O = 0d0


      ALLOCATE( MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MEK' )
      MEK = 0d0
      ! -- end VOCs

      !---------------------------------------------------
      ! Pre-store array for grid box surface area in cm2
      !---------------------------------------------------

      ! Allocate array
      ALLOCATE( A_CM2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      !---------------------------------------------------
      ! Read & Regrid masks for Streets' emissions
      !---------------------------------------------------

      ALLOCATE( MASK_CHINA_1x1( I1x1, J1x1-1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CHINA_1x1' )
      MASK_CHINA_1x1 = 0

      ALLOCATE( MASK_CHINA_05x0666( I05x0666, J05x0666 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CHINA_05x0666' )
      MASK_CHINA_05x0666 = 0

      ALLOCATE( MASK_CHINA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CHINA' )
      MASK_CHINA = 0d0

      ALLOCATE( MASK_SE_ASIA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_SE_ASIA' )
      MASK_SE_ASIA = 0d0

      ! Read China & SE Asia masks from disk
#if   defined( GRID05x0666 )
      CALL READ_STREETS_MASKS_05x0666    ! GEOS-5 nested grids
#else
      CALL READ_STREETS_MASKS            ! Global simulations
#endif

      ! Return to calling program
      END SUBROUTINE INIT_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine CLEANUP_STREETS deallocates all module arrays
!  (bmy, 8/16/06, 9/5/06)
!
!  NOTES:
!  (1 ) Now deallocate MASK_CHINA_1x1 (bmy, 9/5/06)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_STREETS begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2          ) ) DEALLOCATE( A_CM2          )
      IF ( ALLOCATED( MASK_CHINA_1x1 ) ) DEALLOCATE( MASK_CHINA_1x1 )
      IF ( ALLOCATED( MASK_CHINA_05x0666 ) )
     &  DEALLOCATE( MASK_CHINA_05x0666 )                   !(dan)
      IF ( ALLOCATED( MASK_CHINA     ) ) DEALLOCATE( MASK_CHINA     )
      IF ( ALLOCATED( MASK_SE_ASIA   ) ) DEALLOCATE( MASK_SE_ASIA   )
      IF ( ALLOCATED( NOx            ) ) DEALLOCATE( NOx            )
      IF ( ALLOCATED( CO             ) ) DEALLOCATE( CO             )
      IF ( ALLOCATED( SO2            ) ) DEALLOCATE( SO2            )
      IF ( ALLOCATED( NH3            ) ) DEALLOCATE( NH3            )
      IF ( ALLOCATED( CH4            ) ) DEALLOCATE( CH4            )
      IF ( ALLOCATED( CO2            ) ) DEALLOCATE( CO2            )

      ! Now deallocate VOCs (phs, 3/7/08)
      IF ( ALLOCATED( C3H8           ) ) DEALLOCATE( C3H8          )
      IF ( ALLOCATED( C2H6           ) ) DEALLOCATE( C2H6          )
      IF ( ALLOCATED( ALK4           ) ) DEALLOCATE( ALK4          )
      IF ( ALLOCATED( ALD2           ) ) DEALLOCATE( ALD2          )
      IF ( ALLOCATED( PRPE           ) ) DEALLOCATE( PRPE          )
      IF ( ALLOCATED( MEK            ) ) DEALLOCATE( MEK           )
      IF ( ALLOCATED( CH2O           ) ) DEALLOCATE( CH2O          )
      IF ( ALLOCATED( ACET           ) ) DEALLOCATE( ACET          )

      ! Return to calling program
      END SUBROUTINE CLEANUP_STREETS_ANTHRO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE STREETS_ANTHRO_MOD
