! $Id: seasalt_mod.f,v 1.2 2010/03/09 15:03:46 daven Exp $
      MODULE SEASALT_MOD
!
!******************************************************************************
!  Module SEASALT_MOD contains arrays and routines for performing either a
!  coupled chemistry/aerosol run or an offline seasalt aerosol simulation.
!  Original code taken from Mian Chin's GOCART model and modified accordingly.
!  (bec, rjp, bmy, 6/22/00, 7/18/08)
!
!  Seasalt aerosol species: (1) Accumulation mode (usually 0.1 -  0.5 um)
!                           (2) Coarse mode       (usually 0.5 - 10.0 um)
!
!  NOTE: You can change the bin sizes for accumulation mode and coarse
!        mode seasalt in the "input.geos" file in v7-yy-zz and higher.
!
!  Module Variables:
!  ============================================================================
!  (1 ) DRYSALA  (INTEGER) : Drydep index for accumulation mode sea salt
!  (2 ) DRYSALC  (INTEGER) : Drydep index for coarse mode sea salt
!  (3 ) NSALT    (INTEGER) : Number of sea salt tracers
!  (4 ) IDDEP    (INTEGER) : Drydep index array for sea salt tracers
!  (5 ) REDGE    (REAL*8 ) : Array for edges of seasalt radius bins
!  (6 ) RMID     (REAL*8 ) : Array for centers of seasalt radius bins
!  (7 ) SRC      (REAL*8 ) : Array for baseline seasalt emission/bin [kg/m2]
!  (7 ) SRC_N    (REAL*8 ) : Array for baseline seasalt emission/bin [#/m2]
!  (8 ) SS_DEN   (REAL*8 ) : Sea salt density [kg/m3]
!  (9 ) ALK_EMIS (REAL*8 ) : Array for alkalinity [kg]
!  (10) N_DENS   (REAL*8 ) : Number density of seasalt emissions [#/m3]
!  (11) SALT_V   (REAL*8)  : Log-normal volum size distribution for sea salt
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMSEASALT        : Driver routine for sea salt loss processes
!  (2 ) WET SETTLING       : Routine which performs wet settling of sea salt
!  (3 ) DRY_DEPOSITION     : Routine which performs dry deposition of sea salt
!  (4 ) EMISSSEASALT       : Driver routine for sea salt emissions
!  (5 ) SRCSALT            : Updates surface mixing ratio for sea salt
!  (6 ) GET_ALK            : Gets the alkalinity of seasalt emissions
!  (6 ) INIT_SEASALT       : Allocates all module arrays
!  (7 ) CLEANUP_SEASALT    : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by "seasalt_mod.f":
!  ============================================================================
!  (1 ) dao_mod.f          : Module w/ arrays for GMAO met fields
!  (2 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (3 ) drydep_mod.f       : Module w/ GEOS-CHEM drydep routines
!  (4 ) error_mod.f        : Module w/ I/O error and NaN check routines
!  (5 ) grid_mod.f         : Module w/ horizontal grid information
!  (6 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (7 ) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (8 ) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (9 ) time_mod.f         : Module w/ routines to compute date & time
!  (10) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
!
!  NOTES:
!  (1 ) Now references "logical_mod.f" and "tracer_mod.f".  Comment out
!        SS_SIZE, this has been replaced by SALA_REDGE_um and SALC_REDGE_um
!        from "tracer_mod.f".  Increased NR_MAX to 200. (bmy, 7/20/04)
!  (2 ) Added error check in EMISSSEASALT (bmy, 1/20/05)
!  (3 ) Now references "pbl_mix_mod.f" (bmy, 2/22/05)
!  (4 ) Added routine GET_ALK to account for alkalinity. (bec, bmy, 4/13/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Now only call dry deposition routine if LDRYD=T (bec, bmy, 5/23/06)
!  (7 ) Remove unused variables from GET_ALK.  Also fixed variable declaration
!        bug in WET_SETTLING. (bec, bmy, 9/5/06)
!  (8 ) Extra error check for low RH in WET_SETTLING (phs, 6/11/08)
!  (9 ) Bug fix to remove a double-substitution in GET_ALK (bec, bmy, 7/18/08)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "seasalt_mod.f"
      !=================================================================

      ! Make everyting PRIVATE ...
      PRIVATE

      ! ... except these variables (jaegle 5/11/11)
      PUBLIC :: SALT_V
      PUBLIC :: DMID

      ! ... except these routines
      PUBLIC :: CHEMSEASALT
      PUBLIC :: EMISSSEASALT
      PUBLIC :: CLEANUP_SEASALT
      PUBLIC :: GET_ALK

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER, PARAMETER   :: NSALT = 2
      INTEGER, PARAMETER   :: NR_MAX = 200
      INTEGER              :: DRYSALA, DRYSALC

      ! Arrays
      INTEGER              :: IDDEP(NSALT)
      REAL*8,  ALLOCATABLE :: REDGE(:,:)
      REAL*8,  ALLOCATABLE :: RMID(:,:)
      REAL*8,  ALLOCATABLE :: SRC(:,:)
      REAL*8,  ALLOCATABLE :: SRC_N(:,:)
      REAL*8,  ALLOCATABLE :: ALK_EMIS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: N_DENS(:,:,:,:)
      ! Add SALT_V and DMID (jaegle 5/11/11)
      REAL*8,  ALLOCATABLE :: SALT_V(:)
      REAL*8,  ALLOCATABLE :: DMID(:)
      REAL*8               :: SS_DEN(NSALT)    = (/ 2200.d0, 2200.d0 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMSEASALT
!
!******************************************************************************
!  Subroutine CHEMSEASALT is the interface between the GEOS-CHEM main program
!  and the seasalt chemistry routines that mostly calculates seasalt dry
!  deposition (rjp, bmy, 1/24/02, 5/23/06)
!
!  NOTES:
!  (1 ) Now reference STT from "tracer_mod.f".  Now references LPRT from
!        "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Now only call DRY_DEPOSITION if LDRYD=T (bec, bmy, 5/23/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT,    LDRYD
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSALA, IDTSALC

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: N

      !=================================================================
      ! CHEMSEASALT begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Initialize (if necessary)
         CALL INIT_SEASALT

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'SALA' )
                  DRYSALA = N
               CASE ( 'SALC' )
                  DRYSALC = N
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO

         ! Store in IDDEP array
         IDDEP(1) = DRYSALA
         IDDEP(2) = DRYSALC

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Maybe someday we should merge these two separate calculations
      ! into one (rjp, 4/3/04)
      !=================================================================

      !-------------------
      ! Accumulation mode
      !-------------------
      CALL WET_SETTLING( STT(:,:,:,IDTSALA), 1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Accum' )

      IF ( LDRYD ) THEN
      ! If LNLPBL (non local PBL mixing) is turned on, do sea salt
      ! dry deposition in vdiff as for all the other aerosols (jaegle 5/11/11)
      !IF ( LDRYD .AND. .NOT. LNLPBL ) THEN  ! This is not yet supported for adjoint code
         CALL DRY_DEPOSITION( STT(:,:,:,IDTSALA), 1 )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Accum' )
      ENDIF

      !-------------------
      ! Coarse mode
      !-------------------
      CALL WET_SETTLING( STT(:,:,:,IDTSALC), 2 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Coarse' )

      IF ( LDRYD ) THEN
      ! If LNLPBL (non local PBL mixing) is turned on, do sea salt
      ! dry deposition in vdiff as for all the other aerosols (jaegle 5/11/11)
      !IF ( LDRYD .AND. .NOT. LNLPBL ) THEN  ! This is not yet supported for adjoint code
         CALL DRY_DEPOSITION( STT(:,:,:,IDTSALC), 2 )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Coarse')
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMSEASALT

!------------------------------------------------------------------------------

      SUBROUTINE WET_SETTLING( TC, N )
!
!******************************************************************************
!  Subroutine WET_SETTLING performs wet settling of sea salt.
!  (bec, rjp, bmy, 4/20/04, 6/11/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer [kg]
!  (2 ) N  (INTEGER) : N=1 is accum mode; N=2 is coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified tracer
!
!  NOTES:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Bug fix: DTCHEM has to be REAL*8, not integer. (bmy, 9/7/06)
!  (4 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 5/1/08)
!  (5 ) Update sea salt density calculation using Tang et al. (1997) (bec, jaegle 5/11/11)
!  (6 ) Update hygroscopic growth for sea salt using Lewis and Schwartz (2006) and and density
!       calculation based on Tang et al. (1997) (bec, jaegle 5/11/11)
!  (7 ) Itegrate settling velocity over entire size distribution (jaegle 5/11/11)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : T, BXHEIGHT, RH
      USE DIAG_MOD,      ONLY : AD44
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : SALA_REDGE_um, SALC_REDGE_um, XNUMOL
      USE TRACERID_MOD,  ONLY : IDTSALA,       IDTSALC
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      ! add (jaegle 5/11/11)
      USE ERROR_MOD,     ONLY : DEBUG_MSG, ERROR_STOP

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_GCTM"      ! g0
#     include "CMN_DIAG"      ! ND44

      ! Argumetns
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,      J,     L
      REAL*8                 :: DELZ,   DELZ1, REFF,     DEN
      REAL*8                 :: P,      DP,    PDP,      TEMP
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      ! replace RCM with RUM (radis in micron) jaegle 5/11/11
      REAL*8                 :: RUM,    RWET,  RATIO_R,  RHO
      REAL*8                 :: TOT1,   TOT2,  DTCHEM
      REAL*8                 :: VTS(LLPAR)
      REAL*8                 :: TC0(LLPAR)
      ! New variables (jaegle 5/11/11)
      REAL*8                 :: SW
      REAL*8                 :: R0,       R1, NR, DEDGE, SALT_MASS
      REAL*8                 :: SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW
      REAL*8                 :: WTP, RHO1
      INTEGER                :: ID
      LOGICAL, SAVE          :: FIRST = .TRUE.

      ! Parameters
      REAL*8,  PARAMETER     :: C1 =  0.7674d0
      REAL*8,  PARAMETER     :: C2 =  3.079d0
      REAL*8,  PARAMETER     :: C3 =  2.573d-11
      REAL*8,  PARAMETER     :: C4 = -1.424d0
      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) (jaegle 5/11/11)
      REAL*8,  PARAMETER     :: A1 =  7.93d-3
      REAL*8,  PARAMETER     :: A2 = -4.28d-5
      REAL*8,  PARAMETER     :: A3 =  2.52d-6
      REAL*8,  PARAMETER     :: A4 = -2.35d-8
      ! increment of radius for integration of settling velocity (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2
      ! parameter for convergence
      REAL*8,  PARAMETER     :: EPSI = 1.0D-4
      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols (jaegle 5/11/11)
      ! geometric dry mean diameters (microns)
      REAL*8,  PARAMETER     ::   RG_A = 0.085d0
      REAL*8,  PARAMETER     ::   RG_C = 0.4d0
      ! sigma of the size distribution
      REAL*8,  PARAMETER     ::   SIG_A = 1.5d0
      REAL*8,  PARAMETER     ::   SIG_C = 1.8d0


      !=================================================================
      ! WET_SETTLING begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Sea salt density [kg/m3]
      DEN  = SS_DEN( N )

      ! Seasalt effective radius (i.e. midpt of radius bin) [m]
      SELECT CASE ( N )

         ! Accum mode
         ! add R0 and R1 = edges if the sea salt size bins (jaegle 5/11/11)
         CASE( 1 )
            REFF = 0.5d-6 * ( SALA_REDGE_um(1) + SALA_REDGE_um(2) )
            R0 = SALA_REDGE_um(1)
            R1 = SALA_REDGE_um(2)

         ! Coarse mode
         CASE( 2 )
            REFF = 0.5d-6 * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            R0 = SALC_REDGE_um(1)
            R1 = SALC_REDGE_um(2)

      END SELECT

      ! Number of dry radius size bins between lowest radius (accumulation
      ! mode) and largest radii (coarse mode) (jaegle 5/11/11)

      NR = INT( ( ( SALC_REDGE_um(2) - SALA_REDGE_um(1) ) / DR )
     &                  + 0.5d0 )

      ! Error check
      IF ( NR > NR_MAX ) THEN
        CALL ERROR_STOP( 'Too many bins!', 'SRCSALT (seasalt_mod.f)')
      ENDIF

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a coarse mode
      ! and accumulation model log-normal distribution functions (jaegle 5/11/11)
      !=================================================================
      IF ( FIRST) THEN

        ! Lower edge of 0th bin
	DEDGE=SALA_REDGE_um(1) * 2d0

	! Loop over diameters
        DO ID = 1, NR
           ! Diameter of mid-point in microns
           DMID(ID)  = DEDGE + ( DR )

	   ! Calculate the dry volume size distribution as the sum of two log-normal
	   ! size distributions. The parameters for the size distribution are
	   ! based on Reid et al. and Quinn et al.
	   ! The scaling factors 13. and 0.8 for acc and coarse mode aerosols are
	   ! chosen to obtain a realistic distribution
	   ! SALT_V (D) = dV/dln(D) [um3]
	   SALT_V(ID) = PI / 6d0* (DMID(ID)**3) * (
     &         13d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_A*2d0) )**2d0/
     &                   LOG(SIG_A)**2d0 )
     &         /( sqrt(2d0 * PI) * LOG(SIG_A) )  +
     &         0.8d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_C*2d0) )**2d0/
     &                   LOG(SIG_C)**2d0)
     &         /( sqrt(2d0 * PI) * LOG(SIG_C) )  )
	   ! update the next edge
	   DEDGE = DEDGE + DR*2d0
        ENDDO

        ! Reset after the first time
        IF ( FIRST ) FIRST = .FALSE.
      ENDIF


      ! Sea salt radius [cm]
      !RCM  = REFF * 100d0
      ! The radius used in the Gerber formulation for hygroscopic growth
      ! of sea salt should be in microns (RUM) instead of cm (RCM). Replace RCM
      ! with RUM (jaegle 5/11/11)
      !RUM  = REFF * 1d6

      ! Exponential factors
      !FAC1 = C1 * ( RCM**C2 )
      !FAC2 = C3 * ( RCM**C4 )
      ! Replace with RUM (jaegle 5/11/11)
      !FAC1 = C1 * ( RUM**C2 )
      !FAC2 = C3 * ( RUM**C4 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET )
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             )
!$OMP+PRIVATE( SW,      ID,    SALT_MASS_TOTAL,      VTS_WEIGHT       ) !jaegle 5/11/11
!$OMP+PRIVATE( DMIDW,   RHO1,  WTP,  SALT_MASS                        ) !jaegle 5/11/11
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize
         DO L = 1, LLPAR
            VTS(L) = 0d0
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1d0

            ! Temperature [K]
            TEMP    = T(I,J,L)

            ! Cap RH at 0.99
            RHB     = MIN( 0.99d0, RH(I,J,L) * 1d-2 )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m]
            ! (Gerber, 1985)
            !RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHB))+RCM**3.d0)**0.33d0
	    ! Fix bugs in the Gerber formula:  a log10 (instead of ln) should be used and the
	    ! dry radius should be expressed in micrometers (instead of cm) also add more significant
	    ! digits to the exponent (should be 1/3) (jaegle 5/11/11)
            !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.d0)**0.33333d0

            ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth (bec, jaegle 5/11/11)
            RWET = REFF * (4.d0 / 3.7d0) *
     &              ( (2.d0 - RHB)/(1.d0 - RHB) )**(1.d0/3.d0)


            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.d0

            ! Density of the wet aerosol (kg/m3)
            RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0

            ! Above density calculation is chemically unsound because it ignores chemical solvation.
            ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of wet aerosol (kg/m3)
            ! (bec, jaegle 5/11/11)
            RATIO_R = ( REFF / RWET )
            ! Assume an initial density of 1000 kg/m3
            RHO  = 1000.D0
            RHO1 = 0.d0 !initialize (bec, 6/21/10)
            DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
                ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
                WTP    = 100.d0 * DEN/RHO * RATIO_R**3.d0
                ! Then calculate density of wet aerosol using equation 5
                ! in Tang et al., 1997 [kg/m3]
                RHO1   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) +
     $               (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
                ! Now calculate new weight percent using above density calculation
                WTP    = 100.d0 * DEN/RHO1 * RATIO_R**3.d0
                ! Now recalculate new wet density [kg/m3]
                RHO   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) +
     $              (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
            ENDDO

            ! Dp = particle diameter [um]
            DP      = 2.d0 * RWET * 1.d6

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.d0 * RHO * RWET**2 * g0 / 9.d0

            !===========================================================
            ! NOTE: Slip correction factor calculations following
            ! Seinfeld, pp464 which is thought to be more accurate
            ! but more computation required. (rjp, 1/24/02)
            !
            ! # air molecule number density
            ! num = P * 1d3 * 6.023d23 / (8.314 * Temp)
            !
            ! # gas mean free path
            ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
            !
            ! # Slip correction
            ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp
            !     &     / (2. * lamda))) / Dp
            !
            ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
            ! which produces slip correction factore with small error
            ! compared to the above with less computation.
            !===========================================================

            ! Slip correction factor (as function of P*dp)
            Slip = 1.d0+(15.60d0 + 7.0d0 * EXP(-0.059d0 * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature
            VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )

            ! Settling velocity [m/s]
            VTS(L) = CONST * Slip / VISC

            ! This settling velocity is for the mid-point of the size bin.
            ! In the following we derive scaling factors to take into account
	    ! the strong dependence on radius of the settling velocity and the
	    ! mass size distribution:
	    !  VTS_WEIGHTED = total( M(k) x VTS(k)) / total( M(k) )
	    ! The settling velocity is a function of the radius squared (see definition
	    ! of CONST above)
	    ! so VTS(k) = VTS * (RMID(k)/RWET)^2
            ! (jaegle 5/11/11)

	    SALT_MASS_TOTAL = 0d0
	    VTS_WEIGHT      = 0d0
	    DO ID = 1, NR
	       ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diamter):
	       ! dM/dlnDw = dV/dlnDw * RHO, we assume that the density of sea-salt
	       ! doesn't change much over the size range.
	       ! and
	       ! dV/dlnDw = dV/dlnD * dlnD/dlnDw = dV/dlnD * Dw/D = dV/dlnD * Rwet/Rdry
	       ! Further convert to dM/dDw = dM/dln(Dw) * dln(Dw)/Dw = dM/dln(Dw)/Dw
	       ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * RHO /Rw
	       !
	       IF (DMID(ID) .ge. R0*2d0 .and. DMID(ID) .le. R1*2d0 ) THEN
	         DMIDW = DMID(ID) * RWET/REFF  ! wet radius [um]
	         SALT_MASS   = SALT_V(ID) * RWET/REFF * RHO / (DMIDW*0.5d0)
	         VTS_WEIGHT  = VTS_WEIGHT +
     &              SALT_MASS * VTS(L) * (DMIDW/(RWET*1d6*2d0) )**2d0 *
     &                            (2d0 * DR *  RWET/REFF)
	         SALT_MASS_TOTAL=SALT_MASS_TOTAL+SALT_MASS *
     &                            (2d0 * DR *  RWET/REFF)
	       ENDIF

            ENDDO
            ! Calculate the weighted settling velocity:
            VTS(L) = VTS_WEIGHT/SALT_MASS_TOTAL
         ENDDO

         ! Method is to solve bidiagonal matrix which is
         ! implicit and first order accurate in z (rjp, 1/24/02)

         ! Save initial tracer concentration in column
         DO L = 1, LLPAR
            TC0(L) = TC(I,J,L)
         ENDDO

         ! We know the boundary condition at the model top
         L    = LLTROP
         DELZ = BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

         DO L = LLTROP-1, 1, -1
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)
            TC(I,J,L) = 1.d0 / ( 1.d0 + DTCHEM * VTS(L) / DELZ )
     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
     &                *  TC(I,J,L+1) )
         ENDDO

         !==============================================================
         ! ND44 diagnostic: sea salt loss [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 ) THEN

            ! Initialize
            TOT1 = 0d0
            TOT2 = 0d0

            ! Compute column totals of TCO(:) and TC(I,J,:,N)
            DO L = 1, LLPAR
               TOT1 = TOT1 + TC0(L)
               TOT2 = TOT2 + TC(I,J,L)
            ENDDO

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * XNUMOL(IDTSALA) / AREA_CM2

            ! Store in AD44 array
            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE WET_SETTLING

!------------------------------------------------------------------------------

      SUBROUTINE DRY_DEPOSITION( TC, N )
!
!******************************************************************************
!  Subroutine DRY_DEPOSITION computes the loss of sea salt by dry deposition
!  at the surface, using an implicit method. (bec, rjp, bmy, 4/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer [kg]
!  (2 ) N  (INTEGER) : N=1 is accum mode; N=2 is coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified tracer
!
!  NOTES:
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (2 ) Update to calculate the drydep throughout the entire PBL instead of just
!       at the surface. This is more in line with what is done in dry_dep.f. This
!       is only used if LNLPBL is turned off (or for GEOS-4 and prior met fields).
!       (jaegle 5/11/11)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTSALA,   IDTSALC
      USE TIME_MOD,     ONLY : GET_MONTH, GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      ! Add PBL variables (jaegle  5/5/11)
      USE PBL_MIX_MOD, ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L


#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! g0
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,        J,     L,      DTCHEM
      REAL*8                 :: OLD,      NEW,   G,      REFF
      REAL*8                 :: DIAM,     U_TS0, REYNOL, ALPHA
      REAL*8                 :: BETA,     GAMMA, DENS,   FLUX
      REAL*8                 :: AREA_CM2, TOT1,  TOT2
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      ! New variables for applying dry dep thoughout the PBL (jaegle 5/11/11)
      INTEGER                :: PBL_MAX
      REAL*8                 :: F_UNDER_TOP, FREQ

      ! Parameters
      REAL*8,  PARAMETER     :: RHOA = 1.25d-3

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Maximum extent of the PBL [model layers] (jaegle 5/11/11)
      PBL_MAX = GET_PBL_MAX_L()

      ! Zero temporary array for drydep diagnostic
      IF ( ND44 > 0 ) ND44_TMP = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, OLD, NEW, FLUX )
!$OMP+PRIVATE( L, F_UNDER_TOP   , FREQ ) ! (jaegle 5/11/11)
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over levels under PBL
      DO L = 1, PBL_MAX

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Fraction of box (I,J,L) under PBL top [unitless]
            F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

            ! Only apply drydep to boxes w/in the PBL
            IF ( F_UNDER_TOP > 0d0 ) THEN

              ! Sea salt dry deposition frequency [1/s] accounting
              ! for fraction of each grid box located beneath the PBL top
              FREQ = DEPSAV(I,J,IDDEP(N)) * F_UNDER_TOP
              ! Only apply drydep loss if FREQ is nonzero
              IF ( FREQ > 0d0 ) THEN
                  ! Old tracer concentration [kg]
                  OLD  = TC(I,J,L)

                  ! New tracer concentration [kg]
                  NEW  = OLD * EXP( -FREQ * DTCHEM )


            ! Old tracer concentration [kg]
            !OLD  = TC(I,J,1)

            ! New tracer concentration [kg]
            !NEW  = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )

            !===========================================================
            ! ND44 diagnostic: sea salt drydep loss [molec/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Convert drydep loss from [kg/s] to [molec/cm2/s]
               FLUX = ( OLD - NEW ) / DTCHEM
               FLUX = FLUX * XNUMOL(IDTSALA) / AREA_CM2

               ! Store in AD44
               !AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

            ! Update tracer array
            TC(I,J,L) = NEW
                  ENDIF
            ENDIF
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !====================================================================
      ! ND44 diagnostic: save into AD44 array, summing in the vertical
      !====================================================================
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD44(I,J,IDDEP(N),1) = SUM( ND44_TMP(I,J,:) )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE DRY_DEPOSITION

!------------------------------------------------------------------------------

      SUBROUTINE EMISSSEASALT
!
!******************************************************************************
!  Subroutine EMISSSEASALT is the interface between the GEOS-CHEM model
!  and the SEASALT emissions routines in "seasalt_mod.f".
!  (bec, rjp, bmy, 3/24/03, 2/22/05)
!
!  NOTES:
!  (1 ) Now references LPRT from "logical_mod.f" and STT from "tracer_mod.f".
!        (bmy, 7/20/04)
!  (2 ) Now make sure IDTSALA, IDTSALC are nonzero before calling SRCSALT.
!        (bmy, 1/26/05)
!  (3 ) Remove reference to header file "CMN" (bmy, 2/22/05)
!  (4 ) Now call INIT_SEASALT on the first timestep.  Also initialize ALK_EMIS
!        and N_DENS on each timestep. (bec, bmy, 4/13/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSALA, IDTSALC

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I, J, L, N

      !=================================================================
      ! EMISSSEASALT begins here!
      !=================================================================

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### in EMISSEASALT' )

      ! Allocate all module arrays (bec, bmy, 4/13/05)
      IF ( FIRST ) THEN
         CALL INIT_SEASALT
         FIRST = .FALSE.
      ENDIF

      ! Initialize for each timestep (bec, bmy, 4/13/05)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NSALT
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ALK_EMIS(I,J,L,N) = 0d0
         N_DENS(I,J,L,N)   = 0d0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Accumulation mode emissions
      IF ( IDTSALA > 0 ) THEN
         CALL SRCSALT( STT(:,:,:,IDTSALA), 1 )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Accum' )
      ENDIF

      ! Coarse mode emissions
      IF ( IDTSALC > 0 ) THEN
         CALL SRCSALT( STT(:,:,:,IDTSALC), 2 )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Coarse' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE EMISSSEASALT

!-----------------------------------------------------------------------------

      SUBROUTINE SRCSALT( TC, N )
!
!******************************************************************************
!  The new SRCSALT is based on the sea salt source function of Gong (2003) with
!  the empirical sea surface temperature (SST) dependence of Jaegle et al. (2011). This
!  SST dependence was derived based on comparisons to cruise observations of
!  coarse mode sea salt mass concentrations.
!
!  Contact: Lyatt Jaegle (jaegle@uw.edu)
!
! Old:
!!  Subroutine SRCSALT updates the surface mixing ratio of dry sea salt
!!  aerosols for NSALT size bins.  The generation of sea salt aerosols
!!  has been parameterized following Monahan et al. [1986] parameterization
!!  as described by Gong et al. [1997].  (bec, rjp, bmy, 4/20/04, 11/23/09)
!
!!  Contact: Becky Alexander (bec@io.harvard.edu) or
!!           Rokjin Park     (rjp@io.harvard.edu)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer array [v/v]
!  (2 ) N  (INTEGER) : N=1 denotes accumulation mode; N=2 denotes coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified sea salt concentration [v/v]
!
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
!  (3 ) Gong, S. L., "A parameterization of sea-salt aerosol source function
!        for sub- and super-micron particles", Global Biogeochem.  Cy., 17(4),
!        1097, doi:10.1029/2003GB002079, 2003.
!  (4 ) Jaegle, L., P.K. Quinn, T.S. Bates, B. Alexander, J.-T. Lin, "Global
!        distribution of sea salt aerosols: New constraints from in situ and
!        remote sensing observations", Atmos. Chem. Phys., 11, 3137-3157,
!        doi:10.5194/acp-11-3137-2011.
!
!  NOTES:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_TOP_L from "pbl_mix_mod.f".
!        Removed reference to header file CMN.  Removed reference to
!        "pressure_mod.f".  (bmy, 2/22/05)
!  (3 ) Now also compute alkalinity and number density of seasalt emissions.
!        (bec, bmy, 4/13/05)
!  (4 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) The source function is for wet aerosol radius (RH=80%, with a radius
!        twice the size of dry aerosols) so BETHA should be set to 2
!        instead of 1.  Also now use LOG10 instead of LOG in the expressions
!        for the seasalt base source, since we need the logarithm to the base
!        10. (jaegle, bec, bmy, 11/23/09)
!  (6 ) Update to use the Gong (2003) source function (jaegle 5/11/11)
!  (7 ) Apply an empirical sea surface temperature dependence to Gong (2003) (jaegle 5/11/11)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : PBL, AD, IS_WATER, AIRVOL
      ! Add TSKIN (jaegle 5/11/11)
      USE DAO_MOD,       ONLY : TSKIN ! jaegle
      USE DIAG_MOD,      ONLY : AD08
      USE ERROR_MOD,     ONLY : DEBUG_MSG, ERROR_STOP
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE TRACER_MOD,    ONLY : SALA_REDGE_um, SALC_REDGE_um, XNUMOL

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND44, ND08
#     include "CMN_GCTM"      ! PI

      ! Arguments
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,     J,      L
      INTEGER                :: R,     NR,     NTOP
      REAL*8                 :: W10M,  DTEMIS, R0
      REAL*8                 :: R1,    CONST,  CONST_N
      REAL*8                 :: FEMIS, A_M2
      REAL*8                 :: SALT(IIPAR,JJPAR)
      REAL*8                 :: SALT_N(IIPAR,JJPAR)
      ! New variables (jaegle 5/11/11)
      REAL*8                 :: A, B, SST, SCALE

      ! Increment of radius for Emission integration (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2
      REAL*8, PARAMETER      :: BETHA = 2.d0

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCSALT begins here!
      !=================================================================

      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! no longer used (jaegle 5/11/11)
      ! Constant [volume * time * other stuff??]
      !CONST   = 4d0/3d0 * PI * DR * DTEMIS * 1.d-18 * 1.373d0

      !CONST_N = DTEMIS * DR * 1.373d0
      !  Constant for converting from [#/m2/s/um] to [#/m2]
      CONST_N = DTEMIS * (DR * BETHA)

      ! Lower and upper limit of size bin N [um]
      ! Note that these are dry size bins. In order to
      ! get wet (RH=80%) sizes, we need to multiply by
      ! BETHA.
      SELECT CASE( N )

         ! Accum mode
         CASE( 1 )
            R0 = SALA_REDGE_um(1)
            R1 = SALA_REDGE_um(2)

         ! Coarse mode
         CASE( 2 )
            R0 = SALC_REDGE_um(1)
            R1 = SALC_REDGE_um(2)

      END SELECT

      ! Number of radius size bins
      NR = INT( ( ( R1 - R0 ) / DR ) + 0.5d0 )

      ! Error check
      IF ( NR > NR_MAX ) THEN
         CALL ERROR_STOP( 'Too many bins!', 'SRCSALT (seasalt_mod.f)' )
      ENDIF

      ! Initialize source
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SALT(I,J)   = 0d0
         SALT_N(I,J) = 0d0
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Define edges and midpoints of each incrmental radius bin
      ! This only has to be done once per sea salt type
      !=================================================================
      IF ( FIRST ) THEN

         ! Lower edge of 0th bin
         REDGE(0,N) = R0

         ! Loop over the # of radius bins
         DO R = 1, NR

            ! Midpoint of IRth bin
            RMID(R,N)  = REDGE(R-1,N) + ( DR / 2d0 )

            ! Upper edge of IRth bin
            REDGE(R,N) = REDGE(R-1,N) + DR


            ! Sea salt base source [#/m2]. Note that the Gong formulation
            ! is for r80 (radius at 80% RH), so we need to multiply RMID
            ! by the scaling factor BETHA=2.
            A =  4.7*(1.+30.*(BETHA*RMID(R,N)))
     &                       **(-0.017*(BETHA*RMID(R,N))**(-1.44))
            B =  (0.433d0-LOG10(BETHA*RMID(R,N))) / 0.433d0
            SRC_N(R,N) = CONST_N * 1.373 * (1.d0/(BETHA*RMID(R,N))**(A))
     &           * (1.d0+0.057d0*(BETHA*RMID(R,N))**3.45d0)
     &           * 10d0**(1.607d0*EXP(-(B**2)))

            ! Sea salt base source [kg/m2]: multiply the number of particles
            ! by the dry volume multiplied by the dry density of sea-salt.
            SRC(R,N)  =   SRC_N(R,N) * 4d0/3d0 * PI * 1.d-18
     &           *  SS_DEN( N )  *  (RMID(R,N))**3


            !-----------------------------------------------------------
            ! IMPORTANT NOTE!
            !
            ! In mathematics, "LOG" means "log10".
            ! In Fortran,     "LOG" means "ln" (and LOG10 is "log10").
            !
            ! The following equations require log to the base 10, so
            ! we need to use the Fortran function LOG10 instead of LOG.
            ! (jaegle, bmy, 11/23/09)
            !-----------------------------------------------------------

!            ! Old Monahan et al. (1986) formulation
!            ! Sea salt base source [kg/m2]
!            CONST_N = DTEMIS * (DR * BETHA)
!            SRC(R,N)  = CONST * SS_DEN( N )
!     &           * ( 1.d0 + 0.057d0*( BETHA * RMID(R,N) )**1.05d0 )
!     &           * 10d0**( 1.19d0*
!     &             EXP(-((0.38d0-LOG10(BETHA*RMID(R,N)))/0.65d0)**2))
!     &           / BETHA**2

!            ! Sea salt base source [#/m2] (bec, bmy, 4/13/05)
!            SRC_N(R,N) = CONST_N * (1.d0/RMID(R,N)**3)
!     &           * (1.d0+0.057d0*(BETHA*RMID(R,N))**1.05d0)
!     &           * 10d0**(1.19d0*EXP(-((0.38d0-LOG10(BETHA*RMID(R,N)))
!     &           /0.65d0)**2))/ BETHA**2

!### Debug
!###           WRITE( 6, 100 ) R,REDGE(R-1,N),RMID(R,N),REDGE(R,N),SRC(R,N)
!### 100        FORMAT( 'IR, R0, RMID, R1: ', i3, 3f11.4,2x,es13.6 )
         ENDDO

         ! Reset only after N=NSALT
         IF ( FIRST .and.  N == NSALT ) FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Emission is integrated over a given size range for each bin
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, R, A_M2, W10M )
!$OMP+PRIVATE( SST , SCALE         ) !(jaegle 5/11/11)
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         A_M2 = GET_AREA_M2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Test if this is a water box
            IF ( IS_WATER(I,J) ) THEN

               ! Wind speed at 10 m altitude [m/s]
               W10M = SQRT( SFCWINDSQR(I,J) )

               ! Sea surface temperature in Celcius (jaegle 5/11/11)
               SST = TSKIN(I,J) - 273.15d0
               ! Limit SST to 0-30C range
               SST = MAX( SST , 0d0 ) ! limit to  0C
               SST = MIN( SST , 30d0 ) ! limit to 30C
               ! Empirical SST scaling factor (jaegle 5/11/11)
               SCALE = 0.329d0 + 0.0904d0*SST -
     &                 0.00717d0*SST**2d0 + 0.000207d0*SST**3d0

               ! Reset to using original Gong (2003) emissions (jaegle 6/30/11)
	       !SCALE = 1.0d0

! The source function calculated with GEOS-4 2x2.5 wind speeds is too high compared to GEOS-5
! at the same resolution. The 10m winds in GEOS-4 are too rapid. To correct this, apply a global
! scaling factor of 0.72 (jaegle 5/11/11)
#if   defined( GEOS_4 )
                SCALE = SCALE * 0.72d0
#endif

               ! Loop over size bins
               DO R = 1, NR

                  ! Update seasalt source into SALT [kg]
                  SALT(I,J)   = SALT(I,J) +
     &               ( SCALE * SRC(R,N)   * A_M2 * W10M**3.41d0 )

                  ! Update seasalt source into SALT_N [#]
                  ! (bec, bmy, 4/13/05)
                  SALT_N(I,J) = SALT_N(I,J) +
     &               ( SCALE * SRC_N(R,N) * A_M2 * W10M**3.41d0 )

               ENDDO
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Now partition seasalt emissions through boundary layer
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, FEMIS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Layer in which the PBL top occurs
         NTOP = CEILING( GET_PBL_TOP_L( I, J ) )

         ! Loop thru the boundary layer
         DO L = 1, NTOP

            ! Fraction of the PBL spanned by box (I,J,L) [unitless]
            FEMIS             = GET_FRAC_OF_PBL( I, J, L )

            ! Add seasalt emissions into box (I,J,L) [kg]
            TC(I,J,L)         = TC(I,J,L) + ( FEMIS * SALT(I,J) )

	    ! Alkalinity [kg] (bec, bmy, 4/13/05)
            ALK_EMIS(I,J,L,N) = SALT(I,J)

	    ! Number density [#/m3] (bec, bmy, 4/13/05)
	    N_DENS(I,J,L,N)   = SALT_N(I,J) / AIRVOL(I,J,L)

         ENDDO

         ! ND08 diagnostic: sea salt emissions [kg]
         IF ( ND08 > 0 ) THEN
            AD08(I,J,N) = AD08(I,J,N) + SALT(I,J)
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCSALT

!------------------------------------------------------------------------------

      SUBROUTINE GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N )
!
!******************************************************************************
!  Function GET_ALK returns the seasalt alkalinity emitted at each timestep to
!  sulfate_mod.f for chemistry on seasalt aerosols.
!  (bec, 12/7/04, 9/5/06)
!
!  Arguments as Input:
!  ============================================================================
!
!  NOTES:
!  (1 ) Becky Alexander says we can remove AREA1, AREA2 (bec, bmy, 9/5/06)
!  (2 ) Bug fix to remove a double-substitution.  Replace code lines for
!        TERM{123}A, TERM{123}B, TERM{123}AN, TERM{123}BN. (bec, bmy, 7/18/08)
!  (3 ) Updated hygroscopic growth parameters (bec, bmy, 11/23/09)
!******************************************************************************
!
      USE DAO_MOD,      ONLY : AD, RH
      USE ERROR_MOD,    ONLY : IT_IS_NAN
      USE TRACER_MOD,   ONLY : SALA_REDGE_um, SALC_REDGE_um

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L

      ! Return value
      REAL*8, INTENT(OUT)  :: ALK1, ALK2 ! [kg]
      REAL*8, INTENT(OUT)  :: Kt1, Kt2, Kt1N, Kt2N ! [s-1]

      REAL*8,  PARAMETER :: PI = 3.14159265
      REAL*8             :: N1, N2, Kt
      REAL*8             :: HGF, ALK
      REAL*8             :: RAD1, RAD2, RAD3
      REAL*8             :: term1a, term2a, term3a
      REAL*8             :: term1b, term2b, term3b
      REAL*8             :: term1aN, term2aN, term3aN
      REAL*8             :: term1bN, term2bN, term3bN
      REAL*8             :: const1, const2, const1N, const2N
      REAL*8             :: a1, a2, b1, b2, a1N, a2N, b1N, b2N
      REAL*8,  PARAMETER :: MINDAT = 1.d-20
      INTEGER            :: IRH
      REAL*8,  PARAMETER   :: gamma_SO2 = 0.11d0 !from Worsnop et al. (1989)
      REAL*8,  PARAMETER   :: gamma_HNO3 = 0.2d0 !from JPL [2001]
      REAL*8,  PARAMETER   :: Dg = 0.2d0 !gas phase diffusion coeff. [cm2/s]
      REAL*8,  PARAMETER   :: v = 3.0d4 !cm/s

      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! GET_ALK begins here!
      !=================================================================

      ! Zero variables
      ALK1  = 0.D0
      ALK2  = 0.D0
      KT1   = 0.D0
      KT2   = 0.D0
      KT1N  = 0.D0
      KT2N  = 0.D0
      N1    = 0.D0
      N2    = 0.D0

      ! [kg] use this when not transporting alk
      ALK1  = ALK_EMIS(I,J,L,1)
      ALK2  = ALK_EMIS(I,J,L,2)

      !-----------------------------------------------------------------------
      ! NOTE: If you want to transport alkalinity then uncomment this section
      ! (bec, bmy, 4/13/05)
      !
      !! alkalinity [v/v] to [kg] use this when transporting alk
      !! or using Liao et al [2004] assumption of a continuous supply of
      ! alkalinity based on Laskin et al. [2003]
      !ALK1 = STT(I,J,L,IDTSALA) * AD(I,J,L)/TCVV(IDTSALA)
      !ALK2 = STT(I,J,L,IDTSALC) * AD(I,J,L)/TCVV(IDTSALC)
      !-----------------------------------------------------------------------

      ! Conversion from [m-3] --> [cm-3]
      N1 = N_DENS(I,J,L,1) * 1.d-6
      N2 = N_DENS(I,J,L,2) * 1.d-6

      ALK = ALK1 + ALK2

      ! If there is any alkalinity ...
      IF ( ALK > MINDAT ) THEN

         ! set humidity index IRH as a percent
         IRH = RH(I,J,L)
         IRH = MAX(  1, IRH )
         IRH = MIN( 99, IRH )

         !--------------------------------------------------------------------
         ! Prior to 11/23/09:
         ! These hygroscopic growth factors are incorrect (bec, bmy, 11/23/09)
         !! hygroscopic growth factor for sea-salt from Chin et al. (2002)
         !IF ( IRH < 100 ) HGF = 2.2d0
         !IF ( IRH < 99  ) HGF = 1.9d0
         !IF ( IRH < 95  ) HGF = 1.8d0
         !IF ( IRH < 90  ) HGF = 1.6d0
         !IF ( IRH < 80  ) HGF = 1.5d0
         !IF ( IRH < 70  ) HGF = 1.4d0
         !IF ( IRH < 50  ) HGF = 1.0d0
         !--------------------------------------------------------------------

         ! Hygroscopic growth factor for sea-salt from Chin et al. (2002)
         ! Updated (bec, bmy, 11/23/09)
         IF ( IRH < 100 ) HGF = 4.8d0
         IF ( IRH < 99  ) HGF = 2.9d0
         IF ( IRH < 95  ) HGF = 2.4d0
         IF ( IRH < 90  ) HGF = 2.0d0
         IF ( IRH < 80  ) HGF = 1.8d0
         IF ( IRH < 70  ) HGF = 1.6d0
         IF ( IRH < 50  ) HGF = 1.0d0

         ! radius of sea-salt aerosol size bins [cm] accounting for
         ! hygroscopic growth
         RAD1 = SALA_REDGE_um(1) * HGF * 1.d-4
         RAD2 = SALA_REDGE_um(2) * HGF * 1.d-4
         RAD3 = SALC_REDGE_um(2) * HGF * 1.d-4

         !----------------------------------
         ! SO2 uptake onto fine particles
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of
	 ! SO2 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1 = 4.D0/(V*GAMMA_SO2)
         A1     = (RAD1/DG)+CONST1
         B1     = (RAD2/DG)+CONST1
         TERM1A = ((B1**2)/2.0d0) - ((A1**2)/2.0d0)
         TERM2A = 2.D0*CONST1*(B1-A1)
         TERM3A = (CONST1**2)*LOG(B1/A1)
         KT1    = 4.D0*PI*N1*(DG**3)*(TERM1A - TERM2A + TERM3A)

         !----------------------------------
         ! SO2 uptake onto coarse particles
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of
	 ! SO2 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2 = 4.D0/(V*GAMMA_SO2)
         A2     = (RAD2/DG)+CONST2
         B2     = (RAD3/DG)+CONST2
         TERM1B = ((B2**2)/2.0d0) - ((A2**2)/2.0d0)
         TERM2B = 2.D0*CONST2*(B2-A2)
         TERM3B = (CONST2**2)*LOG(B2/A2)
         KT2    = 4.D0*PI*N2*(DG**3)*(TERM1B - TERM2B + TERM3B)
         KT     = KT1 + KT2

         !----------------------------------
         ! HNO3 uptake onto fine particles
         !----------------------------------

         ! calculate gas-to-particle rate constant for uptake of
         ! HNO3 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1N = 4.D0/(V*GAMMA_HNO3)
         A1N     = (RAD1/DG)+CONST1N
         B1N     = (RAD2/DG)+CONST1N
         TERM1AN = ((B1N**2)/2.0d0) - ((A1N**2)/2.0d0)
         TERM2AN = 2.D0*CONST1N*(B1N-A1N)
         TERM3AN = (CONST1N**2)*LOG(B1N/A1N)
         KT1N    = 4.D0*PI*N1*(DG**3)*(TERM1AN - TERM2AN + TERM3AN)

         !----------------------------------
         ! HNO3 uptake onto coarse particles
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of
	 ! HNO3 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2N = 4.D0/(V*GAMMA_HNO3)
         A2N     = (RAD2/DG)+CONST2N
         B2N     = (RAD3/DG)+CONST2N
         TERM1BN = ((B2N**2)/2.0d0) - ((A2N**2)/2.0d0)
         TERM2BN = 2.D0*CONST2N*(B2N-A2N)
         TERM3BN = (CONST2N**2)*LOG(B2N/A2N)
         KT2N    = 4.D0*PI*N2*(DG**3)*(TERM1BN - TERM2BN + TERM3BN)


      ELSE

         ! If no alkalinity, set everything to zero
         KT1  = 0.D0
         KT1N = 0.D0
         KT2  = 0.D0
         KT2N = 0.D0

      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_ALK

!------------------------------------------------------------------------------

      SUBROUTINE INIT_SEASALT
!
!******************************************************************************
!  Subroutine INIT_SEASALT initializes and zeroes all module arrays
!  (bmy, 4/26/04, 4/13/05)
!
!  NOTES:
!  (1 ) Now exit if we have allocated arrays before.  Now also allocate
!        ALK_EMIS & N_DENS.  Now reference CMN_SIZE. (bec, bmy, 4/13/05)
!  (2 ) Added SALT_V and DMID (jaegle 5/11/11)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_SEASALT begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Allocate arrays
      ALLOCATE( REDGE( 0:NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REDGE' )
      REDGE = 0d0

      ALLOCATE( RMID( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RMID' )
      RMID = 0d0

      ALLOCATE( SRC( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRC' )
      SRC = 0d0

      ALLOCATE( SRC_N( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRC_N' )
      SRC_N = 0d0

      ALLOCATE( ALK_EMIS( IIPAR, JJPAR, LLTROP, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK_EMIS' )
      ALK_EMIS = 0d0

      ALLOCATE( N_DENS( IIPAR, JJPAR, LLTROP, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'N_DENS' )
      N_DENS = 0d0

      ALLOCATE( SALT_V( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALT_V' )
      SALT_V = 0d0

      ALLOCATE( DMID( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMID' )
      DMID = 0d0

      ! Reset IS_INIT
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_SEASALT

!----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SEASALT
!
!******************************************************************************
!  Subroutine INIT_SEASALT deallocates all module arrays
!  (bmy, 4/26/04, 4/13/05)
!
!  NOTES:
!  (1 ) Now deallocates ALK_EMIS, N_DENS, SRC_N (bec, bmy, 4/13/05)
!  (2 ) Deallocated SALT_V and DMID (jaegle 5/11/11)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_SEASALT begins here!
      !=================================================================
      IF ( ALLOCATED( REDGE    ) ) DEALLOCATE( REDGE    )
      IF ( ALLOCATED( RMID     ) ) DEALLOCATE( RMID     )
      IF ( ALLOCATED( SRC      ) ) DEALLOCATE( SRC      )
      IF ( ALLOCATED( SRC_N    ) ) DEALLOCATE( SRC_N    )
      IF ( ALLOCATED( ALK_EMIS ) ) DEALLOCATE( ALK_EMIS )
      IF ( ALLOCATED( N_DENS   ) ) DEALLOCATE( N_DENS   )
      IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_SEASALT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE SEASALT_MOD
