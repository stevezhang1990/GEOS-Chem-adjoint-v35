! $Id: diag_oh_mod.f,v 1.2 2012/03/01 22:00:26 daven Exp $
      MODULE DIAG_OH_MOD
!
!******************************************************************************
!  Module DIAG_OH_MOD contains routines and variables to archive OH mass
!  and air mass concentrations.  These are then used to print out the mass-
!  weighted mean OH concentration in 1e5 molec/cm3.  This is a metric of
!  how certain chemisry simulations are performing. (bmy, 7/20/04, 6/24/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AIR_MASS (REAL*8) : Array used to sum mean air mass [molec/cm3]
!  (2 ) OH_MASS  (REAL*8) : Array used to sum mean OH mass  [molec/cm3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_DIAG_OH        : Driver routine for mean OH diagnostic (fullchem)
!  (2 ) DO_DIAG_OH_CH4    : Driver routine for mean OH diagnostic (CH4 sim)
!  (3 ) PRINT_DIAG_OH     : Prints the mean OH concentration [1e5 molec/cm3]
!  (4 ) INIT_DIAG_OH      : Allocates and zeroes module arrays
!  (5 ) CLEANUP_DIAG_OH   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "diag_oh_mod.f":
!  ============================================================================
!  (1 ) comode_mod.f      : Module w/ SMVGEAR allocatable arrays
!  (2 ) error_mod.f       : Module w/ I/O error and NaN check routines
!  (3 ) logical_mod.f     : Module w/ GEOS-CHEM logical switches
!  (4 ) tracer_mod.f      : Module w/ GEOS-CHEM tracer array STT etc.
!  (5 ) tracerid_mod.f    : Module w/ pointers to tracers & emissions
!
!  NOTES:
!  (1 ) Remove code for obsolete CO-OH simulation (bmy, 6/24/05)
!  (2 ) Add OH_LOSS array (kjw, dkh, 02/12/12, adj32_023)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "diag_oh_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_DIAG_OH
      PUBLIC :: DO_DIAG_OH
      PUBLIC :: DO_DIAG_OH_CH4
      PUBLIC :: INIT_DIAG_OH
      PUBLIC :: PRINT_DIAG_OH

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      LOGICAL             :: DO_SAVE_OH

      ! Arrays
      REAL*8, ALLOCATABLE :: OH_MASS(:,:,:)
      REAL*8, ALLOCATABLE :: AIR_MASS(:,:,:)
      REAL*8, ALLOCATABLE :: OH_LOSS(:,:,:)
      REAL*8, ALLOCATABLE :: OHCH4_LOSS(:,:,:)
      REAL*8, ALLOCATABLE :: CH4_MASS(:,:,:)
      REAL*8, ALLOCATABLE :: CH4_TROPMASS(:,:,:)
      REAL*8, ALLOCATABLE :: CH4_EMIS(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_DIAG_OH
!
!******************************************************************************
!  Subroutine DO_DIAG_OH sums the OH and air mass (from SMVGEAR arrays) for
!  the mean OH concentration diagnostic. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! Reference to F90 modules
      USE COMODE_MOD,   ONLY : AIRDENS, CSPEC, JLOP, T3, VOLUME
      USE TRACERID_MOD, ONLY : IDOH

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NPVERT, NLAT, NLONG

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I,     J,       L,       JLOOP
      REAL*8        :: XLOSS, XOHMASS, XAIRMASS

      !=================================================================
      ! DO_DIAG_OH begins here!
      !=================================================================

      ! Safety valve -- avoid seg faults
      IF ( .not. DO_SAVE_OH ) RETURN

      ! Loop over boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, XAIRMASS, XOHMASS, XLOSS )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, NPVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Cycle if this isn't a valid SMVGEAR gridbox
         IF ( JLOOP > 0 ) THEN

            ! Sum air mass term into AIR_MASS array
            XAIRMASS        = AIRDENS(JLOOP)    * VOLUME(JLOOP)
            AIR_MASS(I,J,L) = AIR_MASS(I,J,L)   + XAIRMASS

            ! Sum OH mass term into OH_MASS array
            XOHMASS         = CSPEC(JLOOP,IDOH) * XAIRMASS
            OH_MASS(I,J,L)  = OH_MASS(I,J,L)    + XOHMASS

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DO_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE DO_DIAG_OH_CH4( I, J, L, XOHMASS, XAIRMASS, XLOSS,
     &                   XCH4LOSS, XCH4TROPMASS, XCH4EMIS, XCH4MASS )
!
!******************************************************************************
!  Subroutine DO_DIAG_OH_CH4 passes the OH loss, OH mass, and air mass terms
!  from "global_ch4_mod.f" to "diag_oh_mod.f" (bmy, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, & altitude indices
!  (4  ) XOHMASS  (REAL*8 ) : OH mass term from "global_ch4_mod.f"
!  (5  ) XAIRMASS (REAL*8 ) : air mass term from "global_ch4_mod.f"
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I          ! Longitude index
      INTEGER, INTENT(IN) :: J          ! Latitude index
      INTEGER, INTENT(IN) :: L          ! Level index
      REAL*8,  INTENT(IN) :: XOHMASS    ! OH Mass  (from global_ch4_mod.f)
      REAL*8,  INTENT(IN) :: XAIRMASS   ! Air mass (from global_ch4_mod.f)
      REAL*8,  INTENT(IN) :: XLOSS      ! Loss of ch3ccl3 by OH
      REAL*8,  INTENT(IN) :: XCH4LOSS   ! Loss of ch4 by OH
      REAL*8,  INTENT(IN) :: XCH4MASS   ! CH4 Mass  (from global_ch4_mod.f)
      REAL*8,  INTENT(IN) :: XCH4TROPMASS   ! CH4 Mass  (from global_ch4_mod.f)
      REAL*8,  INTENT(IN) :: XCH4EMIS   ! CH4 emissions

      !=================================================================
      ! DO_DIAG_OH_CH4 begins here!
      !=================================================================

      ! Sum air mass & OH mass into arrays
      AIR_MASS(I,J,L)     = AIR_MASS(I,J,L)     + XAIRMASS
      OH_MASS(I,J,L)      = OH_MASS(I,J,L)      + XOHMASS
      OH_LOSS(I,J,L)      = OH_LOSS(I,J,L)      + XLOSS
      OHCH4_LOSS(I,J,L)   = OHCH4_LOSS(I,J,L)   + XCH4LOSS
      CH4_MASS(I,J,L)     = CH4_MASS(I,J,L)     + XCH4MASS
      CH4_TROPMASS(I,J,L) = CH4_TROPMASS(I,J,L) + XCH4TROPMASS
      CH4_EMIS(I,J,L)     = CH4_EMIS(I,J,L)     + XCH4EMIS

      ! Return to calling program
      END SUBROUTINE DO_DIAG_OH_CH4

!------------------------------------------------------------------------------

      SUBROUTINE PRINT_DIAG_OH
!
!******************************************************************************
!  Subroutine PRINT_DIAG_OH prints the mass-weighted OH concentration at
!  the end of a simulation. (bmy, 10/21/03, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! Reference to F90 modules
      USE TRACER_MOD,   ONLY : ITS_A_CH4_SIM

      ! Local variables
      REAL*8 :: SUM_OHMASS, SUM_MASS, SUM_OHLOSS, OHCONC, LIFETIME
      REAL*8 :: SUM_OHCH4LOSS, SUM_CH4MASS, SUM_CH4EMIS, SUM_CH4TROPMASS

      !=================================================================
      ! PRINT_DIAG_OH begins here!
      !=================================================================

      ! Return if this diagnostic is turned off
      IF ( .not. DO_SAVE_OH ) RETURN

      ! Total Mass-weighted OH [molec OH/cm3] * [molec air]
      SUM_OHMASS = SUM( OH_MASS )

      ! Atmospheric air mass [molec air]
      SUM_MASS   = SUM( AIR_MASS )
      
      ! OH Loss from CH4 + OH [molec / box / s]
      SUM_OHCH4LOSS = SUM( OHCH4_LOSS )

      ! Atmospheric mass of CH4
      SUM_CH4MASS = SUM( CH4_MASS )

      ! Atmospheric mass of tropospheric CH4
      SUM_CH4TROPMASS = SUM( CH4_TROPMASS )

      ! Atmospheric ch4 emissions
      SUM_CH4EMIS = SUM( CH4_EMIS )

      ! Avoid divide-by-zero errors
      IF ( SUM_MASS > 0d0 ) THEN

         ! Divide OH by [molec air] and report as [1e5 molec/cm3]
         OHCONC = ( SUM_OHMASS / SUM_MASS ) / 1d5

         ! Write value to log file
         WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
         WRITE( 6, *       ) 'ND23: Mass-Weighted OH Concentration'
         WRITE( 6, *       ) 'Mean OH = ', OHCONC, ' [1e5 molec/cm3]'
         WRITE( 6, '(  a)' ) REPEAT( '=', 79 )

         ! Avoid divide-by-zero errors
         IF ( ITS_A_CH4_SIM() ) THEN
            IF ( SUM_OHLOSS > 0 ) THEN

               ! Mass weighted lifetimes printed below
               WRITE( 6, * ) 'All lifetimes printed below ' //
     &                       'are mass-weighted'
               WRITE( 6, '(  a)' ) REPEAT( '-', 79 )

               ! Calculate CH4 lifetime w/r/t OH loss [years]
               LIFETIME = ( SUM_MASS / SUM_OHCH4LOSS ) / 
     &                    ( 3600d0*365d0*24d0 )
               ! Write value to log file
               WRITE( 6, *       ) 'Methane (CH4)'
               WRITE( 6, *       ) 'Tropospheric Lifetime w/r/t OH  = ', 
     &                                 LIFETIME, ' [years]'
               WRITE( 6, '(  a)' ) REPEAT( '=', 79 )

               ! Calculate CH4 lifetime  [years]
               LIFETIME = ( SUM_CH4TROPMASS / SUM_CH4EMIS ) / 
     &                    ( 3600d0*365d0*24d0 )
               ! Write value to log file
               WRITE( 6, *       ) 'Methane (CH4)'
               WRITE( 6, *       ) 'Tropospheric Lifetime (total)  = ', 
     &                                 LIFETIME, ' [years]'
               WRITE( 6, '(  a)' ) REPEAT( '=', 79 )

               ! Calculate CH4 lifetime  [years]
               LIFETIME = ( SUM_CH4MASS / SUM_CH4EMIS ) / 
     &                    ( 3600d0*365d0*24d0 )
               ! Write value to log file
               WRITE( 6, *       ) 'Methane (CH4)'
               WRITE( 6, *       ) 'Global Lifetime (total)  = ', 
     &                                 LIFETIME, ' [years]'
               WRITE( 6, '(  a)' ) REPEAT( '=', 79 )

            ENDIF
         ENDIF
      ELSE

         ! Write error msg if SUM_MASS is zero
         WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(  a)' ) 'Could not print mass-weighted OH!'
         WRITE( 6, '(  a)' ) 'Atmospheric air mass is zero!'
         WRITE( 6, '(  a)' ) REPEAT( '=', 79 )

      ENDIF

      ! Return to MAIN program
      END SUBROUTINE PRINT_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG_OH
!
!******************************************************************************
!  Subroutine INIT_DIAG_OH initializes all module arrays.
!  (bmy, 7/20/04, 6/24/05)
!
!  NOTES:
!  (1 ) Remove references to CO-OH simulation and to CMN_DIAG (bmy, 6/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LCHEM
      USE TRACER_MOD,  ONLY : ITS_A_FULLCHEM_SIM, ITS_A_CH4_SIM

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER :: AS, LMAX

      !=================================================================
      ! INIT_DIAG_OH begins here!
      !=================================================================

      ! Initialize
      DO_SAVE_OH = .FALSE.

      ! Return if we are not doing chemistry
      IF ( .not. LCHEM ) RETURN

      ! Set vertical levels and decide whether to print CH3CCl3
      ! lifetime or just mean mass-weighted OH concentration
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Fullchem: tropopshere only
         LMAX       = LLTROP
         DO_SAVE_OH = .TRUE.

      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         ! CH4: all levels
         LMAX       = LLPAR
         DO_SAVE_OH = .TRUE.

      ENDIF

      ! Echo info
      WRITE( 6, 100 ) DO_SAVE_OH
 100  FORMAT( /, 'Turn on Mean OH diagnostic (ND23)? :', L5 )

      ! Return if we aren't saving mean OH
      IF ( .not. DO_SAVE_OH ) RETURN

      !=================================================================
      ! Allocate arrays
      !=================================================================

      ! Air mass array
      ALLOCATE( AIR_MASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIR_MASS' )
      AIR_MASS = 0d0

      ! OH mass array
      ALLOCATE( OH_MASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH_MASS' )
      OH_MASS = 0d0

      ! OH LOSS array
      ALLOCATE( OH_LOSS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH_LOSS' )
      OH_LOSS = 0d0

      ALLOCATE( OHCH4_LOSS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OHCH4_LOSS' )
      OHCH4_LOSS = 0d0

      ALLOCATE( CH4_MASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_MASS' )
      CH4_MASS = 0d0

      ALLOCATE( CH4_TROPMASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_TROPMASS' )
      CH4_TROPMASS = 0d0

      ALLOCATE( CH4_EMIS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_EMIS' )
      CH4_EMIS = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG_OH
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG_OH deallocates all module arrays. (bmy, 7/20/04)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG_OH begins here!
      !=================================================================
      IF ( ALLOCATED( OH_MASS  ) ) DEALLOCATE( OH_MASS  )
      IF ( ALLOCATED( AIR_MASS ) ) DEALLOCATE( AIR_MASS )
      IF ( ALLOCATED( OH_LOSS  ) ) DEALLOCATE( OH_LOSS  )
      IF ( ALLOCATED( OHCH4_LOSS  ) ) DEALLOCATE( OHCH4_LOSS  )
      IF ( ALLOCATED( CH4_MASS  ) ) DEALLOCATE( CH4_MASS  )
      IF ( ALLOCATED( CH4_TROPMASS  ) ) DEALLOCATE( CH4_TROPMASS  )
      IF ( ALLOCATED( CH4_EMIS  ) ) DEALLOCATE( CH4_EMIS  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG_OH

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG_OH_MOD
