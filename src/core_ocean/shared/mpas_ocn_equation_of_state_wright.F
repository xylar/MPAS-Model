! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state_wright
!
!> \brief MPAS ocean Wright (1997) equation of state
!> \author Xylar Asay-Davis
!> \date   June 30, 2019
!> \details
!>  This module contains the routines for computing density from
!>  temperature, salinity and depth using an equation of state by Wright (1997)
!>  that approximates the nonlinear equation of state using a set of
!>  rational polynomials.
!>  See Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!
!-----------------------------------------------------------------------

module ocn_equation_of_state_wright

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use mpas_constants
   use ocn_constants

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_equation_of_state_wright_density, &
             ocn_equation_of_state_pressure_grad_wright, &
             ocn_equation_of_state_wright_init

   !*** generic interface for case of density only or density and
   !*** expansion coeffs
   interface ocn_equation_of_state_wright_density
      module procedure ocn_equation_of_state_wright_density_only
      module procedure ocn_equation_of_state_wright_density_exp
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !*** temporary array to hold modified T, S, p for EOS calculation
   !*** declared here to make it thread-shared - could be replaced
   !***   by a local subroutine temporary in a different thread model

   real (kind=RKIND), dimension(:,:), allocatable :: &
      tracerTemp, tracerSalt, boussinesqPres, lambda0, alpha0, p0, &
      interfacePressure, vertPressureIntegral

   !*** valid range of T,S for Wright (1997) EOS
   !***   allowiwng T to become more negative because of sub-ice-shelf cavities

   real (kind=RKIND), parameter :: &
      ocnEqStateTmin = -3.0_RKIND, &! valid pot. temp. range
      ocnEqStateTmax = 30.0_RKIND, &
      ocnEqStateSmin = 28.0_RKIND, &! valid salinity, in psu
      ocnEqStateSmax = 38.0_RKIND

   !***  Wright (1997) constants from Table 1, last column

   real (kind=RKIND), parameter ::      &
      a0 = 7.057924e-4_RKIND,           &
      a1 = 3.480336e-7_RKIND,           &
      a2 = -1.112733e-7_RKIND

   real (kind=RKIND), parameter ::      &
      b0 = 5.790749e8_RKIND,            &
      b1 = 3.516535e6_RKIND,            &
      b2 = -4.002714e4_RKIND,           &
      b3 = 2.084372e2_RKIND,            &
      b4 = 5.944068e5_RKIND,            &
      b5 = -9.643486e3_RKIND

   real (kind=RKIND), parameter ::      &
      c0 = 1.704853e5_RKIND,            &
      c1 = 7.904722e2_RKIND,            &
      c2 = -7.984422_RKIND,             &
      c3 = 5.140652e-2_RKIND,           &
      c4 = -2.302158e2_RKIND,           &
      c5 = -3.079464_RKIND

   !***  6th-order quadrature points and weights (integral from 0 to 1, not -1
   !***  to 1)
   !*** xq = 0.5*(x + 1); wq = 0.5*w where x and w are the usual 6th-order
   !*** weights

   integer, parameter :: nQuad = 5

   real (kind=RKIND), dimension(nQuad), parameter ::                  &
      xq = (/0.0_RKIND, 0.17267316464601146_RKIND, 0.5_RKIND,     &
             0.8273268353539885_RKIND, 1.0_RKIND/),               &
      wq = (/0.05_RKIND, 0.2722222222222222_RKIND,                &
             0.35555555555555557_RKIND, 0.2722222222222222_RKIND, &
             0.05_RKIND/)


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_density_only
!
!> \brief   Computes equation of state
!> \author  Xylar Asay-Davis
!> \date    June 30, 2019
!> \details
!>  This routine computes the density from model temperature,
!>  salinity, and depth using potential-temperature formulation from
!>  Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!>  and a simple linearization of pressure p = -rho0*g*z with the Boussinesq
!>  reference density
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to zero pressure) can be computed
!>      using displacementType = 'absolute'. kDisplaced is ignored.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_density_only(nVertLevels,        &
                               nCells, kDisplaced, displacementType,        &
                               indexT, indexS, tracers, zMid, maxLevelCell, &
                               density, err)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density         ! computed density

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      integer :: &
         iCell, k ! cell and vertical level loop indices

      real (kind=RKIND) :: &
         T, S, T2, T3,     &! adjusted T, S, T squared, T cubed,
         p,                &! Boussinesq pressure at reference level,
         alpha0, lambda0,  &! functions of theta and S defined by Wright (1997)
         p0

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !$omp master
      allocate(tracerTemp(nVertLevels, nCells))
      allocate(tracerSalt(nVertLevels, nCells))
      allocate(boussinesqPres(nVertLevels, nCells))
      !$omp end master
      !$omp barrier

      ! compute tracerTemp, tracerSalt, boussinesqPres
      call compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, tracers, &
                               tracerTemp, tracerSalt)
      call compute_boussinesq_p(nVertLevels, nCells, kDisplaced, &
                                displacementType, zMid, maxLevelCell, &
                                boussinesqPres)

#ifdef MPAS_OPENACC
      !$omp master
      !$acc enter data copyin(density, boussinesqPres, tracerTemp, tracerSalt)
      !$acc parallel loop gang vector collapse(2) &
      !$acc&   present(density, boussinesqPres, tracerTemp, tracerSalt)
#else
      !$omp  do schedule(runtime)
#endif
      do iCell=1,nCells
         do k=1,nVertLevels

            S = tracerSalt(k, iCell)
            T = tracerTemp(k, iCell)
            p = boussinesqPres(k, iCell)

            T2 = T*T
            T3 = T*T2

            p0 = b0 + b1*T + b2*T2 + b3*T3 + b4*S + b5*T*S
            lambda0 = c0 + c1*T + c2*T2 + c3*T3 + c4*S + c5*T*S
            alpha0 = a0 + a1*T + a2*S

            density(k, iCell) = (p + p0)/(lambda0 + alpha0*(p + p0))

         end do
      end do

#ifdef MPAS_OPENACC
      !$acc update host(density)
      !$acc exit data delete(density, boussinesqPres, tracerTemp, tracerSalt)
      !$omp end master
      !$omp barrier
#else
      !$omp end do
#endif

      !$omp master
      deallocate(tracerTemp)
      deallocate(tracerSalt)
      deallocate(boussinesqPres)
      !$omp end master

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_density_only!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_density_exp
!
!> \brief   Computes equation of state with expansion/contraction coeffs
!> \author  Xylar Asay-Davis
!> \date    July 1, 2019
!> \details
!>  This routine computes the density from model temperature,
!>  salinity, and depth using potential-temperature formulation from
!>  Wright (1997), doi: 10.1175/1520-0426(1997)014<0735:AEOSFU>2.0.CO;2
!>  and a simple linearization of pressure p = -rho0*g*z with the Boussinesq
!>  reference density
!>
!>  Density can be computed in-situ using kDisplaced=0 and
!>      displacementType = 'relative'.
!>
!>  Potential density (referenced to zero pressure) can be computed
!>      using displacementType = 'absolute'. kDisplaced is ignored.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_density_exp(nVertLevels,         &
                               nCells, kDisplaced, displacementType,        &
                               indexT, indexS, tracers, zMid, maxLevelCell, &
                               density, thermalExpansionCoeff,              &
                               salineContractionCoeff, err)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced,         &! target layer for displacement
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density         ! computed density

      ! Thermal expansion coeff, $-1/\rho d\rho/dT$ (note negative sign)
      ! Saline contraction coeff, $1/\rho d\rho/dS$
      real (kind=RKIND), dimension(:,:), intent(out) :: &
         thermalExpansionCoeff,  &! Thermal expansion  coeff (alpha)
         salineContractionCoeff   ! Saline contraction coeff (beta)

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, k ! cell and vertical level loop indices

      real (kind=RKIND) :: &
         T, S, T2, T3,     &! adjusted T, S, T squared, T cubed,
         p,                &! Boussinesq pressure at reference level,
         alpha0, lambda0,  &! functions of theta and S defined by Wright (1997)
         p0,               &
         denom, denom2,    &! denominator of rational polynomial for rho
         dalpha0dT,        &! d(alpha0)/d(temperature)
         dalpha0dS,        &! d(alpha0)/d(salinity)
         dlambda0dT,       &! d(lambda0)/d(temperature)
         dlambda0dS,       &! d(lambda0)/d(salinity)
         dp0dT,            &! d(p0)/d(temperature)
         dp0dS,            &! d(p0)/d(salinity)
         drhodT,           &! derivative of density with respect to temperature
         drhodS             ! derivative of density with respect to salinity

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0


      !$omp master
      allocate(tracerTemp(nVertLevels, nCells))
      allocate(tracerSalt(nVertLevels, nCells))
      allocate(boussinesqPres(nVertLevels, nCells))
      !$omp end master
      !$omp barrier

      ! compute tracerTemp, tracerSalt, boussinesqPres
      call compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, tracers, &
                               tracerTemp, tracerSalt)
      call compute_boussinesq_p(nVertLevels, nCells, kDisplaced, &
                                displacementType, zMid, maxLevelCell, &
                                boussinesqPres)

#ifdef MPAS_OPENACC
      !$omp master
      !$acc enter data copyin(density, boussinesqPres, tracerTemp, tracerSalt, &
      !$acc&                  thermalExpansionCoeff,  &
      !$acc&                  salineContractionCoeff))
      !$acc parallel loop gang vector collapse(2) &
      !$acc&   present(density, boussinesqPres, tracerTemp, tracerSalt, &
      !$acc&           thermalExpansionCoeff,         &
      !$acc&           salineContractionCoeff)
#else
      !$omp  do schedule(runtime)
#endif
      do iCell=1,nCells
         do k=1,nVertLevels

            S = tracerSalt(k, iCell)
            T = tracerTemp(k, iCell)
            p = boussinesqPres(k, iCell)

            T2 = T*T
            T3 = T*T2

            p0 = b0 + b1*T + b2*T2 + b3*T3 + b4*S + b5*T*S
            dp0dT = b1 + 2.0_RKIND*b2*T + 3.0_RKIND*b3*T2 + b5*S
            dp0dS = b4 + b5*T

            lambda0 = c0 + c1*T + c2*T2 + c3*T3 + c4*S + c5*T*S
            dlambda0dT = c1 + 2.0_RKIND*c2*T + 3.0_RKIND*c3*T2 + c5*S
            dlambda0dS = c4 + c5*T

            alpha0 = a0 + a1*T + a2*S
            dalpha0dT = a1
            dalpha0dS = a2

            denom = 1.0_RKIND/(lambda0 + alpha0*(p + p0))
            denom2 = denom*denom

            drhodT = (lambda0*dp0dT - (p + p0)*(dlambda0dT + (p + p0)*dalpha0dT))*denom2
            drhodS = (lambda0*dp0dS - (p + p0)*(dlambda0dS + (p + p0)*dalpha0dS))*denom2

            density(k, iCell) = (p + p0)*denom
            thermalExpansionCoeff(k, iCell) = -drhodT/density(k, iCell)
            salineContractionCoeff(k, iCell) = drhodS/density(k, iCell)

         end do
      end do

#ifdef MPAS_OPENACC
      !$acc update host(density,                &
      !$acc&            thermalExpansionCoeff,  &
      !$acc&            salineContractionCoeff)
      !$acc exit data delete(density, boussinesqPres, tracerTemp, tracerSalt, &
      !$acc&                 thermalExpansionCoeff,  &
      !$acc&                 salineContractionCoeff)
      !$omp end master
      !$omp barrier
#else
      !$omp end do
#endif

      !$omp master
      deallocate(tracerTemp)
      deallocate(tracerSalt)
      deallocate(boussinesqPres)
      !$omp end master

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_density_exp!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_pressure_grad_wright
!
!> \brief   Computes the semi-analytic horizontal pressure-gradient force
!> \author  Xylar Asay-Davis
!> \date    July 2019
!> \details
!>  Performs the contour integral around a control volume surrounding an edge
!>  to get the horizontal pressure-gradient force (HPGF).  The integral is
!>  performed assuming vertically constant T and S in each layer but taking
!>  into account that density (and therefore pressure) varies with depth within
!>  a layer, in a Boussinesq adaptation of Adcroft et al. (2008, doi:
!>  10.1016/j.ocemod.2008.02.001)
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_pressure_grad_wright(meshPool, statePool, &
                                                         diagnosticsPool,     &
                                                         timeLevelIn,         &
                                                         tend_normalVelocity, &
                                                         err)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information

      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state fields

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend_normalVelocity          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! returned error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: surfacePressure, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer, dimension(:), pointer :: nEdgesArray, nCellsArray
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      integer :: nEdges, nCells, nVertLevels, indexT, indexS
      integer, pointer :: indexTptr, indexSptr, nVertLevelsPtr

      real (kind=RKIND) :: &
         T, T2, T3, S, P, alpha, Lambda, deltaPi, Pi_bar, eps, eps2, series

      integer :: &
         iCell, iEdge, k

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

      call mpas_pool_get_array(diagnosticsPool, 'surfacePressure', &
                                                 surfacePressure)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, &
                               timeLevelIn)

      call mpas_pool_get_subpool  (statePool,   'tracers', &
                                                 tracersPool)
      call mpas_pool_get_array    (tracersPool, 'activeTracers', &
                                                 activeTracers, timeLevelIn)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTptr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSptr)
      call mpas_pool_get_dimension(meshPool,    'nVertLevels', &
                                                 nVertLevelsPtr)


      ! Only need pressure gradient on 0 halo edges
      nEdges = nEdgesArray( 1 )
      ! Only need tracers on 0 and 1 halo cells
      nCells = nCellsArray( 2 )

      indexT      = indexTptr
      indexS      = indexSptr
      nVertLevels = nVertLevelsPtr

      !$omp master
      allocate(tracerTemp(nVertLevels, nCells))
      allocate(tracerSalt(nVertLevels, nCells))
      allocate(boussinesqPres(nVertLevels+1, nCells))
      allocate(p0(nVertLevels, nCells))
      allocate(lambda0(nVertLevels, nCells))
      allocate(alpha0(nVertLevels, nCells))
      allocate(interfacePressure(nVertLevels+1, nCells))
      allocate(vertPressureIntegral(nVertLevels, nCells))
      !$omp end master
      !$omp barrier

      ! compute tracerTemp, tracerSalt
      call compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, &
                               activeTracers, tracerTemp, tracerSalt)

      call compute_vert_pressure_integral(nVertLevels, nCells, &
                                          surfacePressure, layerThickness, &
                                          tracerTemp, tracerSalt, &
                                          boussinesqPres, p0, lambda0,  &
                                          alpha0, interfacePressure, &
                                          vertPressureIntegral)

      call compute_total_pressure_integral(nVertLevels, nEdges, cellsOnEdge, &
                                           dcEdge, edgeMask, boussinesqPres, &
                                           p0, lambda0, alpha0, &
                                           interfacePressure, &
                                           vertPressureIntegral, &
                                           tend_normalVelocity)

      !$omp master
      deallocate(tracerTemp)
      deallocate(tracerSalt)
      deallocate(boussinesqPres)
      deallocate(interfacePressure)
      deallocate(vertPressureIntegral)
      deallocate(p0)
      deallocate(lambda0)
      deallocate(alpha0)
      !$omp end master

      !-----------------------------------------------------------------

   end subroutine ocn_equation_of_state_pressure_grad_wright!}}}

!***********************************************************************
!
!  routine ocn_equation_of_state_wright_init
!
!> \brief   Initializes Wright equation of state
!> \author  Xylar Asay-Davis
!> \date    June 30, 2019
!> \details
!>  This routine initializes a variety of quantities related to
!>  the equation of state from Wright (1997).
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_wright_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: &
         domain        ! domain containing all state, mesh info

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err  ! error flag

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_wright_init!}}}

!***********************************************************************
!
!  routine compute_bounded_T_S
!
!> \brief   Compute T, S within valid bounds of Wright EOS
!> \author  Xylar Asay-Davis
!> \date    July 8, 2019
!
!-----------------------------------------------------------------------

   subroutine compute_bounded_T_S(nVertLevels, nCells, indexT, indexS, &
                                  tracers, tracerTemp, tracerSalt)!{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nVertLevels,        &! max number of verical levels
         nCells,             &! number of horizontal cells
         indexT,             &! temperature index in tracer array
         indexS               ! salinity    index in tracer array

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers              ! array of tracers including T,S

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         tracerTemp, tracerSalt

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      integer :: &
         iCell, k

      real (kind=RKIND) :: &
         T, S               ! adjusted T, S

      !$omp do schedule(runtime)
      do iCell=1,nCells
         do k = 1, nVertLevels
            T = min(tracers(indexT, k, iCell), ocnEqStateTmax)
            S = min(tracers(indexS, k, iCell), ocnEqStateSmax)
            tracerTemp(k, iCell) = max(T, ocnEqStateTmin)
            tracerSalt(k, iCell) = max(S, ocnEqStateSmin)
         end do
      end do
      !$omp end do

   !--------------------------------------------------------------------

   end subroutine compute_bounded_T_S!}}}

!***********************************************************************
!
!  routine compute_boussinesq_p
!
!> \brief   Compute Boussinesq pressure, referenced to requested layer
!> \author  Xylar Asay-Davis
!> \date    July 8, 2019
!
!-----------------------------------------------------------------------

   subroutine compute_boussinesq_p(nVertLevels, nCells, kDisplaced,      &
                                   displacementType, zMid, maxLevelCell, &
                                   boussinesqPres)
   !{{{
   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horizontal cells
         nVertLevels,        &! max number of verical levels
         kDisplaced           ! target layer for displacement

      character(len=*), intent(in) :: &
         displacementType     ! choice of displacement

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         zMid

      integer, dimension(:), intent(in) :: &
         maxLevelCell

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         boussinesqPres

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer ::          &
         iCell, k,        &! cell and vertical level loop indices
         kPressure         ! index to determine ref level for pressure

      !-----------------------------------------------------------------
      !  if displacementType = 'relative', density is calculated
      !     with pressure referenced to level k + kDisplaced
      !     If kDisplaced=0, in-situ density is returned (no displacement)
      !     If kDisplaced/=0, a displaced density is returned
      !
      !  if displacementType = 'absolute', potential density is calculated
      !     referenced to zero pressure

      if (displacementType == 'relative') then

         !$omp do schedule(runtime)
         do iCell=1,nCells
            do k = 1, nVertLevels
               kPressure = min(k + kDisplaced, maxLevelCell(iCell))
               kPressure = max(kPressure, 1)
               boussinesqPres(k, iCell) = -rho_sw*gravity*zMid(kPressure, iCell)
            end do
         end do
         !$omp end do

      else ! displacementType == 'absolute'

         ! The reference gauge pressure is hard coded to the surface (zero)
         ! for now but this could be replaced by a namelist option with a
         ! different reference pressure in the future
         boussinesqPres(:, :) = 0.0_RKIND

      endif

   !--------------------------------------------------------------------

   end subroutine compute_boussinesq_p!}}}


!***********************************************************************
!
!  routine compute_vert_pressure_integral
!
!> \brief   Computes pressure-gradient terms on vertical faces
!> \author  Xylar Asay-Davis
!> \date    July 2019
!
!-----------------------------------------------------------------------

   subroutine compute_vert_pressure_integral(nVertLevels, nCells, &
                                             surfacePressure, layerThickness, &
                                             tracerTemp, tracerSalt, &
                                             boussinesqPres, p0, lambda0,  &
                                             alpha0, interfacePressure, &
                                             vertPressureIntegral)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nVertLevels, nCells

      real (kind=RKIND), dimension(:), intent(in) :: surfacePressure
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, tracerTemp, tracerSalt

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         boussinesqPres, p0, lambda0, alpha0, &
         interfacePressure, vertPressureIntegral

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND) :: &
         T, T2, T3, S, P, alpha, Lambda, deltaPi, Pi_bar, eps, eps2, series

      integer :: &
         iCell, iEdge, k

      !$omp  do schedule(runtime)
      do iCell=1,nCells
         boussinesqPres(1, iCell) = surfacePressure(iCell)
         interfacePressure(1, iCell) = surfacePressure(iCell)
         do k=1,nVertLevels

            S = tracerSalt(k, iCell)
            T = tracerTemp(k, iCell)

            T2 = T*T
            T3 = T*T2

            p0(k, iCell) = b0 + b1*T + b2*T2 + b3*T3 + b4*S + b5*T*S
            lambda0(k, iCell) = c0 + c1*T + c2*T2 + c3*T3 + c4*S + c5*T*S
            alpha0(k, iCell) = a0 + a1*T + a2*S

            deltaPi = rho_sw*gravity*layerThickness(k, iCell)
            boussinesqPres(k+1, iCell) = boussinesqPres(k, iCell) + deltaPi
            Pi_bar = 0.5_RKIND*(boussinesqPres(k+1, iCell) &
                                + boussinesqPres(k, iCell))

            alpha = alpha0(k, iCell)*rho_sw
            Lambda = lambda0(k, iCell)/alpha0(k, iCell)
            P = Pi_bar + p0(k, iCell) + Lambda
            eps = deltaPi/(2.0_RKIND*P)
            eps2 = eps*eps
            series = 1.0_RKIND + eps2*(1.0_RKIND/3.0_RKIND &
               + eps2*(1.0_RKIND/5.0_RKIND + eps2*(1.0_RKIND/7.0_RKIND &
               + eps2*(1.0_RKIND/9.0_RKIND + eps2*(1.0_RKIND/11.0_RKIND)))))

            interfacePressure(k+1, iCell) = interfacePressure(k, iCell) &
               + deltaPi/alpha*(1.0_RKIND - Lambda/P * series)

            ! -\rho_0 g \int_{z_{t}}^{z_{b}} p dz
            vertPressureIntegral(k, iCell) = deltaPi &
                * (interfacePressure(k, iCell) + deltaPi/(2.0_RKIND*alpha) &
                   + Lambda/alpha*(1 - (1 + eps) * series))

         end do
      end do

      !$omp end do

      !-----------------------------------------------------------------

   end subroutine compute_vert_pressure_integral!}}}

!***********************************************************************
!
!  routine compute_total_pressure_integral
!
!> \brief   Computes pressure-gradient terms along layer interfaces, and total
!> \author  Xylar Asay-Davis
!> \date    July 2019
!
!-----------------------------------------------------------------------

   subroutine compute_total_pressure_integral(nVertLevels, nEdges, cellsOnEdge, &
                                              dcEdge, edgeMask, boussinesqPres, &
                                              p0, lambda0, alpha0, &
                                              interfacePressure, &
                                              vertPressureIntegral, &
                                              tend_normalVelocity)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nVertLevels, nEdges

      integer, dimension(:,:), intent(in) :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), intent(in) :: dcEdge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         boussinesqPres, p0, lambda0, alpha0, interfacePressure, &
         vertPressureIntegral

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend_normalVelocity          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(nQuad) :: &
         p0x, lambda0x, Px, alphax, Lambdax, deltaPix, Pi_barx, &
         epsx, eps2x, seriesx

      real (kind=RKIND) :: &
         layerIntegral, layerPressureIntegralK, layerPressureIntegralKP1

      integer :: &
         iEdge, k, cell1, cell2, iq

      !$omp  do schedule(runtime)
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         layerIntegral = 0.5_RKIND*(boussinesqPres(1, cell1) &
                                    + boussinesqPres(1, cell2))

         layerPressureIntegralK = (boussinesqPres(1, cell2) &
                                   - boussinesqPres(1, cell1)) &
                                  * layerIntegral

         do k=1,nVertLevels

            deltaPix(1) = boussinesqPres(k+1, cell1) &
               - boussinesqPres(k, cell1)
            deltaPix(nQuad) = boussinesqPres(k+1, cell2) &
               - boussinesqPres(k, cell2)

            Pi_barx(1) = 0.5_RKIND*(boussinesqPres(k+1, cell1) &
                                    + boussinesqPres(k, cell1))
            Pi_barx(nQuad) = 0.5_RKIND*(boussinesqPres(k+1, cell2) &
                                    + boussinesqPres(k, cell2))

            alphax(1) = alpha0(k, cell1)
            alphax(nQuad) = alpha0(k, cell2)

            lambda0x(1) = lambda0(k, cell1)
            lambda0x(nQuad) = lambda0(k, cell2)

            p0x(1) = p0(k, cell1)
            p0x(nQuad) = p0(k, cell2)

            do iq=2,nQuad-1
               deltaPix(iq) = deltaPix(1) + xq(iq)*(deltaPix(nQuad) - deltaPix(1))
               Pi_barx(iq) = Pi_barx(1) + xq(iq)*(Pi_barx(nQuad) - Pi_barx(1))
               alphax(iq) = alphax(1) + xq(iq)*(alphax(nQuad) - alphax(1))
               lambda0x(iq) = lambda0x(1) + xq(iq)*(lambda0x(nQuad) - lambda0x(1))
               p0x(iq) = p0x(1) + xq(iq)*(p0x(nQuad) - p0x(1))
            end do

            Lambdax(:) = lambda0x(:)/alphax(:)
            alphax(:) = alphax(:)*rho_sw
            Px(:) = Pi_barx(:) + p0x(:) + Lambdax(:)
            epsx(:) = deltaPix(:)/(2.0_RKIND*Px(:))
            eps2x(:) = epsx(:)*epsx(:)
            seriesx(:) = 1.0_RKIND + eps2x(:)*(1.0_RKIND/3.0_RKIND &
               + eps2x(:)*(1.0_RKIND/5.0_RKIND + eps2x(:)*(1.0_RKIND/7.0_RKIND &
               + eps2x(:)*(1.0_RKIND/9.0_RKIND + eps2x(:)*(1.0_RKIND/11.0_RKIND)))))


            layerIntegral = layerIntegral + &
               sum(wq(:) * deltaPix(:)/alphax(:) &
                   * (1.0_RKIND - Lambdax(:)/Px(:)*seriesx(:)))

            layerPressureIntegralKp1 = (boussinesqPres(k+1, cell2) &
                                        - boussinesqPres(k+1, cell1)) &
                                       * layerIntegral

            tend_normalVelocity(k, iEdge) = tend_normalVelocity(k, iEdge) &
               + edgeMask(k, iEdge)/(dcEdge(iEdge)*rho_sw**2*gravity) &
               * (vertPressureIntegral(k, cell2) &
                  - vertPressureIntegral(k, cell1) &
                  + layerPressureIntegralK - layerPressureIntegralKp1)

            layerPressureIntegralK = layerPressureIntegralKp1

         end do
      end do

      !$omp end do

      !-----------------------------------------------------------------

   end subroutine compute_total_pressure_integral!}}}

!***********************************************************************

end module ocn_equation_of_state_wright

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
