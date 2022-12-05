!> Calculates and applies diffusive fluxes as a parameterization of horizontal mixing (non-neutral) by
!! mesoscale eddies near the top and bottom (to be implemented) boundary layers of the ocean.

module MOM_hor_bnd_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE
use MOM_checksums,             only : hchksum, uvchksum
use MOM_domains,               only : pass_var, sum_across_PEs
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_diag_vkernels,         only : reintegrate_column
use MOM_error_handler,         only : MOM_error, FATAL, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type, log_param
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, remapping_core_h
use MOM_remapping,             only : remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_io,                    only : stdout, stderr, MOM_read_data, slasher
use MOM_io,                    only : file_exists, field_exists, field_size
use MOM_string_functions,      only : extractWord

implicit none ; private

public near_boundary_unit_tests, hor_bnd_diffusion, hor_bnd_diffusion_init
public boundary_k_range, hor_bnd_diffusion_end

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, public, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface bopundary
integer, public, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary
#include <MOM_memory.h>

!> Sets parameters for horizontal boundary mixing module.
type, public :: hbd_CS ; private
  logical :: debug           !< If true, write verbose checksums for debugging.
  integer :: deg             !< Degree of polynomial reconstruction.
  integer :: hbd_nk          !< Maximum number of levels in the HBD grid [nondim]
  integer :: surface_boundary_scheme !< Which boundary layer scheme to use
                             !! 1. ePBL; 2. KPP
  logical :: limiter         !< Controls whether a flux limiter is applied in the
                             !! native grid (default is true).
  logical :: limiter_remap   !< Controls whether a flux limiter is applied in the
                             !! remapped grid (default is false).
  logical :: transition      !< Controls whether diffusion will transition from neutral to
                             !! horizontal (default is false).
  real    :: delta_h         !< The minimum change in thickness between to adjacent cells used to
                             !! calculate the minimum isopycnal slope (delta_h/dx) that can be
                             !! represented in the transition layer.
  real    :: H_subroundoff   !< A thickness that is so small that it can be added to a thickness of
                             !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
                             !! If Angstrom is 0 or exceedingly small, this is negligible compared to 1e-17 m.
  real :: min_thickness      !< Minimum thickness allowed when building HBD grid [H ~> m or kg m-2].

  ! HBD dynamic grids
  real,    allocatable, dimension(:,:,:,:) :: hbd_grd_u   !< HBD thicknesses at t-points adjecent to
                                                          !! u-points                       [H ~> m or kg m-2]
  real,    allocatable, dimension(:,:,:,:) :: hbd_grd_v   !< HBD thicknesses at t-points adjacent to
                                                          !! v-points (left and right)      [H ~> m or kg m-2]
  integer, allocatable, dimension(:,:)     :: hbd_u_kmax  !< Maximum vertical index in hbd_grd_u      [nondim]
  integer, allocatable, dimension(:,:)     :: hbd_v_kmax  !< Maximum vertical index in hbd_grd_v      [nondim]
  ! neutral slopes at BLD base
  real,    allocatable, dimension(:,:)     :: slp_x       !< Zonal isopycnal slope at the base of the boundary
                                                          !! layer                                    [nondim]
  real,    allocatable, dimension(:,:)     :: slp_y       !< Meridional isopycnal slope at the base of the
                                                          !! boundary layer                           [nondim]
  real,    allocatable, dimension(:,:)     :: hbd_slp_x   !< Zonal isopycnal slope at the base of the boundary
                                                          !! layer used in the HBD module             [nondim]
  real,    allocatable, dimension(:,:)     :: hbd_slp_y   !< Meridional isopycnal slope at the base of the
                                                          !! boundary layer used in the HBD module    [nondim]

  type(remapping_CS)    :: remap_CS !< Control structure to hold remapping configuration.
  type(KPP_CS),           pointer :: KPP_CSp => NULL()  !< KPP control structure needed to get BLD.
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL() !< ePBL control structure needed to get BLD.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_slp_x = -1, id_slp_y = -1
  integer :: id_hbd_slp_x = -1, id_hbd_slp_y = -1
  !>@}

end type hbd_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40) :: mdl = "MOM_hor_bnd_diffusion" !< Name of this module
integer :: id_clock_hbd                            !< CPU clock for hbd

contains

!> Initialization routine that reads runtime parameters and sets up pointers to other control structures that might be
!! needed for horizontal boundary diffusion.
logical function hor_bnd_diffusion_init(Time, G, GV, US, param_file, diag, diabatic_CSp, CS)
  type(time_type), target,          intent(in)    :: Time          !< Time structure
  type(ocean_grid_type),            intent(in)    :: G             !< Grid structure
  type(verticalGrid_type),          intent(in)    :: GV            !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US            !< A dimensional unit scaling type
  type(param_file_type),            intent(in)    :: param_file    !< Parameter file structure
  type(diag_ctrl), target,          intent(inout) :: diag          !< Diagnostics control structure
  type(diabatic_CS),                pointer       :: diabatic_CSp  !< KPP control structure needed to get BLD
  type(hbd_CS),                     pointer       :: CS            !< Horizontal boundary mixing control structure

  ! local variables
  logical :: boundary_extrap   ! controls if boundary extrapolation is used in the HBD code
  character(len=80)  :: string ! Temporary strings

  if (ASSOCIATED(CS)) then
    call MOM_error(FATAL, "hor_bnd_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call get_param(param_file, mdl, "USE_HORIZONTAL_BOUNDARY_DIFFUSION", hor_bnd_diffusion_init, &
                 default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
           "This module implements horizontal diffusion of tracers near boundaries", &
           all_default=.not.hor_bnd_diffusion_init)
  call get_param(param_file, mdl, "USE_HORIZONTAL_BOUNDARY_DIFFUSION", hor_bnd_diffusion_init, &
                 "If true, enables the horizontal boundary tracer's diffusion module.", &
                 default=.false.)
  if (.not. hor_bnd_diffusion_init) return

  allocate(CS)
  CS%diag => diag
  CS%H_subroundoff = GV%H_subroundoff
  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

  call get_param(param_file, mdl, "USE_HORIZONTAL_BOUNDARY_DIFFUSION", hor_bnd_diffusion_init, &
                 "If true, enables the horizontal boundary tracer's diffusion module.", &
                 default=.false.)

  call get_param(param_file, mdl, "HBD_NK", CS%hbd_nk, &
                 "The maximum number of vertical layers in the HBD grid.", units="nondim", &
                 default=(GV%ke*2)+2)

  call get_param(param_file, mdl, "HBD_MIN_THICKNESS", CS%min_thickness, &
                 "When regridding the HBD grid, this is the minimum layer "//&
                 "thickness allowed.", units="m", scale=GV%m_to_H, &
                 default=0.001)

  ! allocate the hbd grids and k_max
  allocate(CS%hbd_grd_u(SZIB_(G),SZJ_(G),CS%hbd_nk,2), source=0.0)
  allocate(CS%hbd_grd_v(SZI_(G),SZJB_(G),CS%hbd_nk,2), source=0.0)
  allocate(CS%hbd_u_kmax(SZIB_(G),SZJ_(G)), source=0)
  allocate(CS%hbd_v_kmax(SZI_(G),SZJB_(G)), source=0)

  CS%surface_boundary_scheme = -1
  !GMM, lines below should be commented out in the idealized experiments
  !if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
  !  call MOM_error(FATAL,"Horizontal boundary diffusion is true, but no valid boundary layer scheme was found")
  !endif

  call get_param(param_file, mdl, "HBD_TRANSITION", CS%transition, &
                   "Controls if diffusion will transition from neutral to \n"//&
                   "horizontal.", default=.false.)
  if (CS%transition) then
    call get_param(param_file, mdl, "HBD_DELTA_H", CS%delta_h, &
                   "The maximum change in thickness between two adjacent cells used to \n"//&
                   "calculate the maximum isopycnal slope (delta_h/dx) that can be \n"    //&
                   "represented in the transition from neutral to horizontal.", units="m",  &
                   scale=US%m_to_Z, default=0.1)
  endif

  call get_param(param_file, mdl, "APPLY_LIMITER", CS%limiter, &
                   "If True, apply a flux limiter in the native grid.", default=.true.)
  call get_param(param_file, mdl, "APPLY_LIMITER_REMAP", CS%limiter_remap, &
                   "If True, apply a flux limiter in the remapped grid.", default=.false.)
  call get_param(param_file, mdl, "HBD_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in HBD code", &
                 default=.false.)
  call get_param(param_file, mdl, "HBD_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  !### Revisit this hard-coded answer_date.
  call initialize_remapping( CS%remap_CS, string, boundary_extrapolation = boundary_extrap ,&
       check_reconstruction=.false., check_remapping=.false., answer_date=20190101)
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
  call get_param(param_file, mdl, "HBD_DEBUG", CS%debug, &
                 "If true, write out verbose debugging data in the HBD module.", &
                 default=.false.)

  ! diagnostics
  CS%id_slp_x =  register_diag_field('ocean_model', 'bld_slp_x', CS%diag%axesCu1, Time, &
           'Zonal slope of neutral surface at the base of the BLD', 'nondim', conversion=US%Z_to_L)
  CS%id_slp_y =  register_diag_field('ocean_model', 'bld_slp_y', CS%diag%axesCv1, Time, &
           'Meridional slope of neutral surface at the base of the BLD', 'nondim', conversion=US%Z_to_L)
  CS%id_hbd_slp_x =  register_diag_field('ocean_model', 'hbd_slp_x', CS%diag%axesCu1, Time, &
           'Zonal slope of neutral surface used in the HBD module', 'nondim', conversion=US%Z_to_L)
  CS%id_hbd_slp_y =  register_diag_field('ocean_model', 'hbd_slp_y', CS%diag%axesCv1, Time, &
           'Meridional slope of neutral surface used in the HBD module', 'nondim', conversion=US%Z_to_L)

  if (CS%id_slp_x > 0)  allocate(CS%slp_x(SZIB_(G),SZJ_(G)), source=0.)
  if (CS%id_slp_y > 0)  allocate(CS%slp_y(SZI_(G),SZJB_(G)), source=0.)
  if (CS%id_hbd_slp_x > 0)  allocate(CS%hbd_slp_x(SZIB_(G),SZJ_(G)), source=0.)
  if (CS%id_hbd_slp_y > 0)  allocate(CS%hbd_slp_y(SZI_(G),SZJB_(G)), source=0.)

  id_clock_hbd = cpu_clock_id('(Ocean HBD)', grain=CLOCK_MODULE)

end function hor_bnd_diffusion_init

!> Driver routine for calculating horizontal diffusive fluxes near the top and bottom boundaries.
!! Diffusion is applied using only information from neighboring cells, as follows:
!! 1) remap tracer to a z* grid (HBD grid)
!! 2) calculate diffusive tracer fluxes (F) in the HBD grid using a layer by layer approach
!! 3) remap fluxes to the native grid
!! 4) update tracer by adding the divergence of F
subroutine hor_bnd_diffusion(G, GV, US, h, Coef_x, Coef_y, dt, VarMix, Reg, CS)
  type(ocean_grid_type),                intent(inout) :: G      !< Grid type
  type(verticalGrid_type),              intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),                intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                        intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),    intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),    intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                 intent(in)    :: dt     !< Tracer time step * I_numitts
                                                                !! (I_numitts in tracer_hordiff) [T ~> s]
  type(VarMix_CS),                      intent(in)    :: VarMix !< Variable mixing type
  type(tracer_registry_type),           pointer       :: Reg    !< Tracer registry
  type(hbd_CS),                         pointer       :: CS     !< Control structure for this module

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))           :: hbl         !< Boundary layer depth [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uFlx        !< Zonal flux of tracer [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vFlx        !< Meridional flux of tracer
                                                            !! [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZIB_(G),SZJ_(G))          :: uwork_2d    !< Layer summed u-flux transport
                                                            !! [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZI_(G),SZJB_(G))          :: vwork_2d    !< Layer summed v-flux transport
                                                            !! [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tendency    !< tendency array for diagnostic [conc T-1 ~> conc s-1]
  real, dimension(SZI_(G),SZJ_(G))           :: tendency_2d !< depth integrated content tendency for diagn
  type(tracer_type), pointer                 :: tracer => NULL() !< Pointer to the current tracer
  real, dimension(SZK_(GV)) :: tracer_1d                    !< 1d-array used to remap tracer change to native grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tracer_old  !< local copy of the initial tracer concentration,
                                                            !! only used to compute tendencies.
  real, dimension(SZI_(G),SZJ_(G))           :: tracer_int  !< integrated tracer before HBD is applied
                                                            !! [conc H L2 ~> conc m3 or conc kg]
  real, dimension(SZI_(G),SZJ_(G))           :: tracer_end  !< integrated tracer after HBD is applied.
                                                            !! [conc H L2 ~> conc m3 or conc kg]
  integer :: i, j, k, m   !< indices to loop over
  real    :: Idt          !< inverse of the time step [T-1 ~> s-1]
  real    :: tmp1, tmp2   !< temporary variables [conc H L2 ~> conc m3 or conc kg]

  call cpu_clock_begin(id_clock_hbd)
  Idt = 1./dt

  !GMM, in the idealized experiments, the following must be set
  hbl(:,:) = 100.
  hbl(4:6,:) = 300.

  !if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G, US, m_to_BLD_units=GV%m_to_H)
  !if (ASSOCIATED(CS%energetic_PBL_CSp)) call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US, &
  !                                                                 m_to_MLD_units=GV%m_to_H)

  call pass_var(hbl,G%Domain)

  ! build HBD grid
  call hbd_grid(SURFACE, G, GV, VarMix, hbl, h, CS)

  do m = 1,Reg%ntr
    ! current tracer
    tracer => Reg%tr(m)

    if (CS%debug) then
      call hchksum(tracer%t, "before HBD "//tracer%name,G%HI)
    endif

    ! for diagnostics
    if (tracer%id_hbdxy_conc > 0 .or. tracer%id_hbdxy_cont > 0 .or. tracer%id_hbdxy_cont_2d > 0 .or. CS%debug) then
      tendency(:,:,:) = 0.0
      tracer_old(:,:,:) = tracer%t(:,:,:)
    endif

    ! Diffusive fluxes in the i- and j-direction
    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    ! HBD layer by layer
    do j=G%jsc,G%jec
      do i=G%isc-1,G%iec
        if (G%mask2dCu(I,j)>0.) then
          if (hbl(I,j) == 0. .or. hbl(I+1,j) == 0.) then
            return
          else
            call fluxes_layer_method(SURFACE, G%ke, CS%hbd_u_kmax(I,j), &
               h(I,j,:), h(I+1,j,:), tracer%t(I,j,:), tracer%t(I+1,j,:),&
               CS%hbd_grd_u(I,j,:,1), CS%hbd_grd_u(I,j,:,2), hbl(I,j),  &
               hbl(I+1,j), Coef_x(I,j), uFlx(I,j,:), G%areaT(I,j),      &
               G%areaT(I+1,j), CS)
          endif
        endif
      enddo
    enddo
    do J=G%jsc-1,G%jec
      do i=G%isc,G%iec
        if (G%mask2dCv(i,J)>0.) then
          if (hbl(i,J) == 0. .or. hbl(i,J+1) == 0.) then
            return
          else
            call fluxes_layer_method(SURFACE, GV%ke, CS%hbd_v_kmax(i,J), &
              h(i,J,:), h(i,J+1,:), tracer%t(i,J,:), tracer%t(i,J+1,:),  &
              CS%hbd_grd_v(i,J,:,1) , CS%hbd_grd_v(i,J,:,2), hbl(i,J),   &
              hbl(i,J+1),  Coef_y(i,J), vFlx(i,J,:), G%areaT(i,J),       &
              G%areaT(i,J+1), CS)
          endif
        endif
      enddo
    enddo

    ! Update the tracer fluxes
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        tracer%t(i,j,k) = tracer%t(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))* &
                          G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff )

        if (tracer%id_hbdxy_conc > 0  .or. tracer%id_hbdxy_cont > 0 .or. tracer%id_hbdxy_cont_2d > 0 ) then
          tendency(i,j,k) = ((uFlx(I-1,j,k)-uFlx(I,j,k)) + (vFlx(i,J-1,k)-vFlx(i,J,k)))  * &
                            G%IareaT(i,j) * Idt
        endif
      endif
    enddo ; enddo ; enddo

    ! Do user controlled underflow of the tracer concentrations.
    if (tracer%conc_underflow > 0.0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        if (abs(tracer%t(i,j,k)) < tracer%conc_underflow) tracer%t(i,j,k) = 0.0
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(tracer%t, "after HBD "//tracer%name,G%HI)
      tracer_int(:,:) = 0.0; tracer_end(:,:) = 0.0
      ! tracer (native grid) before and after HBD
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        do k=1,GV%ke
          tracer_int(i,j) = tracer_int(i,j) + tracer_old(i,j,k) * &
                    (h(i,j,k)*(G%mask2dT(i,j)*G%areaT(i,j)))
          tracer_end(i,j) = tracer_end(i,j) + tracer%t(i,j,k) * &
                      (h(i,j,k)*(G%mask2dT(i,j)*G%areaT(i,j)))
        enddo
      enddo; enddo

      tmp1 = SUM(tracer_int)
      tmp2 = SUM(tracer_end)
      call sum_across_PEs(tmp1)
      call sum_across_PEs(tmp2)
      if (is_root_pe()) write(*,*)'Total '//tracer%name//' before/after HBD:', tmp1, tmp2
    endif

    ! Post the tracer diagnostics
    if (tracer%id_hbd_dfx>0)      call post_data(tracer%id_hbd_dfx, uFlx(:,:,:)*Idt, CS%diag)
    if (tracer%id_hbd_dfy>0)      call post_data(tracer%id_hbd_dfy, vFlx(:,:,:)*Idt, CS%diag)
    if (tracer%id_hbd_dfx_2d>0) then
      uwork_2d(:,:) = 0.
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        uwork_2d(I,j) = uwork_2d(I,j) + (uFlx(I,j,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbd_dfx_2d, uwork_2d, CS%diag)
    endif

    if (tracer%id_hbd_dfy_2d>0) then
      vwork_2d(:,:) = 0.
      do k=1,GV%ke ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        vwork_2d(i,J) = vwork_2d(i,J) + (vFlx(i,J,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbd_dfy_2d, vwork_2d, CS%diag)
    endif

    ! post tendency of tracer content
    if (tracer%id_hbdxy_cont > 0) then
      call post_data(tracer%id_hbdxy_cont, tendency, CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_hbdxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        do k=1,GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_hbdxy_cont_2d, tendency_2d, CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array and its units.
    if (tracer%id_hbdxy_conc > 0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + CS%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_hbdxy_conc, tendency, CS%diag)
    endif

  enddo
  ! diagnostics that are not tracer dependent
  if (CS%id_hbd_slp_x > 0) call post_data(CS%id_hbd_slp_x, CS%hbd_slp_x, CS%diag)
  if (CS%id_hbd_slp_y > 0) call post_data(CS%id_hbd_slp_y, CS%hbd_slp_y, CS%diag)
  if (CS%id_slp_x > 0) call post_data(CS%id_slp_x, CS%slp_x, CS%diag)
  if (CS%id_slp_y > 0) call post_data(CS%id_slp_y, CS%slp_y, CS%diag)

  call cpu_clock_end(id_clock_hbd)

end subroutine hor_bnd_diffusion

!> Calculate the harmonic mean of two quantities
!! See \ref section_harmonic_mean.
real function harmonic_mean(h1,h2)
  real :: h1 !< Scalar quantity
  real :: h2 !< Scalar quantity
  if (h1 + h2 == 0.) then
    harmonic_mean = 0.
  else
    harmonic_mean = 2.*(h1*h2)/(h1+h2)
  endif
end function harmonic_mean

!> Returns the location of the minimum value in a 1D array
!! between indices s and e.
integer function  find_minimum(x, s, e)
  integer, intent(in) :: s              !< start index
  integer, intent(in) :: e              !< end index
  real, dimension(e), intent(in) :: x   !< 1D array to be checked

  ! local variables
  real :: minimum
  integer :: location
  integer :: i

  minimum  = x(s)   ! assume the first is the min
  location = s      ! record its position
  do i = s+1, e     ! start with next elements
    if (x(i) < minimum) then !   if x(i) less than the min?
      minimum  = x(i)   !      Yes, a new minimum found
      location = i                !      record its position
    end if
  enddo
  find_minimum = location          ! return the position
end function  find_minimum

!> Swaps the values of its two formal arguments.
subroutine swap(a, b)
  real, intent(inout) :: a  !< First value to be swaped
  real, intent(inout) :: b  !< Second value to be swaped

  ! local variables
  real :: tmp

  tmp = a
  a = b
  b = tmp
end subroutine swap

!> Receives a 1D array x and sorts it into ascending order.
subroutine sort(x, n)
  real, dimension(n),  intent(inout) :: x        !< 1D array to be sorted
  integer,             intent(in   ) :: n        !< # of pts in the array

  ! local variables
  integer :: i, location

  do i = 1, n-1
    location = find_minimum(x, i, n) ! find min from this to last
    call swap(x(i), x(location))     ! swap this and the minimum
  enddo
end subroutine sort

!> Returns the unique values in a 1D array.
subroutine unique(val, n, val_unique, val_max)
  integer,                         intent(in   ) :: n          !< # of pts in the array.
  real, dimension(n),              intent(in   ) :: val        !< 1D array to be checked.
  real, dimension(:), allocatable, intent(inout) :: val_unique !< Returned 1D array with unique values.
  real,                  optional, intent(in   ) :: val_max    !< sets the maximum value in val_unique to
                                                               !! this value.
  ! local variables
  real, dimension(n) :: tmp
  integer :: i, j, ii
  real :: min_val, max_val
  logical :: limit

  limit = .false.
  if (present(val_max)) then
    limit = .true.
    if (val_max > MAXVAL(val)) then
      if (is_root_pe()) write(*,*)'val_max, MAXVAL(val)',val_max, MAXVAL(val)
      call MOM_error(FATAL,"Houston, we've had a problem in unique (val_max cannot be > MAXVAL(val))")
    endif
  endif

  tmp(:) = 0.
  min_val = MINVAL(val)-1
  max_val = MAXVAL(val)
  i = 0
  do while (min_val<max_val)
      i = i+1
      min_val = MINVAL(val, mask=val>min_val)
      tmp(i) = min_val
  enddo
  ii = i
  if (limit) then
    do j=1,ii
      if (tmp(j) <= val_max) i = j
    enddo
  endif
  allocate(val_unique(i), source=tmp(1:i))
end subroutine unique

!> Given layer thicknesses (and corresponding interfaces) and BLDs in two adjacent columns,
!! return a set of 1-d layer thicknesses whose interfaces cover all interfaces in the left
!! and right columns plus the two BLDs. This can be used to accurately remap tracer tendencies
!! in both columns.
subroutine merge_interfaces_hbd(nk_L, h_L, nk_R, h_R, hbl_L, hbl_R, H_subroundoff, h)
  integer,                         intent(in   ) :: nk_L   !< Number of layers left column            [nondim]
  real, dimension(nk_L),             intent(in   ) :: h_L    !< Layer thicknesses in the left column    [H ~> m or kg m-2]
  integer,                         intent(in   ) :: nk_R   !< Number of layers right column           [nondim]
  real, dimension(nk_R),             intent(in   ) :: h_R    !< Layer thicknesses in the right column   [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_L  !< Thickness of the boundary layer in the left column
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_R  !< Thickness of the boundary layer in the right column
  type(hbd_CS),         pointer       :: CS       !< Horizontal diffusion control structure
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: H_subroundoff !< GV%H_subroundoff                 [H ~> m or kg m-2]
  real, dimension(:), allocatable, intent(inout) :: h     !< Combined thicknesses                     [H ~> m or kg m-2]

  ! Local variables
  integer                         :: n           !< Number of layers in eta_all
  real, dimension(nk_L+1)         :: eta_L       !< Interfaces in the left coloumn
  real, dimension(nk_R+1)         :: eta_R       !< Interfaces in the right coloumn
  real, dimension(:), allocatable :: eta_all     !< Combined interfaces in the left/right columns + hbl_L and hbl_R
  real, dimension(:), allocatable :: eta_unique  !< Combined interfaces (eta_L, eta_R), possibly hbl_L and hbl_R
  real                            :: min_depth   !< Minimum depth
  real                            :: max_depth   !< Maximum depth
  real                            :: max_bld     !< Deepest BLD
  integer :: k, kk, nk1                          !< loop indices (k and kk) and array size (nk1)
  logical :: debug = .false.

  n = (nk_L + nk_R)+4 ! GMM, check this
  allocate(eta_all(n))
  ! compute and merge interfaces
  eta_L(:) = 0.0; eta_R(:) = 0.0; eta_all(:) = 0.0

  ! left
  kk = 0
  do k=2,nk_L+1
    eta_L(k) = eta_L(k-1) + h_L(k-1)
    if (debug .and. is_root_pe()) write(*,*)'L -> k, eta_L(k)',k,eta_L(k)
    kk = kk + 1
    eta_all(kk)   = eta_L(k)
  enddo

  ! right
  do k=2,nk_R+1
    eta_R(k) = eta_R(k-1) + h_R(k-1)
    if (debug .and. is_root_pe()) write(*,*)'R -> k, eta_R(k)',k,eta_R(k)
    kk = kk + 1
    eta_all(kk) = eta_R(k)
  enddo

  ! add hbl_L and hbl_R into eta_all
  eta_all(kk+1) = hbl_L
  eta_all(kk+2) = hbl_R

  ! find maximum depth
  min_depth = MIN(MAXVAL(eta_L), MAXVAL(eta_R))
  max_bld = MAX(hbl_L, hbl_R)
  max_depth = MIN(min_depth, max_bld)

  ! sort eta_all
  call sort(eta_all, n)
  ! remove duplicates from eta_all and sets maximum depth
  call unique(eta_all, n, eta_unique, max_depth)

  nk1 = SIZE(eta_unique)
  allocate(h(nk1-1))
  do k=1,nk1-1
    if (debug .and. is_root_pe()) write(*,*)'merge -> k, eta(k)',k,eta_unique(k)
    h(k) = (eta_unique(k+1) - eta_unique(k)) + H_subroundoff
  enddo
  if (debug .and. is_root_pe()) write(*,*)'merge -> k+1, eta(k+1)',k,eta_unique(k)

end subroutine merge_interfaces_hbd

!> Given layer thicknesses (and corresponding interfaces) and BLDs in two adjacent columns,
!! return a set of 1-d layer thicknesses whose interfaces cover all interfaces in the left
!! and right columns plus the two BLDs. This can be used to accurately remap tracer tendencies
!! in both columns.
subroutine merge_interfaces(nk, h_L, h_R, hbl_L, hbl_R, H_subroundoff, h)
  integer,                         intent(in   ) :: nk     !< Number of layers                        [nondim]
  real, dimension(nk),             intent(in   ) :: h_L    !< Layer thicknesses in the left column    [H ~> m or kg m-2]
  real, dimension(nk),             intent(in   ) :: h_R    !< Layer thicknesses in the right column   [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_L  !< Thickness of the boundary layer in the left column
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: hbl_R  !< Thickness of the boundary layer in the right column
  type(hbd_CS),         pointer       :: CS       !< Horizontal diffusion control structure
                                                           !!                                         [H ~> m or kg m-2]
  real,                            intent(in   ) :: H_subroundoff !< GV%H_subroundoff                 [H ~> m or kg m-2]
  real, dimension(:), allocatable, intent(inout) :: h     !< Combined thicknesses                     [H ~> m or kg m-2]

  ! Local variables
  integer                         :: n           !< Number of layers in eta_all
  real, dimension(nk+1)           :: eta_L, eta_R!< Interfaces in the left and right coloumns
  real, dimension(:), allocatable :: eta_all     !< Combined interfaces in the left/right columns + hbl_L and hbl_R
  real, dimension(:), allocatable :: eta_unique  !< Combined interfaces (eta_L, eta_R), possibly hbl_L and hbl_R
  real                            :: min_depth   !< Minimum depth
  real                            :: max_depth   !< Maximum depth
  real                            :: max_bld     !< Deepest BLD
  integer :: k, kk, nk1                          !< loop indices (k and kk) and array size (nk1)

  n = (2*nk)+3
  allocate(eta_all(n))
  ! compute and merge interfaces
  eta_L(:) = 0.0; eta_R(:) = 0.0; eta_all(:) = 0.0
  kk = 0
  do k=2,nk+1
    eta_L(k) = eta_L(k-1) + h_L(k-1)
    eta_R(k) = eta_R(k-1) + h_R(k-1)
    kk = kk + 2
    eta_all(kk)   = eta_L(k)
    eta_all(kk+1) = eta_R(k)
  enddo

  ! add hbl_L and hbl_R into eta_all
  eta_all(kk+2) = hbl_L
  eta_all(kk+3) = hbl_R

  ! find maximum depth
  min_depth = MIN(MAXVAL(eta_L), MAXVAL(eta_R))
  max_bld = MAX(hbl_L, hbl_R)
  max_depth = MIN(min_depth, max_bld)

  ! sort eta_all
  call sort(eta_all, n)
  ! remove duplicates from eta_all and sets maximum depth
  call unique(eta_all, n, eta_unique, max_depth)

  nk1 = SIZE(eta_unique)
  allocate(h(nk1-1))
  do k=1,nk1-1
    h(k) = (eta_unique(k+1) - eta_unique(k)) + H_subroundoff
  enddo
end subroutine merge_interfaces

!> Calculates the maximum flux that can leave a cell and uses that to apply a
!! limiter to F_layer.
subroutine flux_limiter(F_layer, area_L, area_R, phi_L, phi_R, h_L, h_R)
  real, intent(inout) :: F_layer !< Tracer flux to be checked [H L2 conc ~> m3 conc]
  real, intent(in) :: area_L     !< Area of left cell [L2 ~> m2]
  real, intent(in) :: area_R     !< Area of right cell [L2 ~> m2]
  real, intent(in) :: h_L        !< Thickness of left cell [H ~> m or kg m-2]
  real, intent(in) :: h_R        !< Thickness of right cell [H ~> m or kg m-2]
  real, intent(in) :: phi_L      !< Tracer concentration in the left cell [conc]
  real, intent(in) :: phi_R      !< Tracer concentration in the right cell [conc]

  ! local variables
  real :: F_max !< maximum flux allowed
  ! limit the flux to 0.2 of the tracer *gradient*
  ! Why 0.2?
  !  t=0         t=inf
  !   0           .2
  ! 0 1 0       .2.2.2
  !   0           .2
  !
  F_max = -0.2 * ((area_R*(phi_R*h_R))-(area_L*(phi_L*h_L)))

  if ( SIGN(1.,F_layer) == SIGN(1., F_max)) then
    ! Apply flux limiter calculated above
    if (F_max >= 0.) then
      F_layer = MIN(F_layer,F_max)
    else
      F_layer = MAX(F_layer,F_max)
    endif
  else
    F_layer = 0.0
  endif
end subroutine flux_limiter

!> Find the k-index range corresponding to the layers that are within the boundary-layer region
subroutine boundary_k_range(boundary, nk, h, hbl, k_top, zeta_top, k_bot, zeta_bot)
  integer,             intent(in   ) :: boundary !< SURFACE or BOTTOM                       [nondim]
  integer,             intent(in   ) :: nk       !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the column         [H ~> m or kg m-2]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [H ~> m or kg m-2]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables
  real :: htot ! Summed thickness [H ~> m or kg m-2]
  integer :: k
  ! Surface boundary layer
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    k_bot = 1
    zeta_bot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_bot = nk
      zeta_bot = 1.
      return
    endif
    do k=1,nk
      htot = htot + h(k)
      if ( htot >= hbl) then
        k_bot = k
        zeta_bot = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  ! Bottom boundary layer
  elseif ( boundary == BOTTOM ) then
    k_top = nk
    zeta_top = 1.
    k_bot = nk
    zeta_bot = 0.
    htot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_top = 1
      zeta_top = 1.
      return
    endif
    do k=nk,1,-1
      htot = htot + h(k)
      if (htot >= hbl) then
        k_top = k
        zeta_top = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  else
    call MOM_error(FATAL,"Houston, we've had a problem in boundary_k_range")
  endif

end subroutine boundary_k_range

!> Build the HBD grid where tracers will be rammaped to.
subroutine build_hbd_grid(nk, k_max, slope_bld, dx, hbd1, hbd2, CS)
  integer,             intent(in   ) :: nk     !< number of layers in hbd1 and hbd2 [nondim]
  integer,             intent(in   ) :: k_max  !< vertical index to start neutral-horizontal transition [nondim]
  real,                intent(in   ) :: slope_bld !< Max. isopycnal slope              [nondim]
  real,                intent(in   ) :: dx    !< dx or dy                              [m]
  !real, dimension(k_max), intent(in   ) :: h_vel !< Layer thickness at vel pts               [H ~> m or kg m-2]
  real, dimension(nk), intent(in   ) :: hbd1  !< Layer thickness in the first HBD grid     [H ~> m or kg m-2]
  real, dimension(nk), intent(inout) :: hbd2  !< Layer thickness in the second HBD grid    [H ~> m or kg m-2]
  type(hbd_CS),        pointer       :: CS    !< Horizontal diffusion control structure

  ! local variables
  integer :: k
  real, dimension(k_max+1) :: e1, e2, e_vel
  real    :: dh
  real    :: abs_slope
  real    :: slope, a

  ! slope_bld is always > 0
  abs_slope = abs(slope_bld)

  ! build interfaces for e_vel
  e1(:) = 0.0
  e2(:) = 0.0
  e_vel(:) = 0.0
  do k=2,k_max+1
    e1(k) = e1(k-1) + hbd1(k-1)
  enddo

  ! isopycnal slope decreases linearly
  do k=k_max+1,1,-1
    slope = (abs_slope/e1(k_max+1)) * e1(k)
    dh = (atan(slope)*dx)
    e2(k) = e1(k) - (1.0 * dh)
  enddo

  ! build thicknesses hbd2
  do k=1,k_max
    hbd2(k) = (e2(k+1) - e2(k)) + CS%H_subroundoff
  enddo

  if (CS%debug) then
    if (is_root_pe()) then
      do k=1,nk
        write(*,*)'k, hbd1(k), hbd2(k), h_mean(k), h_vel(k):', k, hbd1(k), hbd2(k), 0.5*(hbd1(k)+hbd2(k)) !, h_vel(k)
      enddo
      write(*,*)'SUM(hbd1(:)), SUM(hbd2(:))',SUM(hbd1(:)), SUM(hbd2(:))
      write(*,*)'###### abs_slope (applied slope) #####', abs_slope
      if (abs(abs_slope)> 0.) then
        do k=k_max+1,1,-1
          slope = TAN((e1(k)-e2(k))/dx)
          write(*,*)'k, slope', k, slope
        enddo
      endif
    endif
  endif

end subroutine build_hbd_grid

!> Build the HBD grid where tracers will be rammaped to.
subroutine hbd_grid(boundary, G, GV, VarMix, hbl, h, CS)
  integer,                 intent(in   ) :: boundary !< Which boundary layer SURFACE or BOTTOM           [nondim]
  type(ocean_grid_type),   intent(inout) :: G    !< Grid type
  type(verticalGrid_type), intent(in)    :: GV   !< ocean vertical grid structure
  type(VarMix_CS),         intent(in)    :: VarMix !< Variable mixing type
  real, dimension(SZI_(G),SZJ_(G))       :: hbl      !< Boundary layer depth                   [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          intent(in)     :: h        !< Layer thickness in the native grid     [H ~> m or kg m-2]
  type(hbd_CS),          pointer         :: CS       !< Horizontal diffusion control structure

  ! Local variables
  real, allocatable :: dz_top(:)     !< temporary HBD grid given by merge_interfaces        [H ~> m or kg m-2]
  real :: slp_x !< Zonal isopycnal slope at the base of the boundary layer                            [nondim]
  real :: slp_y !< Meridional isopycnal slope at the base of the boundary layer                       [nondim]
  real :: hbd_slp_x
  real :: hbd_slp_y
  real :: h_vel(SZK_(GV)) !< Thicknesses at u- and v-points in the native grid
                          !! The harmonic mean is used to avoid zero values      [H ~> m or kg m-2]
  real    :: max_hbl      !< maximum boundary layer depth between left and right columns [H ~> m or kg m-2]
  integer :: k_bot_min, k_bot_max !< k-indices min and max, respectively.
  integer :: k_bot_L, k_bot_R !< k-indices for left and right columns, respectively.
  integer :: k_bot_diff   !< different between left and right k-indices.
  integer :: k_top, k_bot !< indices used to store position of maximum isopycnal slope.
  real    :: zeta_top     !< distance from the top of a layer to the boundary
                          !! layer depth in the native grid                                [nondim]
  real    :: zeta_bot     !< distance from the bottom of a layer to the boundary
                                     !! layer depth in the native grid                                [nondim]
  real    :: max_slope    !< Maximum isopycnal slope that can be represented by the HBD grid          [nondim]
  integer :: nk, i, j, k  !< number of layers in the HBD grid, and integers used in do-loops

  ! reset arrays
  CS%hbd_grd_u(:,:,:,:) = CS%H_subroundoff
  CS%hbd_grd_v(:,:,:,:) = CS%H_subroundoff
  CS%hbd_u_kmax(:,:)  = 0
  CS%hbd_v_kmax(:,:)  = 0

  ! neutral surface slopes at the base of the boundary layer in i- and j-direction
  slp_x = 0.
  slp_y = 0.

  do j=G%jsc,G%jec
    do I=G%isc-1,G%iec
      if (G%mask2dCu(I,j)>0.) then
        call merge_interfaces(GV%ke, h(I,j,:), h(I+1,j,:), hbl(I,j), hbl(I+1,j), &
                              CS%H_subroundoff, dz_top)
        nk = SIZE(dz_top)
        if (nk > CS%hbd_nk) then
          write(*,*)'nk, CS%hbd_nk', nk, CS%hbd_nk
          call MOM_error(FATAL,"Houston, we've had a problem in hbd_grid, u-points (nk cannot be > CS%hbd_nk)")
        endif

        CS%hbd_u_kmax(I,j) = nk

        ! set the HBD grid to dz_top
        do k=1,nk
          CS%hbd_grd_u(I,j,k,1) = dz_top(k)
          CS%hbd_grd_u(I,j,k,2) = dz_top(k)
        enddo

        ! thicknesses at u-points on the native grid
        h_vel(:) = 0.0
        do k=1,GV%ke
          h_vel(k) = harmonic_mean(h(I,j,k), h(I+1,j,k))
        enddo

        ! max hbl
        max_hbl = max(hbl(I,j), hbl(I+1,j))

        ! get zonal slope
        call boundary_k_range(boundary, GV%ke, h_vel(:), max_hbl, k_top, zeta_top, k_bot, zeta_bot)
        ! slp_x is at the interface, hence the + 1
        slp_x = VarMix%slope_x(I,j,k_bot+1)
        if (CS%id_slp_x > 0) CS%slp_x(I,j) = slp_x

        ! Calculate vertical indices containing the boundary layer in dz_top
        call boundary_k_range(boundary, nk, dz_top, hbl(I,j  ), k_top, zeta_top, k_bot_L, zeta_bot)
        call boundary_k_range(boundary, nk, dz_top, hbl(I+1,j), k_top, zeta_top, k_bot_R, zeta_bot)
        k_bot_min = MIN(k_bot_L, k_bot_R)
        k_bot_max = MAX(k_bot_L, k_bot_R)
        k_bot_diff = (k_bot_max - k_bot_min)

        if ((k_bot_diff > 1) .and. (CS%transition)) then
          ! compute the maximum slope allowed
          max_slope = CS%delta_h/G%dxCu(I,j)
          !max_slope = CS%delta_h

          ! slope is > 0, right column is shallower
          if (slp_x >= 0.) then
            ! limit zonal slope
            if (abs(slp_x) > max_slope) slp_x = max_slope
            call build_hbd_grid(nk, k_bot_max, slp_x, G%dxCu(I,j), &
                                CS%hbd_grd_u(I,j,:,1), &
                                CS%hbd_grd_u(I,j,:,2), &
                                CS)

          ! slope is < 0, left column is shallower
          else
            ! limit zonal slope
            if (abs(slp_x) > max_slope) slp_x = -max_slope
            call build_hbd_grid(nk, k_bot_max, slp_x, G%dxCu(I,j), &
                                CS%hbd_grd_u(I,j,:,2), &
                                CS%hbd_grd_u(I,j,:,1), &
                                CS)
          endif
        else
          ! slope used is zero by definition
          slp_x = 0.0
        endif
        if (CS%id_hbd_slp_x > 0) CS%hbd_slp_x(I,j) = slp_x
        deallocate(dz_top)
      endif
    enddo
  enddo

  do J=G%jsc-1,G%jec
    do i=G%isc,G%iec
      if (G%mask2dCv(i,J)>0.) then
        call merge_interfaces(GV%ke, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1), &
                              CS%H_subroundoff, dz_top)

        nk = SIZE(dz_top)
        if (nk > CS%hbd_nk) then
          write(*,*)'nk, CS%hbd_nk', nk, CS%hbd_nk
          call MOM_error(FATAL,"Houston, we've had a problem in hbd_grid, v-points (nk cannot be > CS%hbd_nk)")
        endif

        CS%hbd_v_kmax(i,J) = nk

        ! set the HBD grid to dz_top
        do k=1,nk
          CS%hbd_grd_v(i,J,k,1) = dz_top(k)
          CS%hbd_grd_v(i,J,k,2) = dz_top(k)
        enddo

        ! thicknesses at v-points on the native grid
        h_vel(:) = 0.0
        do k=1,GV%ke
          h_vel(k) = harmonic_mean(h(i,J,k), h(i,J+1,k))
        enddo

        ! max hbl
        max_hbl = max(hbl(i,J), hbl(i,J+1))

        call boundary_k_range(boundary, GV%ke, h_vel(:), max_hbl, k_top, zeta_top, k_bot, zeta_bot)
        ! slp_y is at the interface, hence the + 1
        slp_y = VarMix%slope_y(i,J,k_bot+1)
        if (CS%id_slp_y > 0) CS%slp_y(i,J) = slp_y

        ! Calculate vertical indices containing the boundary layer in dz_top
        call boundary_k_range(boundary, nk, dz_top, hbl(i,J  ), k_top, zeta_top, k_bot_L, zeta_bot)
        call boundary_k_range(boundary, nk, dz_top, hbl(i,J+1), k_top, zeta_top, k_bot_R, zeta_bot)
        k_bot_min = MIN(k_bot_L, k_bot_R)
        k_bot_max = MAX(k_bot_L, k_bot_R)
        k_bot_diff = (k_bot_max - k_bot_min)

        if ((k_bot_diff > 1) .and. (CS%transition)) then
          ! compute the maximum slope allowed
          max_slope = CS%delta_h/G%dyCv(i,J)
          !max_slope = CS%delta_h

          ! slope is > 0, right column is shallower
          if (slp_y >= 0.) then
            ! limit meridional slope
            if (abs(slp_y) > max_slope) slp_y = max_slope
            call build_hbd_grid(nk, k_bot_max, slp_y, G%dyCv(i,J), &
                                CS%hbd_grd_v(i,J,:,1), &
                                CS%hbd_grd_v(i,J,:,2), &
                                CS)
          ! slope is < 0, left column is shallower
          else
            ! limit meridional slope
            if (abs(slp_y) > max_slope) slp_y = -max_slope
            call build_hbd_grid(nk, k_bot_max, slp_y, G%dyCv(i,J), &
                                CS%hbd_grd_v(i,J,:,2), &
                                CS%hbd_grd_v(i,J,:,1), &
                                CS)
          endif
        else
          ! slope used is zero by definition
          slp_y = 0.0
        endif
        if (CS%id_hbd_slp_y > 0) CS%hbd_slp_y(i,J) = slp_y
        deallocate(dz_top)
      endif
    enddo
  enddo

end subroutine hbd_grid


!> Calculate the horizontal boundary diffusive fluxes using the layer by layer method.
!! See \ref section_method
subroutine fluxes_layer_method(boundary, nk, hbd_nk, h_L, h_R, phi_L, phi_R, hbd_L, hbd_R, &
                               hbl_L, hbl_R, khtr_u, F_layer, area_L, area_R, CS)

  integer,              intent(in   ) :: boundary !< Which boundary layer SURFACE or BOTTOM           [nondim]
  integer,              intent(in   ) :: nk       !< Number of layers in the native grid              [nondim]
  integer,              intent(in   ) :: hbd_nk   !< Number of layers in the HBD grid                 [nondim]
  real, dimension(nk),  intent(in   ) :: h_L      !< Thicknesses in the native grid (left)  [H ~> m or kg m-2]
  real, dimension(nk),  intent(in   ) :: h_R      !< Thicknesses in the native grid (right) [H ~> m or kg m-2]
  real, dimension(nk),  intent(in   ) :: phi_L    !< Tracer values in the native grid (left)            [conc]
  real, dimension(nk),  intent(in   ) :: phi_R    !< Tracer values in the native grid (right)           [conc]
  real, dimension(hbd_nk),  intent(in   ) :: hbd_L    !< Thicknesses in the HBD grid (left)  [H ~> m or kg m-2]
  real, dimension(hbd_nk),  intent(in   ) :: hbd_R    !< Thicknesses in the HBD grid (right) [H ~> m or kg m-2]
  real,                 intent(in   ) :: hbl_L    !< Thickness of the boundary boundary
                                                  !! layer (left)                           [H ~> m or kg m-2]
  real,                 intent(in   ) :: hbl_R    !< Thickness of the boundary boundary
                                                  !! layer (right)
  real,                 intent(in   ) :: khtr_u   !< Horizontal diffusivities times the time step
                                                  !! at a velocity point                            [L2 ~> m2]
  real, dimension(nk),  intent(  out) :: F_layer  !< Layerwise diffusive flux at U- or V-point
                                                  !! in the native grid                 [H L2 conc ~> m3 conc]
  real,                 intent(in   ) :: area_L   !< Area of the horizontal grid (left)             [L2 ~> m2]
  real,                 intent(in   ) :: area_R   !< Area of the horizontal grid (right)            [L2 ~> m2]
  type(hbd_CS),         pointer       :: CS       !< Horizontal diffusion control structure

  ! Local variables
  real, allocatable :: phi_L_z(:)    !< Tracer values in the ztop grid (left)                           [conc]
  real, allocatable :: phi_R_z(:)    !< Tracer values in the ztop grid (right)                          [conc]
  real, allocatable :: F_layer_z(:)  !< Diffusive flux at U/V-point in the ztop grid    [H L2 conc ~> m3 conc]
  real, allocatable :: hbd_h_vel(:)  !< Thicknesses at u- and v-points in the HBD grid
                                     !! The harmonic mean is used to avoid zero values      [H ~> m or kg m-2]
  real              :: h_vel(nk)     !< Thicknesses at u- and v-points in the native grid
                                     !! The harmonic mean is used to avoid zero values      [H ~> m or kg m-2]
  real    :: khtr_avg                !< Thickness-weighted diffusivity at the velocity-point [L2 T-1 ~> m2 s-1]
                                     !! This is just to remind developers that khtr_avg should be
                                     !! computed once khtr is 3D.
  real    :: htot                    !< Total column thickness                              [H ~> m or kg m-2]
  integer :: k, k_bot_min, k_top_max !< k-indices, min and max for bottom and top, respectively
  integer :: k_bot_max, k_top_min    !< k-indices, max and min for bottom and top, respectively
  integer :: k_bot_diff, k_top_diff  !< different between left and right k-indices for bottom and top, respectively
  integer :: k_top_L, k_bot_L        !< k-indices left native grid
  integer :: k_top_R, k_bot_R        !< k-indices right native grid
  real    :: zeta_top_L, zeta_top_R  !< distance from the top of a layer to the boundary
                                     !! layer depth in the native grid                                [nondim]
  real    :: zeta_bot_L, zeta_bot_R  !< distance from the bottom of a layer to the boundary
                                     !! layer depth in the native grid                                [nondim]
  real    :: wgt                     !< weight to be used in the linear transition to the interior    [nondim]
  real    :: a                       !< coefficient to be used in the linear transition to the interior [nondim]
  real    :: tmp1, tmp2              !< dummy variables                                     [H ~> m or kg m-2]
  real    :: htot_max                !< depth below which no fluxes should be applied       [H ~> m or kg m-2]

  F_layer(:) = 0.0

  ! allocate arrays
  allocate(phi_L_z(hbd_nk), source=0.0)
  allocate(phi_R_z(hbd_nk), source=0.0)
  allocate(F_layer_z(hbd_nk), source=0.0)
  allocate(hbd_h_vel(hbd_nk), source=0.0)

  ! remap tracer to the hbd grid
  call remapping_core_h(CS%remap_cs, nk, h_L(:), phi_L(:), hbd_nk, hbd_L(:), &
                        phi_L_z(:), CS%H_subroundoff, CS%H_subroundoff)
  call remapping_core_h(CS%remap_cs, nk, h_R(:), phi_R(:), hbd_nk, hbd_R(:), &
                        phi_R_z(:), CS%H_subroundoff, CS%H_subroundoff)

  ! Calculate vertical indices containing the boundary layer in dz_top
  call boundary_k_range(boundary, hbd_nk, hbd_L, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, hbd_nk, hbd_R, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)

  if (boundary == SURFACE) then
    k_bot_min = MIN(k_bot_L, k_bot_R)
    k_bot_max = MAX(k_bot_L, k_bot_R)
    k_bot_diff = (k_bot_max - k_bot_min)

    do k = k_bot_max,1,-1
      !! GMM, this is the newest method, more general. It uses the actual thicknesses for
      !! the left and right cells
      F_layer_z(k) = -khtr_u * ((phi_R_z(k) * hbd_R(k)) - (phi_L_z(k) * hbd_L(k)))
      if (CS%limiter_remap) call flux_limiter(F_layer_z(k), area_L, area_R, phi_L_z(k), &
                                             phi_R_z(k), hbd_L(k), hbd_R(k))
    enddo

    do k = 1, hbd_nk
      hbd_h_vel(k) = harmonic_mean(hbd_L(k), hbd_R(k))
    enddo
  endif

  ! TODO, boundary == BOTTOM
  !  if (boundary == BOTTOM) then
  !  endif

  ! thicknesses at velocity points
  do k = 1,nk
    h_vel(k) = harmonic_mean(h_L(k), h_R(k))
  enddo

  ! remap flux to h_vel (native grid)
  call reintegrate_column(hbd_nk, hbd_h_vel(:), F_layer_z(:), nk, h_vel(:), 0.0, F_layer(:))

  htot_max = SUM(hbd_h_vel)

  tmp1 = 0.0; tmp2 = 0.0
  do k = 1,nk
    ! apply flux_limiter
    if (CS%limiter .and. F_layer(k) /= 0.) then
       call flux_limiter(F_layer(k), area_L, area_R, phi_L(k), phi_R(k), h_L(k), h_R(k))
    endif

    ! if tracer point is below htot_max, set flux to zero
    if (MAX(tmp1+(h_L(k)*0.5), tmp2+(h_R(k)*0.5)) > htot_max) then
      F_layer(k) = 0.
    endif

    tmp1 = tmp1 + h_L(k)
    tmp2 = tmp2 + h_R(k)
  enddo

  ! deallocated arrays
  deallocate(phi_L_z)
  deallocate(phi_R_z)
  deallocate(F_layer_z)
  deallocate(hbd_h_vel)

end subroutine fluxes_layer_method

!> Unit tests for near-boundary horizontal mixing
logical function near_boundary_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables
  integer, parameter    :: nk = 2               ! Number of layers
  real, dimension(nk+1) :: eta1                 ! Updated interfaces with one extra value           [m]
  real, dimension(:), allocatable :: h1         ! Upates layer thicknesses                          [m]
  real, dimension(nk)   :: phi_L, phi_R         ! Tracer values (left and right column)             [conc]
  real, dimension(nk)   :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real                  :: khtr_u               ! Horizontal diffusivities at U-point               [m2 s-1]
  real                  :: hbl_L, hbl_R         ! Depth of the boundary layer (left and right)      [m]
  real, dimension(nk)   :: F_layer              ! Diffusive flux within each layer at U-point       [conc m3 s-1]
  character(len=120)    :: test_name            ! Title of the unit test
  integer               :: k_top                ! Index of cell containing top of boundary
  real                  :: zeta_top             ! Nondimension position                             [nondim]
  integer               :: k_bot                ! Index of cell containing bottom of boundary
  real                  :: zeta_bot             ! Nondimension position                             [nondim]
  type(hbd_CS), pointer :: CS

  allocate(CS)
  ! fill required fields in CS
  call initialize_remapping( CS%remap_CS, 'PLM', boundary_extrapolation = .true. ,&
       check_reconstruction = .true., check_remapping = .true.)
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
  CS%H_subroundoff = 1.0E-20
  CS%debug=.false.
  CS%limiter=.false.
  CS%limiter_remap=.false.

  near_boundary_unit_tests = .false.
  write(stdout,*) '==== MOM_hor_bnd_diffusion ======================='

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 1., 2, 0., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 1., 2, 0., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Surface boundary is deeper than column thickness'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 21.0, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 0., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 0., test_name, verbose)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed boundary_k_range'

  ! unit tests for sorting array and finding unique values
  test_name = 'Sorting array'
  eta1 = (/1., 0., 0.1/)
  call sort(eta1, nk+1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, eta1, (/0., 0.1, 1./) )

  test_name = 'Unique values'
  call unique((/0., 1., 1., 2./), nk+2, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0., 1., 2./) )
  deallocate(h1)

  test_name = 'Unique values with maximum depth'
  call unique((/0., 1., 1., 2., 3./), nk+3, h1, 2.)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0., 1., 2./) )
  deallocate(h1)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed sort and unique'

  ! unit tests for merge_interfaces
  test_name = 'h_L = h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/1., 2./), (/1., 2./), 1.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 0.5/) )
  deallocate(h1)

  test_name = 'h_L = h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/1., 2./), (/1., 2./), 0.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0.5, 0.5, 0.5/) )
  deallocate(h1)

  test_name = 'h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/1., 3./), (/2., 2./), 1.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 0.5/) )
  deallocate(h1)

  test_name = 'h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/1., 3./), (/2., 2./), 0.5, 1.5, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+1, test_name, h1, (/0.5, 0.5, 0.5/) )
  deallocate(h1)

  test_name = 'Left deeper than right, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/2., 3./), (/2., 2./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Left has zero thickness, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/4., 0./), (/2., 2./), 2.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk-1, test_name, h1, (/2./) )
  deallocate(h1)

  test_name = 'Left has zero thickness, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/4., 0./), (/2., 2./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Right has zero thickness, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk, (/2., 2./), (/0., 4./), 2.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk-1, test_name, h1, (/2./) )
  deallocate(h1)

  test_name = 'Right has zero thickness, h_L /= h_R and BLD_L /= BLD_R'
  call merge_interfaces(nk, (/2., 2./), (/0., 4./), 1.0, 2.0, CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/1., 1./) )
  deallocate(h1)

  test_name = 'Right deeper than left, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk+1, (/2., 2., 0./), (/2., 2., 1./), 4., 4., CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, h1, (/2., 2./) )
  deallocate(h1)

  test_name = 'Right and left small values at bottom, h_L /= h_R and BLD_L = BLD_R'
  call merge_interfaces(nk+2, (/2., 2., 1., 1./), (/1., 1., .5, .5/), 3., 3., CS%H_subroundoff, h1)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk+2, test_name, h1, (/1., 1., .5, .5/) )
  deallocate(h1)

  if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed merge interfaces'

  ! All cases in this section have hbl which are equal to the column thicknesses
  test_name = 'Equal hbl and same layer thicknesses (gradient from right to left)'
  hbl_L = 2.; hbl_R = 2.
  h_L = (/2.,2./) ; h_R = (/2.,2./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  khtr_u = 1.
  ! GMM, TODO: fix these calls
  !call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
  !                           khtr_u, F_layer, 1., 1., CS)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-2.0,0.0/) )

  !test_name = 'Equal hbl and same layer thicknesses (gradient from left to right)'
  !hbl_L = 2.; hbl_R = 2.
  !h_L = (/2.,2./) ; h_R = (/2.,2./)
  !phi_L = (/2.,1./) ; phi_R = (/1.,1./)
  !khtr_u = 0.5
  !call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
  !                           khtr_u, F_layer, 1., 1., CS)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/1.0,0.0/) )

  !test_name = 'hbl < column thickness, hbl same, linear profile right, khtr=2'
  !hbl_L = 2; hbl_R = 2
  !h_L = (/1.,2./) ; h_R = (/1.,2./)
  !phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  !khtr_u = 2.
  !call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
  !                           khtr_u, F_layer, 1., 1., CS)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.0,-4.0/) )

  !test_name = 'Different hbl and different column thicknesses (zero gradient)'
  !hbl_L = 12; hbl_R = 20
  !h_L = (/6.,6./) ; h_R = (/10.,10./)
  !phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  !khtr_u = 1.
  !call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
  !                           khtr_u, F_layer, 1., 1., CS)
  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.,0./) )

  !test_name = 'Different hbl and different column thicknesses (gradient from left to right)'

  !hbl_L = 15; hbl_R = 10.
  !h_L = (/10.,5./) ; h_R = (/10.,0./)
  !phi_L = (/1.,1./) ; phi_R = (/0.,0./)
  !khtr_u = 1.
  !call fluxes_layer_method(SURFACE, nk, hbl_L, hbl_R, h_L, h_R, phi_L, phi_R, &
  !                           khtr_u, F_layer, 1., 1., CS)

  !near_boundary_unit_tests = near_boundary_unit_tests .or. &
  !                           test_layer_fluxes( verbose, nk, test_name, F_layer, (/10.,0.0/) )

  !if (.not. near_boundary_unit_tests) write(stdout,*) 'Passed fluxes_layer_method'

end function near_boundary_unit_tests

!> Returns true if output of near-boundary unit tests does not match correct computed values
!! and conditionally writes results to stream
logical function test_layer_fluxes(verbose, nk, test_name, F_calc, F_ans)
  logical,                    intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),          intent(in) :: test_name !< Brief description of the unit test
  integer,                    intent(in) :: nk        !< Number of layers
  real, dimension(nk),        intent(in) :: F_calc    !< Fluxes of the unitless tracer from the algorithm [s^-1]
  real, dimension(nk),        intent(in) :: F_ans     !< Fluxes of the unitless tracer calculated by hand [s^-1]
  ! Local variables
  integer :: k

  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdout,*) "MOM_hor_bnd_diffusion, UNIT TEST FAILED: ", test_name
      write(stdout,10) k, F_calc(k), F_ans(k)
    elseif (verbose) then
      write(stdout,10) k, F_calc(k), F_ans(k)
    endif
  enddo

10 format("Layer=",i3," F_calc=",f20.16," F_ans",f20.16)
end function test_layer_fluxes

!> Return true if output of unit tests for boundary_k_range does not match answers
logical function test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, k_top_ans, zeta_top_ans,&
                                       k_bot_ans, zeta_bot_ans, test_name, verbose)
  integer :: k_top               !< Index of cell containing top of boundary
  real    :: zeta_top            !< Nondimension position
  integer :: k_bot               !< Index of cell containing bottom of boundary
  real    :: zeta_bot            !< Nondimension position
  integer :: k_top_ans           !< Index of cell containing top of boundary
  real    :: zeta_top_ans        !< Nondimension position
  integer :: k_bot_ans           !< Index of cell containing bottom of boundary
  real    :: zeta_bot_ans        !< Nondimension position
  character(len=80) :: test_name !< Name of the unit test
  logical :: verbose             !< If true always print output

  test_boundary_k_range = k_top /= k_top_ans
  test_boundary_k_range = test_boundary_k_range .or. (zeta_top /= zeta_top_ans)
  test_boundary_k_range = test_boundary_k_range .or. (k_bot /= k_bot_ans)
  test_boundary_k_range = test_boundary_k_range .or. (zeta_bot /= zeta_bot_ans)

  if (test_boundary_k_range) write(stdout,*) "UNIT TEST FAILED: ", test_name
  if (test_boundary_k_range .or. verbose) then
    write(stdout,20) "k_top", k_top, "k_top_ans", k_top_ans
    write(stdout,20) "k_bot", k_bot, "k_bot_ans", k_bot_ans
    write(stdout,30) "zeta_top", zeta_top, "zeta_top_ans", zeta_top_ans
    write(stdout,30) "zeta_bot", zeta_bot, "zeta_bot_ans", zeta_bot_ans
  endif

  20 format(A,"=",i3,X,A,"=",i3)
  30 format(A,"=",f20.16,X,A,"=",f20.16)


end function test_boundary_k_range

!> Deallocates hor_bnd_diffusion control structure
subroutine hor_bnd_diffusion_end(CS)
  type(hbd_CS), pointer :: CS  !< Horizontal boundary diffusion control structure

  if (associated(CS)) deallocate(CS)

end subroutine hor_bnd_diffusion_end

!> \namespace mom_hor_bnd_diffusion
!!
!! \section section_HBD The Horizontal Boundary Diffusion (HBD) framework
!!
!! The HBD framework accounts for the effects of diabatic mesoscale fluxes
!! within surface and bottom boundary layers. Unlike the equivalent adiabatic
!! fluxes, which is applied along neutral density surfaces, HBD is purely
!! horizontal. To assure that diffusive fluxes are strictly horizontal
!! regardless of the vertical coordinate system, this method relies on
!! regridding/remapping techniques.
!!
!! The bottom boundary layer fluxes remain to be implemented, although some
!! of the steps needed to do so have already been added and tested.
!!
!! Boundary horizontal diffusion is applied as follows:
!!
!! 1) remap tracer to a z* grid (HBD grid)
!! 2) calculate diffusive tracer fluxes (F) in the HBD grid using a layer by layer approach (@ref section_method)
!! 3) remap fluxes to the native grid
!! 4) update tracer by adding the divergence of F
!!
!! \subsection section_method Along layer approach
!!
!! Here diffusion is applied layer by layer using only information from neighboring cells.
!!
!! Step #1: define vertical grid using interfaces and surface boundary layers from left and right
!! columns (see merge_interfaces).
!!
!! Step #2: compute vertical indices containing boundary layer (boundary_k_range).
!! For the TOP boundary layer, these are:
!!
!! k_top, k_bot, zeta_top, zeta_bot
!!
!! Step #2: calculate the diffusive flux at each layer:
!!
!! \f[ F_{k} = -KHTR \times h_{eff}(k) \times (\phi_R(k) - \phi_L(k)),  \f]
!! where h_eff is the [harmonic mean](@ref section_harmonic_mean) of the layer thickness
!! in the left and right columns.
!!
!! Step #3: option to linearly decay the flux from k_bot_min to k_bot_max:
!!
!! If HBD_LINEAR_TRANSITION = True and k_bot_diff > 1, the diffusive flux will decay
!! linearly between the top interface of the layer containing the minimum boundary
!! layer depth (k_bot_min) and the lower interface of the layer containing the
!! maximum layer depth (k_bot_max).
!!
!! Step #4: remap the fluxes back to the native grid. This is done at velocity points, whose vertical grid
!! is determined using [harmonic mean](@ref section_harmonic_mean). To assure monotonicity,
!! tracer fluxes are limited so that 1) only down-gradient fluxes are applied,
!! and 2) the flux cannot be larger than F_max, which is defined using the tracer
!! gradient:
!!
!! \f[ F_{max} = -0.2 \times [(V_R(k) \times \phi_R(k)) - (V_L(k) \times \phi_L(k))],  \f]
!! where V is the cell volume. Why 0.2?
!!          t=0         t=inf
!!           0           .2
!!         0 1 0       .2.2.2
!!           0           .2
!!
!! \subsection section_harmonic_mean Harmonic Mean
!!
!! The harmonic mean (HM) betwen h1 and h2 is defined as:
!!
!! \f[ HM = \frac{2 \times h1 \times h2}{h1 + h2} \f]
!!
end module MOM_hor_bnd_diffusion
