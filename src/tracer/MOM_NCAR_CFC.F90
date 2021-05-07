!> Simulates CFCs using the OCMIP2 protocols
module MOM_NCAR_CFC

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coupler_types, only : extract_coupler_type_data, set_coupler_type_data
use MOM_coupler_types, only : atmos_ocn_coupler_flux
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_NCAR_CFC, initialize_NCAR_CFC !, flux_init_NCAR_CFC
public NCAR_CFC_column_physics, NCAR_CFC_surface_state, NCAR_CFC_fluxes
public NCAR_CFC_stock, get_surface_CFC, NCAR_CFC_end


integer, parameter :: NTR = 2 !< the number of tracers in this module.

!> The control structure for the  NCAR_CFC tracer package
type, public :: NCAR_CFC_CS ; private
  character(len=200) :: IC_file !< The file in which the CFC initial values can
                                !! be found, or an empty string for internal initilaization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM6 tracer registry
  real, pointer, dimension(:,:,:) :: &
    CFC11 => NULL(), &     !< The CFC11 concentration [mol m-3].
    CFC12 => NULL()        !< The CFC12 concentration [mol m-3].
  ! In the following variables a suffix of _11 refers to CFC11 and _12 to CFC12.
  !>@{ Coefficients used in the CFC11 and CFC12 solubility calculation
  real :: a1_11, a1_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [nondim]
  real :: a2_11, a2_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-1]
  real :: a3_11, a3_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-2]
  real :: a4_11, a4_12   ! Coefficients for calculating CFC11 and CFC12 Schmidt numbers [degC-3]

  real :: d1_11, d1_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [nondim]
  real :: d2_11, d2_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-1]
  real :: d3_11, d3_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [log(hectoKelvin)-1]
  real :: d4_11, d4_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-2]

  real :: e1_11, e1_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1]
  real :: e2_11, e2_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1 hectoKelvin-1]
  real :: e3_11, e3_12   ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-2 hectoKelvin-2]
  !>@}
  real :: CFC11_IC_val = 0.0    !< The initial value assigned to CFC11 [mol m-3].
  real :: CFC12_IC_val = 0.0    !< The initial value assigned to CFC12 [mol m-3].
  real :: CFC11_land_val = -1.0 !< The value of CFC11 used where land is masked out [mol m-3].
  real :: CFC12_land_val = -1.0 !< The value of CFC12 used where land is masked out [mol m-3].
  logical :: tracers_may_reinit !< If true, tracers may be reset via the initialization code
                                !! if they are not found in the restart files.
  character(len=16) :: CFC11_name !< CFC11 variable name
  character(len=16) :: CFC12_name !< CFC12 variable name

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                                             !! the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()  !< Model restart control structure

  ! The following vardesc types contain a package of metadata about each tracer.
  type(vardesc) :: CFC11_desc !< A set of metadata for the CFC11 tracer
  type(vardesc) :: CFC12_desc !< A set of metadata for the CFC12 tracer
end type NCAR_CFC_CS

contains

!> Register the OCMIP2 CFC tracers to be used with MOM and read the parameters
!! that are used with this tracer package
function register_NCAR_CFC(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(NCAR_CFC_CS),       pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module.
  type(tracer_registry_type), &
                           pointer    :: tr_Reg     !< A pointer to the tracer registry.
  type(MOM_restart_CS),    pointer    :: restart_CS !< A pointer to the restart control structure.
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.

  ! Local variables
  character(len=40)  :: mdl = "MOM_NCAR_CFC" ! This module's name.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  real, dimension(:,:,:), pointer :: tr_ptr => NULL()
  real :: a11_dflt(4), a12_dflt(4) ! Default values of the various coefficients
  real :: d11_dflt(4), d12_dflt(4) ! In the expressions for the solubility and
  real :: e11_dflt(3), e12_dflt(3) ! Schmidt numbers.
  character(len=48) :: flux_units ! The units for tracer fluxes.
  logical :: register_NCAR_CFC
  integer :: isd, ied, jsd, jed, nz, m

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_NCAR_CFC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! GMM, the following call might not be necessary

  ! This call sets default properties for the air-sea CFC fluxes and obtains the
  ! indicies for the CFC11 and CFC12 flux coupling.
  !call flux_init_NCAR_CFC(CS, verbosity=3)
  !if ((CS%ind_cfc_11_flux < 0) .or. (CS%ind_cfc_12_flux < 0)) then
  !  ! This is most likely to happen with the dummy version of atmos_ocn_coupler_flux
  !  ! used in ocean-only runs.
  !  call MOM_ERROR(WARNING, "CFCs are currently only set up to be run in " // &
  !                 " coupled model configurations, and will be disabled.")
  !  deallocate(CS)
  !  register_NCAR_CFC = .false.
  !  return
  !endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CFC_IC_FILE", CS%IC_file, &
                 "The file in which the CFC initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/CFC_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "CFC_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, CFC_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)

  !   The following vardesc types contain a package of metadata about each tracer,
  ! including, the name; units; longname; and grid information.
  CS%CFC11_name = "CFC11" ; CS%CFC12_name = "CFC12"
  CS%CFC11_desc = var_desc(CS%CFC11_name,"mol m-3","CFC-11 Concentration", caller=mdl)
  CS%CFC12_desc = var_desc(CS%CFC12_name,"mol m-3","CFC-12 Concentration", caller=mdl)
  if (GV%Boussinesq) then ; flux_units = "mol s-1"
  else ; flux_units = "mol m-3 kg s-1" ; endif

  allocate(CS%CFC11(isd:ied,jsd:jed,nz)) ; CS%CFC11(:,:,:) = 0.0
  allocate(CS%CFC12(isd:ied,jsd:jed,nz)) ; CS%CFC12(:,:,:) = 0.0

  ! This pointer assignment is needed to force the compiler not to do a copy in
  ! the registration calls.  Curses on the designers and implementers of F90.
  tr_ptr => CS%CFC11
  ! Register CFC11 for horizontal advection, diffusion, and restarts.

  ! CFC11 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/42625c97b8fe75124a345962c4430982.html
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC11_desc, registry_diags=.true., &
                       cmor_name='cfc11' ,&
                       cmor_units='mol m-3'  ,&
                       cmor_longname='Mole Concentration of CFC11 in Sea Water' ,&
                       flux_units=flux_units, &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)
  ! Do the same for CFC12
  ! CFC12 cmor conventions: http://clipc-services.ceda.ac.uk/dreq/u/3ab8e10027d7014f18f9391890369235.html
  tr_ptr => CS%CFC12
  call register_tracer(tr_ptr, Tr_Reg, param_file, HI, GV, &
                       tr_desc=CS%CFC12_desc, registry_diags=.true., &
                       cmor_name='cfc12' ,&
                       cmor_units='mol m-3'  ,&
                       cmor_longname='Mole Concentration of CFC12 in Sea Water' ,&
                       flux_units=flux_units, &
                       restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

  ! Set and read the various empirical coefficients.

!-----------------------------------------------------------------------
! Default Schmidt number coefficients for CFC11 (_11) and CFC12 (_12) are given
! by Zheng et al (1998), JGR vol 103, C1.
!-----------------------------------------------------------------------
  a11_dflt(:) = (/ 3501.8, -210.31,  6.1851, -0.07513 /)
  a12_dflt(:) = (/ 3845.4, -228.95,  6.1908, -0.06743 /)
  call get_param(param_file, mdl, "CFC11_A1", CS%a1_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="nondim", default=a11_dflt(1))
  call get_param(param_file, mdl, "CFC11_A2", CS%a2_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-1", default=a11_dflt(2))
  call get_param(param_file, mdl, "CFC11_A3", CS%a3_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-2", default=a11_dflt(3))
  call get_param(param_file, mdl, "CFC11_A4", CS%a4_11, &
                 "A coefficient in the Schmidt number of CFC11.", &
                 units="degC-3", default=a11_dflt(4))

  call get_param(param_file, mdl, "CFC12_A1", CS%a1_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="nondim", default=a12_dflt(1))
  call get_param(param_file, mdl, "CFC12_A2", CS%a2_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-1", default=a12_dflt(2))
  call get_param(param_file, mdl, "CFC12_A3", CS%a3_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-2", default=a12_dflt(3))
  call get_param(param_file, mdl, "CFC12_A4", CS%a4_12, &
                 "A coefficient in the Schmidt number of CFC12.", &
                 units="degC-3", default=a12_dflt(4))

!-----------------------------------------------------------------------
! Solubility coefficients for alpha in mol/l/atm for CFC11 (_11) and CFC12 (_12)
! after Warner and Weiss (1985) DSR, vol 32.
!-----------------------------------------------------------------------
  d11_dflt(:) = (/ -229.9261, 319.6552, 119.4471, -1.39165 /)
  e11_dflt(:) = (/ -0.142382, 0.091459, -0.0157274 /)
  d12_dflt(:) = (/ -218.0971, 298.9702, 113.8049, -1.39165 /)
  e12_dflt(:) = (/ -0.143566, 0.091015, -0.0153924 /)

  call get_param(param_file, mdl, "CFC11_D1", CS%d1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(1))
  call get_param(param_file, mdl, "CFC11_D2", CS%d2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK", default=d11_dflt(2))
  call get_param(param_file, mdl, "CFC11_D3", CS%d3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="none", default=d11_dflt(3))
  call get_param(param_file, mdl, "CFC11_D4", CS%d4_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="hK-2", default=d11_dflt(4))
  call get_param(param_file, mdl, "CFC11_E1", CS%e1_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1", default=e11_dflt(1))
  call get_param(param_file, mdl, "CFC11_E2", CS%e2_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-1", default=e11_dflt(2))
  call get_param(param_file, mdl, "CFC11_E3", CS%e3_11, &
                 "A coefficient in the solubility of CFC11.", &
                 units="PSU-1 hK-2", default=e11_dflt(3))

  call get_param(param_file, mdl, "CFC12_D1", CS%d1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(1))
  call get_param(param_file, mdl, "CFC12_D2", CS%d2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK", default=d12_dflt(2))
  call get_param(param_file, mdl, "CFC12_D3", CS%d3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="none", default=d12_dflt(3))
  call get_param(param_file, mdl, "CFC12_D4", CS%d4_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="hK-2", default=d12_dflt(4))
  call get_param(param_file, mdl, "CFC12_E1", CS%e1_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1", default=e12_dflt(1))
  call get_param(param_file, mdl, "CFC12_E2", CS%e2_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-1", default=e12_dflt(2))
  call get_param(param_file, mdl, "CFC12_E3", CS%e3_12, &
                 "A coefficient in the solubility of CFC12.", &
                 units="PSU-1 hK-2", default=e12_dflt(3))

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS

  register_NCAR_CFC = .true.
end function register_NCAR_CFC

!!! GMM we might not need this subroutine
!> This subroutine initializes the air-sea CFC fluxes, and optionally returns
!! the indicies of these fluxes.  It can safely be called multiple times.
!subroutine flux_init_NCAR_CFC(CS, verbosity)
!  type(NCAR_CFC_CS), optional, pointer :: CS !< An optional pointer to the control structure
!                                               !! for this module; if not present, the flux indicies
!                                               !! are not stored.
!  integer,             optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.
!
!  ! These can be overridden later in via the field manager?
!  character(len=128) :: default_ice_restart_file = 'ice_ocmip2_cfc.res.nc'
!  character(len=128) :: default_ocean_restart_file = 'ocmip2_cfc.res.nc'
!  integer :: ind_flux(2) ! Integer indices of the fluxes
!
!  ! These calls obtain the indices for the CFC11 and CFC12 flux coupling.  They
!  ! can safely be called multiple times.
!!  ind_flux(1) = atmos_ocn_coupler_flux('cfc_11_flux', &
!!       flux_type = 'air_sea_gas_flux', implementation='ocmip2', &
!!       param=(/ 9.36e-07, 9.7561e-06 /), &
!!       ice_restart_file = default_ice_restart_file, &
!!       ocean_restart_file = default_ocean_restart_file, &
!!       caller = "register_NCAR_CFC", verbosity=verbosity)
!!  ind_flux(2) = atmos_ocn_coupler_flux('cfc_12_flux', &
!!       flux_type='air_sea_gas_flux', implementation='ocmip2', &
!!       param = (/ 9.36e-07, 9.7561e-06 /), &
!!       ice_restart_file=default_ice_restart_file, &
!!       ocean_restart_file=default_ocean_restart_file, &
!!       caller="register_NCAR_CFC", verbosity=verbosity)
!
!  ind_flux(1) = 1
!  ind_flux(2) = 2
!  if (present(CS)) then ; if (associated(CS)) then
!    CS%ind_cfc_11_flux = ind_flux(1)
!    CS%ind_cfc_12_flux = ind_flux(2)
!  endif ; endif
!
!end subroutine flux_init_NCAR_CFC

!> Initialize the OCMP2 CFC tracer fields and set up the tracer output.
subroutine initialize_NCAR_CFC(restart, day, G, GV, US, h, diag, OBC, CS, &
                                 sponge_CSp)
  logical,                        intent(in) :: restart    !< .true. if the fields have already been
                                                           !! read from a restart file.
  type(time_type), target,        intent(in) :: day        !< Time of the start of the run.
  type(ocean_grid_type),          intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in) :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(in) :: h          !< Layer thicknesses [H ~> m or kg m-2].
  type(diag_ctrl), target,        intent(in) :: diag       !< A structure that is used to regulate
                                                           !! diagnostic output.
  type(ocean_OBC_type),           pointer    :: OBC        !< This open boundary condition type
                                                           !! specifies whether, where, and what
                                                           !! open boundary conditions are used.
  type(NCAR_CFC_CS),              pointer    :: CS         !< The control structure returned by a
                                                           !! previous call to register_NCAR_CFC.
  type(sponge_CS),                pointer    :: sponge_CSp !< A pointer to the control structure for
                                                           !! the sponges, if they are in use.
                                                           !! Otherwise this may be unassociated.
!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

  logical :: from_file = .false.

  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag


  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC11, CS%CFC11_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC11, CS%CFC11_name, CS%CFC11_land_val, &
                         CS%CFC11_IC_val, G, GV, US, CS)

  if (.not.restart .or. (CS%tracers_may_reinit .and. &
      .not.query_initialized(CS%CFC12, CS%CFC12_name, CS%restart_CSp))) &
    call init_tracer_CFC(h, CS%CFC12, CS%CFC12_name, CS%CFC12_land_val, &
                         CS%CFC12_IC_val, G, GV, US, CS)

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  ! GMM: TODO this must be coded
  endif

end subroutine initialize_NCAR_CFC

!>This subroutine initializes a tracer array.
subroutine init_tracer_CFC(h, tr, name, land_val, IC_val, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: tr   !< The tracer concentration array
  character(len=*),                          intent(in)  :: name !< The tracer name
  real,                                      intent(in)  :: land_val !< A value the tracer takes over land
  real,                                      intent(in)  :: IC_val   !< The initial condition value for the tracer
  type(NCAR_CFC_CS),                         pointer     :: CS   !< The control structure returned by a
                                                                 !! previous call to register_NCAR_CFC.

  ! This subroutine initializes a tracer array.

  logical :: OK
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (len_trim(CS%IC_file) > 0) then
    !  Read the tracer concentrations from a netcdf file.
    if (.not.file_exists(CS%IC_file, G%Domain)) &
      call MOM_error(FATAL, "initialize_NCAR_CFC: Unable to open "//CS%IC_file)
    if (CS%Z_IC_file) then
      OK = tracer_Z_init(tr, h, CS%IC_file, name, G, GV, US)
      if (.not.OK) then
        OK = tracer_Z_init(tr, h, CS%IC_file, trim(name), G, GV, US)
        if (.not.OK) call MOM_error(FATAL,"initialize_NCAR_CFC: "//&
                "Unable to read "//trim(name)//" from "//&
                trim(CS%IC_file)//".")
      endif
    else
      call MOM_read_data(CS%IC_file, trim(name), tr, G%Domain)
    endif
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) < 0.5) then
        tr(i,j,k) = land_val
      else
        tr(i,j,k) = IC_val
      endif
    enddo ; enddo ; enddo
  endif

end subroutine init_tracer_CFC

!>  This subroutine applies diapycnal diffusion, souces and sinks and any other column
!! tracer physics or chemistry to the OCMIP2 CFC tracers.
!! CFCs are relatively simple, as they are passive tracers with only a surface flux as a source.
subroutine NCAR_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(NCAR_CFC_CS),     pointer    :: CS   !< The control structure returned by a
                                              !! previous call to register_NCAR_CFC.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! CFCs are relatively simple, as they are passive tracers. with only a surface
! flux as a source.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real, pointer, dimension(:,:,:) :: CFC11 => NULL(), CFC12 => NULL()
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  integer :: i, j, k, m, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  CFC11 => CS%CFC11 ; CFC12 => CS%CFC12

  ! These two calls unpack the fluxes from the input arrays.
  !   The -GV%Rho0 changes the sign convention of the flux and changes the units
  ! of the flux from [Conc. m s-1] to [Conc. kg m-2 T-1].
  !call extract_coupler_type_data(fluxes%tr_fluxes, CS%ind_cfc_11_flux, CFC11_flux, &
  !                               scale_factor=-G%US%R_to_kg_m3*GV%Rho0*US%T_to_s, idim=idim, jdim=jdim)
  !call extract_coupler_type_data(fluxes%tr_fluxes, CS%ind_cfc_12_flux, CFC12_flux, &
  !                               scale_factor=-G%US%R_to_kg_m3*GV%Rho0*US%T_to_s, idim=idim, jdim=jdim)

  ! Use a tridiagonal solver to determine the concentrations after the
  ! surface source is applied and diapycnal advection and diffusion occurs.
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CFC11, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CFC11, G, GV, sfc_flux=fluxes%CFC11_flux)

    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CFC12, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CFC12, G, GV, sfc_flux=fluxes%CFC12_flux)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CFC11, G, GV, sfc_flux=fluxes%CFC11_flux)
    call tracer_vertdiff(h_old, ea, eb, dt, CFC12, G, GV, sfc_flux=fluxes%CFC12_flux)
  endif

  ! Write out any desired diagnostics from tracer sources & sinks here.

end subroutine NCAR_CFC_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function NCAR_CFC_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),           intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),         intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(:),              intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                           !! tracer, in kg times concentration units [kg conc].
  type(NCAR_CFC_CS),             pointer       :: CS     !< The control structure returned by a
                                                           !! previous call to register_NCAR_CFC.
  character(len=*), dimension(:),  intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),  intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,               intent(in)    :: stock_index !< The coded index of a specific
                                                                !! stock being sought.
  integer                                        :: NCAR_CFC_stock !< The number of stocks calculated here.

  ! Local variables
  real :: mass
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  NCAR_CFC_stock = 0
  if (.not.associated(CS)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  call query_vardesc(CS%CFC11_desc, name=names(1), units=units(1), caller="NCAR_CFC_stock")
  call query_vardesc(CS%CFC12_desc, name=names(2), units=units(2), caller="NCAR_CFC_stock")
  units(1) = trim(units(1))//" kg" ; units(2) = trim(units(2))//" kg"

  stocks(1) = 0.0 ; stocks(2) = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    mass = G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j) * h(i,j,k)
    stocks(1) = stocks(1) + CS%CFC11(i,j,k) * mass
    stocks(2) = stocks(2) + CS%CFC12(i,j,k) * mass
  enddo ; enddo ; enddo
  stocks(1) = GV%H_to_kg_m2 * stocks(1)
  stocks(2) = GV%H_to_kg_m2 * stocks(2)

  NCAR_CFC_stock = 2

end function NCAR_CFC_stock

!> This subroutine extracts the surface CFC concentrations and put them in sfc_state.
subroutine NCAR_CFC_surface_state(sfc_state, G, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  type(NCAR_CFC_CS),       pointer       :: CS !< The control structure returned by a previous
                                               !! call to register_NCAR_CFC.

  ! Local variables
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  do j=js,je ; do i=is,ie
    sfc_state%sfc_cfc11(i,j) = CS%CFC11(i,j,1)
    sfc_state%sfc_cfc12(i,j) = CS%CFC12(i,j,1)
  enddo ; enddo

end subroutine NCAR_CFC_surface_state

!> This subroutine calculate CFC fluxes.
subroutine NCAR_CFC_fluxes(cfc11_atm, cfc12_atm, fluxes, sfc_state, G) !, GV, CS)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: cfc11_atm !< cfc11 atm concentration [??].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: cfc12_atm !< cfc12 atm concentration [??].
  !type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure
  type(surface),           intent(in) :: sfc_state !< A structure containing fields that
                                                   !! describe the surface state of the ocean.
  type(forcing),           intent(inout) :: fluxes !< A structure containing pointers to thermodynamic
                                                !! and tracer forcing fields.  Unused fields have NULL ptrs.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    CFC11_Csurf, &  ! The CFC-11 surface concentrations times the Schmidt number term [mol m-3].
    CFC12_Csurf, &  ! The CFC-12 surface concentrations times the Schmidt number term [mol m-3].
    CFC11_alpha, &  ! The CFC-11 solubility [mol m-3 pptv-1].
    CFC12_alpha, &  ! The CFC-12 solubility [mol m-3 pptv-1].
    kw,          &  ! gas transfer velocity [?????] GMM TODO
    cair            ! The surface gas concentration in equilibrium with the atmosphere (saturation concentration)
  real :: ta        ! Absolute sea surface temperature [hectoKelvin] (Why use such bizzare units?)
  real :: sal       ! Surface salinity [PSU].
  real :: SST       ! Sea surface temperature [degC].
  real :: alpha_11  ! The solubility of CFC 11 [mol m-3 pptv-1].
  real :: alpha_12  ! The solubility of CFC 12 [mol m-3 pptv-1].
  real :: sc_11, sc_12 ! The Schmidt numbers of CFC 11 and CFC 12.
  real :: sc_no_term   ! A term related to the Schmidt number.
  real :: d11_dflt(4), d12_dflt(4) ! values of the various coefficients
  real :: e11_dflt(3), e12_dflt(3) ! in the expressions for the solubility
  real :: d1_11, d1_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [nondim]
  real :: d2_11, d2_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-1]
  real :: d3_11, d3_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [log(hectoKelvin)-1]
  real :: d4_11, d4_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [hectoKelvin-2]
  real :: e1_11, e1_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1]
  real :: e2_11, e2_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-1 hectoKelvin-1]
  real :: e3_11, e3_12 ! Coefficients for calculating CFC11 and CFC12 solubilities [PSU-2 hectoKelvin-2]
  integer :: i, j, m, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  !-----------------------------------------------------------------------
  ! Solubility coefficients for alpha in mol/l/atm for CFC11 (_11) and CFC12 (_12)
  ! after Warner and Weiss (1985) DSR, vol 32.
  !-----------------------------------------------------------------------
  d11_dflt(:) = (/ -229.9261, 319.6552, 119.4471, -1.39165 /)
  e11_dflt(:) = (/ -0.142382, 0.091459, -0.0157274 /)
  d12_dflt(:) = (/ -218.0971, 298.9702, 113.8049, -1.39165 /)
  e12_dflt(:) = (/ -0.143566, 0.091015, -0.0153924 /)

  d1_11 = d11_dflt(1)
  d2_11 = d11_dflt(2)
  d3_11 = d11_dflt(3)
  d4_11 = d11_dflt(4)

  e1_11 = e11_dflt(1)
  e2_11 = e11_dflt(2)
  e3_11 = e11_dflt(3)


  d1_12 = d12_dflt(1)
  d2_12 = d12_dflt(2)
  d3_12 = d12_dflt(3)
  d4_12 = d12_dflt(4)

  e1_12 = e12_dflt(1)
  e2_12 = e12_dflt(2)
  e3_12 = e12_dflt(3)

  ! GMM, TODO: add OMP calls

  do j=js,je ; do i=is,ie
    ! cite a paper that uses this formula?
    ta = max(0.01, (sfc_state%SST(i,j) + 273.15) * 0.01) ! Why is this in hectoKelvin?
    sal = sfc_state%SSS(i,j) ; SST = sfc_state%SST(i,j)
    !    Calculate solubilities using Warner and Weiss (1985) DSR, vol 32.
    ! The final result is in mol/cm3/pptv (1 part per trillion 1e-12)
    ! Use Bullister and Wisegavger for CCl4.
    ! The factor 1.e-09 converts from mol/(l * atm) to mol/(m3 * pptv).
    alpha_11 = exp(d1_11 + d2_11/ta + d3_11*log(ta) + d4_11*ta**2 +&
                   sal * ((e3_11 * ta + e2_11) * ta + e1_11)) * &
               1.0e-09 * G%mask2dT(i,j)
    alpha_12 = exp(d1_12 + d2_12/ta + d3_12*log(ta) + d4_12*ta**2 +&
                   sal * ((e3_12 * ta + e2_12) * ta + e1_12)) * &
               1.0e-09 * G%mask2dT(i,j)

    ! Calculate Schmidt numbers using coefficients given by
    ! Wanninkhof (2014); doi:10.4319/lom.2014.12.351.
    call comp_CFC_schmidt(sfc_state%SST(i,j), sc_11, sc_12, G%mask2dT(i,j))
    !write(*,*)'sfc_state%SST(i,j), sc_11, sc_12, G%mask2dT(i,j)',sfc_state%SST(i,j), sc_11, sc_12, G%mask2dT(i,j)
    ! GMM, the following was stopping the model because of unbounded values. Using the GFDL approach instead.
    !sc_no_term = sqrt(660.0 / sc_11)
    ! The abs here is to avoid NaNs. The model should be failing at this point.
    sc_no_term = sqrt(660.0 / (abs(sc_11) + 1.0e-30))

    CFC11_alpha(i,j) = alpha_11 * sc_no_term
    CFC11_Csurf(i,j) = sfc_state%sfc_CFC11(i,j) * sc_no_term
    sc_no_term = sqrt(660.0 / (abs(sc_12) + 1.0e-30))
    CFC12_alpha(i,j) = alpha_12 * sc_no_term
    CFC12_Csurf(i,j) = sfc_state%sfc_CFC12(i,j) * sc_no_term

    !---------------------------------------------------------------------
    !     Gas exchange/piston velocity parameter
    !---------------------------------------------------------------------
    ! pop uses xkw_coeff = 6.97e-9_r8 ! in s/cm, from a = 0.251 cm/hr s^2/m^2 in Wannikhof 2014
    ! 6.97e-07 m/s
    kw(i,j) = 6.97e-07 *  ((1 - fluxes%ice_fraction(i,j))*fluxes%u10_sqr(i,j))

    ! air concentrations and cfcs BC's fluxes
    cair(i,j) = 9.7561e-06 * CFC11_alpha(i,j) * cfc11_atm(i,j) * fluxes%p_surf_full(i,j)
    fluxes%cfc11_flux(i,j) = kw(i,j) * (cair(i,j) - CFC11_Csurf(i,j))
    cair(i,j) = 9.7561e-06 * CFC12_alpha(i,j) * cfc12_atm(i,j) * fluxes%p_surf_full(i,j)
    fluxes%cfc12_flux(i,j) = kw(i,j) * (cair(i,j) - CFC12_Csurf(i,j))
  enddo ; enddo

  ! GMM, TODO

  !   These calls load these values into the appropriate arrays in the
  ! coupler-type structure.

  ! GMM commented here:
  ! These calls are setting a single 2-d field in a coupler_3d_bc_type from a two-dimensional array

  !call set_coupler_type_data(CFC11_alpha, CS%ind_cfc_11_flux, sfc_state%tr_fields, &
  !                           solubility=.true., idim=idim, jdim=jdim)
  !call set_coupler_type_data(CFC11_Csurf, CS%ind_cfc_11_flux, sfc_state%tr_fields, &
  !                           idim=idim, jdim=jdim)
  !call set_coupler_type_data(CFC12_alpha, CS%ind_cfc_12_flux, sfc_state%tr_fields, &
  !                           solubility=.true., idim=idim, jdim=jdim)
  !call set_coupler_type_data(CFC12_Csurf, CS%ind_cfc_12_flux, sfc_state%tr_fields, &
  !                           idim=idim, jdim=jdim)

end subroutine NCAR_CFC_fluxes

!> Compute Schmidt numbers of CFCs following Wanninkhof (2014); doi:10.4319/lom.2014.12.351
!! Range of validity of fit is -2:40.
subroutine comp_CFC_schmidt(sst_in, cfc11_sc, cfc12_sc, mask)
  real, intent(in)    :: sst_in   !< The sea surface temperature [degC].
  real, intent(inout) :: cfc11_sc !< Schmidt number of CFC11 [nondim].
  real, intent(inout) :: cfc12_sc !< Schmidt number of CFC12 [nondim].
  real, intent(in)    :: mask     !< 0 for land points and 1 for ocean points on the h-grid [nondim].

  !local variables
  real , parameter :: a_11 = 3579.2
  real , parameter :: b_11 = -222.63
  real , parameter :: c_11 = 7.5749
  real , parameter :: d_11 = -0.14595
  real , parameter :: e_11 = 0.0011874
  real , parameter :: a_12 = 3828.1
  real , parameter :: b_12 = -249.86
  real , parameter :: c_12 = 8.7603
  real , parameter :: d_12 = -0.1716
  real , parameter :: e_12 = 0.001408
  real             :: sst


  if (mask>0.) then
    ! clip SST to avoid bad values
    sst = MAX(-2.0, MIN(40.0, sst_in))
    cfc11_sc = a_11 + sst * (b_11 + sst * (c_11 + sst * (d_11 + sst * e_11)))
    cfc12_sc = a_12 + sst * (b_12 + sst * (c_12 + sst * (d_12 + sst * e_12)))
  else
    cfc11_sc = 0.0
    cfc12_sc = 0.0
  endif

end subroutine comp_CFC_schmidt

!> Copies surface CFCs from control structure into CFC11 and CFC12.
subroutine get_surface_CFC(CS, CFC11, CFC12, G, US)
  type(NCAR_CFC_CS),                pointer       :: CS  !< The control structure returned by a previous
                                                         !! call to register_NCAR_CFC.
  type(ocean_grid_type),            intent(in)    :: G   !< Grid structure.
  type(unit_scale_type),            intent(in)    :: US  !< A dimensional unit scaling type.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: CFC11 !< Surface CFC11 concentration [???].
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: CFC12 !< Surface CFC12 concentration [???].

  ! Local variables
  integer :: i,j

  !scale = US%m_to_Z ; if (present(m_to_BLD_units)) scale = m_to_BLD_units

  !$OMP parallel do default(shared)
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    CFC11(i,j) = CS%CFC11(i,j,1)
    CFC12(i,j) = CS%CFC12(i,j,1)
  enddo ; enddo

end subroutine get_surface_CFC

!> Deallocate any memory associated with the OCMIP2 CFC tracer package
subroutine NCAR_CFC_end(CS)
  type(NCAR_CFC_CS), pointer :: CS   !< The control structure returned by a
                                       !! previous call to register_NCAR_CFC.
!   This subroutine deallocates the memory owned by this module.
! Argument: CS - The control structure returned by a previous call to
!                register_NCAR_CFC.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%CFC11)) deallocate(CS%CFC11)
    if (associated(CS%CFC12)) deallocate(CS%CFC12)

    deallocate(CS)
  endif
end subroutine NCAR_CFC_end


!> \namespace mom_ncar_cfc
!!
!!   By Gustavo Marques, 2021
!!
!!     This module contains the code that is needed to set
!!   up and use CFC-11 and CFC-12 in CESM using the OCMIP2 protocols

end module MOM_NCAR_CFC
