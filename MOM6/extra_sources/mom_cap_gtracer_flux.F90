!> Contains routines for handling FMS coupler_bc_type tracer flux structures
!! when using MOM generic tracers

module MOM_cap_gtracer_flux

use MOM_domains,               only: domain2d
use MOM_coupler_types,         only: coupler_1d_bc_type, coupler_2d_bc_type
use MOM_ocean_model_nuopc,     only: ocean_model_flux_init
use coupler_types_mod,         only: coupler_type_register_restarts, coupler_type_restore_state
use atmos_ocean_fluxes_mod,    only: atmos_ocean_type_fluxes_init, atmos_ocean_fluxes_init
use fms2_io_mod,               only: FmsNetcdfDomainFile_t
use fms2_io_mod,               only: fms2_check_if_open => check_if_open, fms2_close_file => close_file
use fms2_io_mod,               only: fms2_write_data => write_data, fms2_read_restart => read_restart, fms2_write_restart => write_restart
use fms2_io_mod,               only: fms2_get_global_io_domain_indices => get_global_io_domain_indices
use field_manager_mod,         only: fm_field_name_len, fm_type_name_len, fm_loop_over_list, fm_change_list
use fm_util_mod,               only: fm_util_get_real_array
use mpp_mod,                   only: mpp_error, FATAL

implicit none; private

! Public member functions
public :: gas_exchange_init
public :: gas_fields_restore
public :: gas_fields_restart
public :: add_gas_fluxes_param
public :: get_coupled_field_name
public :: UNKNOWN_CMEPS_FIELD

character(len=*), parameter      :: mod_name = 'mom_cap_gtracer_flux'
character(len=*), parameter      :: UNKNOWN_CMEPS_FIELD = "UNKNOWN_FIELD"

!> FMS coupler_bc_types for additional tracer fields when using generic tracers
logical :: gas_fluxes_initialized = .false.  ! This is set to true when the following types are initialized.
type(coupler_1d_bc_type), target :: ex_gas_fields_atm ! tracer fields in atm
    !< Structure containing atmospheric surface variables that are used in the
    !! calculation of the atmosphere-ocean gas fluxes, as well as parameters
    !! regulating these fluxes. The fields in this structure are never actually
    !! set, but the structure is used for initialisation of components and to 
    !! spawn other structure whose fields are set.
type(coupler_1d_bc_type), target :: ex_gas_fields_ocn ! tracer fields atop the ocean
    !< Structure containing ocean surface variables that are used in the
    !! calculation of the atmosphere-ocean gas fluxes, as well as parameters
    !! regulating these fluxes. The fields in this structure are never actually
    !! set, but the structure is used for initialisation of components and to 
    !! spawn other structure whose fields are set.
type(coupler_1d_bc_type), target :: ex_gas_fluxes ! tracer fluxes between the atm and ocean
    !< A structure for exchanging gas or tracer fluxes between the atmosphere and
    !! ocean, defined by the field table, as well as a place holder of 
    !! intermediate calculations, such as piston velocities, and parameters that 
    !! impact the fluxes. The fields in this structure are never actually set, 
    !! but the structure is used for initialisation of components and to spawn 
    !! other structure whose fields are set.

contains

!> \brief Gas and tracer field initialization routine for running with MOM generic tracers.
!! Copied and adapted slightly from
!! https://github.com/NOAA-GFDL/FMScoupler/blob/77618869f48507c8629f28457cb701e25e1ea4fc/full/flux_exchange.F90#L626.
subroutine gas_exchange_init (gas_fields_atm, gas_fields_ocn, gas_fluxes)
    type(coupler_1d_bc_type), optional, pointer :: gas_fields_atm ! tracer fields in atm
        !< Pointer to a structure containing atmospheric surface variables that
        !! are used in the calculation of the atmosphere-ocean gas fluxes, as well
        !! as parameters regulating these fluxes.
    type(coupler_1d_bc_type), optional, pointer :: gas_fields_ocn ! tracer fields atop the ocean
        !< Pointer to a structure containing ocean surface variables that are 
        !! used in the calculation of the atmosphere-ocean gas fluxes, as well as
        !! parameters regulating these fluxes.
    type(coupler_1d_bc_type), optional, pointer :: gas_fluxes ! tracer fluxes between the atm and ocean
        !< Pointer to a structure for exchanging gas or tracer fluxes between the
        !! atmosphere and ocean, defined by the field table, as well as a place holder
        !! of intermediate calculations, such as piston velocities, and parameters
        !! that impact the fluxes.

    if (.not.gas_fluxes_initialized) then
        call atmos_ocean_type_fluxes_init( )
        call ocean_model_flux_init( )
        call atmos_ocean_fluxes_init(ex_gas_fluxes, ex_gas_fields_atm, ex_gas_fields_ocn)
        gas_fluxes_initialized = .true.
    endif

    if (present(gas_fields_atm)) gas_fields_atm => ex_gas_fields_atm
    if (present(gas_fields_ocn)) gas_fields_ocn => ex_gas_fields_ocn
    if (present(gas_fluxes)) gas_fluxes => ex_gas_fluxes

end subroutine gas_exchange_init

!> \brief Restore FMS coupler_bc_type state from ocean restart file
! See https://github.com/NOAA-GFDL/FMScoupler/blob/main/full/coupler_main.F90#L1896
subroutine gas_fields_restore(gas_fields, domain, directory)
    type(coupler_2d_bc_type), intent(inout) :: gas_fields !< FMS coupler_bc_type to be registered for restarts
    type(domain2D), intent(in)              :: domain     !< The FMS domain to use for this registration call
    character(len=*), optional, intent(in)  :: directory  !< Directory containing the restart file

    ! local variables
    type(FmsNetcdfDomainFile_t), dimension(:), pointer :: ocn_bc_restart => NULL() !< Structures describing the restart files
    integer                                            :: num_ocn_bc_restart = 0   !< The number of restart files to use
    integer                                            :: n
    
    
    call coupler_type_register_restarts(gas_fields, ocn_bc_restart, num_ocn_bc_restart, &
                domain, to_read=.true., ocean_restart=.true., directory=directory)

    ! Restore the fields from the restart files
    do n = 1, num_ocn_bc_restart
      if (fms2_check_if_open(ocn_bc_restart(n))) then
        call fms2_read_restart(ocn_bc_restart(n))
      endif
    enddo

    ! Check whether the restarts were read successfully.
    call coupler_type_restore_state(gas_fields, use_fms2_io=.true., test_by_field=.true.)

    do n = 1, num_ocn_bc_restart
      if(fms2_check_if_open(ocn_bc_restart(n))) call fms2_close_file(ocn_bc_restart(n))
    enddo

end subroutine gas_fields_restore

!> \brief Write ocean restart file for FMS coupler_bc_type state
! See https://github.com/NOAA-GFDL/FMScoupler/blob/main/full/coupler_main.F90#L2086
subroutine gas_fields_restart(gas_fields, domain, directory)
    type(coupler_2d_bc_type), intent(inout) :: gas_fields !< FMS coupler_bc_type to be registered for restarts
    type(domain2D), intent(in)              :: domain     !< The FMS domain to use for this registration call
    character(len=*), optional, intent(in)  :: directory  !< Directory containing the restart file
    
    ! local variables
    type(FmsNetcdfDomainFile_t), dimension(:), pointer :: ocn_bc_restart => NULL() !< Structures describing the restart files
    integer                                            :: num_ocn_bc_restart = 0   !< The number of restart files to use
    integer                                            :: n
    
    call coupler_type_register_restarts(gas_fields, ocn_bc_restart, num_ocn_bc_restart, &
            domain, to_read=.false., ocean_restart=.true., directory=directory)
    
    do n = 1, num_ocn_bc_restart
        if (fms2_check_if_open(ocn_bc_restart(n))) then
        call fms2_write_restart(ocn_bc_restart(n))
        call add_domain_dimension_data(ocn_bc_restart(n))
        call fms2_close_file(ocn_bc_restart(n))
        endif
    enddo

end subroutine gas_fields_restart

!> Register the axis data as a variable in the netcdf file and add some dummy data.
!! This is needed so the combiner can work correctly when the io_layout is not 1,1
!! Copied from https://github.com/NOAA-GFDL/FMScoupler/blob/77618869f48507c8629f28457cb701e25e1ea4fc/full/coupler_main.F90#L2026
subroutine add_domain_dimension_data(fileobj)
    type(FmsNetcdfDomainFile_t)        :: fileobj !< Fms2io domain decomposed fileobj

    ! local variables
    integer, dimension(:), allocatable :: buffer !< Buffer with axis data
    integer                            :: is, ie !< Starting and Ending indices for data
  
    call fms2_get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call fms2_write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)
  
    call fms2_get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call fms2_write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)
  
end subroutine add_domain_dimension_data

!> Retrieve param array from field_manager and add to FMS coupler_bc_type. This is
!! needed because the coupler_type_spawn routine does not copy the param array into
!! the spawned type. Hopefully we can get rid of this. This routine is based on
!! https://github.com/NOAA-GFDL/FMS/blob/7f585284f1487c0629f8075be350385e6e75dd2f/coupler/atmos_ocean_fluxes.F90#L448
subroutine add_gas_fluxes_param(gas_fluxes)
    type(coupler_2d_bc_type), intent(inout) :: gas_fluxes !< FMS coupler_bc_type to add param to

    ! local variables
    integer                                 :: n
    character(len=fm_field_name_len)        :: name
    character(len=fm_type_name_len)         :: typ
    integer                                 :: ind

    character(len=*), parameter             :: sub_name = 'add_gas_fluxes_param'
    character(len=*), parameter             :: error_header =&
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    n = 0
    do while (fm_loop_over_list('/coupler_mod/fluxes', name, typ, ind))
        if (typ .ne. 'list') then
            call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
        endif

        n = n + 1

        if (.not. fm_change_list('/coupler_mod/fluxes/' // trim(name))) then
            call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
        endif

        if (gas_fluxes%bc(n)%name .eq. name) then
            gas_fluxes%bc(n)%param => fm_util_get_real_array('param')
        else
            call mpp_error(FATAL, trim(error_header) // ' Problem setting param array pointer')
        endif
    enddo
end subroutine add_gas_fluxes_param

!> Return the CMEPS standard_name of the coupled field required for a given coupled
!! generic_tracer flux name.
function get_coupled_field_name(name)
    character(len=64)                  :: get_coupled_field_name !< CMEPS standard_name
    character(len=*), intent(in)       :: name                   !< gtracer flux name

    select case(trim(name))
        case( 'co2_flux' )
            get_coupled_field_name = "Sa_co2prog"
        case( 'o2_flux' )
            get_coupled_field_name = "Sa_o2"
        case default
            get_coupled_field_name = UNKNOWN_CMEPS_FIELD
    end select
end function get_coupled_field_name

end module MOM_cap_gtracer_flux