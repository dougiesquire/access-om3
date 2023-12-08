!> Contains routine for calculating FMS coupler_bc_type ocean tracer
!! fluxes when using MOM generic tracers
!
!! This module was copied from FMScoupler at
!! https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_fluxes_calc.F90
!! and subsequently modified in the following ways:
!! - Operate on 2D inputs, rather than 1D
!! - Add calculation for 'air_sea_deposition' taken from
!!   https://github.com/NOAA-GFDL/FMScoupler/blob/6442d387153064644325c96a5e9e2935139d5e3c/full/atmos_ocean_dep_fluxes_calc.F90
!! - Multiply fluxes by ice_fraction input, rather than masking based on seawater input
!! - Use MOM over FMS modules where easy to do so
!! - Make tsurf input optional, as it is only used by a few implementations
!! - Use ind_runoff rather than ind_deposition in runoff flux calculation (note, their
!!   values are equal)
!! - Rename gas_fields_ice to gas_fields_ocn

module atmos_ocean_fluxes_calc_mod

use MOM_coupler_types,         only: coupler_2d_bc_type
use MOM_coupler_types,         only: ind_flux, ind_deltap, ind_kw, ind_flux0
use MOM_coupler_types,         only: ind_pcair, ind_u10, ind_psurf
use MOM_coupler_types,         only: ind_alpha, ind_csurf, ind_sc_no
use MOM_coupler_types,         only: ind_runoff, ind_deposition
use mpp_mod,                   only: mpp_error, FATAL
use FMSconstants,              only: wtmair, rdgas, vonkarm

implicit none; private

public atmos_ocean_fluxes_calc

character(len=*), parameter :: mod_name = "atmos_ocean_fluxes_calc_mod"
real, parameter             :: epsln=1.0e-30

contains

!> \brief Calculate the FMS coupler_bc_type ocean tracer fluxes. Units should be mol/m^2/s.
!! Upward flux is positive.
subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ocn, gas_fluxes,&
  ice_fraction, isc, iec, jsc, jec, tsurf, ustar, cd_m)
  type(coupler_2d_bc_type), intent(in)     :: gas_fields_atm ! fields in atm
      !< Structure containing atmospheric surface variables that are used in the calculation
      !! of the atmosphere-ocean tracer fluxes.
  type(coupler_2d_bc_type), intent(in)     :: gas_fields_ocn ! fields atop the ocean
      !< Structure containing ocean surface variables that are used in the calculation of the
      !! atmosphere-ocean tracer fluxes.
  type(coupler_2d_bc_type), intent(inout)  :: gas_fluxes ! fluxes between the atm and ocean
      !< Structure containing the gas fluxes between the atmosphere and the ocean and
      !! parameters related to the calculation of these fluxes.
  real, intent(in)                         :: ice_fraction(isc:iec,jsc:jec) !< sea ice fraction
  integer, intent(in)                      :: isc !< The start i-index of cell centers within
                                                  !! the computational domain
  integer, intent(in)                      :: iec !< The end i-index of cell centers within the
                                                  !! computational domain
  integer, intent(in)                      :: jsc !< The start j-index of cell centers within
                                                  !! the computational domain
  integer, intent(in)                      :: jec !< The end j-index of cell centers within the
                                                  !! computational domain
  real, intent(in), optional               :: tsurf(isc:iec,jsc:jec) !< surface temperature
  real, intent(in), optional               :: ustar(isc:iec,jsc:jec) !< friction velocity, not
                                                                     !! used
  real, intent(in), optional               :: cd_m (isc:iec,jsc:jec) !< drag coefficient, not
                                                                     !! used

  character(len=*), parameter   :: sub_name = 'atmos_ocean_fluxes_calc'
  character(len=*), parameter   :: error_header =&
      & '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  integer                                 :: n
  integer                                 :: i
  integer                                 :: j
  real, dimension(:,:), allocatable       :: kw
  real, dimension(:,:), allocatable       :: cair
  character(len=128)                      :: error_string

  real, parameter                         :: permeg=1.0e-6

  ! Return if no fluxes to be calculated
  if (gas_fluxes%num_bcs .le. 0) return

  if (.not. associated(gas_fluxes%bc)) then
    if (gas_fluxes%num_bcs .ne. 0) then
      call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
    else
      return
    endif
  endif

  do n = 1, gas_fluxes%num_bcs
    ! only do calculations if the flux has not been overridden
    if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then
      if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then
        if (.not. allocated(kw)) then
          allocate( kw(isc:iec,jsc:jec) )
          allocate ( cair(isc:iec,jsc:jec) )
        elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
          call mpp_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                  & (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) * &
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
              cair(i,j) = &
                  gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & sqrt(660. / (gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j) + epsln)) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & sqrt(660. / (gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j) + epsln)) *&
                  & gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j)) / &
                (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
          if (.not. present(tsurf)) then
            call mpp_error(FATAL, trim(error_header) // ' Implementation ' //&
                trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                ' requires input tsurf')
          endif
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) = &
                  & (1 - ice_fraction(i,j)) * gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) /&
                  & (770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) *&
                  & 101325./(rdgas*wtmair*1e-3*tsurf(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln))
              !alpha: mol/m3/atm
              cair(i,j) = &
                  gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) * &
                  gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
              cair(i,j) = max(cair(i,j),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                  & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then
          if (.not. present(tsurf)) then
            call mpp_error(FATAL, trim(error_header) // ' Implementation ' //&
                trim(gas_fluxes%bc(n)%implementation) // ' for ' // trim(gas_fluxes%bc(n)%name) //&
                ' requires input tsurf')
          endif
          !f1p: not sure how to pass salinity. For now, just force at 35.
          do j = jsc,jec
            do i = isc,iec
              !calc_kw(tk,p,u10,h,vb,mw,sc_w,ustar,cd_m)
              gas_fluxes%bc(n)%field(ind_kw)%values(i,j) =&
                  & (1 - ice_fraction(i,j)) * calc_kw(tsurf(i,j),&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j),&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j),&
                  & 101325./(rdgas*wtmair*1e-3*tsurf(i,j)*max(gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j),epsln)),&
                  & gas_fluxes%bc(n)%param(2),&
                  & gas_fluxes%bc(n)%param(1),&
                  & gas_fields_ocn%bc(n)%field(ind_sc_no)%values(i,j))
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * 9.86923e-6
              cair(i,j) = max(cair(i,j),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j))
              gas_fluxes%bc(n)%field(ind_flux0)%values(i,j) =&
                  & gas_fluxes%bc(n)%field(ind_kw)%values(i,j) *&
                  & max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i,j) =&
                  & (max(gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j),0.) - cair(i,j)) /&
                  & (gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) * permeg + epsln)
            enddo
          enddo
        else
          call mpp_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then
        if (.not. allocated(kw)) then
          allocate( kw(isc:iec,jsc:jec) )
          allocate ( cair(isc:iec,jsc:jec) )
        elseif ((size(kw(:,:), dim=1) .ne. iec-isc+1) .or. (size(kw(:,:), dim=2) .ne. jec-jsc+1)) then
          call mpp_error(FATAL, trim(error_header) // ' Sizes of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & gas_fields_atm%bc(n)%field(ind_u10)%values(i,j)**2
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
          do j = jsc,jec
            do i = isc,iec
              kw(i,j) = (1 - ice_fraction(i,j)) * gas_fluxes%bc(n)%param(1) *&
                  & max(0.0, gas_fields_atm%bc(n)%field(ind_u10)%values(i,j) - gas_fluxes%bc(n)%param(2))
              cair(i,j) =&
                  & gas_fields_ocn%bc(n)%field(ind_alpha)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_pcair)%values(i,j) *&
                  & gas_fields_atm%bc(n)%field(ind_psurf)%values(i,j) * gas_fluxes%bc(n)%param(3)
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = kw(i,j) *&
                  & (gas_fields_ocn%bc(n)%field(ind_csurf)%values(i,j) - cair(i,j))
            enddo
          enddo
        else
          call mpp_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call mpp_error(FATAL, 'Bad parameter (' // trim(error_string) //&
              & ') for air_sea_deposition for ' // trim(gas_fluxes%bc(n)%name))
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  gas_fields_atm%bc(n)%field(ind_deposition)%values(i,j) / gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        else
          call mpp_error(FATAL, 'Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then
        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call mpp_error(FATAL, ' Bad parameter (' // trim(error_string) //&
              & ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'river') then
          do j = jsc,jec
            do i = isc,iec
              gas_fluxes%bc(n)%field(ind_flux)%values(i,j) = (1 - ice_fraction(i,j)) *&
                  & gas_fields_atm%bc(n)%field(ind_runoff)%values(i,j) /&
                  & gas_fluxes%bc(n)%param(1)
            enddo
          enddo
        else
          call mpp_error(FATAL, ' Unknown implementation (' //&
              & trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      else
        call mpp_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) //&
            & ') for ' // trim(gas_fluxes%bc(n)%name))
      endif
    endif
  enddo

  if (allocated(kw)) then
    deallocate(kw)
    deallocate(cair)
  endif
end subroutine  atmos_ocean_fluxes_calc

!> Calculate \f$k_w\f$
!!
!! Taken from Johnson, Ocean Science, 2010. (http://doi.org/10.5194/os-6-913-2010)
!!
!! Uses equations defined in Liss[1974],
!! \f[
!!  F = K_g(c_g - H C_l) = K_l(c_g/H - C_l)
!! \f]
!! where \f$c_g\f$ and \f$C_l\f$ are the bulk gas and liquid concentrations, \f$H\f$
!! is the Henry's law constant (\f$H = c_{sg}/C_{sl}\f$, where \f$c_{sg}\f$ is the
!! equilibrium concentration in gas phase (\f$g/cm^3\f$ of air) and \f$C_{sl}\f$ is the
!! equilibrium concentration of unionised dissolved gas in liquid phase (\f$g/cm^3\f$
!! of water)),
!! \f[
!!    1/K_g = 1/k_g + H/k_l
!! \f]
!! and
!! \f[
!!    1/K_l = 1/k_l + 1/{Hk_g}
!! \f]
!! where \f$k_g\f$ and \f$k_l\f$ are the exchange constants for the gas and liquid
!! phases, respectively.
real function calc_kw(tk, p, u10, h, vb, mw, sc_w, ustar, cd_m)
  real, intent(in) :: tk !< temperature at surface in kelvin
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
  real, intent(in) :: h !< Henry's law constant (\f$H=c_sg/C_sl\f$) (unitless)
  real, intent(in) :: vb !< Molar volume
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: sc_w
  real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                      !! ustar = \f$u_{10} \sqrt{C_D}\f$.
  real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                      !! ustar is provided.
                                      !! If ustar is not provided,
                                      !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

  real :: ra,rl,tc

  tc = tk-273.15
  ra = 1./max(h*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
  rl = 1./max(calc_kl(tc,u10,sc_w),epsln)
  calc_kw = 1./max(ra+rl,epsln)
end function calc_kw

!> Calculate \f$k_a\f$
!!
!! See calc_kw
real function calc_ka(t, p, mw, vb, u10, ustar, cd_m)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< molar volume
  real, intent(in) :: u10 !< wind speed at 10m above the surface in m/s
  real, intent(in), optional :: ustar !< Friction velocity (m/s).  If not provided,
                                      !! ustar = \f$u_{10} \sqrt{C_D}\f$.
  real, intent(in), optional :: cd_m !< Drag coefficient (\f$C_D\f$).  Used only if
                                      !! ustar is provided.
                                      !! If ustar is not provided,
                                      !! cd_m = \f$6.1 \times 10^{-4} + 0.63 \times 10^{-4} *u_10\f$

  real             :: sc
  real             :: ustar_t, cd_m_t

  if (.not. present(ustar)) then
    !drag coefficient
    cd_m_t = 6.1e-4 +0.63e-4*u10
    !friction velocity
    ustar_t = u10*sqrt(cd_m_t)
  else
    cd_m_t = cd_m
    ustar_t = ustar
  end if
  sc = schmidt_g(t,p,mw,vb)
  calc_ka = 1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))
end function calc_ka

!> Calculate \f$k_l\f$
!!
!! See calc_kw, and Nightingale, Global Biogeochemical Cycles, 2000
!! (https://doi.org/10.1029/1999GB900091)
real function calc_kl(t, v, sc)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: v !< wind speed at surface in m/s
  real, intent(in) :: sc

  calc_kl = (((0.222*v**2)+0.333*v)*(max(sc,epsln)/600.)**(-0.5))/(100.*3600.)
end function calc_kl

!> Schmidt number of the gas in air
real function schmidt_g(t, p, mw, vb)
  real, intent(in) :: t !< temperature at surface in C
  real, intent(in) :: p !< pressure at surface in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< molar volume

  real :: d,v

  d = d_air(t,p,mw,vb)
  v = v_air(t)
  schmidt_g = v / d
end function schmidt_g

!> From Fuller, Industrial & Engineering Chemistry (https://doi.org/10.1021/ie50677a007)
real function d_air(t, p, mw, vb)
  real, intent(in) :: t  !< temperature in c
  real, intent(in) :: p  !< pressure in pa
  real, intent(in) :: mw !< molecular weight (g/mol)
  real, intent(in) :: vb !< diffusion coefficient (\f$cm3/mol\f$)

  real, parameter :: ma = 28.97d0 !< molecular weight air in g/mol
  real, parameter :: va = 20.1d0  !< diffusion volume for air (\f$cm^3/mol\f$)

  real            :: pa

  ! convert p to atm
  pa = 9.8692d-6*p
  d_air = 1d-3 *&
      & (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
  ! d_air is in cm2/s convert to m2/s
  d_air = d_air * 1d-4
end function d_air

!> kinematic viscosity in air
real function p_air(t)
  real, intent(in) :: t

  real, parameter :: sd_0 = 1.293393662d0,&
      & sd_1 = -5.538444326d-3,&
      & sd_2 = 3.860201577d-5,&
      & sd_3 = -5.2536065d-7
  p_air = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
end function p_air

!> Kinematic viscosity in air (\f$m^2/s\f$
real function v_air(t)
  real, intent(in) :: t !< temperature in C
  v_air = n_air(t)/p_air(t)
end function v_air

!> dynamic viscosity in air
real function n_air(t)
  real, intent(in) :: t !< temperature in C

  real, parameter :: sv_0 = 1.715747771d-5,&
      & sv_1 = 4.722402075d-8,&
      & sv_2 = -3.663027156d-10,&
      & sv_3 = 1.873236686d-12,&
      & sv_4 = -8.050218737d-14
  ! in n.s/m^2 (pa.s)
  n_air = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
end function n_air

end module atmos_ocean_fluxes_calc_mod
