module thermochem
  use precision, only: WP
  use messager, only: die
   implicit none
   private

   ! Make public
   public :: Rcst
   public :: lookup_molarmass

   ! Thermodynamic constants
   real(WP), parameter :: Rcst = 8.314472_WP  ! [J/(mol.K)]

 contains
   ! Molar Mass
   subroutine lookup_molarmass(spec_name,mm)
    implicit none

    character(len=*), intent(in) :: spec_name
    real(WP), intent(out) :: mm

    ! Get molar mass [kg/mol]
    select case (trim(spec_name))
    case ('H')
       mm = 1.00794E-3_WP
    case ('H2')
       mm = 2.01588E-3_WP
    case ('N')
       mm = 14.0067E-3_WP
    case ('N2')
       mm = 28.0134E-3_WP
    case ('O')
       mm = 15.9994E-3_WP
    case ('O2')
       mm = 31.9988E-3_WP
    case ('CO2')
       mm = 44.0095E-3_WP
    case('air')
       mm = 28.96443E-3_WP
    case('water', 'H20')
       mm = 18.01528E-3_WP
    case('hexane')
       mm = 0.086177_WP
    case('heptane')
       mm = 0.100204_WP
    case('dodecane')
       mm = 0.170338_WP
    case('acetone')
       mm = 0.05808_WP
    case default
       call die('Unknown species')
    end select
  
    return
  end subroutine lookup_molarmass

   ! Cp
   ! Cv

 end module thermochem
 
