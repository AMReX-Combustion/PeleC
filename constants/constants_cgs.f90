! Fundamental constants taken from NIST's 2010 CODATA recommended values

module fundamental_constants_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module, only: M_PI
  
  implicit none

  ! newton's gravitational constant
  real(amrex_real), parameter :: Gconst = 6.67428e-8_amrex_real      ! cm^3/g/s^2
! new value; if uncommented initial models will need to be re-HSE'ed
!  real(amrex_real), parameter :: Gconst = 6.67384e-8_amrex_real      ! cm^3/g/s^2

  ! boltzmann's constant
  real(amrex_real), parameter :: k_B    = 1.3806488e-16_amrex_real   ! erg/K

  ! planck's constant over 2pi
  real(amrex_real), parameter :: hbar   = 1.054571726e-27_amrex_real ! erg s

  ! planck's constant 
  real(amrex_real), parameter :: hplanck = 6.62606957e-27_amrex_real ! erg s

  ! avogradro's Number
  real(amrex_real), parameter :: n_A    = 6.02214129e23_amrex_real   ! mol^-1

  ! convert eV to erg
  real(amrex_real), parameter :: ev2erg = 1.602176487e-12_amrex_real

  ! convert MeV to eV
  real(amrex_real), parameter :: MeV2eV = 1.0e6_amrex_real

  ! mass of proton
  real(amrex_real), parameter :: m_p     = 1.672621777e-24_amrex_real ! g

  ! mass of neutron
  real(amrex_real), parameter :: m_n      = 1.674927351e-24_amrex_real ! g

  ! mass of electron
  real(amrex_real), parameter :: m_e     = 9.10938291e-28_amrex_real  ! g

  ! speed of light in vacuum
  real(amrex_real), parameter :: c_light = 2.99792458e10_amrex_real   ! cm/s

  ! electron charge
  ! NIST: q_e = 1.602176565e-19 C
  !
  ! C is the SI unit Coulomb; in cgs we have the definition:
  !     1 C = 0.1 * |c_light| * 1 statC
  ! where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
  ! and |c_light| is the speed of light in cgs (but without units)
  real(amrex_real), parameter :: q_e     = 4.80320451e-10_amrex_real  ! erg^1/2 cm^1/2

  ! stefan-boltzmann constant
  real(amrex_real), parameter :: sigma_SB = 5.670373e-5_amrex_real    ! erg/s/cm^2/K^4

  ! radiation constant
  real(amrex_real), parameter :: a_rad = 4.0_amrex_real*sigma_SB/c_light

  ! Number of centimeters in a parsec and an AU.
  ! Note that since the length of an AU is defined exactly
  ! by IAU convention, the length of a parsec is also
  ! defined exactly as (6.48e5 / pi) AU.
  real(amrex_real), parameter :: AU = 1.49597871e13_amrex_real            ! cm
  real(amrex_real), parameter :: parsec = 3.085677587679311e18_amrex_real ! cm

  ! Hubble constant (in s^{-1}, converted from 100 (km/s)/Mpc by dividing by 3.08568025e19km/Mpc)
  real(amrex_real), parameter :: Hubble_const = 32.407764868e-19_amrex_real

  ! solar mass (from http://asa.usno.navy.mil/SecK/Constants.html)
  real(amrex_real), parameter :: M_solar = 1.9884e33_amrex_real

end module fundamental_constants_module
