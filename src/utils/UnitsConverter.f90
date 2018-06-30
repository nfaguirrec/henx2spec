!!**********************************************************************************
!!    Universidad Nacional de Colombia                                             !
!!    Grupo de Química Teórica                                                     !
!!    http://www.gqt-un.unal.edu.co/                                               !
!!                                                                                 !
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2007-2012) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
!!                nfaguirrec@iff.csic.es                                           !
!!                                                                                 !
!!    This program is free software; you can redistribute it and/or modify         !
!!    it under the terms of the GNU General Public License as published by         !
!!    the Free Software Foundation; either version 2 of the License, or            !
!!    (at your option) any later version.                                          !
!!                                                                                 !
!!    this program is distributed in the hope that it will be useful,              !
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of               !
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
!!    GNU General Public License for more details.                                 !
!!                                                                                 !
!!    You should have received a copy of the GNU General Public License            !
!!    along with thisPtr program. If not, write to the Free Software Foundation,   !
!!    Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              !
!!                                                                                 !
!!**********************************************************************************

! module Unit_
! 	implicit none
! 	private
! 	
! 	type, public :: Unit
! 		character(*) :: name
! 		character(*) :: symbol
! 		character(*) :: category
! 		real(8) :: value
! 	end type Unit
! 	
! end module Unit_

! Esto debe ser una clase singleton
! y entre lo que contendrá será
! el sistema de unidades con el que
! actualmente se está trabajando
! module SystemOfUnits_
! 	implicit none
! 	private
! 	
! 	type, public :: SystemOfUnits
! 		real(8) :: mass
! 	end type SystemOfUnits
! 	
! end module SystemOfUnits_

! module CM1Units_
! 	implicit none
! 	private
! 	
! 	!!
! 	! Mass   => amu
! 	! Lenght => Å
! 	! Time   => 9.14294598 10⁻¹¹ s
! 	! Energy => cm⁻¹
! 	!!
! 	real(8), public, parameter :: amu = 1.0_8
! 	real(8), public, parameter :: angs = 1.0_8
! 	real(8), public, parameter :: t_cm1 = 1.0_8
! 	real(8), public, parameter :: cm1 = 1.0_8
! 	real(8), public, parameter :: g = 1.0_8/1.660538782e-24
! 	real(8), public, parameter :: kg = 1.0_8/1.660538782e-27
! 	real(8), public, parameter :: m_au = 1.0_8/1822.88853_8
! 	real(8), public, parameter :: e_au = 1.0_8/4.55633538e-6
! 	real(8), public, parameter :: Eh = e_au
! 	real(8), public, parameter :: Ha = e_au
! end module CM1Units_

!!! Será renombrada UnitsManager
module UnitsConverter_
	implicit none
	private
	
	!**
	! Mass   => amu
	! Lenght => Å
	! Time   => 9.14294598 10⁻¹¹ s
	! Energy => cm⁻¹
	!**
	real(8), public, parameter :: amu = 1822.88853_8
	real(8), public, parameter :: angs = 1.0_8/0.52917726_8
	real(8), public, parameter :: t_cm1 = 1.0_8
	real(8), public, parameter :: cm1 = 4.55633538e-6
	real(8), public, parameter :: bohr = 1.0_8
! 	real(8), public, parameter :: g = 1.0_8/1.660538782e-24
! 	real(8), public, parameter :: kg = 1.0_8/1.660538782e-27
	real(8), public, parameter :: m_au = 1.0_8
	real(8), public, parameter :: e_au = 1.0_8
	real(8), public, parameter :: Eh = e_au
	real(8), public, parameter :: Ha = e_au
	
	type, public :: Unit
		character(:), allocatable :: name
		character(:), allocatable :: symbol
		character(:), allocatable :: category
		real(8) :: value
	end type Unit
	
	type, private :: UnitsConverter
		type(Unit) :: globalUnits(4)
		
		contains
			procedure init
			procedure destroy
			procedure str
			procedure show
! 			procedure setGlobalUnits
	end type UnitsConverter
	
	type(UnitsConverter), target, public :: IntegralManager_instance
	
	contains
	
	!**
	! @brief Constructor
	!**
	subroutine init( this )
		implicit none
		class(UnitsConverter) :: this 
		
	end subroutine init
	
	!**
	! @brief Copy constructor
	!**
	subroutine copy( this, other )
		implicit none
		class(UnitsConverter) :: this
		type(UnitsConverter), intent(in) :: other

! 		this.val = other.val
	end subroutine copy
	
	!**
	! @brief Destructor
	!**
	subroutine destroy( this )
		implicit none
		class(UnitsConverter) :: this
		
	end subroutine destroy
	
	!**
	! @brief Convert to string
	!**
	function str( this ) result( output )
		implicit none
		class(UnitsConverter) :: this 
		character(len=200) :: output
		
		integer :: fmt
		character(len=200) :: strBuffer
		
		output = ""
		
		output = trim(output)//"<UnitsConverter:"
		
! 		output = trim(output)//"min="
! 		fmt = int(log10(this.min+1.0))+1
! 		write(strBuffer, "(f<fmt+7>.6)") this.min
! 		output = trim(output)//trim(strBuffer)
! 		
! 		output = trim(output)//",size="
! 		fmt = int(log10(float(this.size+1)))+1
! 		write(strBuffer, "(i<fmt>)") this.size
! 		output = trim(output)//trim(strBuffer)
		
		output = trim(output)//">"
	end function str
	
	!**
	! @brief Show 
	!**
	subroutine show( this, unit )
		implicit none
		class(UnitsConverter) :: this
		integer, optional, intent(in) :: unit
		
		integer :: effunit
		
		if( present(unit) ) then
			effunit = unit
		else
			effunit = 6
		end if
		
		write(effunit,"(a)") trim(str(this))
	end subroutine show
	
end module UnitsConverter_
