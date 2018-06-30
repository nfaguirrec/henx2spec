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

module NIntegration_
	use Grid_
	use NFunction_
	implicit none
	private
	
	!**
	! @brief Public parameters
	!**
	integer, public, parameter :: SIMPSON  = 0
	integer, public, parameter :: EXTSIMPSON  = 1
	integer, public, parameter :: SIMPSON38  = 2
	integer, public, parameter :: TRAPEZOIDAL = 3
	integer, public, parameter :: FIXED_QUADRATURE = 4
	integer, public, parameter :: QUADRATURE = 5
	integer, public, parameter :: ADAPTIVE_QUADRATURE = 6
	integer, public, parameter :: BOOLE = 7
	
	public :: &
		NIntegration_test
	
	type, public :: NIntegration
		type(NFunction) :: func
		integer :: method
		
		contains
			procedure init
			procedure destroy
			procedure str
			procedure show
			procedure evaluate
	end type NIntegration
	
	contains
	
	subroutine init( this, func, method )
		class(NIntegration) :: this 
		type(NFunction) :: func
		integer, optional :: method
		
		this.func = func
		
		if( present(method) ) then
			this.method = method
		else
			this.method = SIMPSON
		end if
	end subroutine init
	
	subroutine destroy( this )
		class(NIntegration) :: this
		
! 		call this.func.destroy()
	end subroutine destroy
	
	function str( this ) result( output )
		class(NIntegration) :: this 
		character(len=200) :: output
		
		integer :: fmt
		character(len=200) :: strBuffer
		
		output = ""
		
! 		output = trim(output)//"<NIntegration:"
! 		
! 		output = trim(output)//"min="
! 		fmt = int(log10(this.min+1.0))+1
! 		write(strBuffer, "(f<fmt+7>.6)") this.min
! 		output = trim(output)//trim(strBuffer)
! 		
! 		output = trim(output)//",max="
! 		fmt = int(log10(this.max+1.0))+1
! 		write(strBuffer, "(f<fmt+7>.6)") this.max
! 		output = trim(output)//trim(strBuffer)
! 		
! 		output = trim(output)//",h="
! 		fmt = int(log10(this.h+1.0))+1
! 		write(strBuffer, "(f<fmt+7>.6)") this.h
! 		output = trim(output)//trim(strBuffer)
! 		
! 		output = trim(output)//",size="
! 		fmt = int(log10(float(this.size+1)))+1
! 		write(strBuffer, "(i<fmt>)") this.size
! 		output = trim(output)//trim(strBuffer)
! 		
! 		output = trim(output)//">"
	end function str
	
	subroutine show( this, unit )
		class(NIntegration) :: this
		integer, optional, intent(in) :: unit
		
		integer :: effunit
		
		if( present(unit) ) then
			effunit = unit
		else
			effunit = 6
		end if
		
		write(effunit,"(a)") trim(this.str())
	end subroutine show
	
	function evaluate( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
		select case( this.method )
			case( SIMPSON )
				output = simpsonRule( this )
			case( EXTSIMPSON )
				output = extendedSimpsonRule( this )
			case( SIMPSON38 )
				output = simpson38Rule( this )
			case( TRAPEZOIDAL )
				output = trapezoidalRule( this )
			case( FIXED_QUADRATURE )
				output = fixedQuadratureRule( this )
			case( QUADRATURE )
				output = quadratureRule( this )
			case( ADAPTIVE_QUADRATURE )
				output = adaptiveQuadratureRule( this )
			case( BOOLE )
				output = booleRule( this )
		end select
	end function evaluate
	
	!**
	! \int_a^b f(x)dx = \frac{1}{3}\left[ f(x_0) + 2\sum_{i=1}^{n/2-1}f(x_{2i}) + 4\sum_{i=1}^{n/2-1}f(x_{2i-1}) + f(x_n)\right]
	!**
	function simpsonRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
		real(8) :: sum
		integer :: size
		integer :: i
		
		size = this.func.xGrid.size
		
		sum = this.func.yArray(1)
		
		do i=1,size/2-1
			sum = sum + 2.0_8*this.func.yArray(2*i) + 4.0_8*this.func.yArray(2*i+1)
		end do
		
		sum = sum + this.func.yArray(size)
		
		output = sum*( this.func.xGrid.stepSize/3.0_8 )
	end function simpsonRule
	
	function extendedSimpsonRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
		real(8) :: sum
		integer :: size
		integer :: i
		
		size = this.func.xGrid.size
		
		sum = 0.0
		do i=5,size-4
			sum = sum + this.func.yArray(i)
		end do
		sum = 48.0_8*sum
		
		sum = sum + 17.0_8*this.func.yArray(1) + 59.0_8*this.func.yArray(2) + 43.0_8*this.func.yArray(3) + 49.0_8*this.func.yArray(4)
		sum = sum + 49.0_8*this.func.yArray(size-3) + 43.0_8*this.func.yArray(size-2) + 59.0_8*this.func.yArray(size-1) + 17.0_8*this.func.yArray(size)
		
		output = sum*( this.func.xGrid.stepSize/48.0_8 )
	end function extendedSimpsonRule
	
	function simpson38Rule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
	end function simpson38Rule

	function trapezoidalRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
	end function trapezoidalRule

	function fixedQuadratureRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
	end function fixedQuadratureRule

	function quadratureRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
	end function quadratureRule

	function adaptiveQuadratureRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
	end function adaptiveQuadratureRule
	
	function booleRule( this ) result( output )
		class(NIntegration) :: this
		real(8) :: output
		
		real(8) :: sum
		integer :: size
		integer :: i
		
		size = this.func.xGrid.size
		
		sum = 0.0_8
		
		do i=1,min(size,size-4),4
			sum = sum + 14.0_8*this.func.yArray(i) + 64.0_8*this.func.yArray(i+1) &
			        + 24.0_8*this.func.yArray(i+2) + 64.0_8*this.func.yArray(i+3) + 14.0_8*this.func.yArray(i+4)
		end do
		
! 		sum = this.func.yArray(1)
! 		do i=1,(size-2)/4
! 			sum = sum + 14.0_8*this.func.yArray(4*i-2) + 64.0_8*this.func.yArray(4*i-1) &
! 			        + 24.0_8*this.func.yArray(4*i) + 64.0_8*this.func.yArray(4*i+1) + 14.0_8*this.func.yArray(4*i+2)
! 		end do
		
		output = sum*( this.func.xGrid.stepSize/45.0_8 )
	end function booleRule
	
	!**
	! This is neccesary only for NFunction_test()
	!**
	function funcTest( x ) result( output )
		real(8), intent(in) :: x
		real(8) :: output
		
		output = exp(-0.44*x)*sin(x)**2.0_8
	end function funcTest
	
	subroutine NIntegration_test()
		type(Grid) :: xGrid
		type(NFunction) :: nFunc
		type(NIntegration) :: integrator
		real(8) :: exactValue
		real(8) :: value
		integer :: i
		
		call xGrid.init( 0.0_8, 50.0_8, 102 )
		call xGrid.show()
		
		call nFunc.fromFunction( xGrid, funcTest )
		call nFunc.show()
		
		exactValue = 1.083902743303941_8
		
		write(*,*) "SIMPSON"
		write(*,*) "======="
		call integrator.init( nFunc, SIMPSON )
		write(*,'(A,F30.8)') "Exact     = ", exactValue
		write(*,'(A,F30.8)') "Numerical = ", integrator.evaluate()
		write(*,'(A,F27.5,A3)') "Error(%)  = ", 100.0_8*( integrator.evaluate()-exactValue )/exactValue, "%"
		write(*,*) ""
		
		write(*,*) "BOOLE"
		write(*,*) "====="
		call integrator.init( nFunc, BOOLE )
		write(*,'(A,F30.8)') "Exact     = ", exactValue
		write(*,'(A,F30.8)') "Numerical = ", integrator.evaluate()
		write(*,'(A,F27.5,A3)') "Error(%)  = ", 100.0_8*( integrator.evaluate()-exactValue )/exactValue, "%"
		write(*,*) ""
		
		call integrator.destroy()
		call nFunc.destroy()
		call xGrid.destroy()
	end subroutine NIntegration_test
	
end module NIntegration_
