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

module NFunction_
	use Grid_
	implicit none
	private
	
	public :: &
		NFunction_test
	
	type, public:: NFunction
		integer :: size
		type(Grid) :: xGrid
		real(8), allocatable :: yArray(:)
		
		contains
			procedure :: fromArray
			procedure :: fromFunction
			procedure :: copy
			procedure :: destroy
			
			procedure :: addition
			procedure :: subtraction
			procedure :: multiplication
			procedure :: division
			procedure :: exponentiation
			
			procedure :: additionFC
			procedure :: subtractionFC
			procedure :: multiplicationFC
			procedure :: divisionFC
			procedure :: exponentiationFC
			
			generic :: assignment(=) => copy
			generic :: operator(+) => addition, additionFC
			generic :: operator(-) => subtraction, subtractionFC
			generic :: operator(*) => multiplication, multiplicationFC
			generic :: operator(/) => division, divisionFC
			generic :: operator(**) => exponentiation, exponentiationFC
			
			procedure :: str
			procedure :: show
			procedure :: setUnits
	end type NFunction
	
	interface
		function prototypeFunction( x ) result( output )
			real(8), intent(in) :: x
			real(8) :: output
		end function prototypeFunction
	end interface
	
	contains
	
	!**
	! @brief Constructor
	!**
	subroutine fromArray( this, xGrid, yArray, units )
		class(NFunction) :: this 
		type(Grid), intent(in) :: xGrid
		real(8), intent(in) :: yArray(:)
		real(8), optional, intent(in) :: units(2)
		
		integer :: i
		real(8) :: unitsEff(2)
		
		if( allocated(this.yArray) ) then
			call this.destroy()
		end if
		
		if( present(units) ) then
			unitsEff = units
		else
			unitsEff = [1.0_8, 1.0_8]
		end if
		
		this.size = xGrid.size
		call this.xGrid.copy( xGrid )
		call this.xGrid.setUnits( unitsEff(1) )
		
		allocate( this.yArray(xGrid.size) )
		this.yArray = yArray*unitsEff(2)
	end subroutine fromArray
	
	!**
	! @brief Constructor
	!**
	subroutine fromFunction( this, xGrid, func, units )
		class(NFunction) :: this 
		type(Grid), intent(in) :: xGrid
		procedure(prototypeFunction) :: func
		real(8), optional, intent(in) :: units(2)
		
		integer :: i
		real(8) :: unitsEff(2)
		
		if( allocated(this.yArray) ) then
			call this.destroy()
		end if
		
		if( present(units) ) then
			unitsEff = units
		else
			unitsEff = [1.0_8, 1.0_8]
		end if
		
		this.size = xGrid.size
		call this.xGrid.copy( xGrid )
		call this.xGrid.setUnits( unitsEff(1) )
		
		allocate( this.yArray(this.size) )
		
		do i=1,this.size
			this.yArray(i) = func( xGrid.data(i) )*unitsEff(2)
		end do
	end subroutine fromFunction
	
	!**
	! @brief Copy constructor
	!**
	subroutine copy( this, other )
		class(NFunction), intent(inout) :: this
		type(NFunction), intent(in) :: other
		
		if( allocated(this.yArray) ) then
			call this.destroy()
		end if
		
		this.size = other.size
		call this.xGrid.copy( other.xGrid )
		allocate( this.yArray(this.size) )
		this.yArray = other.yArray
	end subroutine copy
	
	!**
	! @brief Destructor
	!**
	subroutine destroy( this )
		class(NFunction) :: this
		
		this.size = -1
		call this.xGrid.destroy()
		deallocate( this.yArray )
	end subroutine destroy
	
	!**
	! @brief
	!**
	function addition( this, other ) result( output )
		class(NFunction), intent(in) :: this
		class(NFunction), intent(in) :: other
		type(NFunction) :: output
		
		if( this.size /= other.size ) then
			write(*,*) "## ERROR ## the Numerical Functions have not the same size"
			stop
		end if
		
		call output.copy( this )
		output.yArray = this.yArray + other.yArray
	end function addition
	
	!**
	! @brief
	!**
	function additionFC( this, constant ) result( output )
		class(NFunction), intent(in) :: this
		real(8), intent(in) :: constant
		type(NFunction) :: output
		
		call output.copy( this )
		output.yArray = this.yArray+constant
	end function additionFC
	
	!**
	! @brief
	!**
	function subtraction( this, other ) result( output )
		class(NFunction), intent(in) :: this
		type(NFunction), intent(in) :: other
		type(NFunction) :: output
		
		if( this.size /= other.size ) then
			write(*,*) "## ERROR ## the Numerical Functions have not the same size"
			stop
		end if
		
		call output.copy( this )
		output.yArray = this.yArray - other.yArray
	end function subtraction
	
	!**
	! @brief
	!**
	function subtractionFC( this, constant ) result( output )
		class(NFunction), intent(in) :: this
		real(8), intent(in) :: constant
		type(NFunction) :: output
		
		call output.copy( this )
		output.yArray = this.yArray-constant
	end function subtractionFC
	
	!**
	! @brief
	!**
	function multiplication( this, other ) result( output )
		class(NFunction), intent(in) :: this
		type(NFunction), intent(in) :: other
		type(NFunction) :: output
		
		if( this.size /= other.size ) then
			write(*,*) "## ERROR ## the Numerical Functions have not the same size"
			stop
		end if
		
		call output.copy( this )
		output.yArray = this.yArray*other.yArray
	end function multiplication
	
	!**
	! @brief
	!**
	function multiplicationFC( this, constant ) result( output )
		class(NFunction), intent(in) :: this
		real(8), intent(in) :: constant
		type(NFunction) :: output
		
		call output.copy( this )
		output.yArray = this.yArray*constant
	end function multiplicationFC
	
	!**
	! @brief
	!**
	function division( this, other ) result( output )
		class(NFunction), intent(in) :: this
		type(NFunction), intent(in) :: other
		type(NFunction) :: output
		
		if( this.size /= other.size ) then
			write(*,*) "## ERROR ## the Numerical Functions have not the same size"
			stop
		end if
		
		call output.copy( this )
		output.yArray = this.yArray/other.yArray
	end function division
	
	!**
	! @brief
	!**
	function divisionFC( this, constant ) result( output )
		class(NFunction), intent(in) :: this
		real(8), intent(in) :: constant
		type(NFunction) :: output
		
		call output.copy( this )
		output.yArray = this.yArray/constant
	end function divisionFC
	
	!**
	! @brief
	!**
	function exponentiation( this, other ) result( output )
		class(NFunction), intent(in) :: this
		type(NFunction), intent(in) :: other
		type(NFunction) :: output
		
		if( this.size /= other.size ) then
			write(*,*) "## ERROR ## the Numerical Functions have not the same size"
			stop
		end if
		
		call output.copy( this )
		output.yArray = this.yArray**other.yArray
	end function exponentiation
	
	!**
	! @brief
	!**
	function exponentiationFC( this, constant ) result( output )
		class(NFunction), intent(in) :: this
		real(8), intent(in) :: constant
		type(NFunction) :: output
		
		call output.copy( this )
		output.yArray = this.yArray**constant
	end function exponentiationFC
	
	!**
	! @brief String representation of the object
	!**
	function str( this ) result( output )
		class(NFunction) :: this 
		character(len=200) :: output
		
		integer :: fmt
		character(len=200) :: strBuffer
		
		output = ""
		
		output = trim(output)//"<NFunction:"
		
		output = trim(output)//this.xGrid.str()
		
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
		
		output = trim(output)//">"
	end function str
	
	!**
	! Write the string representation of the object
	! in a selected unit
	!**
	subroutine show( this, unit )
		class(NFunction) :: this
		integer, optional, intent(in) :: unit
		
		integer :: unitEff
		
		if( present(unit) ) then
			unitEff = unit
		else
			unitEff = 6
		end if
		
		write(unitEff,"(a)") trim(str(this))
	end subroutine show
	
	subroutine setUnits( this, units )
		class(NFunction) :: this
		real(8), intent(in) :: units(2)
		
		call this.xGrid.setUnits( units(1) )
		this.yArray = this.yArray*units(2)
	end subroutine setUnits
	
	!**
	! This is neccesary only for NFunction_test()
	!**
	function funcTest( x ) result( output )
		real(8), intent(in) :: x
		real(8) :: output
		
		output = 5.0_8*( exp(2.0_8*(2.0_8-x))-2.0_8*exp(2.0_8-x) )
	end function funcTest
	
	!**
	! This is neccesary only for NFunction_test()
	!**
	function funcTest2( x ) result( output )
		real(8), intent(in) :: x
		real(8) :: output
		
		output = sin(x)
	end function funcTest2
	
	!**
	! Test method of this class
	!**
	subroutine NFunction_test()
		type(Grid) :: xGrid
		type(NFunction) :: nFunc
		type(NFunction) :: nFunc2
		type(NFunction) :: nFunc3
		real(8) :: value
		real(8), allocatable :: data(:)
		integer :: i
		
		call xGrid.init( 1.0_8, 10.0_8, 100 )
		call xGrid.show()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Test from function
		write(*,*) "---"
		write(*,*) "Testing from function"
		write(*,*) "---"
		
		call nFunc.fromFunction( xGrid, func=funcTest )
		call nFunc.show()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Test for copy constructor
		write(*,*) "---"
		write(*,*) "Testing copy constructor"
		write(*,*) "---"
		
		call nFunc2.copy( nFunc )
		call nFunc2.show()
		
		call nFunc.destroy()
		call nFunc2.destroy()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Test from array
		write(*,*) "---"
		write(*,*) "Testing from array"
		write(*,*) "---"
		
		allocate( data(xGrid.size) )
		do i=1,xGrid.size
			data(i) = funcTest( xGrid.data(i) )
		end do
		
		call nFunc.fromArray( xGrid, yArray=data )
		call nFunc.show()
		
		call xGrid.destroy()
		call nFunc.destroy()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Test from IFStream
! 		write(*,*) "---"
! 		write(*,*) "Testing from IFStream"
! 		write(*,*) "---"
! 		
! 		call ifile.init( "morse.dat" )
! 		call nFunc.fromFStream( ifile )
! 		call nFunc.show()
! 		call ifile.destroy()
! 		call nFunc.save( "salida4" )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Test operators
		write(*,*) "---"
		write(*,*) "Testing operators"
		write(*,*) "---"
		
		call xGrid.init( 1.0_8, 10.0_8, 100 )
		call nFunc.fromFunction( xGrid, func=funcTest )
		call nFunc2.fromFunction( xGrid, func=funcTest2 )
		
		nFunc3 = nFunc+nFunc2
		
		nFunc3 = nFunc+3.0_8
		
		nFunc3 = nFunc-nFunc2
		
		nFunc3 = nFunc-3.0_8
		
		nFunc3 = nFunc*nFunc2
		
		nFunc3 = nFunc*3.0_8
		
		nFunc3 = nFunc/nFunc2
		
		nFunc3 = nFunc/3.0_8
		
		nFunc3 = nFunc**nFunc2
		
		nFunc3 = nFunc**2.0_8
		
		call xGrid.destroy()
		call nFunc.destroy()
		call nFunc2.destroy()
		call nFunc3.destroy()
	end subroutine NFunction_test
	
end module NFunction_
