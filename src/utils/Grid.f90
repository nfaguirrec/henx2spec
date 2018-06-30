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

!**
! This class represents a one-dimensional grid 
!**
module Grid_
	implicit none
	private
	
	public :: &
		Grid_test
	
	type, public :: Grid
		real(8) :: min
		real(8) :: max
		real(8) :: stepSize
		integer :: size
		logical :: isEquallyspaced
		real(8), allocatable :: data(:)
		
		contains
			procedure :: init
			procedure :: fromArray
			procedure :: copy
			procedure :: destroy
			procedure :: str
			procedure :: show
			procedure :: setUnits
	end type Grid
	
	contains
	
	subroutine init( this, min, max, size, stepSize )
		class(Grid) :: this 
		real(8) :: min
		real(8) :: max
		integer, optional :: size
		real(8), optional :: stepSize
		
		integer :: i
		
		this.min = min
		this.max = max
		this.isEquallyspaced = .true.
		
		if( present(stepSize) ) then
			this.stepSize = stepSize
			this.size = int( abs(max-min)/stepSize )+1
		else if( present(size) ) then
			this.size = size
			this.stepSize = abs(max-min)/float(size-1)
		end if
		
		allocate( this.data(this.size) )
		do i=1,this.size
			this.data( i ) = this.min+float(i-1)*this.stepSize
		end do
	end subroutine init
	
	!**
	! @brief Constructor
	!**
	subroutine fromArray( this, array )
		class(Grid) :: this 
		real(8), intent(in) :: array(:)
		
		integer :: i
		
		this.min = minval(array)
		this.max = maxval(array)
		this.size = size(array)
		this.stepSize = 0.0_8
		this.isEquallyspaced = .false.
		
		allocate( this.data(this.size) )
		do i=1,this.size
			this.data( i ) = array( i )
		end do
	end subroutine fromArray
	
	!**
	! @brief Copy constructor
	!**
	subroutine copy( this, other )
		class(Grid) :: this
		type(Grid) :: other

		this.min = other.min
		this.max = other.max
		this.stepSize = other.stepSize
		this.size = other.size
		this.isEquallyspaced = other.isEquallyspaced
		
		allocate( this.data(this.size) )
		this.data = other.data
	end subroutine copy
	
	subroutine destroy( this )
		class(Grid) :: this
		
		this.min = 0.0_8
		this.max = 0.0_8
		this.size = 0
		this.stepSize = 0.0_8
		this.isEquallyspaced = .false.
		
		if( allocated( this.data ) ) then
			deallocate( this.data )
		end if
	end subroutine destroy
	
	function str( this ) result( output )
		class(Grid) :: this 
		character(len=200) :: output
		
		integer :: fmt
		character(len=200) :: strBuffer
		
		output = ""
		
		output = trim(output)//"<Grid:"
		
		output = trim(output)//"min="
		fmt = int(log10(this.min+1.0))+1
		write(strBuffer, "(f<fmt+7>.6)") this.min
		output = trim(output)//trim(strBuffer)
		
		output = trim(output)//",max="
		fmt = int(log10(this.max+1.0))+1
		write(strBuffer, "(f<fmt+7>.6)") this.max
		output = trim(output)//trim(strBuffer)
		
		if ( this.isEquallyspaced ) then
			output = trim(output)//",stepSize="
			fmt = int(log10(this.stepSize+1.0))+1
			write(strBuffer, "(f<fmt+7>.6)") this.stepSize
			output = trim(output)//trim(strBuffer)
		end if
		
		output = trim(output)//",size="
		fmt = int(log10(float(this.size+1)))+1
		write(strBuffer, "(i<fmt>)") this.size
		output = trim(output)//trim(strBuffer)
		
		output = trim(output)//">"
	end function str
	
	subroutine show( this, unit )
		class(Grid) :: this
		integer, optional, intent(in) :: unit
		
		integer :: effunit
		
		if( present(unit) ) then
			effunit = unit
		else
			effunit = 6
		end if
		
		write(effunit,"(a)") trim(this.str())
	end subroutine show
	
	subroutine setUnits( this, unit )
		class(Grid) :: this
		real(8), intent(in) :: unit
		
		this.min = this.min*unit
		this.max = this.max*unit
		this.stepSize = this.stepSize*unit
		this.data = this.data*unit
	end subroutine setUnits
	
	subroutine Grid_test()
		real(8) :: rMin, rMax
		integer :: gridSize
		type(Grid) :: rGrid
		type(Grid) :: rGrid2
		
		integer :: i
			
		rMin=0.0_8
		rMax=5.0_8
		gridSize=10
		
		call rGrid.init( rMin, rMax, gridSize )
		call rGrid.show()
		
		do i=1,rGrid.size
			write(*,"(i5,f10.5)") i, rGrid.data(i)
		end do
		
		write(*,*) "---"
		write(*,*) "Testing copy constructor"
		write(*,*) "---"
		call rGrid2.copy( rGrid )
		call rGrid2.show()
		do i=1,rGrid2.size
			write(*,"(i5,f10.5)") i, rGrid2.data(i)
		end do
		
		call rGrid.destroy()
		call rGrid2.destroy()
	end subroutine Grid_test
	
end module Grid_
