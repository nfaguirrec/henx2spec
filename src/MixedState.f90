!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2011-2012) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
!!                nfaguirrec@iff.csic.es                                           !
!!    (2011-2012) María Pilar de Lara-Castells                                     !
!!                delara@iff.csic.es                                               !
!!                                                                                 !
!!    This file is part of the HenX2Spec program, which concept is based on        !
!!    a F77 code written by:                                                       !
!!                                                                                 !
!!                Pablo Villarreal Herrán                                          !
!!                p.villarreal@csic.es                                             !
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

module MixedState_
	use PureState_
	implicit none
	private
	
	public :: &
		MixedState_boltzmanWeights, &
		MixedState_saveEnergies
	
	type, public :: MixedState
		real(8) :: E
		real(8) :: weight
		
		integer :: nStates
		type(PureState), allocatable :: pureStates(:)
		real(8), allocatable :: coeff(:)
		
		contains
			procedure :: init
			procedure :: destroy
			procedure :: copy
                        generic :: assignment(=) => copy
			procedure :: str
			procedure :: J
	end type MixedState
	
	contains
	
	subroutine init( this, nPureStates )
		class(MixedState) :: this
		integer, intent(in) :: nPureStates
		
		integer :: i
		character(10) :: cID
		
		cID = "X"
		
		this.nStates = nPureStates
		allocate( this.pureStates(nPureStates) )
		
		do i=1,size(this.pureStates)
			call this.pureStates(i).init( -1, -1.0_8, -1.0_8, -1.0_8, -1, -1.0_8, -1.0_8, -100.0_8, 0.0_8, cID)
		end do
		
		allocate( this.coeff(nPureStates) )
	end subroutine init
	
	subroutine destroy( this )
		class(MixedState) :: this
		
		integer :: i
		
		do i=1,size(this.pureStates)
			call this.pureStates(i).destroy()
		end do
		
		this.E = -100.0_8
		this.weight = 0.0_8
		this.nStates = 0
		
		deallocate( this.pureStates )
		deallocate( this.coeff )
	end subroutine destroy
	
	subroutine copy( this, other )
		class(MixedState), intent(inout) :: this
		type(MixedState), intent(in) :: other
		
		integer :: i
		
		if( allocated(this.pureStates) ) call this.destroy()
		
		call this.init( other.nStates )
		
		do i=1,other.nStates
			this.pureStates(i) = other.pureStates(i)
		end do
		
		this.E = other.E
		this.weight = other.weight
		this.nStates = other.nStates
		
		this.coeff = other.coeff
	end subroutine copy
	
	function str( this, t ) result( output )
		class(MixedState) :: this
		integer, intent(in), optional :: t
		character(len=2000) :: output
		
		integer :: tEff
		character(len=200) :: strBuffer
		integer :: i
		
		if ( present(t) ) then
			tEff = t
		else
			tEff = BRA
		end if
		
		output = ""
		do i=1,this.nStates
			if( this.coeff(i) > 0.0_8 ) then
				if( i /= 1 ) then
					write(strBuffer, "(F6.4)") this.coeff(i)
					strBuffer=" +"//trim(adjustl(strBuffer))
				else
					write(strBuffer, "(F6.4)") this.coeff(i)
					strBuffer="  "//trim(adjustl(strBuffer))
				end if
			else
				write(strBuffer, "(F6.4)") abs(this.coeff(i))
				strBuffer=" -"//trim(adjustl(strBuffer))
			end if
			output = trim(output)//trim(strBuffer)//" "//trim(this.pureStates(i).str( tEff ))
		end do
		
	end function str
	
	function J( this ) result( output )
		class(MixedState) :: this
		real(8) :: output
		
		output = this.pureStates(1).J
	end function J
	
	subroutine MixedState_boltzmanWeights( states, T )
		use UnitsConverter_
		type(MixedState), intent(inout) :: states(:)
		real(8), intent(in) :: T ! K
		
		integer :: nstates
		real(8) :: lowest
		real(8) :: kBoltzman, kT
		
		nstates = size(states)
		
		kBoltzman = 3.16682964d-6 ! 1K = 3.16..d-6 au
		
		kT = kBoltzman*T ! a.u.
		
		lowest = minval( states.E )
		
		states.weight = exp( -(states.E-lowest)*cm1/kT )
		states.weight = states.weight/sum(states.weight)
	end subroutine MixedState_boltzmanWeights
	
	subroutine MixedState_saveEnergies( states, oFile )
		type(MixedState), intent(in) :: states(:)
		character(*), intent(in) :: oFile
		
		integer :: i, j
		integer :: nstates
		
		nstates = size(states)
		
		open( 10, file=oFile )
		
		do i=1,nstates
			write(10,'(A,A)') "#", trim(states(i).str())
			write(10,'(F15.2,F15.7)') states(i).J()+0.0, states(i).E
			write(10,'(F15.2,F15.7)') states(i).J()+0.25, states(i).E
			write(10,*) ""
		end do
		
		write(10,*) ""
		
		do i=1,nstates
			do j=1,states(i).nStates
				if( states(i).coeff(j)**2 > 0.01_8 ) then
					write(10,'(F15.2,F15.7)') states(i).J()+0.25, states(i).E
					write(10,'(F15.2,F15.7)') states(i).J()+0.0, states(i).pureStates(j).E
					write(10,*) ""
				end if
			end do
		end do
		
		close( 10 )
	end subroutine MixedState_saveEnergies
	
end module MixedState_

