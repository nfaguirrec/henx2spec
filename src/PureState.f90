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

module PureState_
	implicit none
	private
	
	integer, parameter, public :: BRA = 0
	integer, parameter, public :: KET = 1
	
	public :: &
		PureState_boltzmanWeights, &
		PureState_saveEnergies
	
	type, public :: PureState
		integer :: Nu
		real(8) :: J
		real(8) :: Omega
		real(8) :: S
		integer :: Lambda
		real(8) :: Sigma
		real(8) :: L2
		real(8) :: E
		real(8) :: weight
		real(8) :: Gamma
		character(2) :: pID
		character(2) :: sID
		
		contains
			procedure :: init
			procedure :: destroy
			procedure :: copy
                        generic :: assignment(=) => copy
			procedure :: str
	end type PureState
	
	contains
	
	subroutine init( this, Nu, J, Omega, S, Lambda, Sigma, L2, E, Gamma, pID, sID )
		class(PureState) :: this
		integer, intent(in) :: Nu
		real(8), intent(in) :: J
		real(8), intent(in) :: Omega
		real(8), intent(in) :: S
		integer, intent(in) :: Lambda
		real(8), intent(in) :: Sigma
		real(8), intent(in) :: L2
		real(8), intent(in) :: E
		real(8), intent(in) :: Gamma
		character(2), intent(in) :: pID
		character(2), optional, intent(in) :: sID
		
		this.Nu = Nu
		this.J = J
		this.Omega = Omega
		this.S = S
		this.Lambda = Lambda
		this.Sigma = Sigma
		this.L2 = L2
		this.E = E
		this.Gamma = Gamma
		this.pID = pID
		
		if( present(sID) ) then
			this.sID = sID
		else
			this.sID = ""
		end if
	end subroutine init
	
	subroutine destroy( this )
		class(PureState) :: this
		
	end subroutine destroy
	
	subroutine copy( this, other )
		class(PureState), intent(inout) :: this
		type(PureState), intent(in) :: other
		
		this.Nu = other.Nu
		this.J = other.J
		this.Omega = other.Omega
		this.S = other.S
		this.Lambda = other.Lambda
		this.Sigma = other.Sigma
		this.L2 = other.L2
		this.E = other.E
		this.Gamma = other.Gamma
		this.weight = other.weight
		this.pID = other.pID
		this.sID = other.sID
	end subroutine copy
	
	function str( this, t ) result( output )
		class(PureState) :: this
		integer, intent(in), optional :: t
		character(len=200) :: output
		
		integer :: tEff
		integer :: pot10
		character(len=200) :: strBuffer
		
		if ( present(t) ) then
			tEff = t
		else
			tEff = BRA
		end if
		
		output = ""
		
		if( tEff == BRA ) then
			output = trim(output)//"<"
		else
			output = trim(output)//"|"
		end if
		
		if( trim(this.pID) /= "" ) then
			if( len_trim(this.pID) == 1 ) then
				write(strBuffer, "(A1)") this.pID
			else
				write(strBuffer, "(A2)") this.pID
			end if
			output = trim(output)//trim(strBuffer)//"("
		end if
		
		write(strBuffer, "(I1)") this.Nu
		output = trim(output)//trim(strBuffer)//","
		
		if( this.J < 10.0 ) then
			write(strBuffer, "(F3.1)") this.J
		else if( this.J < 100.0 ) then
			write(strBuffer, "(F4.1)") this.J
		end if
		output = trim(output)//trim(strBuffer)//")"
		
		write(strBuffer, "(I1)") int(2.0*this.S+1.0)
		output = trim(output)//"^"//trim(strBuffer)
		
		select case( abs(this.Lambda) )
			case( 0 )
				output = trim(output)//"\Sigma"
			case( 1 )
				output = trim(output)//"\Pi"
			case( 2 )
				output = trim(output)//"\Delta"
			case( 3 )
				output = trim(output)//"\Phi"
			case( 4 )
				output = trim(output)//"\Gamma"
			case default
				write(strBuffer, "(A1,I3)") "\",abs(this.Lambda)
				output = trim(output)//strBuffer
		end select
		
		if( this.Lambda < 0 ) then
			write(strBuffer, "(A1)") "-"
		else
			write(strBuffer, "(A1)") "+"
		end if
		output = trim(output)//trim(strBuffer)
		
		if( this.Omega < 0 ) then
			write(strBuffer, "(F4.1)") this.Omega
		else
			write(strBuffer, "(F3.1)") this.Omega
		end if
		output = trim(output)//"("//trim(strBuffer)//","
		
		if( this.Sigma < 0 ) then
			write(strBuffer, "(F4.1)") this.Sigma
		else
			write(strBuffer, "(F3.1)") this.Sigma
		end if
		output = trim(output)//trim(strBuffer)//")"
		
		if( trim(this.sID) /= "" ) then
			if( len_trim(this.sID) == 1 ) then
				write(strBuffer, "(A1)") this.sID
			else
				write(strBuffer, "(A2)") this.sID
			end if
			output = trim(output)//trim(strBuffer)
		end if
		
		if( tEff == BRA ) then
			output = trim(output)//"|"
		else
			output = trim(output)//">"
		end if
		
! 		if( tEff == BRA ) then
! 			output = trim(output)//"<"
! 		else
! 			output = trim(output)//"|"
! 		end if
! 		
! 		write(strBuffer, "(I2)") this.Nu
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(I3)") this.J
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(I3)") this.Omega
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(F4.1)") this.S
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(I2)") this.Lambda
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(F5.1)") this.Sigma
! 		output = trim(output)//trim(strBuffer)//","
! 		
! 		write(strBuffer, "(A3)") trim(this.pID)
! 		output = trim(output)//trim(strBuffer)
! 		
! 		if( tEff == BRA ) then
! 			output = trim(output)//" |"
! 		else
! 			output = trim(output)//" >"
! 		end if
	end function str
	
	subroutine PureState_boltzmanWeights( states, T )
		use UnitsConverter_
		type(PureState), intent(inout) :: states(:)
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
		
	end subroutine PureState_boltzmanWeights
	
	subroutine PureState_saveEnergies( states, oFile )
		type(PureState), intent(in) :: states(:)
		character(*), intent(in) :: oFile
		
		integer :: i
		integer :: nstates
		
		nstates = size(states)
		
		open( 10, file=oFile )
		
		do i=1,nstates
			write(10,'(A,A)') "# ", trim(states(i).str())
			write(10,'(F15.2,F15.7)') states(i).J+0.0, states(i).E
			write(10,'(F15.2,F15.7)') states(i).J+0.25, states(i).E
			write(10,*) ""
		end do
		
		close( 10 )
	end subroutine PureState_saveEnergies
	
end module PureState_

