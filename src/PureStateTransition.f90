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

module PureStateTransition_
	use Math_
	use PureState_
	implicit none
	private
	
	public :: &
		PureStateTransition_saveSpectra
		
	integer, parameter :: O_BRACH = -2
	integer, parameter :: P_BRACH = -1
	integer, parameter :: Q_BRACH =  0
	integer, parameter :: R_BRACH =  1
	integer, parameter :: S_BRACH =  2
	integer, parameter :: MAX_BRACH = S_BRACH
	character(1), parameter :: BRANCH_NAME(-MAX_BRACH:MAX_BRACH) = [ "O", "P", "Q", "R", "S" ]
!	real(8), parameter :: MIN_INTENSITY = 1e-5
	real(8), parameter :: MIN_INTENSITY = 1e-12
		
	type, public :: PureStateTransition
		real(8) :: deltaE
		real(8) :: intensity
		real(8) :: Gamma
		type(PureState) :: sBegin
		type(PureState) :: sEnd
		
		contains
			procedure :: init
			procedure :: destroy
			procedure :: str
	end type PureStateTransition
	
	contains
	
	subroutine init( this, deltaE, intensity, sBegin, sEnd )
		class(PureStateTransition) :: this
		real(8), intent(in) :: deltaE
		real(8), intent(in) :: intensity
		type(PureState), intent(in) :: sBegin
		type(PureState), intent(in) :: sEnd
		
		this.deltaE = deltaE
		this.intensity = intensity
		this.Gamma = sEnd.Gamma
		call this.sBegin.copy( sBegin )
		call this.sEnd.copy( sEnd )
	end subroutine init
	
	subroutine destroy( this )
		class(PureStateTransition) :: this
		
		this.deltaE = -100.0_8
		this.intensity = 0.0_8
		call this.sBegin.destroy()
		call this.sEnd.destroy()
	end subroutine destroy
	
	function str( this ) result( output )
		class(PureStateTransition) :: this
		character(len=200) :: output
		
		!! Fue necesario, porque parece que no le gusta hacer cosas como .sBegin.
		!! probablemente los confunda con operadores
		type(PureState) :: sBegin
		type(PureState) :: sEnd
		
		sBegin = this.sBegin
		sEnd = this.sEnd
		
		output = trim(sBegin.str())//"  -->  "//trim(sEnd.str())
	end function str
	
	subroutine PureStateTransition_saveSpectra( statesNu0, statesNu1, oPrefix, alpha, beta, EBegin, EEnd, dE )
		use Math_
		use UnitsConverter_
		type(PureState), intent(in) :: statesNu0(:)
		type(PureState), intent(in) :: statesNu1(:)
		character(*), intent(in) :: oPrefix
		real(8), intent(in) :: alpha
		real(8), intent(in) :: beta
		real(8), intent(in) :: EBegin, EEnd, dE
		
		type(PureStateTransition), allocatable :: transitions(:,:)
		integer, allocatable :: nTransitions(:)
		
		integer :: i, j, k
		integer :: nstates
		type(PureState) :: s1, s2
		
		integer :: deltaJ
		integer :: deltaOmega
		real(8) :: deltaE
		integer :: dn
		real(8) :: minJi
		real(8) :: Mi
		real(8) :: alphai, betai
		real(8) :: intensity
		real(8) :: dip
		real(8) :: Gamma
		real(8) :: maxInt
		real(8) :: normFactor
		real(8) :: E, sumIntensity
		real(8) :: integ
		real(8) :: ssum
		
		logical :: firstIn
		
! 		EBegin = -2.0_8
! 		EEnd = 8.0_8
! 		dE = 0.00001_8
		
		maxInt = -1.0_8
		
		nstates = size(statesNu0)
		
		allocate( transitions(-MAX_BRACH:MAX_BRACH,nstates**2) ) ! 5 por que solo se consideran O,P,Q,R,S
		allocate( nTransitions(-MAX_BRACH:MAX_BRACH) )
		
		write(*,*) "Transitions among pure states"
		write(*,*) "============================="
		write(*,*) ""
		
		write(*,"(A40,A40,A15,A10,A11,A5)") "Initial State", "End State", "dE", "I", "Gamma", "Type"
		write(*,"(A40,A40,A15,A10,A11,A5)") "-------------", "---------", "--", "-", "-----", "----"
		write(*,*) ""
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! La transición analizada es:
		!
		!   s1(nu=0) ----------> s2(nu=0)
		!    w1,E1     dipInd      E2
		!             deltaE1    
		!
		! La  transición está condicionada a la probabilidad
		! de transición por momento de dipolo inducido. El
		! único peso estadístico que afecta el cálculo es el
		! asociado a s1. Básicamente es un espectro de
		! absorción
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		nTransitions = 1
		firstIn = .true.
		do i=1,nstates
			s1 = statesNu0(i)
			
			do j=1,nstates
				s2 = statesNu1(j)
				
				deltaJ = ceiling(s2.J-s1.J)
				deltaOmega = ceiling(s2.Omega-s1.Omega)
				deltaE = s2.E-s1.E
				
				dn = 0 ! Asegura que S y Lambda sean iguales
				if( trim(s1.pID) == trim(s2.pID) ) dn = 1 ! Es el mismo estado heliónico
				
				dip = 0.0_8
				if( deltaOmega == 0 .and. dn == 1 ) then
					minJi = min( s1.J, s2.J )
					
					do Mi=-minJi,minJi
						alphai = alpha*Math_wigner3j( real(s1.J,8), 0.0_8, real(s2.J,8), -real(Mi,8), 0.0_8, real(Mi,8) ) &
							*Math_wigner3j( real(s1.J,8), 0.0_8, real(s2.J,8), -s1.Omega, 0.0_8, s2.Omega )
						
						betai = (2.0_8*beta/3.0_8)*Math_wigner3j( real(s1.J,8), 2.0_8, real(s2.J,8), -real(Mi,8), 0.0_8, real(Mi,8) ) &
							*Math_wigner3j( real(s1.J,8), 2.0_8, real(s2.J,8), -s1.Omega, 0.0_8, s2.Omega )
							
						dip = dip + ( alphai+betai )**2
					end do
				end if
				
				intensity = s1.weight*(2.0_8*s2.J+1.0_8)*dip
				
				if( intensity > MIN_INTENSITY .and. deltaJ <= MAX_BRACH ) then
					
					if( firstIn ) then
						write(*,"(A40)") trim(s1.str())
						firstIn = .false.
					end if
					
					k = nTransitions(deltaJ)
					
					call transitions( deltaJ, k ).init( deltaE, intensity, s1, s2 )
					
					write(*,"(40X,A40,F15.6,F10.6,E11.2,A5)") trim(s2.str()), &
						deltaE, intensity, s2.Gamma, BRANCH_NAME(deltaJ)
					
					nTransitions(deltaJ) = nTransitions(deltaJ) + 1
					
					if( maxInt < intensity ) maxInt = intensity
				end if
			end do
			
			if( .not. firstIn ) then
				write(*,*) ""
				firstIn = .true.
			end if
		end do
		nTransitions = nTransitions - 1
		
		write(*,*) ""
		
		ssum = 0.0
		do i=O_BRACH,S_BRACH
			open( 10+i, file=trim(oPrefix)//"Spec"//BRANCH_NAME(i)//".dat" )
			ssum = ssum + sum( transitions(i,1:nTransitions(i)).intensity*transitions(i,1:nTransitions(i)).Gamma )
		end do
		
		normFactor = 1.0_8/( Math_PI*ssum )
		
		do i=O_BRACH,S_BRACH
			do j=1,nTransitions(i)
				deltaE = transitions(i,j).deltaE
				intensity = normFactor*transitions(i,j).intensity
				
				write(10+i,'(A,2F15.7,A,A)') "#  ", deltaE, intensity, "   ", trim(transitions(i,j).str())
				write(10+i,'(2F15.7)') deltaE, 0.0_8
				write(10+i,'(2F15.7)') deltaE, intensity
				write(10+i,*) ""
			end do
		end do
		
		do i=O_BRACH,S_BRACH
			write(10+i,*) ""
		end do
		
		ssum = 0.0_8
		do E=EBegin,EEnd,dE
		
			integ = 0.0_8
			do i=O_BRACH,S_BRACH
				sumIntensity = 0.0
				do j=1,nTransitions(i)
					deltaE = transitions(i,j).deltaE
					intensity = transitions(i,j).intensity
					Gamma = transitions(i,j).Gamma
					
					sumIntensity = sumIntensity + normFactor*intensity*lorentzian( E, deltaE, Gamma )
				end do
				write(10+i,*) E, sumIntensity
				integ = integ + sumIntensity
			end do
			
! 			write(*,*) "Pure IntegP(", BRANCH_NAME(i), ") = ", integ*dE
			ssum = ssum + integ*dE
		end do
		write(*,*) "Total Integ     = ", ssum
		
		do i=O_BRACH,S_BRACH
			close( 10+i )
		end do
		
		deallocate( transitions )
	end subroutine PureStateTransition_saveSpectra
	
	function lorentzian( x, x0, gamma ) result ( F )
		real(8), intent(in) :: x
		real(8), intent(in) :: x0
		real(8), intent(in) :: gamma
		real(8) :: F
		
		F = gamma**2/( ( x-x0 )**2 + gamma**2 )
	end function lorentzian

end module PureStateTransition_
