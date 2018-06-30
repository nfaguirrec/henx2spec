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

program HenX2Spec
	use InputParser_
	use Matrix_
	use PureState_
	use PureStateTransition_
	use MixedState_
	use MixedStateTransition_
	implicit none
	
	real(8), parameter :: ZERO = 1e-10
	
	character(100) :: iFile
	integer :: Jmax
	logical :: spinHalfInteger
	real(8) :: T = 0.5_8
	logical :: enableJS = .false.
	logical :: enableJL = .false.
	logical :: enableLS = .false.
	
	type(InputParser) :: iparser
	
	real(8) :: hw0 ! Energy for transition nu=0 <-- nu=1 (a.u.) in Cl2
	real(8) :: alpha, beta ! Polarizability components
	real(8) :: Bv(2) ! rotational constants for v=0 and v=1
	
	character(2) :: idUpperBase ! Necesario para verificar el valor de máximo de temperatura permitido
	integer :: nstates
	type(PureState), allocatable :: pureStates_v0(:)
	type(PureState), allocatable :: pureStates_v1(:)
	integer, allocatable :: Jvec(:) ! Contiene el número de estados con un valor de J, el i-ésimo elemento tiene el valor de J=0.5*(i-1)
	integer, allocatable :: Jvec2StatesIDs(:,:) ! El elemento Jvec2StatesIDs( i, 1:Jvec(i) ), retorna la lista de los IDs con J=0.5*(i-1)
	real(8), allocatable :: H_v0(:,:,:)
	real(8), allocatable :: H_v1(:,:,:)
	type(MixedState), allocatable :: mixedStates_v0(:)
	type(MixedState), allocatable :: mixedStates_v1(:)
	
	call main()
	
	contains
	
	!*
	! Convierte del valor de J a identificador
	! dentro de los arreglos
	!*
	function Ji2id( Ji ) result( id )
		real(8), intent(in) :: Ji
		integer :: id
		
		if( spinHalfInteger ) then
			id = int( Ji+0.5 )
		else
			id = int( Ji+1.0 )
		end if
	end function Ji2id
	
	!*
	! Convierte del valor del identificador
	! dentro de los arreglos a J
	!*
	function id2Ji( id ) result( Ji )
		integer, intent(in) :: id
		real(8) :: Ji
		
		if( spinHalfInteger ) then
			Ji = real(id,8)-0.5_8
		else
			Ji = real(id,8)-1.0_8
		end if
	end function id2Ji
	
	!*
	! Retorna el número máximo de valores de J
	! permitidos 
	!*
	function nJ() result( n )
		integer :: n
		
		if( spinHalfInteger ) then
			n = Jmax
		else
			n = Jmax+1
		end if
	end function nJ
	
	!*
	! El programa principal
	!*
	subroutine main()
		integer :: i, k
		integer :: Lambdai, signi
		real(8) :: Ji, Omegai, Si, L2i, Sigmai, Eibase
		character(2) :: pIDi
		real(8) :: GammaEff
		
		call loadCommandLineParameters()
		call iparser.init( trim(iFile)//".inp" )
		
		call MorseCl2()
		
		allocate( Jvec( nJ() ) )
		Jvec = 0
		
		i=0
		do k=1,iparser.nstates
			read( iparser.states(k), * ) pIDi, Si, Lambdai, L2i, Eibase
			
			if( int(Si+0.5) /= int(Si) ) spinHalfInteger = .true.
			
			do Sigmai=-Si,Si
				do signi=-1,1,2
					Omegai=signi*Lambdai+Sigmai
					do Ji=abs(Omegai),Jmax
						Jvec( Ji2id(Ji) ) = Jvec( Ji2id(Ji) ) + 1
						i = i+1
					end do
					
					if( Lambdai == 0 ) exit
				end do
			end do
		end do
		nstates = i
		
		! en el peor de los casos todos los estados son del mismo J
		allocate( Jvec2StatesIDs( nJ(), nstates ) )
		
		write(*,*) ""
		write(*,*) "J distribution"
		write(*,*) "=============="
		write(*,"(A10,A10)") "J", "nstates"
		write(*,"(A10,A10)") "---", "-------"
		do i=1,nJ()
			write(*,"(F10.1,I10)") id2Ji(i), Jvec(i)
		end do
		write(*,*) ""
		
		allocate( pureStates_v0(nstates) )
		allocate( pureStates_v1(nstates) )
		
		call buildStates()
		
		call PureState_saveEnergies( pureStates_v0, trim(iFile)//"PureStates.dat" )
		call PureStateTransition_saveSpectra( pureStates_v0, pureStates_v1, &
			trim(iFile)//"Pure", alpha, beta, iparser.Emin, iparser.Emax, iparser.Estep )
		
		allocate( H_v0( nJ(), maxval(Jvec), maxval(Jvec) ) )
		allocate( H_v1( nJ(), maxval(Jvec), maxval(Jvec) ) )
		
		call buildHMatrix()
		
		allocate( mixedStates_v0(nstates) )
		do i=1,nstates
			call mixedStates_v0(i).init(nstates)
! 			call mixedStates_v0(i).init(20)
		end do
		
		allocate( mixedStates_v1(nstates) )
		do i=1,nstates
			call mixedStates_v1(i).init(nstates)
! 			call mixedStates_v1(i).init(20)
		end do
		
		call diagHMatrix()
		
		deallocate( H_v0 )
		deallocate( H_v1 )
		
		call MixedState_saveEnergies( mixedStates_v0, trim(iFile)//"MixedStates.dat" )
		call MixedStateTransition_saveSpectra( mixedStates_v0, mixedStates_v1, &
			trim(iFile)//"Mixed", alpha, beta, iparser.Emin, iparser.Emax, iparser.Estep )
		
		do i=1,nstates
			call mixedStates_v0(i).destroy()
		end do
		deallocate( mixedStates_v0 )
		
		do i=1,nstates
			call mixedStates_v1(i).destroy()
		end do
		deallocate( mixedStates_v1 )
		
		deallocate( Jvec2StatesIDs )
		call iparser.destroy()
	end subroutine main
	
	!*
	! Interpreta los parámetros desde entrada estandar
	!*
	subroutine loadCommandLineParameters()
		character(50) :: buffer
		integer :: i
		
		if( command_argument_count() < 3 ) then
			write(*,*) "Usage:"
			write(*,*) "    $ HenX2_spec iFile T Jmax [JS] [JL] [LS] [ALL]"
			stop
		else
			call get_command_argument( 1, buffer )
			read( buffer, * ) iFile
			
			call get_command_argument( 2, buffer )
			read( buffer, * ) T
			
			call get_command_argument( 3, buffer )
			read( buffer, * ) Jmax
			
			do i=4,command_argument_count()
				call get_command_argument( i, buffer )
				
				select case( trim(buffer) )
					case( "JS" )
						enableJS = .true.
					case( "JL" )
						enableJL = .true.
					case( "LS" )
						enableLS = .true.
					case( "ALL" )
						enableJS = .true.
						enableJL = .true.
						enableLS = .true.
				end select
			end do
		end if
		
		i = scan( trim(iFile), ".", .true. )
		iFile = iFile( 1:i-1 )
		
		write(*,*) "input file  = ", trim(iFile)//".inp"
		write(*,"(A15,F5.2)") "temperature = ", T
		write(*,  "(A15,I2)") "      Jmax  = ", Jmax
		write(*,*) ""
		write(*,*) "enable JS = ", enableJS
		write(*,*) "enable JL = ", enableJL
		write(*,*) "enable LS = ", enableLS
		write(*,*) ""
	end subroutine loadCommandLineParameters
	
	!*
	! Calcula todas las cosas necesarias de la molecula diatomica aislada
	!*
	subroutine MorseCl2()
		use Morse_
		use Grid_
		use ThrularNumerovMethod_
		use NFunction_
		use NIntegration_
		use UnitsConverter_
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Cl2 X¹Σg+ 
		!! NIST Standard Reference Data Program
		!! -------------------------------------
		!! re    = 1.9879 Å
		!! ωe    = 559.72 cm⁻¹
		!! ωexe  = 2.675 cm⁻¹
		!! m(Cl) = 35.4257 amu
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		real(8) :: Re = 1.9879_8*angs
! ! 		real(8) :: Re = 3.756776_8   <-- Desde Pablo
! 		real(8) :: we = 559.72_8*cm1
! 		real(8) :: wexe = 2.675_8*cm1
! 		real(8) :: rMass = 0.5_8*35.4257_8*amu

		real(8) :: Re
		real(8) :: we
		real(8) :: wexe
		real(8) :: rMass

		type(Morse) :: potential
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Cl2 Polarizabilidad esfericas y anisotropica
		!! G. Maroulis, Mol. Phys. 77 (1992) 1085-1094
		!!--------------------------------------------------------------
		!! alpha = 30.91 +  7.05 (r-re) + 1.73 (r-re)² - 2.17 (r-re)³
		!! beta  = 16.83 + 12.52 (r-re) + 5.49 (r-re)² - 4.50 (r-re)³
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		real(8) :: alphaCi(4) = [ 30.91_8,  7.05_8, 1.73_8, -2.17_8 ]
		real(8) :: alphaCi(4)
! 		real(8) :: betaCi(4)  = [ 16.83_8, 12.52_8, 5.49_8, -4.50_8 ]
		real(8) :: betaCi(4)
		real(8) :: rRe(4)
		
		type(Grid) :: rGrid
		type(ThrularNumerovMethod) :: solver
		type(NFunction) :: nState1, nState2, nFunc
		type(NIntegration) :: integrator
		integer :: i
		
		Re = iparser.Re*angs
		we = iparser.we*cm1
		wexe = iparser.wexe*cm1
		rMass = iparser.rMass*amu
		alphaCi = iparser.polAlpha
		betaCi  = iparser.polBeta
		
		call rGrid.init( 0.5_8, 20.0_8, 10000 )
		call potential.fromExp( rGrid, Re, we, wexe, rMass )
! 		call potential.show()
		
		write(*,*) "Cl2 spectrum"
		write(*,*) "============"
		write(*,*) ""
		
		call solver.init( potential.parent(), rMass=rMass )
		call solver.run()
		
		write(*,"(a5,a20,a20)") "nu", "exact(cm-1)", "numeric(cm-1)"
		write(*,"(a5,a20,a20)") "--", "-----------", "-------------"
		do i=1,solver.nStates
			if ( solver.eigenValue(i) < 0.0_8 ) then
					write(*,"(i5,f20.10,f20.10)") i-1, potential.exactEigenValue(i-1, rMass)/cm1, solver.eigenValue(i)/cm1
			end if
		end do
		
		hw0 = potential.exactEigenValue(2, rMass)/cm1-potential.exactEigenValue(1, rMass)/cm1
		
		write(*,*) ""
		write(*,'(A,F15.10,A)') " Transition energy ( 0 --> 1 ) = ", hw0, " cm-1"
		
		call nState1.copy( solver.eigenFunction(1) )
		call nState2.copy( solver.eigenFunction(2) )
		
		write(*,*) ""
		write(*,*) "Polarizability"
		write(*,*) "=============="
		write(*,*) ""
		write(*,'(A5,A30,A30,A30)') "i", "< 0 | (r-Re)^i | 1 >", "alpha", "beta"
		write(*,'(A5,A30,A30,A30)') "-", "--------------------", "-----", "----"
		
		do i=0,3
			call nFunc.copy( nState1 )
			
			nFunc.yArray = nState1.yArray*nState2.yArray*( nFunc.xGrid.data - Re )**i
			
			call integrator.init( nFunc, BOOLE )
			rRe(i+1) = integrator.evaluate()
			
			write(*,'(I5,F30.10,F30.10,F30.10)') i, rRe(i+1), alphaCi(i+1), betaCi(i+1)
			
			call nFunc.destroy()
		end do
		
		alpha = sum(alphaCi*rRe)
		beta  = sum(betaCi*rRe)
		
		write(*,'(35X,2A30)') "------------", "------------"
		write(*,'(35X,2F30.10)') alpha, beta
		
		write(*,*) ""
		write(*,*) "Rotational constants"
		write(*,*) "===================="
		write(*,*) ""
		
		call nFunc.copy( nState1 )
		
		nFunc.yArray = nState1.yArray*nState1.yArray/( 2.0_8*rMass*nFunc.xGrid.data**2 ) ! < f1 | 1/2mR**2 | f1 >
		call integrator.init( nFunc, BOOLE )
		Bv(1) = integrator.evaluate()/cm1
		
		nFunc.yArray = nState2.yArray*nState2.yArray/( 2.0_8*rMass*nFunc.xGrid.data**2 ) ! < f2 | 1/2mR**2 | f2 >
		call integrator.init( nFunc, BOOLE )
		Bv(2) = integrator.evaluate()/cm1
! 		
		call nFunc.destroy()
		
		write(*,'(A30,F30.10,A10)') "< 0 | 1/2mr**2 | 0 > = ", Bv(1), "cm-1"
		write(*,'(A30,F30.10,A10)') "< 1 | 1/2mr**2 | 1 > = ", Bv(2), "cm-1"
		write(*,*) ""
		
		write(*,*) ""
		write(*,*) "Internuclear distances"
		write(*,*) "======================"
		write(*,*) ""
		
		call nFunc.copy( nState1 )
		
		nFunc.yArray = nState1.yArray*nState1.yArray*nFunc.xGrid.data   ! < f1 | r | f1 >
		call integrator.init( nFunc, BOOLE )
		
		write(*,'(A30,F30.10,A10)') "< 0 | r | 0 > = ", integrator.evaluate()/Angs, "A"
		
		nFunc.yArray = nState2.yArray*nState2.yArray*nFunc.xGrid.data   ! < f2 | r | f2 >
		call integrator.init( nFunc, BOOLE )
! 		
		write(*,'(A30,F30.10,A10)') "< 1 | r | 1 > = ", integrator.evaluate()/Angs, "A"
		
		call nFunc.destroy()
		
		write(*,*) ""
		
		call potential.destroy()
		call rGrid.destroy()
		call solver.destroy()
		call integrator.destroy()
	end subroutine MorseCl2
	
	!*
	! Construye los estados y acoplamientos desde el archivo de entrada
	!*
	subroutine buildStates()
		use UnitsConverter_
		use Math_
		character(50) :: buffer
		type(PureState) :: s1, s2
		
		integer :: Lambdai, signi
		real(8) :: Ji, Omegai, Si, Sigmai
		real(8) :: L2i, Ei, Lij, Eibase, Gammai
		character(2) :: pIDi, pIDj
		character(2) :: sIDi
		integer :: id
		
		integer :: i, j, k
		integer, allocatable :: iterJ(:)
		
		! Necesario para verificar el valor de máximo de temperatura permitido
		real(8) :: energyUpperBase
		real(8) :: testMaxT
		
		 ! Vectores temporales para ordenar los estados por energía
		integer, allocatable :: sIndex(:)
		real(8), allocatable :: energiesToSort(:)
		
		allocate( iterJ( nJ() ) )
		
		idUpperBase = ""
		energyUpperBase = -1e8
		
		iterJ = 1
		i=0
		do k=1,iparser.nstates
			
			sIDi = ""
			select case( InputParser_nItems( trim(iparser.states(k)) ) )
				case( 7 )
					read( iparser.states(k), * ) pIDi, Si, Lambdai, L2i, Eibase, Gammai, sIDi
				case( 6 )
					read( iparser.states(k), * ) pIDi, Si, Lambdai, L2i, Eibase, Gammai
				case default
					write(*,*) "## ERROR ## Bad description in state,"
					write(*,*) "        ", iparser.states(k)
					stop
			end select
			
			if( Eibase > energyUpperBase ) then
				idUpperBase = pIDi
			end if
			
			do Sigmai=-Si,Si
				do signi=-1,1,2
					Omegai=signi*Lambdai+Sigmai
					
					do Ji=abs(Omegai),Jmax
						Ei = Eibase + Bv(1)*( Ji*(Ji+1.0_8) + Si*(Si+1.0_8) + L2i &
							- Omegai**2 - Sigmai**2 - Lambdai**2 )
							
						call pureStates_v0( i+1 ).init( 0, Ji, Omegai, Si, signi*Lambdai, &
							Sigmai, L2i, Ei, Gammai, pIDi, sIDi )
						
						Ei = Eibase + Bv(2)*( Ji*(Ji+1.0_8) + Si*(Si+1.0_8) + L2i &
							- Omegai**2 - Sigmai**2 - Lambdai**2 )
							
						call pureStates_v1( i+1 ).init( 1, Ji, Omegai, Si, signi*Lambdai, &
							Sigmai, L2i, Ei, Gammai, pIDi, sIDi )
						
						Jvec2StatesIDs( Ji2id(Ji), iterJ( Ji2id(Ji) ) ) = i+1
						iterJ( Ji2id(Ji) ) = iterJ( Ji2id(Ji) ) + 1 ! incrementa el contador para el iterador con J=Ji
						
						i = i+1
					end do
					
					if( Lambdai == 0 ) exit
				end do
			end do
		end do
		
		deallocate( iterJ )
		
		call PureState_boltzmanWeights( pureStates_v0, T )
		call PureState_boltzmanWeights( pureStates_v1, T )
		
		write(*,*) "PURE STATES ( nu = 0 )"
		write(*,*) "======================"
		write(*,*) ""
		write(*,*) "Upper primitive state = ", idUpperBase
		write(*,*) ""
		write(*,'(12X,2A10)') "Energy", "Weight"
		write(*,'(12X,2A10)') "------", "------"
			
		do i=1,nJ()
			write(*,"(A,F4.1)") "J = ", id2Ji(i)
			
			allocate( sIndex(Jvec(i)) )
			allocate( energiesToSort(Jvec(i)) )
			
			! organiza los estados por su valor de energía
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				energiesToSort(j) = pureStates_v0(id).E
			end do
			
			call Math_sort( energiesToSort, sIndex )
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, sIndex(j) )
				
				write(*,'(2I5,2X,2F10.5,E11.2,A40)', advance="no") j, sIndex(j), pureStates_v0(id).E, pureStates_v0(id).weight, &
					pureStates_v0(id).Gamma, trim(pureStates_v0(id).str())
				
				if( trim(pureStates_v0(id).pID) == trim(idUpperBase) ) then
					testMaxT = pureStates_v0(id).weight/maxval(pureStates_v0.weight)/(2.0_8*pureStates_v0(id).J+1.0_8)
					write(*,'(F10.5)') testMaxT
					
					if( testMaxT > 0.1_8 ) then
						write(0,*) "## ERROR ## Maximum temperature available have been exceeded"
! 						stop
					end if
				else
					write(*,'(A)') ""
				end if
			end do
			
			deallocate( sIndex )
			deallocate( energiesToSort )
			
			write(*,*) ""
		end do
		
		write(*,*) "PURE STATES ( nu = 1 )"
		write(*,*) "======================"
		write(*,*) ""
		write(*,'(12X,2A10)') "Energy", "Weight"
		write(*,'(12X,2A10)') "------", "------"
			
		do i=1,nJ()
			write(*,"(A,F4.1)") "J = ", id2Ji(i)
			
			allocate( sIndex(Jvec(i)) )
			allocate( energiesToSort(Jvec(i)) )
			
			! organiza los estados por su valor de energía
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				energiesToSort(j) = pureStates_v1(id).E
			end do
			
			call Math_sort( energiesToSort, sIndex )
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, sIndex(j) )
				
				write(*,'(2I5,2X,2F10.5,E11.2,A40)') j, sIndex(j), pureStates_v0(id).E, pureStates_v1(id).weight, &
					pureStates_v1(id).Gamma, trim(pureStates_v1(id).str())
			end do
			
			deallocate( sIndex )
			deallocate( energiesToSort )
			
			write(*,*) ""
		end do
		
		open( 16, file="gEvsJ-Pure.dat" )
		
		do i=1,nJ()
			allocate( sIndex(Jvec(i)) )
			allocate( energiesToSort(Jvec(i)) )
			
			! organiza los estados por su valor de energía
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				energiesToSort(j) = pureStates_v1(id).E
			end do
			
			call Math_sort( energiesToSort, sIndex )
			
			write(16,"(F4.1)",advance='no') id2Ji(i)
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, sIndex(j) )
				
! 				write(16,"(F4.1,F12.5,5X,A)") id2Ji(i), pureStates_v0( id ).E, trim(pureStates_v0( id ).str())
				write(16,"(F15.5,F12.5)",advance='no') pureStates_v0( id ).E, pureStates_v1( id ).E
			end do
			write(16,"(A)") ""
			
			deallocate( sIndex )
			deallocate( energiesToSort )
		end do
		
		write(16,*) ""
		write(16,*) ""
		
		do i=1,nJ()
			allocate( sIndex(Jvec(i)) )
			allocate( energiesToSort(Jvec(i)) )
			
			! organiza los estados por su valor de energía
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				energiesToSort(j) = pureStates_v1(id).E
			end do
			
			call Math_sort( energiesToSort, sIndex )
			
			write(16,"(A,F4.1,5X)",advance='no') "#  ", id2Ji(i)
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, sIndex(j) )
				
				if( j==1 ) then
					write(16,"(A)") trim(pureStates_v0( id ).str())
				else
					write(16,"(A,9X,A)") "#  ", trim(pureStates_v0( id ).str())
				end if
			end do
			
			deallocate( sIndex )
			deallocate( energiesToSort )
		end do
		
		close(16)

		
	end subroutine buildStates
	
	!*
	! Construye la matriz H
	!*
	subroutine buildHMatrix()
		integer :: i, j, k, l
		real(8) :: Ji, id, idj, idk
		type(PureState) :: s1_v0, s2_v0
		type(PureState) :: s1_v1, s2_v1
		real(8) :: Lij, Lpij, Lmij, coup
		character(2) :: pIDi, pIDj
		logical :: thereIsCouplings
		integer :: dOmega, dLambda, dS, dSigma, dn
		logical :: couplingOn
		
		H_v0 = 0.0_8
		H_v1 = 0.0_8
		
		! Diagonal part
		do i=1,nJ()
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				H_v0( i, j, j ) = pureStates_v0(id).E
				H_v1( i, j, j ) = pureStates_v1(id).E
			end do
		end do
		
		write(*,*) ""
		write(*,*) "Couplings"
		write(*,*) "========="
		write(*,*) ""
			
		do i=1,nJ()
			write(*,"(A,F4.1)") "J = ", id2Ji(i)
			
			do j=1,Jvec(i)
				idj = Jvec2StatesIDs( i, j )
				
				s1_v0 = pureStates_v0( idj )
				s1_v1 = pureStates_v1( idj )
				
! 				do k=j+1,Jvec(i)
				do k=1,Jvec(i)
					idk = Jvec2StatesIDs( i, k )
					
					s2_v0 = pureStates_v0( idk )
					s2_v1 = pureStates_v1( idk )
					
					Lpij = 0.0_8
					Lmij = 0.0_8
					
					do l=1,iparser.nCouplingL
						read( iparser.couplingL(l), * ) pIDi, pIDj, Lij
						
						if( ( trim(s1_v0.pID) == trim(pIDi) .and. trim(s2_v0.pID) == trim(pIDj) ) .or. &
						    ( trim(s1_v0.pID) == trim(pIDj) .and. trim(s2_v0.pID) == trim(pIDi) ) ) then
							if( s1_v0.Lambda > s2_v0.Lambda ) Lpij = Lij
							if( s1_v0.Lambda < s2_v0.Lambda ) Lmij = Lij
						end if
					end do
					
					dOmega = int(s2_v0.Omega - s1_v0.Omega)
					dLambda = int(s2_v0.Lambda - s1_v0.Lambda)
					dS = int(s2_v0.S - s1_v0.S)
					dSigma = int(s2_v0.Sigma - s1_v0.Sigma)
					
					dn = 0 ! Asegura que S y Lambda sean iguales
					if( trim(s1_v0.pID) == trim(s2_v0.pID) ) dn = 1 ! Es el mismo estado heliotrónico
					
					coup = 0.0_8
					
					if( enableJS ) then
						if( ( ( dOmega == 1 .and. dSigma == 1 ) &
							.or. ( dOmega == -1 .and. dSigma == -1 ) ) .and. dn == 1 ) then
							
							coup = sqrt( s1_v0.J*(s1_v0.J+1.0) - s1_v0.Omega*s2_v0.Omega ) &
								*sqrt( s1_v0.S*(s1_v0.S+1.0) - s1_v0.Sigma*s2_v0.Sigma )
								
							H_v0( i, j, k ) = H_v0( i, j, k ) - Bv(1)*coup
							H_v1( i, j, k ) = H_v1( i, j, k ) - Bv(2)*coup
							
							write(*,"(2I5,5X,A35,A,A35,A,F8.4)") j, k, &
								trim(s1_v0.str(BRA)), " -B(J+S- + J-S+) ", adjustl(s2_v0.str(KET)), "  = ", -Bv(1)*coup
								
							write(*,"(10X,5X,A35,A,A35,A,F8.4)") &
								trim(s1_v1.str(BRA)), " -B(J+S- + J-S+) ", adjustl(s2_v1.str(KET)), "  = ", -Bv(2)*coup
						end if
					end if
					
					if( enableJL ) then
						couplingOn = .false.
						
						if( dOmega == 1 .and. dLambda == 1 .and. dS == 0 &
							.and. dSigma == 0 .and. abs( Lmij ) > ZERO ) then
						
							coup = Lmij*sqrt( s1_v0.J*(s1_v0.J+1.0) - s1_v0.Omega*s2_v0.Omega )
							
							couplingOn = .true.
								
						else if( dOmega == -1 .and. dLambda == -1 .and. dS == 0 &
							.and. dSigma == 0 .and. abs( Lpij ) > ZERO ) then
						
							coup = Lpij*sqrt( s1_v0.J*(s1_v0.J+1.0) - s1_v0.Omega*s2_v0.Omega )
							
							couplingOn = .true.
						end if
					
						if( abs(dOmega) == 1 .and. abs(dLambda) == 1 .and. abs(dS) == 0 &
							.and. abs(dSigma) == 0 .and. abs( Lmij ) > ZERO ) then
						
							coup = Lmij*sqrt( s1_v0.J*(s1_v0.J+1.0) - s1_v0.Omega*s2_v0.Omega )
							
							couplingOn = .true.
						end if
						
						if( couplingOn ) then
							H_v0( i, j, k ) = H_v0( i, j, k ) - Bv(1)*coup
							H_v1( i, j, k ) = H_v1( i, j, k ) - Bv(2)*coup
							
							write(*,"(2I5,5X,A35,A,A35,A,F8.4)") j, k, &
								trim(s1_v0.str(BRA)), " -B(J+L- + J-L+) ", &
									adjustl(s2_v0.str(KET)), "  = ", -Bv(1)*coup
								
							write(*,"(10X,5X,A35,A,A35,A,F8.4)") &
								trim(s1_v1.str(BRA)), " -B(J+L- + J-L+) ", &
									adjustl(s2_v1.str(KET)), "  = ", -Bv(2)*coup
						end if
					end if
					
					if( enableLS ) then
						couplingOn = .false.
						
						if( dLambda == -1 .and. dSigma == 1 .and. dOmega == 0 &
							.and. dS == 0 .and. abs( Lpij ) > ZERO ) then
							
							coup = Lpij*sqrt( s1_v0.S*(s1_v0.S+1.0) - s1_v0.Sigma*(s1_v0.Sigma+1.0) )
								
							couplingOn = .true.
							
						else if( dLambda == 1 .and. dSigma == -1 .and. dOmega == 0 &
							.and. dS == 0 .and. abs( Lmij ) > ZERO ) then
							
							coup = Lmij*sqrt( s1_v0.S*(s1_v0.S+1.0) - s1_v0.Sigma*(s1_v0.Sigma-1.0) )
							
							couplingOn = .true.
						end if
						
						if( couplingOn ) then
							H_v0( i, j, k ) = H_v0( i, j, k ) + Bv(1)*coup
							H_v1( i, j, k ) = H_v1( i, j, k ) + Bv(2)*coup
							
							write(*,"(2I5,5X,A35,A,A35,A,F8.4)") j, k, &
								trim(s1_v0.str(BRA)), "  B(L+S- + L-S+) ", &
									adjustl(s2_v0.str(KET)), "  = ", Bv(1)*coup
								
							write(*,"(10X,5X,A35,A,A35,A,F8.4)") &
								trim(s1_v1.str(BRA)), "  B(L+S- + L-S+) ", &
									adjustl(s2_v1.str(KET)), "  = ", Bv(2)*coup
						end if
					end if
				end do
			end do
		end do
		
		! Se simetrizan las matrices H
! 		do i=1,nJ()
! 			do j=1,Jvec(i)
! 				do k=j+1,Jvec(i)
! 					H_v0( i, k, j ) = H_v0( i, j, k )
! 					H_v1( i, k, j ) = H_v1( i, j, k )
! 				end do
! 			end do
! 		end do
		
		do i=1,nJ()
			write(*,*) ""
			write(*,"(1X,A,F4.1,2A)") "Hv0 Matrix ( J = ", id2Ji(i), " ) "
			write(*,*) "======================="
			write(*,*) ""
			call Matrix_show( H_v0(i,1:Jvec(i),1:Jvec(i)) )
		end do
		
		do i=1,nJ()
			write(*,*) ""
			write(*,"(1X,A,F4.1,2A)") "Hv1 Matrix ( J = ", id2Ji(i), " )"
			write(*,*) "======================="
			write(*,*) ""
			call Matrix_show( H_v1(i,1:Jvec(i),1:Jvec(i)) )
		end do
		write(*,*) ""
		
	end subroutine buildHMatrix
	
	!*
	! Diagonaliza la matriz H
	!*
	subroutine diagHMatrix()
		real(8), allocatable :: eVecs(:,:)
		real(8), allocatable :: eVals(:)
		real(8) :: norm
		integer :: i, j, k
		integer :: id, id1, id2
		integer :: n, m
		integer :: info
		real(8), allocatable :: H(:,:)
		real(8) :: maxC
		real(8) :: testMaxT
		
		integer :: Lambdai
		real(8) :: Sigmai
		
		do i=1,nJ()
			allocate( H( Jvec(i), Jvec(i) ) )
			allocate( eVecs( Jvec(i), Jvec(i) ) )
			allocate( eVals( Jvec(i) ) )
			
			H = H_v0( i, 1:Jvec(i), 1:Jvec(i) )
			eVecs = 0.0_8
			eVals = 0.0_8
			
			call Matrix_eigen( H, eVecs, eVals )
				
			do j=1,Jvec(i)
				id1 = Jvec2StatesIDs( i, j )
				
				mixedStates_v0( id1 ).E = eVals(j)
				
				maxC = 1e-6
				norm = 0.0_8
				n=1
				do k=1,Jvec(i)
					id2 = Jvec2StatesIDs( i, k )
					
					if( eVecs(k,j)**2 > 1e-4 ) then
						mixedStates_v0( id1 ).coeff(n) = eVecs(k,j)
						mixedStates_v0( id1 ).pureStates(n) = pureStates_v0( id2 )
						
						norm = norm + eVecs(k,j)**2
						
						Lambdai = pureStates_v0( id2 ).Lambda
						Sigmai = pureStates_v1(id2).Sigma
						if( sign( eVecs(k,j)**2, Lambdai+Sigmai ) > maxC**2 ) maxC = eVecs(k,j)
						
						n = n+1
					end if
				end do
				
				mixedStates_v0( id1 ).coeff = (maxC/abs(maxC))*mixedStates_v0( id1 ).coeff
				
				mixedStates_v0( id1 ).nStates = n-1
			end do
			
			deallocate(H)
			deallocate( eVecs )
			deallocate( eVals )
		end do
		
		do i=1,nJ()
			allocate( H( Jvec(i), Jvec(i) ) )
			allocate( eVecs( Jvec(i), Jvec(i) ) )
			allocate( eVals( Jvec(i) ) )
			
			H = H_v1( i, 1:Jvec(i), 1:Jvec(i) )
			eVecs = 0.0_8
			eVals = 0.0_8
			
			call Matrix_eigen( H, eVecs, eVals )
				
			do j=1,Jvec(i)
				id1 = Jvec2StatesIDs( i, j )
				
				mixedStates_v1( id1 ).E = eVals(j)
				
				maxC = 1e-6
				norm = 0.0_8
				n=1
				do k=1,Jvec(i)
					id2 = Jvec2StatesIDs( i, k )
					
					if( eVecs(k,j)**2 > 1e-4 ) then
						mixedStates_v1( id1 ).coeff(n) = eVecs(k,j)
						mixedStates_v1( id1 ).pureStates(n) = pureStates_v1( id2 )
						
						norm = norm + eVecs(k,j)**2
						
						Lambdai = pureStates_v1(id2).Lambda
						Sigmai = pureStates_v1(id2).Sigma
						if( sign( eVecs(k,j)**2, Lambdai+Sigmai ) > maxC**2 ) maxC = eVecs(k,j)
						
						n = n+1
					end if
				end do
				
				mixedStates_v1( id1 ).coeff = (maxC/abs(maxC))*mixedStates_v1( id1 ).coeff
				
				mixedStates_v1( id1 ).nStates = n-1
			end do
			
			deallocate(H)
			deallocate( eVecs )
			deallocate( eVals )
		end do
		
		call MixedState_boltzmanWeights( mixedStates_v0, T )
		call MixedState_boltzmanWeights( mixedStates_v1, T )
		
		open( 16, file="gEvsJ-Mixed.dat" )
		do i=1,nJ()
			write(16,"(F4.1)",advance='no') id2Ji(i)
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				
! 				write(16,"(F4.1,F10.5,A)") id2Ji(i), mixedStates_v0( id ).E, trim(mixedStates_v0(id).str())
				write(16,"(F15.5,F12.5)",advance='no') mixedStates_v0( id ).E, mixedStates_v1( id ).E
			end do
			
			write(16,"(A)") ""
		end do
		
		write(16,"(A)") ""
		write(16,"(A)") ""
		
		do i=1,nJ()
			write(16,"(A,F4.1)",advance='no') "#  ", id2Ji(i)
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				
				if( j==1 ) then
					write(16,"(5X,A)",advance='no') trim(mixedStates_v0(id).str())
				else
					write(16,"(A,9X,A)",advance='no') "#  ", trim(mixedStates_v0(id).str())
				end if
				
				write(16,"(A)") ""
			end do
			
! 			write(16,"(A)") ""
		end do
		
		close(16)
		
		
		write(*,*) "MIXED STATES ( nu = 0 )"
		write(*,*) "========================"
		write(*,*) ""
		write(*,'(7X,2A10)') "Energy", "Weight"
		write(*,'(7X,2A10)') "------", "------"
		do i=1,nJ()
			write(*,"(A,F4.1)") "J = ", id2Ji(i)
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				
				write(*,'(I5,2X,2F10.5,4X)', advance="no") j, mixedStates_v0(id).E, mixedStates_v0(id).weight
				
				do k=1,mixedStates_v0(id).nStates
					if( k == 1 ) then
						write(*,'(F10.4,2X,40A)', advance="no") mixedStates_v0(id).coeff(k), trim(mixedStates_v0(id).pureStates(k).str())
					else
						write(*,'(31X,F10.4,2X,40A)', advance="no") mixedStates_v0(id).coeff(k), trim(mixedStates_v0(id).pureStates(k).str())
					end if
					
					if( trim(mixedStates_v0(id).pureStates(k).pID) == trim(idUpperBase) .and. k == mixedStates_v0(id).nStates ) then
						testMaxT = mixedStates_v0(id).weight/maxval(mixedStates_v0.weight)/(2.0_8*mixedStates_v0(id).J()+1.0_8)
						write(*,'(F20.5)') testMaxT
						
						if( testMaxT > 0.1_8 ) then
							write(0,*) "## ERROR ## Maximum temperature available have been exceeded"
! 							stop
						end if
					else
						write(*,*) ""
					end if
				end do
			end do
			
			write(*,*) ""
		end do
		
		write(*,*) "MIXED STATES ( nu = 1 )"
		write(*,*) "========================"
		write(*,*) ""
		write(*,'(7X,2A10)') "Energy", "Weight"
		write(*,'(7X,2A10)') "------", "------"
		do i=1,nJ()
			write(*,"(A,F4.1)") "J = ", id2Ji(i)
			
			do j=1,Jvec(i)
				id = Jvec2StatesIDs( i, j )
				
				write(*,'(I5,2X,2F10.5,4X)', advance="no") j, mixedStates_v1(id).E, mixedStates_v1(id).weight
				
				do k=1,mixedStates_v1(id).nStates
					if( k == 1 ) then
						write(*,'(F10.4,2X,40A)', advance="no") mixedStates_v1(id).coeff(k), trim(mixedStates_v1(id).pureStates(k).str())
					else
						write(*,'(31X,F10.4,2X,40A)', advance="no") mixedStates_v1(id).coeff(k), trim(mixedStates_v1(id).pureStates(k).str())
					end if
					
					if( trim(mixedStates_v1(id).pureStates(k).pID) == trim(idUpperBase) .and. k == mixedStates_v1(id).nStates ) then
						testMaxT = mixedStates_v1(id).weight/maxval(mixedStates_v1.weight)/(2.0_8*mixedStates_v1(id).J()+1.0_8)
						write(*,'(F20.5)') testMaxT
						
						if( testMaxT > 0.1_8 ) then
							write(0,*) "## ERROR ## Maximum temperature available have been exceeded"
! 							stop
						end if
					else
						write(*,*) ""
					end if
				end do
			end do
			
			write(*,*) ""
		end do
		
	end subroutine diagHMatrix
	
end program HenX2Spec
