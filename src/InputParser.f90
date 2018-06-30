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

module InputParser_
	implicit none
	private
	
	public :: &
		InputParser_test, &
		InputParser_nItems
	
	type, public :: InputParser
		character(:), private, allocatable :: iFile
		integer, private :: numberOfLines
		character(1000), allocatable, private :: lines(:)
		
		integer :: nstates
		character(1000), allocatable :: states(:)
		
		integer :: nCouplingL
		character(1000), allocatable :: couplingL(:)
		
		real(8) :: Emin
		real(8) :: Emax
		real(8) :: Estep
		
		real(8) :: Re
		real(8) :: we
		real(8) :: wexe
		real(8) :: rMass
		real(8) :: polAlpha(4)
		real(8) :: polBeta(4)
		
		contains
			procedure :: init
			procedure :: destroy
			procedure, private :: load
			procedure, private :: parse
	end type InputParser
	
	contains
	
	subroutine init( this, iFile )
		class(InputParser) :: this
		character(*), intent(in) :: iFile
		
		this.iFile = iFile
		this.nstates = 0
		
		call this.load( "#" )
		call this.parse()
	end subroutine init
	
	subroutine destroy( this )
		class(InputParser) :: this
		
		deallocate( this.states )
	end subroutine destroy
	
	subroutine load( this, cComments )
		class(InputParser) :: this
		character(*), optional, intent(in) :: cComments
		
		integer :: iostat
		integer :: posComment
		integer :: i
		character(999) :: buffer
		character(:), allocatable :: line
		character(999) :: fileContent(1000)
		
		open( 10, file=this.iFile, status="old", iostat=iostat )
		
		if( iostat /= 0 ) then
			write(*, *) "### Error ###: The file ( ", this.iFile, " ) cannot be open"
			stop
		else
			this.numberOfLines = 1
			iostat = 1
			do while( iostat /= -1 )
				read(10,'(A)', iostat=iostat) buffer
				
				if( present(cComments) ) then
					posComment = index( buffer, cComments )
					
					if( posComment /= 0 ) then
						line = trim(buffer(1:posComment-1))
					else
						line = trim(buffer)
					end if
				else
					line = trim(buffer)
				end if
				
				if( len_trim(line) /= 0 ) then
					fileContent( this.numberOfLines ) = line
					this.numberOfLines = this.numberOfLines + 1
				end if
				
			end do
		end if
		
		this.numberOfLines = this.numberOfLines - 1
		allocate( this.lines( this.numberOfLines ) )
		
		do i=1,this.numberOfLines
			this.lines(i) = trim( fileContent(i) )
		end do
		
		close( 10 )
	end subroutine load
	
	subroutine parse( this )
		class(InputParser) :: this
		
		integer :: i, j
		character(1000) :: buffer
		real(8) :: rbuffer
		integer :: ssize
		integer :: iostat
		
		write(*,*) "Input file data"
		write(*,*) "==============="
		write(*,*) ""
		
		i=1
		do while( i <= this.numberOfLines )
			if( index( this.lines(i), "STATES" ) /= 0 ) then
				read( this.lines(i), * ) buffer, ssize
				i=i+1
				
				this.nstates = ssize
				allocate( this.states(ssize) )
				
				write(*,*) "STATES from input file ( ", trim(this.iFile), " )"
				write(*,*) ""
				
				do j=1,ssize
					this.states(j) = trim(this.lines(i))
					write(*,*) j, "      ", trim(this.states(j))
					i=i+1
				end do
				
				write(*,*) ""
			else if( index( this.lines(i), "LCOUPLINGS" ) /= 0 ) then
				read( this.lines(i), * ) buffer, ssize
				i=i+1
				
				this.nCouplingL = ssize
				allocate( this.couplingL(ssize) )
				
				write(*,*) "L couplings from input file ( ", trim(this.iFile), " )"
				write(*,*) ""
				
				do j=1,ssize
					this.couplingL(j) = trim(this.lines(i))
					write(*,*) j, "      ", trim(this.couplingL(j))
					i=i+1
				end do
				
				write(*,*) ""
			else if( index( this.lines(i), "SPECTRUM" ) /= 0 ) then
				write(*,*) "Spectrum parameters from input file ( ", trim(this.iFile), " )"
				write(*,*) ""
				i = i+1
				
				do j=1,3
					read( this.lines(i), * ) buffer, rbuffer
					i=i+1
					
					if( index( buffer, "EMIN" ) /= 0 ) then
						this.Emin = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  Emin", this.Emin
					else if( index( buffer, "EMAX" ) /= 0 ) then
						this.Emax = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  Emax", this.Emax
					else if( index( buffer, "STEP" ) /= 0 ) then
						this.Estep = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  Estep", this.Estep
					end if
				end do
				
				write(*,*) ""
			else if( index( this.lines(i), "DIATOMIC" ) /= 0 ) then
				write(*,*) "Diatomic molecule parameters from input file ( ", trim(this.iFile), " )"
				write(*,*) ""
				i = i+1
				
				do j=1,6
					read( this.lines(i), * ) buffer, rbuffer
					i=i+1
					
					if( index( buffer, "RE" ) /= 0 ) then
						this.Re = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  Re", this.Re
					else if( index( buffer, "WE " ) /= 0 ) then
						this.we = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  we", this.we
					else if( index( buffer, "WEXE" ) /= 0 ) then
						this.wexe = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  wexe", this.wexe
					else if( index( buffer, "RMASS" ) /= 0 ) then
						this.rMass = rbuffer
						write(*,"(I12,A12,F16.6)") j, "  rmass", this.rMass
					else if( index( buffer, "POLA" ) /= 0 ) then
						read( this.lines(i-1), * ) buffer, this.polAlpha
						write(*,"(I12,A12,4F16.6)") j, "  polAlpha", this.polAlpha
					else if( index( buffer, "POLB" ) /= 0 ) then
						read( this.lines(i-1), * ) buffer, this.polBeta
						write(*,"(I12,A12,4F16.6)") j, "  polBeta", this.polBeta
					end if
				end do
				
				write(*,*) ""
			else
				i = i+1
			end if
		end do
	end subroutine parse
	
	function InputParser_nItems( line ) result( n )
		character(*) :: line
		integer :: n
		
		integer :: i
		integer :: maxItems = 1000
		character(1000) :: buffer
		
		read (line, *, end=999) ( buffer, i = 0, maxItems )
		999 continue
		
		n = i
	end function InputParser_nItems
	
	subroutine InputParser_test()
		type(InputParser) :: iparser
		integer :: i
		
		call iparser.init( "fermionesN3.inp" )
		
		do i=1,iparser.nstates
			write(*,*) trim(iparser.states(i))
		end do
		write(*,*) ""
		
		do i=1,iparser.nCouplingL
			write(*,*) trim(iparser.couplingL(i))
		end do
		write(*,*) ""
		
		call iparser.destroy()
		
	end subroutine InputParser_test
end module InputParser_