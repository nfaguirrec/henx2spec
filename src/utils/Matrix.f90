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

module Matrix_
	implicit none
	private
	
	public :: &
		Matrix_show, &
		Matrix_isSymmetric, &
		Matrix_eigen, &
		Matrix_inver
	
	contains
	
	subroutine Matrix_show( matrix, ncol )
		real(8), intent(in) :: matrix(:,:)
		integer, intent(in), optional :: ncol
		
		integer :: ncolEff
		
		integer :: i, j, k
		integer :: ssize
		integer :: upper
		
		if ( present(ncol) ) then
			ncolEff = ncol
		else
			ncolEff = 15
		end if
		
		ssize = size(matrix,dim=1)
		
		do i=1,ssize
			if( ssize /= ncolEff ) then
				do j=1,ceiling(1.0*ssize/ncolEff)
					
					if( j /= ceiling(1.0*ssize/ncolEff) ) then
						upper = ncolEff*j
					else
						upper = ncolEff*(j-1)+mod(ssize,ncolEff)
					end if
					
					write (6,"(I5,<upper>F10.4)") i, ( matrix(i,k), k=ncolEff*(j-1)+1,upper )
				end do
			else
				write (6,"(I5,<ssize>F10.4)") i, ( matrix(i,k), k=1,ssize )
			end if
		end do
	end subroutine Matrix_show
	
	character(1) function Matrix_isSymmetric( matrix )
		real(8), intent(in) :: matrix(:,:)
		
		integer :: i, j
		integer :: ssize
		
		ssize = size(matrix,dim=1)
		
		do i=1,ssize
			do j=i+1,ssize
				if( abs( matrix(i,j) - matrix(j,i) ) >= 1e-8 ) then
					Matrix_isSymmetric = 'F'
					return
				end if
			end do
		end do
		
		Matrix_isSymmetric = 'T'
	end function Matrix_isSymmetric
	
	subroutine Matrix_eigen( matrix, eVecs, eVals )
		real(8), intent(in) :: matrix(:,:)
		real(8), intent(inout) :: eVecs(:,:)
		real(8), intent(inout) :: eVals(:)

		real(8), allocatable :: workSpace(:)
		integer :: ssize, info
		
		ssize = size( matrix, dim=2 )
		eVecs = matrix
		
		allocate( workSpace( 3*ssize-1 ) )
		
		call dsyev( 'V', 'U', ssize, eVecs, ssize, eVals, workSpace, 3*ssize-1, info )
		
		if ( info /= 0 ) stop "Get Matrix_eigen values Matrix failed"
		
		deallocate( workSpace )
	end subroutine Matrix_eigen
	
	subroutine Matrix_inver( matrix, inv )
		real(8), intent(in) :: matrix(:,:)
		real(8), intent(inout) :: inv(:,:)
		
                integer, allocatable :: pivotInd(:)
                real(8), allocatable :: workSpace(:)
		real(8), allocatable :: test(:,:)
		integer :: i, ssize, info
		real(8) :: trace
                
                ssize = size( matrix, dim=1 )
		inv = matrix
                
                allocate( pivotInd( ssize ) )
		
		!! Factorizacion LU
		call dgetrf( ssize, ssize, inv, ssize, pivotInd, info )
                if ( info /= 0 ) stop "Get Matrix LU factorization failed"
		
                allocate( workSpace( ssize ) )
		
                !! Invierte la matriz
                call dgetri( ssize, inv, ssize, pivotInd, workSpace, ssize, info )
                if ( info /= 0 ) stop "Get Inverse Matrix failed"
		
                allocate( test( ssize, ssize ) )
		test = matmul(matrix,inv)
		
		trace = 0.0_8
		do i=1,ssize
			trace = trace + test(i,i)
		end do
		
		if( abs( sum(test)-trace ) > 1e-10 ) stop "Get Inverse Matrix failed 2 A*A^-1"
		
                !! libera memoria 
                deallocate(workSpace)
                deallocate(pivotInd)
		deallocate(test)
	end subroutine Matrix_inver
	
end module Matrix_

