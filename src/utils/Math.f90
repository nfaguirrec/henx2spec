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

module Math_
	implicit none
	
	real(8), public, parameter :: Math_PI = acos(-1.0_8)
	
	public :: &
		Math_wigner3j, &
		Math_sort, &
		Math_test
	
	contains
	
	!**
	! This function calculates the value of a Wigner-3j symbol.
	! The equation used is equation 1.5 in "The 3j and 6j symbols" by
	! Manuel Rotenber, R. Bivins, N.Metropolis  and John K. Wooten, JR
	! 1959, The Technology Press, Massachusett.
	!
	! The function has been tested for all sets, where j1d=j2d=5/2
	! and j3d={0, 2, 4}, and m1d+m2d+m3d=0
	! It has been tested to give zero for some vanishing terms
	! If you should want to use it for higher js than the ones I have tested
	! you might have to do something to avoid overflow in the expression with
	! the many factorials
	!**
	function Math_wigner3j( j1, j2, j3, m1, m2, m3 ) result (w3j)
		real(8), intent(in) :: j1,j2,j3,m1,m2,m3
		real(8) :: w3j
		
		integer:: i
		real(8) :: j1d,j2d,j3d,m1d,m2d,m3d
		real(8) :: k,kmin,kmax,ksum,divisor,sign
		real(8), dimension(0:200) ::lnfac
		
		j1d = real(j1,8)
		j2d = real(j2,8)
		j3d = real(j3,8)
		m1d = real(m1,8)
		m2d = real(m2,8)
		m3d = real(m3,8)
		
		lnfac(0)=0.d0
		lnfac(1)=0.d0
		do i=2,200
			lnfac(i)=lnfac(i-1)+log(REAL(i))
		end do
		
		kmax = min(j1d+j2d-j3d, j1d-m1d,j2d+m2d)
		kmin = max(0.d0, j2d-j3d-m1d, j1d-j3d+m2d)
		
		ksum=0.
		do k = kmin, kmax
			divisor=exp(lnfac(nint(k))+lnfac(nint(j1d+j2d-j3d-k))+ &
				lnfac(nint(j1d-m1d-k))+lnfac(nint(j2d+m2d-k))+ &
				lnfac(nint(j3d-j2d+m1d+k))+lnfac(nint(j3d-j1d-m2d+k)))
			ksum=ksum+minusonetothe(int(k))/divisor
		end do !k
		
		w3j= minusonetothe(nint(j1d-j2d-m3d))*sqrt(  &
				exp(lnfac(nint(j1d+j2d-j3d))+lnfac(nint(j1d-j2d+j3d))+&
				lnfac(nint(-j1d+j2d+j3d))+&
				lnfac(nint(j1d+m1d))+lnfac(nint(j1d-m1d))+&
				lnfac(nint(j2d+m2d))+lnfac(nint(j2d-m2d))+&
				lnfac(nint(j3d+m3d))+lnfac(nint(j3d-m3d))-&
				lnfac(nint(j1d+j2d+j3d+1))))*ksum
	end function Math_wigner3j
	
	function minusonetothe(j) result (sgn)
		integer :: i,j
		real*8 :: sgn
		
		i=j
		if (j<=0) i=-j
		if (mod(i,2)==0) then
			sgn=1
		else
			sgn=-1
		endif
	end function minusonetothe
	
	!**
	! indexing array so that array(indexes(j)), j=1..n is in
	! ascending numerical order.
	! method is heapsort, see also subroutine hpsort.
	! taken from numerical recipies, p 233.
	!**
	subroutine Math_sort( array, indexes )
		real(8), allocatable, intent(in) :: array(:)
		integer, allocatable, intent(inout) :: indexes(:)
		
		integer :: i, j, l
		integer :: n
		integer :: id, ir
		real(8) :: value
		
		if( .not. allocated(array) ) then
			write(*,*) "Error in Math_sort, array not allocated"
			stop
		end if
		
		if( .not. allocated(indexes) ) then
			write(*,*) "Error in Math_sort, indexes not allocated"
			stop
		end if
		
		if( size(array) /= size(indexes) ) then
			write(*,*) "Error in Math_sort, array and indexes have different size"
			stop
		end if
		
		n = size(array)
		
		do j=1,n
			indexes(j)=j
		end do
		
		if( n == 1 ) return
		
		l=n/2+1
		ir=n
		
		do while( .true. )
			if( l > 1 ) then
				l = l-1
				id = indexes(l)
				value = array(id)
			else
				id=indexes(ir)
				value=array(id)
				indexes(ir)=indexes(1)
				ir=ir-1
				
				if(ir == 1) then
					indexes(1)=id
					return
				end if
			end if
			
			i = l
			j = 2*l
			
			do while( j <= ir )
				if( j < ir ) then
					if( array(indexes(j)) < array(indexes(j+1)) ) then
						j=j+1
					end if
				end if
				
				if( value < array(indexes(j)) ) then
					indexes(i)=indexes(j)
					i=j
					j=2*j
				else
					j=ir+1
				end if
			end do
			
			indexes(i)=id
		end do
	end subroutine Math_sort
	
	subroutine Math_test()
		real(8), allocatable :: array(:)
		integer, allocatable :: indexes(:)
		
		integer :: i, j, l
		
		write(*,*) ""
		write(*,*) "Sorting vectors"
		write(*,*) "==============="
		allocate( array(19) )
		allocate( indexes(19) )
		
		array = [ -33.89007, -33.89007, -35.21007, -35.42677, -35.90699, &
			  -33.67374, -34.15396, -33.27902, -33.27902, -34.04338, &
			  -34.04338, -31.81726, -31.81726, -33.49329, -33.97351, &
			  -32.91849, -33.39871, -33.24224, -33.24224 ]
			  
		call Math_sort( array, indexes )
		
		write(*,"(5X,2A15)") "original", "sorted"
		write(*,"(5X,2A15)") "--------", "------"
		do i=1,size(array)
			write(*,"(I5,2F15.5)") i, array(i), array( indexes(i) )
		end do
		
		write(*,*) ""
		write(*,*) "3j symbols"
		write(*,*) "=========="
		
		write(*,*) "(3/2  0  3/2;    0  0    0) = ", Math_wigner3j( 1.5_8, 0.0_8, 1.5_8, 0.0_8, 0.0_8, 0.0_8 )
		write(*,*) "(3/2  0  3/2;  1/2  0 -1/2) = ", Math_wigner3j( 1.5_8, 0.0_8, 1.5_8, 0.5_8, 0.0_8,-0.5_8 )
		write(*,*) "(3/2  0  3/2; -1/2  0  1/2) = ", Math_wigner3j( 1.5_8, 0.0_8, 1.5_8,-0.5_8, 0.0_8, 0.5_8 )
		write(*,*) "(3/2  0  3/2;  3/2  0 -3/2) = ", Math_wigner3j( 1.5_8, 0.0_8, 1.5_8, 1.5_8, 0.0_8,-1.5_8 )
		write(*,*) "(3/2  0  3/2; -3/2  0  3/2) = ", Math_wigner3j( 1.5_8, 0.0_8, 1.5_8,-1.5_8, 0.0_8, 1.5_8 )
	end subroutine Math_test
	
end module Math_

