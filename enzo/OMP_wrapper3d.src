      subroutine wrapper3d(x, rank, n1, n2, n3, dir, method)

!     All stride-1, cache-friendly permuted wrapper
!     Written by: Robert Harkness
!     Date:       November, 2003

      implicit none

!     Arguments

      integer :: rank, n1, n2, n3, dir
      complex :: x(n1,n2,n3)
      external :: method

!     Locals

      complex, allocatable :: y(:,:,:), z(:,:,:)
      integer :: n(3)
      integer :: i,j,k

      if( rank /= 3 ) then
        write(0,*) '3D wrapper rank != 3'
        call stop_all_cpus
      end if

      n(1) = n1
      n(2) = 1
      n(3) = 1

!$omp parallel
!$omp-  shared(x, n, dir, n1, n2, n3)
!$omp-  private(j, k)
!$omp-  default(none)
!$omp do
      do k=1,n3
      do j=1,n2
      call fftwrap3d( x(1,j,k), n, dir, method )
      end do
      end do
!$omp end do
!$omp end parallel

      allocate( y(n2,n3,n1) )

      call rotate3d(x,n1,n2,n3,y)

      n(1) = n2

!$omp parallel
!$omp-  shared(y, n, dir, n1, n2, n3)
!$omp-  private(i, k)
!$omp-  default(none)
!$omp do
      do i=1,n1
      do k=1,n3
      call fftwrap3d( y(1,k,i), n, dir, method )
      end do
      end do
!$omp end do
!$omp end parallel

      allocate( z(n3,n1,n2) )

      call rotate3d(y,n2,n3,n1,z)

      deallocate( y)

      n(1) = n3

!$omp parallel
!$omp-  shared(z, n, dir, n1, n2, n3)
!$omp-  private(i, j)
!$omp-  default(none)
!$omp do
      do j=1,n2
      do i=1,n1
      call fftwrap3d( z(1,i,j), n, dir, method )
      end do
      end do
!$omp end do
!$omp end parallel

      call rotate3d(z,n3,n1,n2,x)

      deallocate( z )

      return

      end


      subroutine fftwrap3d( a, n, dir, method )

      implicit none

      complex :: a(*)
      integer :: n(3)
      integer :: dir
      external :: method

      call method(a, n(1), dir)

      return

      end
