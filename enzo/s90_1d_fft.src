      subroutine s90_st1( a, n, dir )

      use singleton

      implicit none

      integer :: n, dir
      COMPLEX(fftkind) :: a(n)

      integer :: i
      integer :: sn(7)

      sn(1) = n

      if( dir == -1 ) then
        call fftn(a, shape(sn), inv=.false. )
        do i = 1,n
        a(i) = a(i) * sqrt(real(n))
        end do
      else
        call fftn(a, shape(sn), inv=.true. )
        do i = 1,n
        a(i) = a(i) / sqrt(real(n))
        end do
      end if

      return
      end
