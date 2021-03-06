      SUBROUTINE INTERP_AP( ya, na, xa, yb, nb, xb )

      ! Linear interpolation xa (ya, na) >> xb (yb, nb)
      ! ya, yb == pressure levels
      ! NOTE: yb:bottom->top ; ya: top->bottom

      integer :: na, nb
      integer :: i, j
      real*4, dimension(na) :: xa, ya
      real*4, dimension(nb) :: xb, yb
      real*4 :: slope, biais

      do j = 1, nb
         do i = 1, na-1
            if (( yb(j) .ge. ya(i)) .and. ( yb(j) .lt. ya(i+1))) then
               if ( (xa(i) -xa(i+1)) .ne. 0.) then
                  slope = ( ya(i) - ya (i+1) ) / (xa(i) -xa(i+1))
                  biais = ya(i) - slope * xa(i)
                  xb(j) = ( yb(j) - biais) / slope
               else
                  xb(j) = xa(i)
               endif
            endif
         enddo
      enddo


      !Return to calling program
      END SUBROUTINE INTERP_AP
