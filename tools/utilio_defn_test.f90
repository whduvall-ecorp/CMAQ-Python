! Test adapter only. Minimal UTILIO_DEFN providing the symbols unchanged
! PREPLM references: the external I/O API POLY (via interface) and M3MSG2.
! No physics is stubbed; POLY is the unchanged reference/ioapi/poly.f.
module UTILIO_DEFN
    implicit none
    interface
        real function POLY(XPT, XPTS, YPTS, NDEG)
            integer, intent(in) :: NDEG
            real, intent(in) :: XPT
            real, intent(in) :: XPTS(NDEG + 1)
            real, intent(in) :: YPTS(NDEG + 1)
        end function POLY
    end interface
contains
    subroutine M3MSG2(MSG)
        character(len=*), intent(in) :: MSG
        write(*, *) trim(MSG)
    end subroutine M3MSG2
end module UTILIO_DEFN
