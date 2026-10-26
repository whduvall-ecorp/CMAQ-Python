! Test adapter only; default REAL matches the unchanged I/O API POLY interface.
program poly_driver
    implicit none
    real, external :: POLY
    character(len=12) :: operation
    integer :: status, ndeg, npts
    real :: xpt, result
    real, allocatable :: xpts(:), ypts(:)
    do
        read(*, *, iostat=status) operation
        if (status < 0) exit
        if (status /= 0) stop 1
        select case (trim(operation))
        case ('POLY')
            read(*, *) ndeg, xpt
            npts = ndeg + 1
            allocate(xpts(npts), ypts(npts))
            read(*, *) xpts
            read(*, *) ypts
            result = POLY(xpt, xpts, ypts, ndeg)
            write(*, '(ES26.17E3)') result
            deallocate(xpts, ypts)
        case default
            stop 2
        end select
    end do
end program poly_driver
