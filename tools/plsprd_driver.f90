! Test adapter only; default REAL matches the unchanged CMAQ interfaces.
program plsprd_driver
    use RUNTIME_VARS, only: LOGDEV
    implicit none
    interface
        subroutine PLSPRD(dthdz, zf, kz, cefstk, pltop, plbot)
            real, intent(in) :: dthdz(:), zf(0:), cefstk
            integer, intent(in) :: kz
            real, intent(out) :: pltop, plbot
        end subroutine PLSPRD
        subroutine DELTA_ZS(emlays, nsrc, source_map, stack, zf, zstk, ddzf)
            integer, intent(in) :: emlays, nsrc, source_map(:)
            real, intent(in) :: stack(:), zf(:,:)
            real, intent(out) :: zstk(:,:), ddzf(:,:)
        end subroutine DELTA_ZS
        subroutine PLMRIS(emlays, lstk, hfx, hmix, stkdm, stkht, stktk, stkve, &
                tstk, ustar, dthdz, ta, wspd, zf, zh, zstk, wstk, zplm)
            integer, intent(in) :: emlays, lstk
            real, intent(in) :: hfx, hmix, stkdm, stkht, stktk, stkve, tstk, ustar
            real, intent(in) :: dthdz(:), ta(:), wspd(:), zf(0:), zh(:), zstk(:)
            real, intent(inout) :: wstk
            real, intent(out) :: zplm
        end subroutine PLMRIS
        subroutine FIRE_PLMRIS(emlays, lstk, hfx, hmix, bflx, tstk, ustar, &
                dthdz, ta, wspd, zf, wstk, zplm)
            integer, intent(in) :: emlays, lstk
            real, intent(in) :: hfx, hmix, bflx, tstk, ustar
            real, intent(in) :: dthdz(:), ta(:), wspd(:), zf(0:)
            real, intent(inout) :: wstk
            real, intent(out) :: zplm
        end subroutine FIRE_PLMRIS
    end interface
    real :: gradients(3), heights(0:3), center, top, bottom
    integer :: status, nlayers, nsrc, nstack, lstk, before, after
    real :: hfx, hmix, dm, ht, tk, ve, ts, us, ws
    character(len=12) :: operation
    integer, allocatable :: source_map(:)
    real, allocatable :: stack(:), zf(:,:), zstk(:,:), ddzf(:,:)
    real, allocatable :: lapse(:), ta(:), wind(:), faces(:), zh(:)
    open(newunit=LOGDEV, file='plmris.log', status='replace', action='write', access='stream', form='formatted')
    heights = [0.0, 100.0, 500.0, 1000.0]
    do
        read(*, *, iostat=status) operation
        if (status < 0) exit
        if (status /= 0) stop 1
        select case (trim(operation))
        case ('PLSPRD')
            read(*, *) gradients, center
            call PLSPRD(gradients, heights, 3, center, top, bottom)
            write(*, '(2(ES26.17E3,1X))') top, bottom
        case ('DELTA_ZS')
            read(*, *) nlayers, nsrc, nstack
            allocate(source_map(nsrc), stack(nstack), zf(nlayers,nsrc), &
                zstk(nlayers,nsrc), ddzf(nlayers,nsrc))
            read(*, *) source_map
            read(*, *) stack
            read(*, *) zf
            call DELTA_ZS(nlayers, nsrc, source_map, stack, zf, zstk, ddzf)
            write(*, '(*(ES26.17E3,1X))') zstk, ddzf
            deallocate(source_map, stack, zf, zstk, ddzf)
        case ('PLMRIS', 'CONNECTED')
            read(*, *) nlayers, lstk
            read(*, *) hfx, hmix, dm, ht, tk, ve, ts, us, ws
            allocate(lapse(nlayers), ta(nlayers), wind(nlayers), faces(0:nlayers), zh(nlayers), &
                zstk(nlayers,1), ddzf(nlayers,1), zf(nlayers,1))
            read(*, *) lapse
            read(*, *) ta
            read(*, *) wind
            read(*, *) faces
            read(*, *) zh
            if (trim(operation) == 'CONNECTED') then
                zf(:,1) = faces(1:)
                call DELTA_ZS(nlayers, 1, [1], [ht], zf, zstk, ddzf)
            else
                read(*, *) zstk(:,1)
            end if
            inquire(unit=LOGDEV, pos=before)
            call PLMRIS(nlayers, lstk, hfx, hmix, dm, ht, tk, ve, ts, us, &
                lapse, ta, wind, faces, zh, zstk(:,1), ws, center)
            flush(LOGDEV)
            inquire(unit=LOGDEV, pos=after)
            if (trim(operation) == 'CONNECTED') then
                call PLSPRD(lapse, faces, nlayers, center, top, bottom)
                write(*, '(*(ES26.17E3,1X))') center, ws, real(merge(1,0,after>before)), &
                    zstk, ddzf, top, bottom
            else
                write(*, '(*(ES26.17E3,1X))') center, ws, real(merge(1,0,after>before))
            end if
            deallocate(lapse, ta, wind, faces, zh, zstk, ddzf, zf)
        case ('FIRE_PLMRIS')
            read(*, *) nlayers, lstk
            read(*, *) hfx, hmix, dm, tk, us, ws
            allocate(lapse(nlayers), ta(nlayers), wind(nlayers), faces(0:nlayers))
            read(*, *) lapse
            read(*, *) ta
            read(*, *) wind
            read(*, *) faces
            inquire(unit=LOGDEV, pos=before)
            call FIRE_PLMRIS(nlayers, lstk, hfx, hmix, dm, tk, us, &
                lapse, ta, wind, faces, ws, center)
            flush(LOGDEV)
            inquire(unit=LOGDEV, pos=after)
            write(*, '(*(ES26.17E3,1X))') center, ws, real(merge(1,0,after>before))
            deallocate(lapse, ta, wind, faces)
        case default
            stop 2
        end select
    end do
    close(LOGDEV)
end program plsprd_driver
