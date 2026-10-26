! Test adapter only; default REAL matches the unchanged PREPLM interface.
program preplm_driver
    implicit none
    interface
        subroutine PREPLM(FIREFLG, EMLAYS, HMIX, HTS, PSFC, TS, DDZF, QV, &
                TA, UW, VW, ZH, ZF, PRES, LSTK, LPBL, TSTK, WSTK, DTHDZ, WSPD)
            logical, intent(in) :: FIREFLG
            integer, intent(in) :: EMLAYS
            real, intent(in) :: HMIX, HTS, PSFC, TS
            real, intent(in) :: DDZF(:), QV(:), TA(:), UW(:), VW(:), ZH(:), ZF(:)
            real, intent(in) :: PRES(0:)
            integer, intent(out) :: LSTK, LPBL
            real, intent(out) :: TSTK, WSTK
            real, intent(out) :: DTHDZ(:), WSPD(:)
        end subroutine PREPLM
    end interface
    character(len=12) :: operation
    integer :: status, emlays, lstk, lpbl, fire_flag, l
    real :: hmix, hts, psfc, ts, tstk, wstk
    real, allocatable :: ddzf(:), qv(:), ta(:), uw(:), vw(:), zh(:), zf(:)
    real, allocatable :: pres(:), dthdz(:), wspd(:)
    logical :: fireflg
    do
        read(*, *, iostat=status) operation
        if (status < 0) exit
        if (status /= 0) stop 1
        select case (trim(operation))
        case ('PREPLM')
            read(*, *) fire_flag, emlays, hmix, hts, psfc, ts
            allocate(ddzf(emlays), qv(emlays), ta(emlays), uw(emlays), vw(emlays), &
                zh(emlays), zf(emlays), pres(0:emlays), dthdz(emlays), wspd(emlays))
            read(*, *) ddzf
            read(*, *) qv
            read(*, *) ta
            read(*, *) uw
            read(*, *) vw
            read(*, *) zh
            read(*, *) zf
            read(*, *) (pres(l), l = 0, emlays)
            fireflg = fire_flag == 1
            call PREPLM(fireflg, emlays, hmix, hts, psfc, ts, ddzf, qv, ta, uw, vw, &
                zh, zf, pres, lstk, lpbl, tstk, wstk, dthdz, wspd)
            write(*, '(2(I8,1X))') lstk, lpbl
            write(*, '(2(ES26.17E3,1X))') tstk, wstk
            write(*, '(*(ES26.17E3,1X))') dthdz
            write(*, '(*(ES26.17E3,1X))') wspd
            deallocate(ddzf, qv, ta, uw, vw, zh, zf, pres, dthdz, wspd)
        case default
            stop 2
        end select
    end do
end program preplm_driver
