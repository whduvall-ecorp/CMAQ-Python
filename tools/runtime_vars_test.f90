! TEST-ONLY adapter: PLMRIS imports LOGDEV and no other runtime state.
! No production runtime module or physics is replaced by this module.
module RUNTIME_VARS
    implicit none
    integer :: LOGDEV
end module RUNTIME_VARS
