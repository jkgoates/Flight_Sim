program test
use connection_m, only: connection
use jsonx_m
! use udp_m, only: udp_initialize, udp_finalize         ! for windows users
implicit none

type(connection) :: ctrls, graphics
real :: s(14), t, vals(7), current, start
type(json_value), pointer :: p1, p2

! call udp_initialize()         ! for windows users

call jsonx_load('connect.json', p1)

call jsonx_get(p1, 'controls', p2)
call ctrls%init(p2, 7)

call jsonx_get(p1, 'graphics', p2)
call graphics%init(p2, 14)

s = 0.
t = 0.
call cpu_time(start)
do while (t < 20.)
    
    vals = ctrls%recv([t])
    
    s(1) = t
    s(8:9) = vals(:2)
    s(10) = -vals(3)
    s(11:) = vals(4:)
    
    call graphics%send(s)
    
    write(*,*) t
    call cpu_time(current)
    t = current - start
end do

! call udp_finalize()           ! for windows users

end program test
