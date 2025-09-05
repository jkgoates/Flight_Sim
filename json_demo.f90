program main
    use jsonx_m
    implicit none

    character(100) :: filename
    type(json_value), pointer :: j_main
    real :: dt

    call get_command_argument(1,filename)

    call jsonx_load(filename,j_main)

    call jsonx_get(j_main, 'simulation.time_step[sec]', dt, 0.0)
    write(*,*) 'time step = ',dt

    ! call jsonx_get(j_main, 'simulation', j_sim)
    ! call jsonx_get(j_sim, 'time_step[sec]', dt, 0.0)
    ! write(*,*) 'time step = ',dt



end program main
