module sim_m

    use vehicle_m
    use jsonx_m
    use micro_time_m
    use connection_m

    implicit none
    
    type(aircraft), dimension(:), allocatable :: vehicles

    type(json_value), pointer :: j_main

    type(connection) :: graphics, controls_conn
    real :: s(13)

    ! NOTES FROM DR. HUNSAKERS CODE
    

contains


    subroutine init(filename)

        implicit none
        
        character(len=100), intent(in) :: filename
        type(json_value), pointer :: j_aircraft, j_initial, p1, p2
        real, dimension(:), allocatable :: euler
        real :: alpha, beta, p, q, r
        real :: da, de, dr, throttle
        integer :: N, cnt, i
        real :: x(2)

        ! Load JSON file
        call jsonx_load(filename, j_main)

        ! Initialize vehicles
        call jsonx_get(j_main, "vehicles", p1)
        N = json_value_count(p1)
        cnt= 0
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            cnt = cnt+1
        end do
        allocate(vehicles(cnt))
        cnt = 1
        do i = 1, N
            call json_value_get(p1, i, p2)
            if (p2%name(1:1) == 'x') cycle
            call vehicles(i)%init(p2)
            cnt = cnt+1
        end do

        call jsonx_get(j_main, "simulation.rk4_verbose", verbose, default_value=.false.)
        call jsonx_get(j_main, "simulation.save_states", save_states, default_value=.false.)

        call udp_initialize()         ! for windows users

        x = newtons_method(2, (/0.0, 0.0/), 0.01, 1.0, 1.0e-12)

        write(*,*) "SOLVED FOR X: ", x
        
    end subroutine init

    subroutine run()

        implicit none

        real :: t, y(13), y_temp(13)
        integer :: io_unit, i
        real :: dt, tf
        logical :: t_R, crashed
        real :: t_p, t_c, start_time, end_time
        type(json_value), pointer :: j_graphics, j_controls
        character(:), allocatable :: temp


        call jsonx_get(j_main, "simulation.time_step[s]", dt, default_value=0.01)
        call jsonx_get(j_main, "simulation.end_time[s]", tf)

        !call jsonx_get(j_main, "graphics", j_graphics)
        !call graphics%init(j_graphics, 13)

        !call jsonx_get(j_main, "controls", j_controls)
        !call controls_conn%init(j_controls, 4)


        !do i = 1, size(vehicles)
            !call vehicles(i)%print_aero_table(norm2(y(1:3)), -y(9))
        !end do

        !open(newunit=io_unit, file='sim_output.csv', status='replace', action='write')
        !write(io_unit,*) 'time[s], u[ft/s], v[ft/s], w[ft/s], p[rad/s], q[rad/s], r[rad/s], xf[ft], yf[ft], zf[ft], e0, ex, ey, ez'
        !write(io_unit,'(ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12, &
                        !A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') &
                            !t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            !,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

        !call quat_norm(y(10:13))

        ! Check for real time sim
        if (dt == 0.0) then
            t_R = .true.
        else
            t_R = .false.
        end if

        if (t_R) then
            y_temp = y
            t_p = get_time()
            do i = 1, size(vehicles)
                y = vehicles(i)%tick_states(dt)
            end do
            t_c = get_time()
            dt = t_c - t_p
            y = y_temp
            t_p = t_c
        end if



        start_time = get_time()
        ! Run simulation
        do while (t < tf)

            if (verbose) then
                write(*,*) 
                write(*,*) 
                write(*,*) "------------------------"
                write(*,*) "time: ", t, " s"
                write(*,*) "------------------------"
                write(*,*) 
            end if

            do i = 1, size(vehicles)
                y = vehicles(i)%tick_states(dt)
            end do
            !y = runge_kutta(t, y, dt)
            t = t + dt
            if (t_R) then
                t_c = get_time()
                dt = t_c - t_p
                t_p = t_c
            end if
            if (verbose) then
                write(*,*) "----------------------"
                write(*,*) "----------------------"
            end if


            write(io_unit,'(ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,&
                            ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') &
                            t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',' &
                            ,y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13)

            s(1) = t
            s(2:10) = y(1:9)
            s(11:13) = quat_to_euler(y(10:13))

            !crashed = vehicle%check_collision(t, y)
            !if (crashed) exit

            !call graphics%send(s)
            !controls = controls_conn%recv()

        end do
        end_time = get_time()

        close(io_unit)

        write(*,*) "Simulation Finished in ", end_time - start_time, " seconds."
        
        if (t_R) then
            write(*,*) "Real time simulation error: ", tf - (end_time - start_time), " seconds."
        end if

        call udp_finalize()

    end subroutine run 

end module sim_m