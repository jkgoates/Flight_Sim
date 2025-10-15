program main
    use simulation_m

    implicit none

    character(len=100) :: input_file
    real :: alt(11)
    real :: Z, rho, P, T, a
    integer :: i


    ! Get input file from command line
    call get_command_argument(1, input_file)

    call init_5(input_file)

    call run()

end program main