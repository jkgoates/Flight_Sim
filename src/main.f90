program main
    use simulation_m

    implicit none

    character(len=100) :: input_file

    ! Get input file from command line
    call get_command_argument(1, input_file)

    call init(input_file)

    call run()

end program main