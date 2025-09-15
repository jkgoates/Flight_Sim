program main
    use simulation_m

    implicit none

    character(len=100) :: input_file
    type(json_value), pointer :: j_main
    real :: dt, V, H, theta


    real :: y(13)

    ! Get input file from command line
    call get_command_argument(1, input_file)

    ! Load JSON file
    call jsonx_load(input_file, j_main)

    call run(j_main)

    
end program main