program main
    use coupling_module
    implicit none
    type(coupling_calculator) :: calculator
    integer :: argc
    character(len=128) :: argv0
    character(len=128) :: argv(3)
    integer :: i

    argc = command_argument_count()

    call get_command_argument(0, argv0)
    if (argc == 3) then
        do i = 1, 3
            call get_command_argument(i, argv(i))
        end do
        calculator  = coupling_calculator(argv(1), argv(2), argv(3))
    else if (argc == 0) then
        calculator = coupling_calculator()
    else
        error stop "Usage: " // trim(argv0) // "prefix_dimer prefix_monomer1 prefix_monomer2"
    end if

    call calculator%read()
    call calculator%calculate()
    call calculator%show()

    stop
end program main

