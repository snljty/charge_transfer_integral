program main
    implicit none

    interface
        subroutine pause_program()
            implicit none
        end subroutine pause_program

        subroutine read_matrix(fl_unit, mat, mat_size)
            implicit none
            integer, intent(in) :: fl_unit
            double precision, intent(out) :: mat(:, :)
            integer, intent(in) :: mat_size
        end subroutine read_matrix

        subroutine read_lt_matrix(fl_unit, mat, mat_size)
            implicit none
            integer, intent(in) :: fl_unit
            double precision, intent(out) :: mat(:, :)
            integer, intent(in) :: mat_size
        end subroutine read_lt_matrix

        subroutine calc_coupling(index1, index2, mat_size, mat_block1_size, F, S, C_sep, &
            J_eff12, e_eff1, e_eff2)
            implicit none
            integer, intent(in) :: index1
            integer, intent(in) :: index2
            integer, intent(in) :: mat_size
            integer, intent(in) :: mat_block1_size
            double precision, intent(in) :: F(mat_size, mat_size)
            double precision, intent(in) :: S(mat_size, mat_size)
            double precision, intent(in) :: C_sep(mat_size, mat_size)
            double precision, intent(out) :: J_eff12
            double precision, intent(out) :: e_eff1
            double precision, intent(out) :: e_eff2
        end subroutine calc_coupling
    end interface

    integer :: argc
    integer :: iarg
    ! integer :: arg_status
    character(kind=1,len=128), allocatable :: argv(:)

    integer, parameter :: stdin_unit = 5
    integer, parameter :: stdout_unit = 6
    integer, parameter :: stderr_unit = 0

    integer, parameter :: exit_success = 0
    integer, parameter :: exit_failure = 1

    double precision, parameter :: h_Planck = 6.62607015D-34
    double precision, parameter :: epsilon_0 = 8.854187817D-12
    double precision, parameter :: q_e = 1.602176634D-19
    double precision, parameter :: m_e = 9.10938215D-31
    double precision, parameter :: Hartree_to_eV = m_e * q_e ** 3 / (4.0D0 * epsilon_0 ** 2 * h_Planck ** 2)

    character(kind=1,len=256) :: buf
    integer :: buf_pos

    character(kind=1,len=128) :: fl_MO_dimer
    character(kind=1,len=128) :: fl_MO_monomer1
    character(kind=1,len=128) :: fl_MO_monomer2
    character(kind=1,len=128) :: fl_CT_out
    integer, parameter :: ifl_unit = 10
    integer, parameter :: ofl_unit = 11

    integer :: num_ele_dimer
    integer :: num_ele_monomer1
    integer :: num_ele_monomer2
    integer :: num_orb_dimer
    integer :: num_orb_monomer1
    integer :: num_orb_monomer2

    integer :: num_ele_alpha ! temporary
    integer :: num_ele_beta  ! temporary

    integer :: index_homo1
    integer :: index_lumo1
    integer :: index_homo2
    integer :: index_lumo2

    double precision, allocatable :: S(:, :)
    double precision, allocatable :: F(:, :)
    double precision, allocatable :: C_sep(:, :)

    double precision :: J_eff12
    double precision :: e_eff1
    double precision :: e_eff2

    integer :: i
    integer :: j

    integer :: io_status

    ! get command arguments
    argc = command_argument_count()
    allocate(argv(0:argc))

    do iarg = 0, argc
        call get_command_argument(iarg, argv(iarg))
    end do

    ! pharse command arguments
    if (argc == 0 .or. argc == 1 .or. argc == 2 .or. argc > 4 ) then
        write(stderr_unit, "(a,4(1x,a))") "Usage:", trim(argv(0)), "dimer.out", "monomer1.out", "monomer2.out"
        write(stderr_unit, "(a)") "If a 4th argument is provided, then the results will be saved to this file."
        if (argc /= 0) then
            stop exit_failure
        else
            stop exit_success
        end if
    else if (argc == 3) then
        fl_CT_out = ""
    else if (argc == 4) then
        fl_CT_out = trim(argv(4))
    end if

    fl_MO_dimer = trim(argv(1))
    fl_MO_monomer1 = trim(argv(2))
    fl_MO_monomer2 = trim(argv(3))

    deallocate(argv)

    ! get number of electrons and orbitals
    ! note that the inner-shell electrons and corresponding nucleal charges are presented in the pseudopotential.

    ! monomer1
    open(ifl_unit, file = fl_MO_monomer1, status = "old", action = "read")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "alpha electrons") /= 0) then
            exit
        end if
    end do
    read(buf, *) num_ele_alpha
    buf = buf(index(buf, "alpha electrons") + len_trim("alpha electrons"):)
    read(buf, *) num_ele_beta
    if (num_ele_alpha /= num_ele_beta) then
        write(stderr_unit, "(a,i0,a,i0,a)") "Error! Number of alpha electrons (", num_ele_alpha, &
            ") and beta electrons (", num_ele_beta, ")in monomer1 mismatch, "
        write(stderr_unit, "(a)") "this does not seem to be a closed-shell system and this program cannot handle this."
        stop exit_failure
    end if
    num_ele_monomer1 = num_ele_alpha + num_ele_beta
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "NBsUse=") /= 0) then
            exit
        end if
    end do
    buf = buf(index(buf, "NBsUse=") + len_trim("NBsUse="):)
    read(buf, *) num_orb_monomer1
    close(ifl_unit)

    ! monomer2
    open(ifl_unit, file = fl_MO_monomer2, status = "old", action = "read")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "alpha electrons") /= 0) then
            exit
        end if
    end do
    read(buf, *) num_ele_alpha
    buf = buf(index(buf, "alpha electrons") + len_trim("alpha electrons"):)
    read(buf, *) num_ele_beta
    if (num_ele_alpha /= num_ele_beta) then
        write(stderr_unit, "(a,i0,a,i0,a)") "Error! Number of alpha electrons (", num_ele_alpha, &
            ") and beta electrons (", num_ele_beta, ")in monomer2 mismatch, "
        write(stderr_unit, "(a)") "this does not seem to be a closed-shell system and this program cannot handle this."
        stop exit_failure
    end if
    num_ele_monomer2 = num_ele_alpha + num_ele_beta
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "NBsUse=") /= 0) then
            exit
        end if
    end do
    buf = buf(index(buf, "NBsUse=") + len_trim("NBsUse="):)
    read(buf, *) num_orb_monomer2
    close(ifl_unit)

    ! dimer
    open(ifl_unit, file = fl_MO_dimer, status = "old", action = "read")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "alpha electrons") /= 0) then
            exit
        end if
    end do
    read(buf, *) num_ele_alpha
    buf = buf(index(buf, "alpha electrons") + len_trim("alpha electrons"):)
    read(buf, *) num_ele_beta
    if (num_ele_alpha /= num_ele_beta) then
        write(stderr_unit, "(a,i0,a,i0,a)") "Error! Number of alpha electrons (", num_ele_alpha, &
            ") and beta electrons (", num_ele_beta, ")in dimer mismatch, "
        write(stderr_unit, "(a)") "this does not seem to be a closed-shell system and this program cannot handle this."
        stop exit_failure
    end if
    num_ele_dimer = num_ele_alpha + num_ele_beta
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "NBsUse=") /= 0) then
            exit
        end if
    end do
    buf = buf(index(buf, "NBsUse=") + len_trim("NBsUse="):)
    read(buf, *) num_orb_dimer
    close(ifl_unit)
    num_ele_alpha = 0
    num_ele_beta = 0

    ! check if the sum of electrons and orbitals of monomers correspond to those in dimer
    if (num_ele_dimer /= num_ele_monomer1 + num_ele_monomer2) then
        write(stderr_unit, "(a,i0,a,i0,a,i0,a)") "Error! Sum of the number of electrons in monomers (", num_ele_monomer1, "+", &
            num_ele_monomer2, ") and the number of electrons in dimer (", num_ele_dimer, ") mismatch."
        stop exit_failure
    end if
    if (num_orb_dimer /= num_orb_monomer1 + num_orb_monomer2) then
        write(stderr_unit, "(a,i0,a,i0,a,i0,a)") "Error! Sum of the number of orbitals in monomers (", num_orb_monomer1, "+", &
            num_orb_monomer2, ") and the number of orbitals in dimer (", num_orb_dimer, ") mismatch."
        stop exit_failure
    end if

    ! allocate memory for matrices
    allocate(F(num_orb_dimer, num_orb_dimer), S(num_orb_dimer, num_orb_dimer), C_sep(num_orb_dimer, num_orb_dimer))

    ! read S (overlap integral matrix) of dimer in passing
    open(ifl_unit, file = fl_MO_dimer, status = "old", action = "read", position = "rewind")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        if (index(buf, "*** Overlap ***") /= 0) then
            exit
        end if
    end do
    call read_lt_matrix(ifl_unit, S, num_orb_dimer)
    close(ifl_unit)

    ! get the indices of HOMOs and LUMOs in monomers.
    index_homo1 = num_ele_monomer2 / 2
    index_homo2 = num_ele_monomer2 / 2
    index_lumo1 = index_homo1 + 1
    index_lumo2 = index_homo2 + 1

    ! read F (Fock matrix) of dimer
    open(ifl_unit, file = fl_MO_dimer, status = "old", action = "read", position = "append")
    do while (.true.)
        backspace(ifl_unit)
        read(ifl_unit, "(a)") buf
        if (index(buf, "Fock matrix (alpha)") /= 0) then
            exit
        end if
        backspace(ifl_unit)
    end do
    call read_lt_matrix(ifl_unit, F, num_orb_dimer)
    close(ifl_unit)

    ! initialize C_sep
    C_sep = 0d0

    ! read C of monomer1
    open(ifl_unit, file = fl_MO_monomer1, status = "old", action = "read", position = "append")
    do while (.true.)
        backspace(ifl_unit)
        read(ifl_unit, "(a)") buf
        if (index(buf, "Alpha MO coefficients") /= 0) then
            exit
        end if
        backspace(ifl_unit)
    end do
    call read_matrix(ifl_unit, C_sep(:num_orb_monomer1, :num_orb_monomer1), num_orb_monomer1)
    close(ifl_unit)

    ! read C of monomer2
    open(ifl_unit, file = fl_MO_monomer2, status = "old", action = "read", position = "append")
    do while (.true.)
        backspace(ifl_unit)
        read(ifl_unit, "(a)") buf
        if (index(buf, "Alpha MO coefficients") /= 0) then
            exit
        end if
        backspace(ifl_unit)
    end do
    call read_matrix(ifl_unit, C_sep(num_orb_monomer1 + 1:, num_orb_monomer1 + 1:), num_orb_monomer2)
    close(ifl_unit)

    ! calculate charge tranfer integrals
    if (trim(fl_CT_out) == "") then
        write(*, "(a)") "# Calculating transfer integrals ..."
        write(*, "(a,i4,a,i4)") "# Monomer 1    HOMO: ", index_homo1, "        LUMO: ", index_lumo1
        write(*, "(a,i4,a,i4)") "# Monomer 2    HOMO: ", index_homo2, "        LUMO: ", index_lumo2
        write(*, "(a)") "# Transfer Integrals between HOMO and LUMO of monomers:"
        write(*, "(a)") "# monomer_1    monomer_2    J_eff_12/meV"
        ! do i = index_homo1, index_lumo1
        !     do j = index_homo2, index_lumo2
        !         call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
        !             J_eff12, e_eff1, e_eff2)
        !         write(*, "(3x,i4,9x,i4,7x,f13.7)") i, j, J_eff12 * Hartree_to_eV * 1.0D3
        !     end do
        ! end do
        write(*, "(a,1x,a)") "#", "HOMO-HOMO:"
        call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
            J_eff12, e_eff1, e_eff2)
        write(*, "(3x,i4,9x,i4,7x,f13.7)") i, j, J_eff12 * Hartree_to_eV * 1.0D3
        i = index_lumo1
        j = index_lumo2
        write(*, "(a,1x,a)") "#", "LUMO-LUMO:"
        call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
            J_eff12, e_eff1, e_eff2)
        write(*, "(3x,i4,9x,i4,7x,f13.7)") i, j, J_eff12 * Hartree_to_eV * 1.0D3      
    else
        write(*, "(a)") "Calculating transfer integrals ..."
        open(ofl_unit, file = fl_CT_out, action = "write", status = "replace")
        write(ofl_unit, "(a,i4,a,i4)") "# Monomer 1    HOMO: ", index_homo1, "        LUMO: ", index_lumo1
        write(ofl_unit, "(a,i4,a,i4)") "# Monomer 2    HOMO: ", index_homo2, "        LUMO: ", index_lumo2
        write(ofl_unit, "(a)") "# monomer_1    monomer_2    J_eff_12/meV"
        do i = 1, num_orb_monomer1
            do j = 1, num_orb_monomer2
                call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
                    J_eff12, e_eff1, e_eff2)
                write(ofl_unit, "(3x,i4,9x,i4,7x,f13.7)") i, j, J_eff12 * Hartree_to_eV * 1.0D3
            end do
        end do
        close(ofl_unit)
        write(*, "(a,a,a)") "Intermolecule charge transfer integrals have been saved to """, &
            trim(fl_CT_out), """."
    end if

    ! deallocate memory for matrices
    deallocate(F, S, C_sep)

    if (argc == 0) then
        call pause_program()
    end if

    stop
end program main

subroutine pause_program()
    implicit none

    write(*, "(a)") "Press <Enter> to exit."
    read(*, *)

    return
end subroutine pause_program

function get_num_items_in_line(line, len_line)
    implicit none
    integer, intent(in) :: len_line
    character(kind=1,len=len_line), intent(in) :: line
    integer :: get_num_items_in_line
    integer :: pos_read_in_line

    pos_read_in_line = len_trim(line)
    get_num_items_in_line = 0
    do while (pos_read_in_line /= 0)
        get_num_items_in_line = get_num_items_in_line + 1
        pos_read_in_line = scan(line(:pos_read_in_line), " ", back = .true.)
        if (pos_read_in_line == 0) exit
        pos_read_in_line = len_trim(line(:pos_read_in_line))
    end do

    return
end function get_num_items_in_line

! Read a ofmated matrix
! This subroutine changes the current position of the file.
! The matrix is allocated before, and
! mat_size should be provided before here. There is way to overcome this, 
! but not provided in this code.
subroutine read_matrix(fl_unit, mat, mat_size)
    implicit none
    integer, intent(in) :: fl_unit
    double precision, intent(out) :: mat(:, :)
    integer, intent(in) :: mat_size
    integer, external :: get_num_items_in_line

    character(kind=1,len=256) :: buf
    integer :: col_in_this_group
    integer :: pos_read_in_line
    integer :: i
    integer :: j

    j = 1
    do while (.true.)
        read(fl_unit, "(a)") buf
        col_in_this_group = get_num_items_in_line(buf, len_trim(buf))
        do i = 1, mat_size
            read(fl_unit, "(a)") buf
            pos_read_in_line = verify(buf, " ")
            pos_read_in_line = pos_read_in_line + scan(buf(pos_read_in_line:), " ") - 1
            read(buf(pos_read_in_line:), *) mat(i, j: j + col_in_this_group - 1)
        end do
        j = j + col_in_this_group
        if (j > mat_size) exit
    end do

    return
end subroutine read_matrix

! Read a lower-triangle ofmat matrix
! This subroutine changes the current position of the file.
! The matrix is allocated before, and
! mat_size should be provided before here. There is way to overcome this, 
! but not provided in this code.
subroutine read_lt_matrix(fl_unit, mat, mat_size)
    implicit none
    integer, intent(in) :: fl_unit
    double precision, intent(out) :: mat(:, :)
    integer, intent(in) :: mat_size
    integer, external :: get_num_items_in_line
    integer :: max_col_use
    character(kind=1,len=256) :: buf
    integer :: i
    integer :: j
    integer :: k
    integer :: pos_read_in_line
    integer :: pos_next
    integer :: start_index
    integer :: num_groups

    read(fl_unit, "(a)") buf
    backspace(fl_unit)
    max_col_use = get_num_items_in_line(buf, len_trim(buf))

    start_index = 1
    num_groups = (mat_size - 1) / max_col_use + 1
    mat = 0

    do k = 1, num_groups
        read(fl_unit, "(a)") buf
        do i = start_index, mat_size
            read(fl_unit, "(a)") buf
            j = start_index
            pos_read_in_line = verify(buf, " ")
            do while (.true.)
                pos_read_in_line = pos_read_in_line + scan(buf(pos_read_in_line:), " ") - 1
                pos_next = verify(buf(pos_read_in_line:), " ")
                if (pos_next == 0) exit
                pos_read_in_line = pos_read_in_line + pos_next - 1
                read(buf(pos_read_in_line:), *) mat(i, j)
                j = j + 1
            end do
        end do
        start_index = start_index + max_col_use
    end do

    do i = 1, mat_size
        do j = 1, i - 1
            mat(j, i) = mat(i, j)
        end do
    end do

    return
end subroutine read_lt_matrix

subroutine calc_coupling(index1, index2, mat_size, mat_block1_size, F, S, C_sep, J_eff12, e_eff1, e_eff2)
    implicit none
    integer, intent(in) :: index1
    integer, intent(in) :: index2
    integer, intent(in) :: mat_size
    integer, intent(in) :: mat_block1_size
    double precision, intent(in) :: F(mat_size, mat_size)
    double precision, intent(in) :: S(mat_size, mat_size)
    double precision, intent(in) :: C_sep(mat_size, mat_size)
    double precision, intent(out) :: J_eff12
    double precision, intent(out) :: e_eff1
    double precision, intent(out) :: e_eff2
    double precision :: e1
    double precision :: e2
    double precision :: J12
    double precision :: S12

    e1  = sum(matmul(matmul(reshape(C_sep(:,                   index1), (/1, mat_size/)), F), &
                            reshape(C_sep(:,                   index1), (/mat_size, 1/))))
    e2  = sum(matmul(matmul(reshape(C_sep(:, mat_block1_size + index2), (/1, mat_size/)), F), &
                            reshape(C_sep(:, mat_block1_size + index2), (/mat_size, 1/))))
    J12 = sum(matmul(matmul(reshape(C_sep(:,                   index1), (/1, mat_size/)), F), &
                            reshape(C_sep(:, mat_block1_size + index2), (/mat_size, 1/))))
    S12 = sum(matmul(matmul(reshape(C_sep(:,                   index1), (/1, mat_size/)), S), &
                            reshape(C_sep(:, mat_block1_size + index2), (/mat_size, 1/))))
    e_eff1 = ((e1 + e2) - 2 * J12 * S12 + (e1 - e2) * sqrt(1 - S12 ** 2)) / (2 * (1 - S12 ** 2))
    e_eff2 = ((e1 + e2) - 2 * J12 * S12 + (e2 - e1) * sqrt(1 - S12 ** 2)) / (2 * (1 - S12 ** 2))
    J_eff12 = (J12 - (e1 + e2) * S12 / 2) / (1 - S12 ** 2)

    return
end subroutine calc_coupling

