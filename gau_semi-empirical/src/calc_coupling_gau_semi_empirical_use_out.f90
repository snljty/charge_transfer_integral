program main
    implicit none

    interface
        subroutine pause_program()
            implicit none
        end subroutine pause_program

        subroutine read_matrix(fl_unit, mat, mat_size)
            implicit none
            integer(kind=4), intent(in) :: fl_unit
            real(kind=8), intent(out) :: mat(:, :)
            integer(kind=4), intent(in) :: mat_size
        end subroutine read_matrix

        subroutine read_lt_matrix(fl_unit, mat, mat_size)
            implicit none
            integer(kind=4), intent(in) :: fl_unit
            real(kind=8), intent(out) :: mat(:, :)
            integer(kind=4), intent(in) :: mat_size
        end subroutine read_lt_matrix

        subroutine calc_coupling(index1, index2, mat_size, mat_block1_size, F, S, C_sep, &
            J_eff12, e_eff1, e_eff2)
            implicit none
            integer(kind=4), intent(in) :: index1
            integer(kind=4), intent(in) :: index2
            integer(kind=4), intent(in) :: mat_size
            integer(kind=4), intent(in) :: mat_block1_size
            real(kind=8), intent(in) :: F(mat_size, mat_size)
            real(kind=8), intent(in) :: S(mat_size, mat_size)
            real(kind=8), intent(in) :: C_sep(mat_size, mat_size)
            real(kind=8), intent(out) :: J_eff12
            real(kind=8), intent(out) :: e_eff1
            real(kind=8), intent(out) :: e_eff2
        end subroutine calc_coupling
    end interface

    integer(kind=4) :: argc
    integer(kind=4) :: iarg
    ! integer(kind=4) :: arg_status
    character(kind=1,len=128), allocatable :: argv(:)

    integer(kind=4), parameter :: stdin_unit = 5
    integer(kind=4), parameter :: stdout_unit = 6
    integer(kind=4), parameter :: stderr_unit = 0

    integer(kind=4), parameter :: exit_success = 0
    integer(kind=4), parameter :: exit_failure = 1

    real(kind=8), parameter :: h_Planck = 6.62607015D-34
    real(kind=8), parameter :: epsilon_0 = 8.854187817D-12
    real(kind=8), parameter :: q_e = 1.602176634D-19
    real(kind=8), parameter :: m_e = 9.10938215D-31
    real(kind=8), parameter :: Hartree_to_eV = m_e * q_e ** 3 / (4.0D0 * epsilon_0 ** 2 * h_Planck ** 2)

    character(kind=1,len=256) :: buf
    integer(kind=4) :: buf_pos

    character(kind=1,len=128) :: fl_MO_dimer
    character(kind=1,len=128) :: fl_MO_monomer1
    character(kind=1,len=128) :: fl_MO_monomer2
    character(kind=1,len=128) :: fl_CT_out
    integer(kind=4), parameter :: ifl_unit = 10
    integer(kind=4), parameter :: ofl_unit = 11

    integer(kind=4) :: num_ele_dimer
    integer(kind=4) :: num_ele_monomer1
    integer(kind=4) :: num_ele_monomer2
    integer(kind=4) :: num_orb_dimer
    integer(kind=4) :: num_orb_monomer1
    integer(kind=4) :: num_orb_monomer2

    integer(kind=4) :: num_ele_alpha ! temporary
    integer(kind=4) :: num_ele_beta  ! temporary

    integer(kind=4) :: index_homo1
    integer(kind=4) :: index_lumo1
    integer(kind=4) :: index_homo2
    integer(kind=4) :: index_lumo2

    real(kind=8), allocatable :: S(:, :)
    real(kind=8), allocatable :: F(:, :)
    real(kind=8), allocatable :: C_sep(:, :)

    real(kind=8) :: J_eff12
    real(kind=8) :: e_eff1
    real(kind=8) :: e_eff2

    integer(kind=4) :: i
    integer(kind=4) :: j

    integer(kind=4) :: io_status

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
        write(*, "(a)") "# monomer_1    monomer_2    J_eff_12/meV      e_eff_1/eV      e_eff_2/eV   "
        do i = index_homo1, index_lumo1
            do j = index_homo2, index_lumo2
                call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
                    J_eff12, e_eff1, e_eff2)
                write(*, "(3x,i4,9x,i4,7x,f13.7,4x,f12.7,4x,f12.7)") i, j, &
                    J_eff12 * Hartree_to_eV * 1.0D3, e_eff1 * Hartree_to_eV, e_eff2 * Hartree_to_eV
            end do
        end do
    else
        write(*, "(a)") "Calculating transfer integrals ..."
        open(ofl_unit, file = fl_CT_out, action = "write", status = "replace")
        write(ofl_unit, "(a,i4,a,i4)") "# Monomer 1    HOMO: ", index_homo1, "        LUMO: ", index_lumo1
        write(ofl_unit, "(a,i4,a,i4)") "# Monomer 2    HOMO: ", index_homo2, "        LUMO: ", index_lumo2
        write(ofl_unit, "(a)") "# monomer_1    monomer_2    J_eff_12/meV      e_eff_1/eV      e_eff_2/eV   "
        do i = 1, num_orb_monomer1
            do j = 1, num_orb_monomer2
                call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
                    J_eff12, e_eff1, e_eff2)
                write(ofl_unit, "(3x,i4,9x,i4,7x,f13.7,4x,f12.7,4x,f12.7)") i, j, &
                    J_eff12 * Hartree_to_eV * 1.0D3, e_eff1 * Hartree_to_eV, e_eff2 * Hartree_to_eV
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
    integer(kind=4), intent(in) :: len_line
    character(kind=1,len=len_line), intent(in) :: line
    integer(kind=4) :: get_num_items_in_line
    integer(kind=4) :: pos_read_in_line

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
    integer(kind=4), intent(in) :: fl_unit
    real(kind=8), intent(out) :: mat(:, :)
    integer(kind=4), intent(in) :: mat_size
    integer(kind=4), external :: get_num_items_in_line

    character(kind=1,len=256) :: buf
    integer(kind=4) :: col_in_this_group
    integer(kind=4) :: pos_read_in_line
    integer(kind=4) :: i
    integer(kind=4) :: j

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
    integer(kind=4), intent(in) :: fl_unit
    real(kind=8), intent(out) :: mat(:, :)
    integer(kind=4), intent(in) :: mat_size
    integer(kind=4), external :: get_num_items_in_line
    integer(kind=4) :: max_col_use
    character(kind=1,len=256) :: buf
    integer(kind=4) :: i
    integer(kind=4) :: j
    integer(kind=4) :: k
    integer(kind=4) :: pos_read_in_line
    integer(kind=4) :: pos_next
    integer(kind=4) :: start_index
    integer(kind=4) :: num_groups

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
    integer(kind=4), intent(in) :: index1
    integer(kind=4), intent(in) :: index2
    integer(kind=4), intent(in) :: mat_size
    integer(kind=4), intent(in) :: mat_block1_size
    real(kind=8), intent(in) :: F(mat_size, mat_size)
    real(kind=8), intent(in) :: S(mat_size, mat_size)
    real(kind=8), intent(in) :: C_sep(mat_size, mat_size)
    real(kind=8), intent(out) :: J_eff12
    real(kind=8), intent(out) :: e_eff1
    real(kind=8), intent(out) :: e_eff2
    real(kind=8) :: e1
    real(kind=8) :: e2
    real(kind=8) :: J12
    real(kind=8) :: S12

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

