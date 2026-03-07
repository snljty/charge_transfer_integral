submodule (coupling_module) matrix_reader
    implicit none

contains

    module subroutine read_gau_vector(ret, ifile, rows, label)
        implicit none
        double precision, allocatable, intent(out) :: ret(:)
        integer, intent(in) :: ifile, rows
        character(len=*), intent(in) :: label

        integer :: status
        character(len=128) :: line
        logical :: found

        found = .false.

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, label) /= 0) then
                found = .true.
                exit
            end if
        end do

        if (.not. found) error stop "Cannot find label."

        allocate(ret(rows))
        read(ifile, *) ret

    end subroutine read_gau_vector

    module subroutine read_gau_matrix(ret, ifile, rows, cols, label)
        implicit none
        double precision, allocatable, intent(out) :: ret(:, :)
        integer, intent(in) :: ifile, rows, cols
        character(len=*), intent(in) :: label

        integer :: status
        character(len=128) :: line
        logical :: found

        found = .false.

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, label) /= 0) then
                found = .true.
                exit
            end if
        end do

        if (.not. found) error stop "Cannot find label."

        allocate(ret(rows, cols))
        read(ifile, *) ret

    end subroutine read_gau_matrix

    module subroutine read_gau_lower_triangle_matrix(ret, ifile, rows, label, use_last_if_exist)
        ! assume square matrix
        implicit none
        double precision, allocatable, intent(out) :: ret(:, :)
        integer, intent(in) :: ifile, rows
        character(len=*), intent(in) :: label
        logical, intent(in), optional :: use_last_if_exist

        integer :: status
        character(len=128) :: line
        logical :: found, use_last_if_exist_actual
        integer :: ibeg, jbeg, jend, i, j
        integer :: nouse

        integer :: nread

        found = .false.
        use_last_if_exist_actual = .false.
        if (present(use_last_if_exist)) use_last_if_exist_actual = use_last_if_exist

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, label) /= 0) then
                found = .true.
                nread = 0
                if (.not. allocated(ret)) allocate(ret(rows, rows))
                ibeg = 1
                do while (ibeg <= rows)
                    read(ifile, '(a)') line
                    jbeg = ibeg
                    do i = ibeg, rows
                        read(ifile, '(a)') line
                        jend = jbeg + 5 - 1
                        if (jend > i) jend = i
                        read(line, *) nouse, ret(i, jbeg:jend)
                        nread = nread + jend + 1 - jbeg
                    end do
                    ibeg = ibeg + 5
                end do
                if (nread * 2 /= rows * (rows + 1)) error stop "Cannot read matrix."
                do i = 1, rows
                    do j = 1, i - 1
                        ret(j, i) = ret(i, j)
                    end do
                end do
                if (.not. use_last_if_exist_actual) exit
            end if
        end do

        if (.not. found .and. .not. use_last_if_exist_actual) error stop "Cannot find label."
        
    end subroutine read_gau_lower_triangle_matrix

    module subroutine read_dimer_ene(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_vector(this%dimer_ene, this%ifile_dimer, this%nbasis_dimer, "Alpha Orbital Energies")

    end subroutine read_dimer_ene

    module subroutine read_monomer1_ene(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_vector(this%monomer1_ene, this%ifile_monomer1, this%nbasis_monomer1, "Alpha Orbital Energies")

    end subroutine read_monomer1_ene

    module subroutine read_monomer2_ene(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_vector(this%monomer2_ene, this%ifile_monomer2, this%nbasis_monomer2, "Alpha Orbital Energies")

    end subroutine read_monomer2_ene

    module subroutine read_dimer_coeff(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_matrix(this%dimer_coeff, this%ifile_dimer, &
            this%nbasis_dimer, this%nbasis_dimer, "Alpha MO coefficients")

    end subroutine read_dimer_coeff

    module subroutine read_monomer1_coeff(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_matrix(this%monomer1_coeff, this%ifile_monomer1, &
            this%nbasis_monomer1, this%nbasis_monomer1, "Alpha MO coefficients")

    end subroutine read_monomer1_coeff

    module subroutine read_monomer2_coeff(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_matrix(this%monomer2_coeff, this%ifile_monomer2, &
            this%nbasis_monomer2, this%nbasis_monomer2, "Alpha MO coefficients")

    end subroutine read_monomer2_coeff

    module subroutine read_dimer_overlap(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_gau_lower_triangle_matrix(this%dimer_overlap, this%ifile_dimer_out, this%nbasis_dimer, &
            "*** Overlap ***")

    end subroutine read_dimer_overlap

    module subroutine read_dimer_Fock(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        double precision, allocatable :: tmp(:, :)
        integer, allocatable :: ipiv(:)
        integer :: i, status
        integer :: n

        n = this%nbasis_dimer
        call read_gau_lower_triangle_matrix(this%dimer_Fock, this%ifile_dimer_out, this%nbasis_dimer, &
            "Fock matrix (alpha)", .true.)

        if (.not. allocated(this%dimer_Fock)) then
            ! write(0, '(a)') "solving Fock ..."
            allocate(this%dimer_Fock(n, n))
            allocate(tmp(n, n))
            allocate(ipiv(n))

            call dgemm('T', 'T', n, n, n, &
                1.D0, this%dimer_coeff, n, this%dimer_overlap, n, &
                0.D0, this%dimer_Fock, n) ! "F" = C.T @ S.T as a temporary matrix
            forall (i = 1:n) tmp(i, :) = this%dimer_ene(i) * this%dimer_Fock(i, :) ! tmp = diag(e) * (C.T S.T)
            this%dimer_Fock = transpose(this%dimer_coeff) ! "F" = C.T
            ! solve "F" @ "tmp_solved" = tmp, a.k.a. (C.T) @ (? == F.T) = (S @ C @ E).T
            call dgesv(n, n, this%dimer_Fock, n, &
                ipiv, tmp, n, status)
            if (status /= 0) error stop "Solving linear equation failed."
            this%dimer_Fock = transpose(tmp) ! F = (F.T).T

            deallocate(tmp)
            deallocate(ipiv)
        end if

    end subroutine read_dimer_Fock

end submodule matrix_reader