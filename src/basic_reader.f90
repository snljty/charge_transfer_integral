submodule (coupling_module) basic_reader
    implicit none

contains

    module subroutine prepare_files(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        logical :: exist

        inquire(file=trim(this%dimer_prefix)//".fch", exist=exist)
        if (exist) then
            open(newunit=this%ifile_dimer, file=trim(this%dimer_prefix)//".fch", status="old", action="read")
        else
            inquire(file=trim(this%dimer_prefix)//".fchk", exist=exist)
            if (.not. exist) error stop "Cannot open dimer fch file for reading."
            open(newunit=this%ifile_dimer, file=trim(this%dimer_prefix)//".fchk", status="old", action="read")
        end if

        inquire(file=trim(this%monomer1_prefix)//".fch", exist=exist)
        if (exist) then
            open(newunit=this%ifile_monomer1, file=trim(this%monomer1_prefix)//".fch", status="old", action="read")
        else
            inquire(file=trim(this%monomer1_prefix)//".fchk", exist=exist)
            if (.not. exist) error stop "Cannot open monomer1 fch file for reading."
            open(newunit=this%ifile_monomer1, file=trim(this%monomer1_prefix)//".fchk", status="old", action="read")
        end if

        inquire(file=trim(this%monomer2_prefix)//".fch", exist=exist)
        if (exist) then
            open(newunit=this%ifile_monomer2, file=trim(this%monomer2_prefix)//".fch", status="old", action="read")
        else
            inquire(file=trim(this%monomer2_prefix)//".fchk", exist=exist)
            if (.not. exist) error stop "Cannot open monomer2 fch file for reading."
            open(newunit=this%ifile_monomer2, file=trim(this%monomer2_prefix)//".fchk", status="old", action="read")
        end if

        inquire(file=trim(this%dimer_prefix)//".out", exist=exist)
        if (exist) then
            open(newunit=this%ifile_dimer_out, file=trim(this%dimer_prefix)//".out", status="old", action="read")
        else
            inquire(file=trim(this%dimer_prefix)//".log", exist=exist)
            if (.not. exist) error stop "Cannot open dimer fch file for reading."
            open(newunit=this%ifile_dimer_out, file=trim(this%dimer_prefix)//".log", status="old", action="read")
        end if

    end subroutine prepare_files

    module subroutine close_files(this)
        implicit none
        class(coupling_calculator), intent(in) :: this

        close(this%ifile_dimer)
        close(this%ifile_monomer1)
        close(this%ifile_monomer2)
        close(this%ifile_dimer_out)
    end subroutine close_files

    module subroutine read_basic_info(ifile, natoms, nelectrons, nbasis)
        implicit none
        integer, intent(in) :: ifile
        integer, intent(out) :: natoms, nelectrons, nbasis

        integer :: status
        character(len=128) :: line
        logical :: found_natoms, found_nelectrons, found_nbasis

        found_natoms = .false.
        found_nelectrons = .false.
        found_nbasis = .false.

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, "Number of atoms") /= 0) then
                found_natoms = .true.
                read(line(45:), *) natoms
                exit
            end if
        end do

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, "Number of electrons") /= 0) then
                found_nelectrons = .true.
                read(line(45:), *) nelectrons
                exit
            end if
        end do

        do while (.true.)
            read(ifile, '(a)', iostat=status) line
            if (status /= 0) exit
            if (index(line, "Number of independent functions") /= 0 .or. &
                index(line, "Number of independant functions") /= 0) then
                found_nbasis = .true.
                read(line(45:), *) nbasis
                exit
            end if
        end do

        if (.not. (found_natoms .and. found_nelectrons .and. found_nbasis)) then
            error stop "Cannot read natoms, nelectrons or nbasis."
        end if

    end subroutine read_basic_info

    module subroutine read_dimer_info(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_basic_info(this%ifile_dimer, this%natoms_dimer, this%nelectrons_dimer, this%nbasis_dimer)
        if (mod(this%nelectrons_dimer, 2) /= 0) error stop "Only closed-shell wavefunction is supported."

    end subroutine read_dimer_info

    module subroutine read_monomer1_info(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_basic_info(this%ifile_monomer1, this%natoms_monomer1, this%nelectrons_monomer1, this%nbasis_monomer1)
        if (mod(this%nelectrons_monomer1, 2) /= 0) error stop "Only closed-shell wavefunction is supported."

        this%iHOMO_monomer1 = this%nelectrons_monomer1 / 2
        this%iLUMO_monomer1 = this%iHOMO_monomer1 + 1

    end subroutine read_monomer1_info

    module subroutine read_monomer2_info(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call read_basic_info(this%ifile_monomer2, this%natoms_monomer2, this%nelectrons_monomer2, this%nbasis_monomer2)
        if (mod(this%nelectrons_monomer2, 2) /= 0) error stop "Only closed-shell wavefunction is supported."

        this%iHOMO_monomer2 = this%nelectrons_monomer2 / 2
        this%iLUMO_monomer2 = this%iHOMO_monomer2 + 1

    end subroutine read_monomer2_info

end submodule basic_reader