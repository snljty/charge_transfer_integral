! calculates charge transfer integral from files contain wavefunction information

program main
    implicit none

    interface
        subroutine pause_program()
            implicit none
        end subroutine pause_program

        subroutine inverse_matrix_inplace(a, n, lda, ipiv, tmp, nb)
            implicit none
            integer, intent(in) :: n
            double precision, intent(inout) :: a(n ** 2)
            integer, intent(in) :: lda
            integer, intent(inout) :: ipiv(n)
            integer, intent(in) :: nb
            double precision, intent(inout) :: tmp(n * nb)
        end subroutine inverse_matrix_inplace

        function get_num_items_in_line(line, len_line)
            implicit none
            integer, intent(in) :: len_line
            character(kind=1,len=len_line), intent(in) :: line
            integer :: get_num_items_in_line
        end function get_num_items_in_line

        subroutine read_matrix(fl_unit, mat, mat_size)
            implicit none
            integer, intent(in) :: fl_unit
            double precision, intent(out) :: mat(:, :)
            integer, intent(in) :: mat_size
        end subroutine read_matrix

        subroutine read_triangle_matrix_compact(fl_unit, mat, mat_size)
            implicit none
            integer, intent(in) :: fl_unit
            double precision, intent(out) :: mat(:, :)
            integer, intent(in) :: mat_size
        end subroutine read_triangle_matrix_compact

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

    double precision, parameter :: h_Planck = 6.62607015D-34
    double precision, parameter :: epsilon_0 = 8.854187817D-12
    double precision, parameter :: q_e = 1.602176634D-19
    double precision, parameter :: m_e = 9.10938215D-31
    double precision, parameter :: Hartree_to_eV = m_e * q_e ** 3 / (4.0D0 * epsilon_0 ** 2 * h_Planck ** 2)

    integer :: argc
    integer :: iarg
    integer :: arg_status
    character(kind=1,len=256) :: argv0
    character(kind=1,len=128), pointer :: argv(:)

    character(kind=1,len=256) :: buf
    integer :: buf_pos
    double precision :: tmp_real

    character(kind=1,len=256) :: fl_MO_dimer
    character(kind=1,len=256) :: fl_MO_monomer1
    character(kind=1,len=256) :: fl_MO_monomer2
    character(kind=1,len=256) :: fl_CT_out
    character(kind=1,len=256) :: fl_Fock_dimer
    character(kind=1,len=256) :: fl_47
    logical(kind=1) :: fl_exist
    integer, parameter :: ifl_unit = 10
    integer, parameter :: ofl_unit = 11
    character(kind=1,len=15), parameter :: fl_Multiwfn_in = "Multiwfn_in.txt"
    character(kind=1,len=16), parameter :: fl_Multiwfn_out = "Multiwfn_out.txt"
    character(kind=1,len=10), parameter :: fl_overlap = "intmat.txt"
    character(kind=1,len=8), parameter :: fl_coefficient = "Cmat.txt"

    integer :: num_ele_dimer
    integer :: num_ele_monomer1
    integer :: num_ele_monomer2
    integer :: num_orb_dimer
    integer :: num_orb_monomer1
    integer :: num_orb_monomer2

    integer :: index_homo1
    integer :: index_lumo1
    integer :: index_homo2
    integer :: index_lumo2

    double precision, pointer :: S(:, :)
    double precision, pointer :: C(:, :)
    double precision, pointer :: C_inv(:, :)
    double precision, pointer :: E(:, :)
    double precision, pointer :: F(:, :)
    double precision, pointer :: C_sep(:, :)

    integer, external :: ilaenv
    integer, pointer :: ipiv(:)
    double precision, pointer :: tmp_arr(:)
    integer :: nb

    double precision :: J_eff12
    double precision :: e_eff1
    double precision :: e_eff2

    integer :: i
    integer :: j

    integer :: io_status
    integer :: sys_status

    ! show basic infomation
    write(*, "(a)") "# This program calculates intermolecule transfer integrals, "
    write(*, "(a)") "# especially between HOMOs and LUMOs."
    write(*, "(a)") "#"
    write(*, "(a)") "# Currently only restricted single-determination method is available."
    write(*, "(a)") "# Phase-match method of wavefunctions is not implied yet."
    write(*, "(a)") "#"
    write(*, "(a)") "# sobereva's Multiwfn (http://sobereva.com/multiwfn/) is required "
    write(*, "(a)") "# to be set properly before running this program."
    write(*, "(a)") "# Libraries lapack and blas are required to be linked to this program, "
    write(*, "(a)") "# and to be used during runtime if they are shared libraries."
    write(*, "(a)") "#"
    write(*, "(a)") "# If you are using Gaussian fch/fchk files, "
    write(*, "(a)") "# ""IOp(3/32=2) NoSymmetry"" should appear in your input files."
    write(*, "(a)") "# Generally, all wavefunction files which contains basis functions and "
    write(*, "(a)") "# MO coefficients as well as energies, and supported by Multiwfn, "
    write(*, "(a)") "# including but not limitted to fch/fchk files of Gaussian, molden files "
    write(*, "(a)") "# of ORCA or xTB, should also be supported by this program."
    write(*, "(a)") "# Please note that fch/fchk files from traditional semi-empirical methods "
    write(*, "(a)") "# such as PM7 in Gaussian is not supported by Multiwfn, and hence is not "
    write(*, "(a)") "# supported by this program."
    write(*, "(a)") "#"
    write(*, "(a)") "# You can also provide a Fock matrix yourself, so we do not need to inversely "
    write(*, "(a)") "# solve the Fock matrix. If you are using Gaussian, the first method is save your "
    write(*, "(a)") "# chk file of dimer, then rerun Gaussian job of dimer with keywords "
    write(*, "(a)") "# ""Guess=Read IOp(5/33=3)"", then copy the last ""Fock matrix (alpha)"" or "
    write(*, "(a)") "# ""Fock matrix (beta)"" and the matrix below in the Gaussian output file "
    write(*, "(a)") "# to dimerFock.txt, and then use ""-dF dimerFock.txt"" to read Fock matrix "
    write(*, "(a)") "# directly. The second method is add ""Guess=Read Population=NBORead"" "
    write(*, "(a)") "# in keywords and add a blank line followed by ""$NBO ARCHIVE FILE=MYDIMER $END"" "
    write(*, "(a)") "# at the end of your Gaussian input file and rerun it, the you'll get a "
    write(*, "(a)") "# MYDIMER.47 file, which is an NBO input file containing Fock matrix of dimer."
    write(*, "(a)") "# Then you can use ""-d47 MYDIMER.47"" to assign this file."
    write(*, "(a)") "#"
    write(*, "(a)") "# For other details, please read ""README.md""."
    write(*, "(a)") "#"

    ! get command arguments
    argc = command_argument_count()
    call get_command_argument(0, argv0, arg_status)
    if (argc > 0) then
        allocate(argv(argc))
        do iarg = 1, argc
            call get_command_argument(iarg, argv(iarg), arg_status)
        end do
    end if

    ! Parse auguments
    fl_MO_dimer = ""
    fl_MO_monomer1 = ""
    fl_MO_monomer2 = ""
    fl_Fock_dimer = ""
    fl_47 = ""
    fl_CT_out = ""
    iarg = 1
    do while (iarg <= argc)
        if ((trim(argv(iarg)) == "-h") .or. (trim(argv(iarg)) == "--help")) then
            if (argc /= 1) then
                write(*, "(a)") "When ""-h"" or ""--help"" is given, no other arguments should be given, "
                write(*, "(a,i0,a)") "buf we got ", (argc - 1), " extra argument(s)."
                stop 1
            end if
            write(*, "(a)") "Usage: "
            write(*, "(7x,a,1x,a)")  trim(argv0), "[-h, --help]"
            write(*, "(7x,a,1x,a)") trim(argv0), "[-d, --dimer DIMER]"
            write(*, "(7x,a,1x,a)") repeat(" ", len_trim(argv0)), "[-dF, --dimerFock DIMERFOCK]"
            write(*, "(7x,a,1x,a)") repeat(" ", len_trim(argv0)), "[-d47, --dimer47 MYDIMER.47]"
            write(*, "(7x,a,1x,a)") repeat(" ", len_trim(argv0)), "[-m1, --monomer1 MONOMER1]"
            write(*, "(7x,a,1x,a)") repeat(" ", len_trim(argv0)), "[-m2, --monomer2 MONOMER2]"
            write(*, "(7x,a,1x,a)") repeat(" ", len_trim(argv0)), "[-f, --full OUTFILE]"
            write(*, "()")
            write(*, "(a)") "--help     : print this help message and stop."
            write(*, "(a)") "--dimer    : file contains necessary infomation of the dimer."
            write(*, "(a)") "--monomer1 : file contains necessary infomation of the 1st monomer."
            write(*, "(a)") "--monomer2 : file contains necessary infomation of the 2nd monomer."
            write(*, "(a)") "--dimerFock: file contains the Fock Matrix of the dimer, lower-triangle format, "          
            write(*, "(a)") "             with the first line as comment. This is optional."
            write(*, "(a)") "--dimer47  : NBO .47 file of dimer, with Fock matrix instead. This is optional."
            write(*, "(a)") "--full     : output all orbital transfer integrals, "
            write(*, "(a)") "             instead of HOMO and LUMO only, to a file. This is optional."
            write(*, "()")
            write(*, "(a)") "Exiting normally."
            stop
        else if ((trim(argv(iarg)) == "-d") .or. (trim(argv(iarg)) == "--dimer")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            fl_MO_dimer = trim(argv(iarg))
            inquire(file = trim(fl_MO_dimer), exist = fl_exist)
            if (.not. fl_exist) then
                write(*, "(a,a,a)") "# File """, trim(fl_MO_dimer), """ not found!"
                write(*, "(a)") "# Please input it in the interactive mode later."
                fl_MO_dimer = ""
            end if
        else if ((trim(argv(iarg)) == "-m1") .or. (trim(argv(iarg)) == "--monomer1")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            fl_MO_monomer1 = trim(argv(iarg))
            inquire(file = trim(fl_MO_monomer1), exist = fl_exist)
            if (.not. fl_exist) then
                write(*, "(a,a,a)") "# File """, trim(fl_MO_monomer1), """ not found!"
                write(*, "(a)") "# Please input it in the interactive mode later."
                fl_MO_monomer1 = ""
            end if
        else if ((trim(argv(iarg)) == "-m2") .or. (trim(argv(iarg)) == "--monomer2")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            fl_MO_monomer2 = trim(argv(iarg))
            inquire(file = trim(fl_MO_monomer2), exist = fl_exist)
            if (.not. fl_exist) then
                write(*, "(a,a,a)") "# File """, trim(fl_MO_monomer2), """ not found!"
                write(*, "(a)") "# Please input it in the interactive mode later."
                fl_MO_monomer2 = ""
            end if
        else if ((trim(argv(iarg)) == "-dF") .or. (trim(argv(iarg)) == "--dimerFock")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            fl_Fock_dimer = trim(argv(iarg))
            inquire(file = trim(fl_Fock_dimer), exist = fl_exist)
            if (.not. fl_exist) then
                write(0, "(a,a,a)") "# File """, trim(fl_Fock_dimer), """ not found!"
                write(0, "(a)") "# Please check your file."
                stop 1
            end if
        else if ((trim(argv(iarg)) == "-d47") .or. (trim(argv(iarg)) == "--dimer47")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            if (trim(fl_Fock_dimer) /= "") then
                write(0, "(a)") "# Triangle format Fock matrix and NBO .47 file cannot be provided together!"
                stop 1
            end if
            fl_47 = trim(argv(iarg))
            do while (.true.)
                if (len_trim(fl_47) >= len_trim(".47")) then
                    if (fl_47(len_trim(fl_47) - len_trim(".47") + 1:) == ".47") exit
                end if
                write(0, "(a)") "# The file extension of the NBO .47 file must be .47!"
                stop 1
                exit
            end do
            inquire(file = trim(fl_47), exist = fl_exist)
            if (.not. fl_exist) then
                write(0, "(a,a,a)") "# File """, trim(fl_Fock_dimer), """ not found!"
                write(0, "(a)") "# Please check your file."
                stop 1
            end if
        else if ((trim(argv(iarg)) == "-f") .or. (trim(argv(iarg)) == "--full")) then
            iarg = iarg + 1
            if (iarg > argc) then
                write(0, "(a,a,a)") "# Missing argument for """, trim(argv(iarg - 1)), """!"
                stop 1
            end if
            fl_CT_out = trim(argv(iarg))
        else
            write(0, "(a,1x,a)") "# Error! Unknown argument:", argv(iarg)
            stop 1
        end if
        iarg = iarg + 1
    end do

    ! release memory of command arguments
    if (argc > 0) deallocate(argv)
    argv => null()

    ! get the files' names if not obtained from the command arguments
    write(*, "(a)") "#"
    if (trim(fl_MO_dimer) == "") then
        write(*, "(a)") "# Input the name of file contains MO infomation of the dimer."
        write(*, "(a)") "# If press <Enter> directly, ""dimer.fchk"" will be used."
        do while (.true.)
            read(*, "(a)") fl_MO_dimer
            if (trim(fl_MO_dimer) == "") then
                fl_MO_dimer = "dimer.fchk"
                inquire(file = trim(fl_MO_dimer), exist = fl_exist)
                if (fl_exist) exit
                fl_MO_dimer = "dimer.fch"
                inquire(file = trim(fl_MO_dimer), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a)") "# Neither ""dimer.fchk"" nor ""dimer.fch"" found!"
            else
                inquire(file = trim(fl_MO_dimer), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a,a,a)") "# File """, trim(fl_MO_dimer), """ not found!"
            end if
            write(*, "(a)") "# Input again."
        end do
    end if
    if (trim(fl_MO_monomer1) == "") then
        write(*, "(a)") "# Input the name of file contains MO infomation of the 1st monomer."
        write(*, "(a)") "# If press <Enter> directly, ""monomer1.fchk"" will be used."
        do while (.true.)
            read(*, "(a)") fl_MO_monomer1
            if (trim(fl_MO_monomer1) == "") then
                fl_MO_monomer1 = "monomer1.fchk"
                inquire(file = trim(fl_MO_monomer1), exist = fl_exist)
                if (fl_exist) exit
                fl_MO_monomer1 = "monomer1.fch"
                inquire(file = trim(fl_MO_monomer1), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a)") "# Neither ""monomer1.fchk"" nor ""monomer1.fch"" found!"
            else
                inquire(file = trim(fl_MO_monomer1), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a,a,a)") "# File """, trim(fl_MO_monomer1), """ not found!"
            end if
            write(*, "(a)") "# Input again."
        end do
    end if
    if (trim(fl_MO_monomer2) == "") then
        write(*, "(a)") "# Input the name of file contains MO infomation of the 2nd monomer."
        write(*, "(a)") "# If press <Enter> directly, ""monomer2.fchk"" will be used."
        do while (.true.)
            read(*, "(a)") fl_MO_monomer2
            if (trim(fl_MO_monomer2) == "") then
                fl_MO_monomer2 = "monomer2.fchk"
                inquire(file = trim(fl_MO_monomer2), exist = fl_exist)
                if (fl_exist) exit
                fl_MO_monomer2 = "monomer2.fch"
                inquire(file = trim(fl_MO_monomer2), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a)") "# Neither ""monomer2.fchk"" nor ""monomer2.fch"" found!"
            else
                inquire(file = trim(fl_MO_monomer2), exist = fl_exist)
                if (fl_exist) exit
                write(*, "(a,a,a)") "# File """, trim(fl_MO_monomer2), """ not found!"
            end if
            write(*, "(a)") "# Input again."
        end do
    end if

    ! show names of files to be used
    write(*, "(a,a,a)") "# Using """, trim(fl_MO_dimer), """ for dimer."
    write(*, "(a,a,a)") "# Using """, trim(fl_MO_monomer1), """ for monomer1."
    write(*, "(a,a,a)") "# Using """, trim(fl_MO_monomer2), """ for monomer2."
    if (trim(fl_Fock_dimer) /= "") write(*, "(a,a,a)") &
        "# Using """, trim(fl_Fock_dimer), """ for Fock matrix of dimer."
    if (trim(fl_47) /= "") write(*, "(a,a,a)") &
        "# Using """, trim(fl_47), """ for Fock matrix of dimer."
    if (trim(fl_CT_out) /= "") write(*, "(a,a,a)") "# Using """, trim(fl_CT_out), """ for output."

    ! first, get amount of electrons and orbitals, to see if the sum of monomers equals dimer

    ! prepare input file used by Multiwfn of basic inofmtation
    open(ofl_unit, file = trim(fl_Multiwfn_in), action = "write", status = "replace")
    write(ofl_unit, "(i3)") -10
    close(ofl_unit)

    ! call Multiwfn of monomer1, get amount of electrons and orbitals
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of monomer1 ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_monomer1) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of monomer1!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_monomer1), """."
        stop 2
    end if

    ! read amound of orbitals and electrons of monomer1 from Multiwfn output
    write(*, "(a)") "# Reading amount of orbitals and electrons of monomer1 ..."
    open(ifl_unit, file = trim(fl_Multiwfn_out), action = "read", status = "old")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Total/Alpha/Beta electrons:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Total/Alpha/Beta electrons:")
    read(buf(buf_pos:), *) tmp_real
    num_ele_monomer1 = nint(tmp_real)
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Basis functions:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Basis functions:")
    read(buf(buf_pos:), *) num_orb_monomer1
    close(ifl_unit, status = "delete")

    ! call Multiwfn of monomer2, get amount of electrons and orbitals
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of monomer2 ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_monomer2) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of monomer2!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_monomer2), """."
        stop 2
    end if

    ! read amound of orbitals and electrons of monomer2 from Multiwfn output
    write(*, "(a)") "# Reading amount of orbitals and electrons of monomer2 ..."
    open(ifl_unit, file = trim(fl_Multiwfn_out), action = "read", status = "old")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Total/Alpha/Beta electrons:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Total/Alpha/Beta electrons:")
    read(buf(buf_pos:), *) tmp_real
    num_ele_monomer2 = nint(tmp_real)
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Basis functions:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Basis functions:")
    read(buf(buf_pos:), *) num_orb_monomer2
    close(ifl_unit, status = "delete")

    ! call Multiwfn of dimer, get amount of electrons and orbitals
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of dimer ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_dimer) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of dimer!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_dimer), """."
        stop 2
    end if

    ! read amound of orbitals and electrons of dimer from Multiwfn output
    write(*, "(a)") "# Reading amount of orbitals and electrons of dimer ..."
    open(ifl_unit, file = trim(fl_Multiwfn_out), action = "read", status = "old")
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Total/Alpha/Beta electrons:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Total/Alpha/Beta electrons:")
    read(buf(buf_pos:), *) tmp_real
    num_ele_dimer = nint(tmp_real)
    do while (.true.)
        read(ifl_unit, "(a)") buf
        buf_pos = index(buf, "Basis functions:")
        if (buf_pos /= 0) exit
    end do
    buf_pos = buf_pos + len_trim("Basis functions:")
    read(buf(buf_pos:), *) num_orb_dimer
    close(ifl_unit, status = "delete")

    ! remove Multiwfn input file
    open(ifl_unit, file = trim(fl_Multiwfn_in))
    close(ifl_unit, status = "delete")

    ! check the sum of electrons and orbitals
    write(*, "(a)") "# number of electrons of monomer1/monomer2/dimer:"
    write(*, "(a1,27x,i5,3x,i5,1x,i5)") "#", num_ele_monomer1, num_ele_monomer2, num_ele_dimer
    if (num_ele_dimer /= num_ele_monomer1 + num_ele_monomer2) then
        write(*, "(a)") "Error! number of electrons mismatches in "
        write(*, "(a,i5,a,i5,a,i5,a)") "dimer (", num_ele_dimer, ") and monomers (", &
            num_ele_monomer1, "+", num_ele_monomer2, ")!"
        stop 3
    end if
    write(*, "(a)") "# number of  orbitals of monomer1/monomer2/dimer:"
    write(*, "(a1,27x,i5,3x,i5,1x,i5)") "#", num_orb_monomer1, num_orb_monomer2, num_orb_dimer
    if (num_orb_dimer /= num_orb_monomer1 + num_orb_monomer2) then
        write(*, "(a)") "Error! number of orbitals mismatches in "
        write(*, "(a,i5,a,i5,a,i5,a)") "dimer (", num_orb_dimer, ") and monomers (", &
            num_orb_monomer1, "+", num_orb_monomer2, ")!"
        stop 3
    end if

    ! allocate memory fos S, F and C_sep.
    allocate(S(num_orb_dimer, num_orb_dimer))
    allocate(F(num_orb_dimer, num_orb_dimer))
    allocate(C_sep(num_orb_dimer, num_orb_dimer))
    nb = ilaenv(1, "dgetri", "", num_orb_dimer, -1, -1, -1)

    ! second, read coefficient matrices of monomers

    ! prepare input file to obtain coefficient matrix
    open(ofl_unit, file = trim(fl_Multiwfn_in), action = "write", status = "replace")
    write(ofl_unit, "(i3)") (/6, 5, 2, -1, -10/) ! Multiwfn commands
    close(ofl_unit)

    ! call Multiwfn to obtain coefficient matrix of monomer1
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of monomer1 ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_monomer1) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of monomer1!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_monomer1), """."
        stop 2
    end if

    ! read coefficient matrix of monomer1
    open(ifl_unit, file = trim(fl_Multiwfn_out))
    close(ifl_unit, status = "delete")
    write(*, "(a)") "# Reading coefficient matrix of monomer1 ..."
    open(ifl_unit, file = trim(fl_coefficient), action = "read", status = "old")
    read(ifl_unit, "(a)") buf
    read(ifl_unit, "(a)") buf
    call read_matrix(ifl_unit, C_sep(:num_orb_monomer1, :num_orb_monomer1), num_orb_monomer1)
    close(ifl_unit, status = "delete")

    ! call Multiwfn to obtain coefficient matrix of monomer2
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of monomer2 ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_monomer2) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of monomer2!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_monomer2), """."
        stop 2
    end if

    ! read coefficient matrix of monomer2
    open(ifl_unit, file = trim(fl_Multiwfn_out))
    close(ifl_unit, status = "delete")
    write(*, "(a)") "# Reading coefficient matrix of monomer2 ..."
    open(ifl_unit, file = trim(fl_coefficient), action = "read", status = "old")
    read(ifl_unit, "(a)") buf
    read(ifl_unit, "(a)") buf
    call read_matrix(ifl_unit, C_sep(num_orb_monomer1 + 1:, num_orb_monomer1 + 1:), num_orb_monomer2)
    close(ifl_unit, status = "delete")

    ! remove Multiwfn input file
    open(ifl_unit, file = trim(fl_Multiwfn_in))
    close(ifl_unit, status = "delete")

    ! third, read overlap matrix, coefficient matrix and orbital energies of dimer.

    ! prepare input file to obtain overlap matrix, &
    ! coefficient matrix (optional) and orbital energy matrix (optional).
    open(ofl_unit, file = trim(fl_Multiwfn_in), action = "write", status = "replace")
    if (trim(fl_Fock_dimer) == "" .and. trim(fl_47) == "") then
        write(ofl_unit, "(i3)") (/6, 3, 5, 2, 7, 1, 2, -1, -10/) ! Multiwfn commands
    else
        write(ofl_unit, "(i3)") (/6, 2, 7, 1, 2, -1, -10/) ! Multiwfn commands
    end if
    close(ofl_unit)

    ! call Multiwfn to obtain overlap matrix, &
    ! coefficient matrix (optional) and orbital energies (optional) of dimer.
    write(*, "(a)") "# Calling Multiwfn to obtain amount of orbitals and electrons of dimer ..."
    call execute_command_line("Multiwfn " // trim(fl_MO_dimer) // " < " // trim(fl_Multiwfn_in) // &
        " > " // trim(fl_Multiwfn_out), wait = .true., exitstat = sys_status)
    if (sys_status /= 0) then
        write(*, "(a)") "Error calling Multiwfn of dimer!"
        write(*, "(a,a,a)") "Please check Multiwfn relative settings and """, trim(fl_MO_dimer), """."
        stop 2
    end if

    ! remove Multiwfn input file
    open(ifl_unit, file = trim(fl_Multiwfn_in))
    close(ifl_unit, status = "delete")

    ! read S from Multiwfn output file
    write(*, "(a)") "# Reading overlap matrix of dimer ..."
    open(ifl_unit, file = trim(fl_overlap), action = "read", status = "old")
    read(ifl_unit, "(a)") buf
    call read_lt_matrix(ifl_unit, S, num_orb_dimer)
    close(ifl_unit, status = "delete")

    ! if not provided Fock matrix of dimer, read coefficient matrix of dimer and get its reverse
    if (trim(fl_Fock_dimer) == "" .and. trim(fl_47) == "") then
        ! allocate memory for tmp_arr, ipiv, C, C_inv and E
        allocate(C(num_orb_dimer, num_orb_dimer))
        allocate(C_inv(num_orb_dimer, num_orb_dimer))
        allocate(E(num_orb_dimer, num_orb_dimer))
        allocate(ipiv(num_orb_dimer))
        allocate(tmp_arr(num_orb_dimer * nb))
        ! read coefficient matrix of dimer
        write(*, "(a)") "# Reading coefficient matrix of dimer ..."
        open(ifl_unit, file = trim(fl_coefficient), action = "read", status = "old")
        read(ifl_unit, "(a)") buf
        read(ifl_unit, "(a)") buf
        call read_matrix(ifl_unit, C, num_orb_dimer)
        close(ifl_unit, status = "delete")
        ! read orbital energies of dimer
        E = 0.0D0
        open(ifl_unit, file = trim(fl_Multiwfn_out), action = "read", status = "old")
        do while (.true.)
            read(ifl_unit, "(a)") buf
            buf_pos = index(buf, "Basic information of all orbitals:")
            if (buf_pos /= 0) exit
        end do
        do i = 1, num_orb_dimer
            read(ifl_unit, "(a)") buf
            buf_pos = index(buf, "Ene(au/eV):")
            buf_pos = buf_pos + len_trim("Ene(au/eV):")
            read(buf(buf_pos:), *) E(i, i)
        end do
        close(ifl_unit, status = "delete")
        ! calculate C_inv from C
        write(*, "(a)") "# Inversing coefficient matrix of dimer ..."
        C_inv = C
        call inverse_matrix_inplace(C_inv, num_orb_dimer, num_orb_dimer, ipiv, tmp_arr, nb)
    else
        ! delete Multiwfn output file
        open(ifl_unit, file = trim(fl_Multiwfn_out), action = "read", status = "old")
        close(ifl_unit, status = "delete")
    end if

    ! calculate or read F
    if (trim(fl_Fock_dimer) == "" .and. trim(fl_47) == "") then
        write(*, "(a)") "# Calculating Fock matrix ..."
        F = matmul(matmul(matmul(S, C), E), C_inv)
        ! release memory for tmp_arr, ipiv, C, C_inv and E
        deallocate(tmp_arr)
        tmp_arr => null()
        deallocate(ipiv)
        ipiv => null()
        deallocate(C)
        C => null()
        deallocate(C_inv)
        C_inv => null()
        deallocate(E)
        E => null()
    else if (trim(fl_Fock_dimer) /= "") then
        write(*, "(a)") "# Reading Fock matrix of dimer from a plain text file ..."
        open(ifl_unit, file = trim(fl_Fock_dimer), action = "read", status = "old")
        read(ifl_unit, "(a)") buf ! title comment line of file for Fock matrix of dimer
        call read_lt_matrix(ifl_unit, F, num_orb_dimer)
        close(ifl_unit)
    else
        write(*, "(a)") "# Reading Fock matrix of dimer from an NBO .47 file ..."
        open(ifl_unit, file = trim(fl_47), action = "read", status = "old")
        do while (.true.)
            read(ifl_unit, "(a)", iostat = io_status) buf
            if (io_status /= 0) then
                write(0, "(a,a,a)") "Error! Cannot find ""$FOCK"" field in """, trim(fl_47), """."
                stop 1
            end if
            if (index(buf, "$FOCK") /= 0) exit
        end do
        call read_triangle_matrix_compact(ifl_unit, F, num_orb_dimer)
        close(ifl_unit)
    end if

    ! calculate charge transfer integral
    write(*, "(a)") "# Calculating transfer integrals ..."
    index_homo1 = num_ele_monomer1 / 2
    index_lumo1 = index_homo1 + 1
    index_homo2 = num_ele_monomer2 / 2
    index_lumo2 = index_homo2 + 1
    if (trim(fl_CT_out) == "") then
        write(*, "(a,i4,a,i4)") "# Monomer 1    HOMO: ", index_homo1, "        LUMO: ", index_lumo1
        write(*, "(a,i4,a,i4)") "# Monomer 2    HOMO: ", index_homo2, "        LUMO: ", index_lumo2
        write(*, "(a)") "# Transfer integrals between HOMO and LUMO of monomers:"
        write(*, "(a)") "# monomer_1    monomer_2    J_eff_12/meV"
        ! do i = index_homo1, index_lumo1
        !     do j = index_homo2, index_lumo2
        !         call calc_coupling(i, j, num_orb_dimer, num_orb_monomer1, F, S, C_sep, &
        !             J_eff12, e_eff1, e_eff2)
        !         write(*, "(3x,i4,9x,i4,7x,f13.7)") i, j, J_eff12 * Hartree_to_eV * 1.0D3
        !     end do
        ! end do
        i = index_homo1
        j = index_homo2
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
        write(*, "(a,a,a)") "# Intermolecule charge transfer integrals have been saved to """, &
            trim(fl_CT_out), """."
    end if

    ! release memory of S, F and C_sep
    deallocate(S)
    S => null()
    deallocate(F)
    F => null()
    deallocate(C_sep)
    C_sep => null()

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

subroutine inverse_matrix_inplace(a, n, lda, ipiv, tmp, nb)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: a(n ** 2)
    integer, intent(in) :: lda
    integer, intent(inout) :: ipiv(n)
    integer, intent(in) :: nb
    double precision, intent(inout) :: tmp(n * nb)
    integer :: info

    ! LU decompostion
    call dgetrf(n, n, a, lda, ipiv, info)
    if (info < 0) then
        write(*, "(a,i1,a)") "The argument with index ", - info, &
            " in subroutine dgetrf has illegal value!"
        stop 4
    else if (info > 0) then
        write(*, "(a,i5,a,i5,a)") "Warning: U(", info, ", ", info, &
            ") is exactly zero!"
    end if

    ! Inverse matrix
    call dgetri(n, a, lda, ipiv, tmp, n ** 2, info)
    if (info < 0) then
        write(*, "(a,i1,a)") "The argument with index ", - info, &
            " in subroutine dgetri has illegal value!"
        stop 4
    else if (info > 0) then
        write(*, "(a,i5,a,i5,a)") "Error: U(", info, ", ", info, ") is exactly zero!"
        write(*, "(a)") "The matrix cannot be inversed!"
        stop 5
    end if

    return
end subroutine inverse_matrix_inplace

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

! If assume this is the upper triangle, then it is column-major.
subroutine read_triangle_matrix_compact(fl_unit, mat, mat_size)
    implicit none
    integer, intent(in) :: fl_unit
    double precision, intent(out) :: mat(:, :)
    integer, intent(in) :: mat_size
    character(kind=1,len=256) :: buf
    integer :: i, j
    integer :: read_pos
    integer :: pos_step

    buf = ""
    i = 1
    j = 0
    read_pos = 1
    do while (.true.)
        j = j + 1
        if (j > i) then
            i = i + 1
            j = 1
            if (i > mat_size) exit
        end if
        pos_step = verify(buf(read_pos:), " ")
        if (pos_step == 0) then
            read(fl_unit, "(a)") buf
            read_pos = 1
            pos_step = verify(buf(read_pos:), " ")
        end if
        read_pos = read_pos + pos_step - 1
        read(buf(read_pos:), *) mat(i, j)
        mat(j, i) = mat(i, j)
        pos_step = index(buf(read_pos:), " ")
        if (pos_step == 0) then
            read(fl_unit, "(a)") buf
            read_pos = 1
        else
            read_pos = read_pos + pos_step - 1
        end if
    end do

    return
end subroutine read_triangle_matrix_compact

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
        do j = 1, i
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

