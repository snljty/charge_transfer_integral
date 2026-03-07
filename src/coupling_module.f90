module coupling_module
    implicit none
    type :: coupling_calculator
        private
        character(len=128) :: dimer_prefix, monomer1_prefix, monomer2_prefix
        integer :: ifile_dimer, ifile_monomer1, ifile_monomer2, ifile_dimer_out

        integer :: natoms_dimer, natoms_monomer1, natoms_monomer2
        integer :: nelectrons_dimer, nelectrons_monomer1, nelectrons_monomer2
        integer :: nbasis_dimer, nbasis_monomer1, nbasis_monomer2
        integer :: iHOMO_monomer1, iLUMO_monomer1, iHOMO_monomer2, iLUMO_monomer2
        double precision, dimension(:), allocatable :: dimer_ene
        double precision, dimension(:, :), allocatable :: dimer_Fock, dimer_overlap, dimer_coeff
        double precision, dimension(:), allocatable :: monomer1_ene, monomer2_ene
        double precision, dimension(:, :), allocatable :: monomer1_coeff, monomer2_coeff
        double precision, dimension(:, :), allocatable :: fragments_coeff
        double precision, dimension(:, :), allocatable :: inter_Hamilton_non_ortho, inter_overlap_non_ortho
        double precision, dimension(:, :), allocatable :: effective_coupling
        double precision, dimension(:, :), allocatable :: effective_ene1, effective_ene2

    contains
        private
        procedure :: prepare_files
        procedure :: read_dimer_info
        procedure :: read_monomer1_info
        procedure :: read_monomer2_info
        procedure :: read_dimer_ene
        procedure :: read_dimer_coeff
        procedure :: read_monomer1_ene
        procedure :: read_monomer1_coeff
        procedure :: read_monomer2_ene
        procedure :: read_monomer2_coeff
        procedure :: generate_fragments_coeff
        procedure :: read_dimer_overlap
        procedure :: read_dimer_Fock
        procedure :: close_files
        procedure :: check

        procedure, public :: set_prefixes
        procedure, public :: read
        procedure, public :: calculate
        procedure, public :: show

        final :: coupling_calculator_destructor

    end type coupling_calculator

    interface coupling_calculator
        module procedure :: coupling_calculator_constructor_def
        module procedure :: coupling_calculator_constructor_args
    end interface coupling_calculator

    interface
        module function coupling_calculator_constructor_def() result(this)
            type(coupling_calculator) :: this
        end function coupling_calculator_constructor_def

        module function coupling_calculator_constructor_args(dimer_prefix, monomer1_prefix, monomer2_prefix) result(this)
            character(len=128), intent(in) :: dimer_prefix, monomer1_prefix, monomer2_prefix
            type(coupling_calculator) :: this
        end function coupling_calculator_constructor_args

        module subroutine coupling_calculator_destructor(this)
            implicit none
            type(coupling_calculator), intent(inout) :: this
        end subroutine coupling_calculator_destructor

        module subroutine set_prefixes(this, dimer_prefix, monomer1_prefix, monomer2_prefix)
            implicit none
            class(coupling_calculator), intent(inout) :: this
            character(len=*), intent(in) :: dimer_prefix, monomer1_prefix, monomer2_prefix
        end subroutine set_prefixes

        module subroutine check(this)
            implicit none
            class(coupling_calculator), intent(in) :: this
        end subroutine check

        module subroutine generate_fragments_coeff(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine generate_fragments_coeff

        module subroutine read(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read

        module subroutine show(this)
            implicit none
            class(coupling_calculator), intent(in) :: this
        end subroutine show

        module subroutine prepare_files(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine prepare_files

        module subroutine close_files(this)
            implicit none
            class(coupling_calculator), intent(in) :: this
        end subroutine close_files

        module subroutine read_basic_info(ifile, natoms, nelectrons, nbasis)
            implicit none
            integer, intent(in) :: ifile
            integer, intent(out) :: natoms, nelectrons, nbasis
        end subroutine read_basic_info

        module subroutine read_dimer_info(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_dimer_info

        module subroutine read_monomer1_info(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer1_info

        module subroutine read_monomer2_info(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer2_info

        module subroutine read_gau_vector(ret, ifile, rows, label)
            implicit none
            double precision, allocatable, intent(out) :: ret(:)
            integer, intent(in) :: ifile, rows
            character(len=*), intent(in) :: label
        end subroutine read_gau_vector

        module subroutine read_gau_matrix(ret, ifile, rows, cols, label)
            implicit none
            double precision, allocatable, intent(out) :: ret(:, :)
            integer, intent(in) :: ifile, rows, cols
            character(len=*), intent(in) :: label
        end subroutine read_gau_matrix

        module subroutine read_gau_lower_triangle_matrix(ret, ifile, rows, label, use_last_if_exist)
            ! assume square matrix
            implicit none
            double precision, allocatable, intent(out) :: ret(:, :)
            integer, intent(in) :: ifile, rows
            character(len=*), intent(in) :: label
            logical, intent(in), optional :: use_last_if_exist
        end subroutine read_gau_lower_triangle_matrix

        module subroutine read_dimer_ene(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_dimer_ene

        module subroutine read_monomer1_ene(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer1_ene

        module subroutine read_monomer2_ene(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer2_ene

        module subroutine read_dimer_coeff(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_dimer_coeff

        module subroutine read_monomer1_coeff(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer1_coeff

        module subroutine read_monomer2_coeff(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_monomer2_coeff

        module subroutine read_dimer_overlap(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_dimer_overlap

        module subroutine read_dimer_Fock(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine read_dimer_Fock

        module subroutine calculate(this)
            implicit none
            class(coupling_calculator), intent(inout) :: this
        end subroutine calculate

    end interface

end module coupling_module

