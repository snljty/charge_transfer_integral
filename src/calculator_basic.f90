submodule (coupling_module) calculator_basic
    implicit none

contains

    module function coupling_calculator_constructor_def() result(this)
        type(coupling_calculator) :: this
        call this%set_prefixes("dimer", "monomer1", "monomer2")

    end function coupling_calculator_constructor_def

    module function coupling_calculator_constructor_args(dimer_prefix, monomer1_prefix, monomer2_prefix) result(this)
        character(len=128), intent(in) :: dimer_prefix, monomer1_prefix, monomer2_prefix
        type(coupling_calculator) :: this
        call this%set_prefixes(dimer_prefix, monomer1_prefix, monomer2_prefix)

    end function coupling_calculator_constructor_args

    module subroutine set_prefixes(this, dimer_prefix, monomer1_prefix, monomer2_prefix)
        implicit none
        class(coupling_calculator), intent(inout) :: this
        character(len=*), intent(in) :: dimer_prefix, monomer1_prefix, monomer2_prefix

        this%dimer_prefix = dimer_prefix
        this%monomer1_prefix = monomer1_prefix
        this%monomer2_prefix = monomer2_prefix

    end subroutine set_prefixes

    module subroutine check(this)
        implicit none
        class(coupling_calculator), intent(in) :: this

        if (this%natoms_dimer /= this%natoms_monomer1 + this%natoms_monomer2) &
            error stop "amount of atoms in dimer mismatches that in monomers."
        if (this%nelectrons_dimer /= this%nelectrons_monomer1 + this%nelectrons_monomer2) &
            error stop "amount of electrons in dimer mismatches that in monomers."
        if (this%nbasis_dimer /= this%nbasis_monomer1 + this%nbasis_monomer2) &
            error stop "amount of basis functions in dimer mismatches that in monomers."

    end subroutine check

    module subroutine generate_fragments_coeff(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        allocate(this%fragments_coeff(this%nbasis_dimer, this%nbasis_dimer))
        this%fragments_coeff = 0.D0
        this%fragments_coeff(:this%nbasis_monomer1, :this%nbasis_monomer1) = this%monomer1_coeff
        this%fragments_coeff(this%nbasis_monomer1 + 1:, this%nbasis_monomer1 + 1:) = this%monomer2_coeff

    end subroutine generate_fragments_coeff

    module subroutine read(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        call this%prepare_files()
        call this%read_dimer_info()
        call this%read_monomer1_info()
        call this%read_monomer2_info()
        call this%check()
        call this%read_monomer1_ene()
        call this%read_monomer1_coeff()
        call this%read_monomer2_ene()
        call this%read_monomer2_coeff()
        call this%generate_fragments_coeff()
        call this%read_dimer_ene()
        call this%read_dimer_coeff()
        call this%read_dimer_overlap()
        call this%read_dimer_Fock()
        call this%close_files()

    end subroutine read

    module subroutine show(this)
        implicit none
        class(coupling_calculator), intent(in) :: this

        double precision, parameter :: scaler = 27.211407953D3
        write(*, '(a, f8.3, a)') "HOMO-HOMO transfer integral: ", &
            this%effective_coupling(this%iHOMO_monomer1, this%iHOMO_monomer2 + this%nbasis_monomer1) * scaler, " meV"
        write(*, '(a, f8.3, a)') "LUMO-LUMO transfer integral: ", &
            this%effective_coupling(this%iLUMO_monomer1, this%iLUMO_monomer2 + this%nbasis_monomer1) * scaler, " meV"

    end subroutine show

    module subroutine coupling_calculator_destructor(this)
        implicit none
        type(coupling_calculator), intent(inout) :: this

        ! will not actually be called unless used on dynamic allocation
        if (allocated(this%dimer_ene)) deallocate(this%dimer_ene)
        if (allocated(this%dimer_Fock)) deallocate(this%dimer_Fock)
        if (allocated(this%dimer_overlap)) deallocate(this%dimer_overlap)
        if (allocated(this%dimer_coeff)) deallocate(this%dimer_coeff)
        if (allocated(this%monomer1_ene)) deallocate(this%monomer1_ene)
        if (allocated(this%monomer2_ene)) deallocate(this%monomer2_ene)
        if (allocated(this%monomer1_coeff)) deallocate(this%monomer1_coeff)
        if (allocated(this%monomer2_coeff)) deallocate(this%monomer2_coeff)
        if (allocated(this%fragments_coeff)) deallocate(this%fragments_coeff)
        if (allocated(this%inter_Hamilton_non_ortho)) deallocate(this%inter_Hamilton_non_ortho)
        if (allocated(this%inter_overlap_non_ortho)) deallocate(this%inter_overlap_non_ortho)
        if (allocated(this%effective_coupling)) deallocate(this%effective_coupling)
        if (allocated(this%effective_ene1)) deallocate(this%effective_ene1)
        if (allocated(this%effective_ene2)) deallocate(this%effective_ene2)
    end subroutine coupling_calculator_destructor

end submodule calculator_basic

