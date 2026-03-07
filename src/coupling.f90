submodule (coupling_module) coupling
    implicit none

contains

    module subroutine calculate(this)
        implicit none
        class(coupling_calculator), intent(inout) :: this

        double precision, allocatable :: Hamilton_diag(:, :)
        integer :: i, n

        n = this%nbasis_dimer

        allocate(this%inter_Hamilton_non_ortho(n, n))
        allocate(this%inter_overlap_non_ortho(n, n))

        allocate(Hamilton_diag(n, n))

        allocate(this%effective_coupling(n, n))
        allocate(this%effective_ene1(n, n))
        allocate(this%effective_ene2(n, n))

        ! first use Hamilton_diag as a temporary matrix for matmul
        ! inter_Hamilton_non_ortho = fragments_coeff.T @ dimer_Fock @ fragments_coeff
        ! inter_overlap_non_ortho = fragments_coeff.T @ dimer_overlap @ fragments_coeff
        call dgemm('T', 'N', n, n, n, &
            1.D0, this%fragments_coeff, n, this%dimer_Fock, n, &
            0.D0, Hamilton_diag, n)
        call dgemm('N', 'N', n, n, n, &
            1.D0, Hamilton_diag, n, this%fragments_coeff, n, &
            0.D0, this%inter_Hamilton_non_ortho, n)
        call dgemm('T', 'N', n, n, n, &
            1.D0, this%fragments_coeff, n, this%dimer_overlap, n, &
            0.D0, Hamilton_diag, n)
        call dgemm('N', 'N', n, n, n, &
            1.D0, Hamilton_diag, n, this%fragments_coeff, n, &
            0.D0, this%inter_overlap_non_ortho, n)

        forall (i = 1:n) Hamilton_diag(i, :) = this%inter_Hamilton_non_ortho(i, i)

        ! e_{i,j}^\text{eff}=\frac{1}{2}\frac{(e_i+e_j)-2J_{ij}S_{ij}\pm(e_i-e_j)\sqrt{1-{S_{ij}}^2}}{1-{S_{ij}}^2}
        ! J_{ij}^\text{eff}=\frac{J_{ij}-\frac{1}{2}(e_i+e_j)S_{ij}}{1-{S_{ij}}^2}
        this%effective_coupling = (this%inter_Hamilton_non_ortho - &
            (Hamilton_diag + transpose(Hamilton_diag)) * this%inter_overlap_non_ortho / 2.D0) / &
            (1. - this%inter_overlap_non_ortho ** 2)
        this%effective_ene1 = ((Hamilton_diag + transpose(Hamilton_diag)) - &
            2.D0 * this%inter_Hamilton_non_ortho * this%inter_overlap_non_ortho + &
            (Hamilton_diag - transpose(Hamilton_diag)) * &
            (1. - this%inter_overlap_non_ortho ** 2) ** .5D0) / &
            (1. - this%inter_overlap_non_ortho ** 2) / 2.D0
        this%effective_ene2 = ((Hamilton_diag + transpose(Hamilton_diag)) - &
            2.D0 * this%inter_Hamilton_non_ortho * this%inter_overlap_non_ortho - &
            (Hamilton_diag - transpose(Hamilton_diag)) * &
            (1. - this%inter_overlap_non_ortho ** 2) ** .5D0) / &
            (1. - this%inter_overlap_non_ortho ** 2) / 2.D0

        deallocate(Hamilton_diag)

    end subroutine calculate

end submodule coupling
