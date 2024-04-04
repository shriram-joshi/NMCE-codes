program main
    
    use kind
    use data
    use library
    use meshing
    use assembly
    use ui
    use omp_lib

    implicit none

    call init_problem()

    if(bft == 1) then
        call lin_basis_phi1d()
    else if (bft == 2) then
        call quad_basis_phi1d()
    end if
    
    call generate_varmesh()

    call global_assembly()

    call print_problem_setup()

    call FullGaussSolverp(JG, RG, nVar)

    print*, "Solution, c(x) ="
    call print_mtrx(RG, 1, size(RG))

    call deallocate_mem()

contains

    subroutine deallocate_mem()
        deallocate(jL)
        deallocate(rL)
        deallocate(RG)
        deallocate(JG)
        deallocate(xMesh)
        deallocate(ph1)
        deallocate(dph1)
    end subroutine deallocate_mem

end program
