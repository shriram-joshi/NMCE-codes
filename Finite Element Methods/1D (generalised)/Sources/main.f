program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: k

    call init_problem()
    
    call Gauss_points()

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
    do k = 1, nVar
        write(*, fmt='(F15.4)', advance='no') RG(k)
    end do
    print*
        
    ! Output data into a file 
    ! open(1, file = 'SSCDVar-Pe_50.dat', status='new')  
    ! do k = 1,nVar  
    !    write(1,*) R(k)   
    ! end do  
    ! close(1)

    call deallocate_mem()

contains

    subroutine print_problem_setup()

        integer:: i, j

        print*, "Global J ="
        do i = 1, nVar
            do j = 1, nVar
                write(*, fmt='(F15.2)', advance='no') JG(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do

        print*, "R ="
        do i = 1, nVar
            write(*, fmt='(F15.2)', advance='no') RG(i)
        end do
        print*

        ! Uncomment to print out mesh points
        print*, "Mesh points, x ="
        do i = 1, nVar
            write(*, fmt='(F15.4)', advance='no') xMesh(i)
        end do
        print*
    end subroutine print_problem_setup

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
