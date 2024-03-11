program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: i, l, k
    Real*8, allocatable :: u(:), Jtemp(:,:)
    real(kind=rk) :: h

    call Gauss_points()
    call basis_phi1d()
    call Init_problem()
    ! Creates a array of mesh points
    call generate_mesh()

    allocate(u(n+1))
    allocate(Jtemp(n+1, n+2))

    do k = 1,n,1

        h = xMesh(k+1)-xMesh(k)

        ! The local j matrix
        jl = reshape((/ 1.0_rk/h - Pe/2.0_rk, -1.0_rk/h - Pe/2.0_rk, -1.0_rk/h + Pe/2.0_rk, 1.0_rk/h + Pe/2.0_rk /), shape(jl))

        ! Assemble the global J matrix [Ju = R]
        J(k:k+nl-1,k:k+nl-1) = J(k:k+nl-1,k:k+nl-1) + jl

        ! You can assemble the R matrix here as well if needed
    end do

    ! Changing J and R according to the BCs
    R(1) = ubc(1)
    R(n+1) = ubc(2)

    J(1,1:n+1) = 0.0_rk
    J(n+1,1:n+1) = 0.0_rk
    J(1,1) = 1.0_rk
    J(n+1,n+1) = 1.0_rk

    print*, "Global J ="
    do i = 1, n+1
        do l = 1, n+1
            write(*, fmt='(F15.2)', advance='no') J(i, l)
        end do
        print *  ! Move to the next line after printing each row
    end do

    print*, "R ="
    do i = 1, n+1
        write(*, fmt='(F15.2)', advance='no') R(i)
    end do
    print *

    u = real(R)
    Jtemp = real(J)

    call FullGaussSolverp(Jtemp, u, n+1)

    print*, "u ="
    do i = 1, n+1
        write(*, fmt='(F15.4)', advance='no') u(i)
    end do
    print *
        
    ! ! Output data into a file 
    ! open(1, file = 'SSCDVar.dat', status='new')  
    ! do i = 1,n+1  
    !    write(1,*) u(i)   
    ! end do  
    ! close(1)

    deallocate(jl)
    deallocate(R)
    deallocate(J)
    deallocate(Jtemp)
    deallocate(u)
end program
