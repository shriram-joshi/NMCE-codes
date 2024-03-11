program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: i, j, k, l
    real(kind=rk), dimension(3) :: dxdxi
    Real*8, allocatable :: u(:), Jtemp(:,:)

    call Init_problem()
    call Gauss_points()
    call lin_basis_phi1d()
    ! Creates a array of mesh points
    call generate_uniform_mesh()

    allocate(u(n+1))
    allocate(Jtemp(n+1, n+2))

    do k = 1,n ! Change according to basis function used

        ! Calculate dx/dxi
        dxdxi = 0
        do j=1,3
            do i=1,nl
                dxdxi(j) = dxdxi(j) + xMesh(k+i-1)*dph1(i,j) 
            end do
        end do

        ! Uncomment to print the dx/dxi term
        ! print*, "dx/dxi = "
        ! print*, dxdxi

        ! Calculate the local j matrix
        jL = 0
        do i = 1,nl
            do j = 1,nl
                do l = 1,3 ! Sum for different quadrature points
                    jL(i,j) = jL(i,j) + wei(l)*(dph1(i,l)*dph1(j,l)/dxdxi(l) - Pe*ph1(i,l)*dph1(j,l))
                end do
            end do
        end do

        print*, "Local j ="
        do i = 1, nl
            do j = 1, nl
                write(*, fmt='(F15.2)', advance='no') jL(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do

        ! Add local j to the global J matrix [Ju = R]
        JG(k:k+nl-1,k:k+nl-1) = JG(k:k+nl-1,k:k+nl-1) + jL

        ! You can assemble the R matrix here as well if needed
    end do

    ! Changing J and R according to the BCs
    R(1) = uBC(1)
    R(n+1) = uBC(2)

    JG(1,1:n+1) = 0.0_rk
    JG(n+1,1:n+1) = 0.0_rk
    JG(1,1) = 1.0_rk
    JG(n+1,n+1) = 1.0_rk

    print*, "Global J ="
    do i = 1, n+1
        do l = 1, n+1
            write(*, fmt='(F15.2)', advance='no') JG(i, l)
        end do
        print *  ! Move to the next line after printing each row
    end do

    ! print*, "R ="
    ! do i = 1, n+1
    !     write(*, fmt='(F15.2)', advance='no') R(i)
    ! end do

    u = real(R)
    Jtemp = real(JG)

    call FullGaussSolverp(Jtemp, u, n+1)

    print*, "u ="
    do i = 1, n+1
        write(*, fmt='(F15.4)', advance='no') u(i)
    end do
        
    ! Output data into a file 
    ! open(1, file = 'SSCD-Pe_x.dat', status='new')  
    ! do i = 1,n+1  
    !    write(1,*) u(i)   
    ! end do  
    ! close(1)

    deallocate(jL)
    deallocate(R)
    deallocate(JG)
    deallocate(Jtemp)
    deallocate(u)
end program
