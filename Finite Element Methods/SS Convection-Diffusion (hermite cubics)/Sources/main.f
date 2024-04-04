program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: i, j, k, l, t
    real(kind=rk), dimension(3) :: dxdxi
    real(kind=rk), dimension(4) :: psi, dpsi

    call Init_problem()
    call Gauss_points()

    call hermite_basis_phi1d() 
    call generate_unifmesh()

    do k = 1,nVar-2,2 ! Change according to basis function used

        ! Calculate dx/dxi
        dxdxi = 0
        do l=1,3
            do i=1,2 ! Hermite cubics only need two points per element
                dxdxi(l) = dxdxi(l) + xMesh(k+i-1)*dph(2*i-1,l) 
            end do
        end do

        print*, "dxdxi ="
        do t = 1, 3
            write(*, fmt='(F15.2)', advance='no') dxdxi(t)
        end do
        print *  ! Move to the next line after printing each row

        ! Calculate the local j matrix
        jL = 0
        do i = 1,nl
            do j = 1,nl
                do l = 1,3 ! Sum for different quadrature points

                    psi = ph(:,l)
                    psi(2) = dxdxi(l)*ph(2,l)
                    psi(4) = dxdxi(l)*ph(4,l)

                    dpsi = dph(:,l)
                    dpsi(2) = dxdxi(l)*dph(2,l)
                    dpsi(4) = dxdxi(l)*dph(4,l)

                    print*, "dpsi ="
                    do t = 1, 4
                        write(*, fmt='(F15.2)', advance='no') dpsi(t)
                    end do
                    print *  ! Move to the next line after printing each row

                    jL(i,j) = jL(i,j) + wei(l)*(dph(i,l)*dpsi(j)/dxdxi(l) + Pe*ph(i,l)*dpsi(j))   
                end do
            end do
        end do

        ! print*, "jL ="
        ! do i = 1, nl
        !     do l = 1, nl
        !         write(*, fmt='(F15.2)', advance='no') jL(i, l)
        !     end do
        !     print *  ! Move to the next line after printing each row
        ! end do

        ! Add local j to the global J matrix [Ju = R]
        JG(k:k+nl-1,k:k+nl-1) = JG(k:k+nl-1,k:k+nl-1) + jL

        ! You can assemble the R matrix here as well if needed
    end do

    ! Changing J and R according to the BCs
    R(1) = uBC(1)
    R(nVar-1) = uBC(2)

    JG(1,1:nVar) = 0.0_rk
    JG(1,1) = 1.0_rk
    JG(nVar-1,1:nVar) = 0.0_rk
    JG(nVar-1,nVar-1) = 1.0_rk

    call printProblemSetup()

    call FullGaussSolverp(JG, R, nVar)

    print*, "Solution, c(x) ="
    do i = 1, nVar
        write(*, fmt='(F15.4)', advance='no') R(i)
    end do
    print*
        
    ! Output data into a file 
    ! open(1, file = 'SSCDVar-Pe_50.dat', status='new')  
    ! do i = 1,n+1  
    !    write(1,*) u(i)   
    ! end do  
    ! close(1)

    deallocate(jL)
    deallocate(R)
    deallocate(JG)

contains

    subroutine printProblemSetup()

        print*, "Global J ="
        do i = 1, nVar
            do l = 1, nVar
                write(*, fmt='(F15.2)', advance='no') JG(i, l)
            end do
            print *  ! Move to the next line after printing each row
        end do

        print*, "R ="
        do i = 1, nVar
            write(*, fmt='(F15.2)', advance='no') R(i)
        end do
        print*

        ! Uncomment to print out mesh points
        print*, "Mesh points, x ="
        do i = 1, n+1
            write(*, fmt='(F15.4)', advance='no') xMesh(i)
        end do
        print*
    end subroutine printProblemSetup
end program
