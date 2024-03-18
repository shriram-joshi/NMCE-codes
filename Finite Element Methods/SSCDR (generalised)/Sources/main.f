program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: i, j, k, l
    real(kind=rk), dimension(3) :: dxdxi

    call Init_problem()
    call Gauss_points()

    if(nl == 2) then
        call lin_basis_phi1d()
    else if (nl == 3) then
        call quad_basis_phi1d()
    end if
    
    call generate_unifmesh()

    do k = 1,nVar-1,nl-1 ! Change according to basis function used

        ! Calculate dx/dxi
        dxdxi = 0
        do l=1,3
            do i=1,nl
                dxdxi(l) = dxdxi(l) + xMesh(k+i-1)*dph1(i,l) 
            end do
        end do

        ! Calculate the local j matrix
        jL = 0
        do i = 1,nl
            do j = 1,nl
                do l = 1,3 ! Sum for different quadrature points
                    jL(i,j) = jL(i,j) + wei(l)*(dph1(i,l)*dph1(j,l)/dxdxi(l) + Pe*ph1(i,l)*dph1(j,l))
                    ! print*, "jL =", jL(i,j)    
                end do
            end do
        end do

        ! Add local j to the global J matrix [Ju = R]
        JG(k:k+nl-1,k:k+nl-1) = JG(k:k+nl-1,k:k+nl-1) + jL

        ! You can assemble the R matrix here as well if needed
    end do

    ! Changing J and R according to the BCs
    R(1) = uBC(1)
    R(nVar) = uBC(2)

    JG(1,1:nVar) = 0.0_rk
    JG(nVar,1:nVar) = 0.0_rk
    JG(1,1) = 1.0_rk
    JG(nVar,nVar) = 1.0_rk

    call printProblemSetup()

    call FullGaussSolverp(JG, R, nVar)

    print*, "Solution, c(x) ="
    do i = 1, nVar
        write(*, fmt='(F15.4)', advance='no') R(i)
    end do
    print*
        
    ! Output data into a file 
    ! open(1, file = 'SSCDVar-Pe_50.dat', status='new')  
    ! do i = 1,nVar  
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
