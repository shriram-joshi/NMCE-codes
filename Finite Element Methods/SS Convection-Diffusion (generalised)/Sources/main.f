program main
    use kind
    use library
    use Data_
    use omp_lib

    implicit none

    integer :: i, l, k
    Real*8, allocatable :: u(:), Jtemp(:,:)

    call Gauss_points()
    call basis_phi1d()
    call Init_problem()
    ! Creates a array of mesh points
    call generate_mesh()

    ! allocate(u(n+1))
    ! allocate(Jtemp(n+1, n+2))

    ! ! Calculate fixed step size
    ! h = (xSpan(2) - xSpan(1))/real(n, rk);
    ! print*, "h = ", h

    ! ! The local j matrix
    ! jL = reshape((/ 1.0_rk/h-Pe/2.0_rk, -1.0_rk/h-Pe/2.0_rk, -1.0_rk/h+Pe/2.0_rk, 1.0_rk/h+Pe/2.0_rk /), shape(jL))

    ! ! Checks if jl matrix is of correct size
    ! if (nl /= size(jL)/size(jL(1,:))) then
    !     print*, "The local j matrix was not initialized correctly. The size of jl is not the same as the specified nl."
    !     stop
    ! end if

    ! print*, "Local j ="
    ! do i = 1, nl
    !     do l = 1, nl
    !         write(*, fmt='(F15.2)', advance='no') jL(i, l)
    !     end do
    !     print *  ! Move to the next line after printing each row
    ! end do

    ! do k = 1,n,1
    !     ! Assemble the global J matrix [Ju = R]
    !     JG(k:k+nl-1,k:k+nl-1) = JG(k:k+nl-1,k:k+nl-1) + jL

    !     ! You can assemble the R matrix here as well if needed
    ! end do

    ! ! Changing J and R according to the BCs
    ! R(1) = uBC(1)
    ! R(n+1) = uBC(2)

    ! JG(1,1:n+1) = 0.0_rk
    ! JG(n+1,1:n+1) = 0.0_rk
    ! JG(1,1) = 1.0_rk
    ! JG(n+1,n+1) = 1.0_rk

    ! print*, "Global J ="
    ! do i = 1, n+1
    !     do l = 1, n+1
    !         write(*, fmt='(F15.2)', advance='no') JG(i, l)
    !     end do
    !     print *  ! Move to the next line after printing each row
    ! end do

    ! print*, "R ="
    ! do i = 1, n+1
    !     write(*, fmt='(F15.2)', advance='no') R(i)
    ! end do

    ! u = real(R)
    ! Jtemp = real(JG)

    ! call FullGaussSolverp(Jtemp, u, n+1)

    ! print*, "u ="
    ! do i = 1, n+1
    !     write(*, fmt='(F15.4)', advance='no') u(i)
    ! end do
        
    ! ! Output data into a file 
    ! ! open(1, file = 'SSCD-Pe_x.dat', status='new')  
    ! ! do i = 1,n+1  
    ! !    write(1,*) u(i)   
    ! ! end do  
    ! ! close(1)

    ! deallocate(jL)
    ! deallocate(R)
    ! deallocate(JG)
    ! deallocate(Jtemp)
    ! deallocate(u)
end program
