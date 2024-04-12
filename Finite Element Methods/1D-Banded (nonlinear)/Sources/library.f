module library

    use kind
    use data
    use omp_lib
    implicit none
contains
!***********************************************************************************************************************
! Linear Basis functions
!***********************************************************************************************************************
subroutine lin_basis_phi1d()
    integer :: i

    ph=0.0_rk
    dph=0.0_rk
    do i=1,3
        ph(1,i) = 1.0_rk - xiarr(i)
        dph(1,i)= -1.0_rk
        ph(2,i) = xiarr(i)
        dph(2,i)= 1.0_rk
    end do
end subroutine lin_basis_phi1d
!***********************************************************************************************************************
! Linear Basis functions End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Quadratic Basis functions
!***********************************************************************************************************************
subroutine quad_basis_phi1d()
    implicit none
    integer :: i

    ph=0.0_rk
    dph=0.0_rk
    do i=1,3
        ph(1,i) = 1.0_rk - 3.0_rk*xiarr(i) + 2.0_rk*xiarr(i)*xiarr(i)
        dph(1,i)= -3.0_rk + 4.0_rk*xiarr(i)
        ph(2,i) = 4.0_rk*(xiarr(i) - xiarr(i)*xiarr(i))
        dph(2,i)= 4.0_rk - 8.0_rk*xiarr(i)
        ph(3,i) = - xiarr(i) + 2.0_rk*xiarr(i)*xiarr(i)
        dph(3,i)= -1.0_rk+4.0_rk*xiarr(i)
    end do
end subroutine quad_basis_phi1d
!***********************************************************************************************************************
! Quadratic Basis functions End
!***********************************************************************************************************************

integer function bgbl(o,p,ij)
    implicit none
    integer, intent(in) :: o, ij
    integer, intent(in) :: p
        ! ij = 1 implies o, p are for equation number (i.e. row # in square coefficient matrix)
        ! ij = 2 implies o, p are for variable number (i.e. col # in square coefficient matrix)
        if (ij==1) then
            bgbl = o
        elseif (ij==2) then
            bgbl = p - o + (bWidth - 1)/2 + 1
        end if
    return
end function bgbl
end module library
