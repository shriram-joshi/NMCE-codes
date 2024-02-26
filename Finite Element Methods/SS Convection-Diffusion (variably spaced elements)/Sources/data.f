Module Data_

    use kind
    use omp_lib
    implicit none

    ! Declare Global variables that will be frequently used during simulation
    ! Examples are:
    real(kind=rk) :: C1, C2
    real(kind=rk) :: Oh, Gbond
    real(kind=rk) :: ph1(3,6),dph1(3,6)
    real(kind=rk), dimension(6) :: xiarr
    real(kind=rk), dimension(3) :: wei
contains

    ! You can make subroutines here
    subroutine Gauss_points()
        implicit none

        xiarr(1)=(1.0_rk-0.774596669241483_rk)*0.5_rk
        xiarr(2)=0.5_rk
        xiarr(3)=(1.0_rk+0.774596669241483_rk)*0.5_rk
        xiarr(4)=0.0_rk
        xiarr(5)=0.5_rk
        xiarr(6)=1.0_rk
        wei(1)=(0.555555555555556_rk)*0.5_rk
        wei(2)=(0.888888888888889_rk)*0.5_rk
        wei(3)=(0.555555555555556_rk)*0.5_rk

    end subroutine Gauss_points
end module Data_
