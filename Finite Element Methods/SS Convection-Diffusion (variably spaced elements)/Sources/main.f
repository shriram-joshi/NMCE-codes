program main
    use kind
    use ns1d_slender
    use Data_
    use omp_lib

    implicit none

    call Gauss_points()
    call basis_phi1d()
    print*, ph1(2,:)
    
end program
