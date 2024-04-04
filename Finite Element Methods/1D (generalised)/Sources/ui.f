module ui

    use kind
    use data
    implicit none
    
contains
    
    ! Prints any general array with size n x m
    subroutine print_mtrx(mtrx, n, m)
    
        implicit none
        integer :: n, m, i, j
        real(kind=rk), dimension(n,m) :: mtrx

        do i = 1, n
            do j = 1, m
                ! Prints upto 4 decimal places
                write(*, fmt='(F15.4)', advance='no') mtrx(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do
            
    end subroutine print_mtrx

    subroutine output_data()
        
        implicit none
        integer :: i
        
        ! Output data into a file 
        open(1, file = 'SSCDVar-Pe_50.dat', status='new')  
        do i = 1,nVar  
           write(1,*) RG(i)   
        end do  
        close(1)
            
    end subroutine output_data

    subroutine print_problem_setup()

        implicit none
        
        print*, "nop ="
        call print_mtrx(real(nop, kind=rk), size(nop(:,1)), size(nop(1,:)))

        print*, "Global J ="
        call print_mtrx(JG, size(JG(:,1)), size(JG(1,:)))

        print*, "Global R ="
        call print_mtrx(RG, 1, size(RG))

        print*, "Mesh ="
        call print_mtrx(xMesh, 1, size(xMesh))

    end subroutine print_problem_setup

end module ui