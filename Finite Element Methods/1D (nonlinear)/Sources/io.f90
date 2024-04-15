module io

    use kind
    use, intrinsic :: ISO_Fortran_env
    implicit none
    
contains
    
    ! Prints any general array with size n x m
    subroutine print_mtrx(mtrx)
    
        implicit none
        integer :: rows, cols, i, j
        real(kind=rk), intent(in) :: mtrx(:,:)

        rows = size(mtrx, 1)
        cols = size(mtrx, 2)

        do i = 1, rows
            do j = 1, cols
                ! Prints upto 4 decimal places
                write(*, fmt='(F15.4)', advance='no') mtrx(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do
            
    end subroutine print_mtrx

    subroutine output_data(mtrx, fname)
        implicit none
        real(kind=rk), intent(in) :: mtrx(:,:)
        character(len=*), intent(in) :: fname
        integer :: rows, cols, i, j
        character(len=100) :: fname_
        logical :: file_exists
    
        ! Determine the dimensions of the matrix
        rows = size(mtrx(:,1))

        ! Concatenate the file name with ".dat" extension
        fname_ = trim(fname) // '.dat'
    
        ! Check if file already exists
        inquire(file=fname_, exist=file_exists)
        if (file_exists) then
            print *, 'File already exists. Overwriting ', fname_
            ! Delete the existing file
            call system('rm ' // fname_)
        endif

    
        ! Open the file for writing
        open(unit=10, file=fname_, status='unknown')
    
        ! Write the data to the file
        do i = 1, rows
            write(10, '(*(E10.5,1x))') mtrx(i, :)
        end do
    
        ! Close the file
        close(unit=10)
    
        print *, 'Data has been written to ', fname_
    end subroutine output_data    

end module io