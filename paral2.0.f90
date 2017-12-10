module Homework
    implicit none
    include "mpif.h"
    
    contains
        subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
        implicit none
	include "mpif.h"
        integer :: mpisize, mpiErr,rank,comm,j,i,i1,i2,request,i4
	integer(4):: mpiRank
        real(8),  dimension(:,:) :: A
        integer :: x1, y1, x2, y2
        integer(4) :: n, L, R, Up, Down, m, tmp
        real(8), allocatable :: current_column(:), B(:,:)
        real(8) :: current_sum, max_sum
        logical :: transpos 
        real(8),allocatable,dimension(:)::max_sumar
	integer(4), allocatable, dimension(:):: X_1, X_2, Y_1, Y_2, coords
	integer(4), dimension(MPI_STATUS_SIZE) :: status


        m = size(A, dim=1) 
        n = size(A, dim=2) 
        transpos = .FALSE.
	allocate(coords(4))


        if (m < n) then 
            transpos = .TRUE.   
            B = transpose(A)
            m = size(B, dim=1) 
            n = size(B, dim=2) 
        else
            B = A     
            endif

        allocate(current_column(m))

 
        max_sum=-huge(0d0)
        x1=1
        y1=1
        x2=1
        y2=1

!	call mpi_init(mpiErr)
	call mpi_comm_size(MPI_COMM_WORLD, mpisize, mpiErr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
	rank=mpiRank+1
	do  L=rank,n,mpisize
		current_column = B(:, L) 
              do R=L,n
 
                     if (R > L) then 
                    	current_column = current_column + B(:, R)
                     endif
                
		     call FindMaxInArray(current_column, current_sum, Up, Down) 
                     if (current_sum > max_sum) then
                    	max_sum = current_sum
                    	x1=Up
                    	x2=Down
                    	y1=L
                    	y2=R
                     endif 
 	      enddo
   	enddo  
     if(mpiSize>1) then
	if (mpiRank/=0) then	
		do i1=1,mpisize
		  call mpi_send(max_sum,1,MPI_REAL8,0,5,MPI_COMM_WORLD,mpiErr)
		enddo
	endif 

	call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
	if(mpiRank==0)then
		allocate(max_sumar(mpisize))
		max_sumar(1)=max_sum
		do i=1,mpisize-1
		   call mpi_recv(max_sum,1,MPI_REAL8,i,MPI_ANY_TAG,MPI_COMM_WORLD,status,mpiErr)
		   max_sumar(i+1)=max_sum
		enddo

	

		max_sum=maxval(max_sumar)
		do j=1,mpisize
	  		if(max_sumar(j)==max_sum) then
				exit
	  		endif
		enddo
		deallocate(max_sumar) 
	endif

	
	call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
	if(mpiRank==0)then
		do i4=1,mpisize-1
		  call MPI_SEND(j,1,MPI_INTEGER4,i4,571,MPI_COMM_WORLD,mpiErr)
		enddo
	else
		call MPI_RECV(j,1,MPI_INTEGER4,0,571,MPI_COMM_WORLD,status,mpiErr)
	endif

	call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
	i2=0
	 if(mpiRank==j-1) then
		do i2=0,mpisize-1
		   if(i2/=j-1) then
		call MPI_SEND(x1,1,MPI_INTEGER4,i2,17,MPI_COMM_WORLD,mpiErr)
		call MPI_SEND(x2,1,MPI_INTEGER4,i2,18,MPI_COMM_WORLD,mpiErr)
		call MPI_SEND(y1,1,MPI_INTEGER4,i2,19,MPI_COMM_WORLD,mpiErr)
		call MPI_SEND(y2,1,MPI_INTEGER4,i2,20,MPI_COMM_WORLD,mpiErr)
		   endif
		enddo
	else
		call MPI_RECV(x1,1,MPI_INTEGER4,j-1,17,MPI_COMM_WORLD,status,mpiErr)
		call MPI_RECV(x2,1,MPI_INTEGER4,j-1,18,MPI_COMM_WORLD,status,mpiErr)
		call MPI_RECV(y1,1,MPI_INTEGER4,j-1,19,MPI_COMM_WORLD,status,mpiErr)
		call MPI_RECV(y2,1,MPI_INTEGER4,j-1,20,MPI_COMM_WORLD,status,mpiErr)
	endif
    endif 
		


	deallocate(current_column)

       	 if (transpos) then  
            tmp = x1
            x1 = y1
            y1 = tmp
    
           tmp = y2
            y2 = x2
            x2 = tmp
         endif
	
!		call MPI_FINALIZE(mpiErr)
        end subroutine


        subroutine FindMaxInArray(a, Sum, Up, Down)
            real(8), intent(in), dimension(:) :: a
            integer(4), intent(out) :: Up, Down
            real(8), intent(out) :: Sum
            real(8) :: cur_sum
            integer(4) :: minus_pos, i

            Sum = a(1)
            Up = 1
            Down = 1
            cur_sum = 0
            minus_pos = 0


   
            do i=1, size(a)
                cur_sum = cur_sum + a(i)
            	if (cur_sum > Sum) then
                	Sum = cur_sum
                	Up = minus_pos + 1
               		 Down = i
                endif
         
            	if (cur_sum < 0) then
            	    cur_sum = 0
            	    minus_pos = i
                endif

            enddo
             

        end subroutine FindMaxInArray


end module Homework



