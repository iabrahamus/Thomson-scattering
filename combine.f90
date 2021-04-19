        PROGRAM combine
        implicit none

        character(len=100) :: arg,file_list,output_file
        character(len=100) :: line,mode,myname
        integer :: num_elements,num_columns
        integer :: const_start,const_stop
        
        call getarg(0,arg)
        read(arg,'(A100)') myname
      
        if (iargc() /= 4) then
                print*,"Usage is '",trim(myname)," mode num_elements file_list &
                         output_file'"
                print*,"Valid modes are 'coherent','incoherent','trajectory'."
                STOP
        endif

        call getarg(1,arg)
        read(arg,*) mode

        call getarg(2,arg)     
        read(arg,*) num_elements

!        call getarg(3,arg)
!        read(arg,*) num_columns

        call getarg(3,arg)
        read(arg,*) file_list

        call getarg(4,arg)
        read(arg,*) output_file

!        print*,"number of elements: ",num_elements
!        print*,"number of columns: ", num_columns
!        print*,"file list: ",file_list
!        print*,"output file: ",output_file

        select case (mode)
                case ("coherent")
                        num_columns = 9
                        const_start = 1
                        const_stop = 3
                        call combine_coh_inc(num_elements,num_columns,&
                                const_start,const_stop,file_list,output_file)
                case ("incoherent") 
                        num_columns = 4
                        const_start = 1
                        const_stop = 3
                        call combine_coh_inc(num_elements,num_columns,&
                                const_start,const_stop,file_list,output_file)
                case ("trajectory")
                        num_columns = 6
                        call combine_trajectory(file_list,output_file,&
                                num_elements,num_columns)
                                
                case default
                        print*, "Invalid mode selection"
                        STOP
        end select 

        END
        
        subroutine combine_trajectory(file_list,output_file,num_elements,num_columns)
                implicit none
                character(len=100), intent(in) :: file_list,output_file
                integer,intent(in) :: num_elements,num_columns
                real(8), ALLOCATABLE :: temparray(:,:)
                integer :: ierr,i,j,ostat=0
                character(len=100) :: var_fmt,current_file
                
                open(unit=11,file=file_list,status='old',action='read')
                open(unit=12,file=output_file,status='new',action='write')

                allocate(temparray(1,num_columns),STAT=ierr)
                if (ierr /= 0) then
                        print*, "Allocating temparray failed!"
                        stop
                endif 

                temparray(:,:) = 0d0

                write(var_fmt,'( "(",I1,"e20.10)")' ) num_columns

                read(11,*,iostat=ostat) current_file
                do while (ostat==0)
!                        print*, "Reading from ",trim(current_file)
                        open(unit=13,file=trim(current_file), &
                        action='read',status='old')
                        do i=1,num_elements
                                read(13,var_fmt,end=10) &
                                (temparray(1,j),j=1,num_columns)
                                do j=1,num_columns
                                        write(12,'(e20.10)',advance='no')&
                                        temparray(1,j)
                                end do 
                                write(12,*)
                        end do
10                      close(13)
                        read(11,*,iostat=ostat) current_file
                end do

        close(11)
        close(12)

        end subroutine combine_trajectory


        subroutine combine_coh_inc(num_elements,num_columns,const_start,&
                                        const_stop,file_list,output_file)
                implicit none
                integer, intent(in) :: num_elements,num_columns,const_start,const_stop
                real(8), ALLOCATABLE :: sumarray(:,:),temparray(:,:)
                integer :: ierr,ostat=0,istat=0,i,j,const_col
                character(len=100), intent(in) :: file_list,output_file
                character(len=100) :: current_file,var_fmt
 
                open(unit=11,file=file_list,status='old',action='read')
                open(unit=12,file=output_file,status='new',action='write')


                allocate(sumarray(num_elements,num_columns),STAT=ierr)
                if (ierr /= 0) then
                        print*, "Allocating sumarray failed!"
                        stop
                endif
        
                allocate(temparray(num_elements,num_columns),STAT=ierr)
                if (ierr /= 0) then
                        print*, "Allocating temparray failed!"
                        stop
                endif

                temparray(:,:) = 0d0

                write(var_fmt,'( "(",I1,"e20.10)")' ) num_columns

                read(11,*,iostat=ostat) current_file
                do while (ostat==0)
        !                print*,"Reading from ",trim(current_file)
                        open(unit=13,file=trim(current_file), &
                             action='read',status='old')
                        do i = 1,num_elements
                                read(13,var_fmt)   &
                                        (temparray(i,j),j=1,num_columns)
	        		do const_col = const_start,const_stop
                                        sumarray(i,const_col) = temparray(i,const_col)
                                end do                         

                                do j = const_stop+1,num_columns
                                sumarray(i,j) = sumarray(i,j) + temparray(i,j)  
                                end do
                        end do
                        close(13)
                        read(11,*,iostat=ostat) current_file
                end do

        
        !        print*, "Input file contents are:"
                do i = 1,num_elements
                        do j=1,num_columns
                                write(12,'(e20.10)',advance='no') sumarray(i,j)
                       end do
                       write(12,*)
                end do

                close(11)
                close(12)
        

                end subroutine combine_coh_inc


