program main

  use functions
  implicit none

  integer, parameter:: N=300, time_steps=400, Density_range=5, beta_range=5, Maximum_Hop=300, s_range=4
  character(*), parameter:: location="C:\Users\Zuriel\Desktop\Graphs\", Row_Format="(500(2X, I3))" 
  double precision, parameter:: initial_density=0.1, density_tic=0.12, beta_tic=0.08, initial_beta=0,initial_s=1.1, s_tic=0.4, fix_density=0.4

  
double precision:: s, entropy, pn, entropy_summand, mean_entropy, min_entropy,max_entropy, r(2), beta, density,&
                   scaled_density, scaled_beta


  integer, dimension(:,:), allocatable:: Matrix, Table
  double precision, dimension(:,:), allocatable:: Entropy_Matrix, Entropy_Line_Matrix
  integer:: counter, entropy_counter, cells, step, i, density_counter, beta_counter
  
 character(len=200)::max_density_char, min_density_char, max_beta_char, min_beta_char, max_entropy_char,density_char, beta_char, min_entropy_char, s_char


!--------------------------------Executable part-----------------------------------------------------------------------------
allocate(Matrix(N,2))
allocate(Table(N,time_steps+1))
allocate(Entropy_Matrix(density_range,beta_range))
allocate(Entropy_Line_Matrix(beta_range,2))



s=1.001   
do density_counter=1, density_range
   density=initial_density+(density_counter-1)*density_tic
   
   do beta_counter=1, beta_range
       beta=initial_beta+Beta_Tic*(beta_counter-1)
       
       !Just to follow the run
       print*, density_counter,beta_counter
       
       !-----------------------------Start the mean entropy computations--------------------------------------------------------
          Matrix=0; table=0
           
           !initial random fill
          do counter=1, N
               call random_number(r)
               if(r(1) < density) matrix(counter, 1)=1
               if(r(2) < density) matrix(counter, 2)=1
          enddo
            
         !Saving initial condition
            Table(:,1)=(/(matrix(i,1)+matrix(i,2), i=1, N)/)
         
         !How many cells are there?
            cells=sum(Table(:,1))   
            
        !Applying migration and reorientation   
         do step=1, time_steps
             call reorientation(matrix, N, dble(beta)) 
             call migration(matrix, N, Maximum_Hop, s)
             Table(:,step+1)=(/(matrix(i,1)+matrix(i,2), i=1, N)/)
         enddo
                      
 
                !--------------------------------------------------Set graph data------------------------------------------------------------------------------------------------------------------
                  write(density_char, "(F10.4)") density; write(beta_char, "(F10.4)") beta;
                  write(max_entropy_char, "(F10.4)") max_entropy; write(min_entropy_char, "(F10.4)") min_entropy
                  write(min_beta_char, "(F10.4)") initial_beta; write(min_density_char, "(F10.4)") initial_density
                  write(max_beta_char, "(F10.4)") (initial_beta+(beta_range-1)*beta_tic)
                  write(max_density_char, "(F10.4)") (initial_density+(density_range-1)*density_tic)
                  beta_char=trim(adjustl(beta_char)); density_char=trim(adjustl(density_char))
               	  write(s_char, "(F10.4)") s; s_char=trim(adjustl(s_char))
			      
        !--------------------------------------------Saving the matrices in a text file------------------------------------------------
		     

				print*, Entropy_Matrix
                open(density_counter+beta_counter, file=location//"Accumulations"//trim(density_char)//trim(beta_char)//".dat", action="write", status="replace", RECL=7000)
                    do counter=1, time_steps+1
					  
	                      write(density_counter+beta_counter, FMT=Row_Format) (/(Table(i,counter), i=1, N)/)

					  
					enddo
                close(1)




   
                 !---------------------------------------------------Graph------------------------------------------------------------------------------------------------------------------------
                         open(3+density_counter+beta_counter, file=location//"Accumulations"//"rho="//trim(density_char)//"beta="//trim(beta_char)//".plt", status="replace", action="write")
                    
                    write(3+density_counter+beta_counter,*) "set title 'LCT"//","//" "//'Density='//trim(density_char)//","//" "//'Beta='//trim(beta_char)//","//" "//'s='//trim(s_char)//"'"       
                    write(3+density_counter+beta_counter,*) "set xlabel 'Node'"
                    write(3+density_counter+beta_counter,*) "set ylabel 'Time'"
                    write(3+density_counter+beta_counter,*) "set autoscale xfix"
                    write(3+density_counter+beta_counter,*) "set autoscale yfix"
                 	!write(3,*)"plot 'Hyperbola.dat' wl" 
                    write(3+density_counter+beta_counter,*) "set xtics auto"
                    write(3+density_counter+beta_counter,*) "set ytics auto"
					! write(3+density_counter+beta_counter,*)"set format x "" "
					  !write(3+density_counter+beta_counter,*)"set format y "" "

                    !write(3,*) "set xrange ["//trim(min_density_char)//":"//trim(max_density_char)//"]"
                    !write(3,*) "set yrange ["//trim(min_beta_char)//":"//trim(max_beta_char)//"]"
                    
                  
                   ! write(3,*)&
                    !"set palette defined"//"("// &
                    !"0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)"
                    !write(3+density_counter+beta_counter,*)"set cbrange ["// trim(min_entropy_char)//":"// trim(max_entropy_char) //"]"
                    write(3+density_counter+beta_counter,*)"set pm3d map"
					write(3+density_counter+beta_counter,*) "set cblabel 'Number of cells'"
                    write(3+density_counter+beta_counter,*)"splot"//"'Accumulations"//trim(density_char)//trim(beta_char)//".dat'"//"matrix with image" 
			


                  
              close(3+density_counter+beta_counter)
            

    
  
   enddo
   
   
enddo

       

 
end program main