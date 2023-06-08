module functions

 implicit none 

contains

 


  subroutine reorientation (X, N, beta)
  	  

    integer, intent(in):: N
    integer, dimension(N,2), intent(inout):: X
    double precision, intent(in):: beta


    integer, allocatable:: copy_X(:,:)
    integer:: i, j1, previous_node, next_node, cells, cells_counter
    double precision:: P1, P2, P_max, Z 
    double precision, allocatable:: R(:) 


    double precision:: gradient

	 
	
	cells=sum(X(:,1))+sum(X(:,2))
    allocate(R(cells))
	allocate(copy_X(N,2))
    

    Copy_X=X
	X=0	
    cells_counter=0
    call random_number(R)
	

	
	do i=1, N
	   
	
	   if(i==1) then
	      previous_node=N
	   else
	      previous_node=i-1 
	   endif

	   if(i==N) then
	      next_node=1
	   else
	      next_node=i+1 
	   endif



	   gradient=-(copy_X(previous_node,1)+copy_X(previous_node,2))+(copy_X(next_node,1)+copy_X(next_node,2))
	   P1=-beta*gradient
	   P2=beta*gradient


	   P_max=max(P1,P2)
	   P1= P1-P_max
	   P2=P2-P_max
       

	   P1=exp(P1)
	   P2=exp(P2)
       

       Z=P2+P1
	   P1=P1/Z
	   P2=P2/Z
	   



	   if(copy_X(i,1)+copy_X(i,2)>0) then

	     do j1=1, copy_X(i,1)+copy_X(i,2)
           cells_counter=cells_counter+1
		   if(R(cells_counter)< P1) then

		     X(i,1)=X(i,1)+1
		     
		   else
		      X(i,2)=X(i,2)+1

		   endif
		 enddo

	   endif
	enddo


   	deallocate(copy_X)
    deallocate(R)

   end subroutine reorientation 

  


   subroutine migration(X, N, max_jumps, s)

   
	  integer:: i, j
	  double precision:: r
	  integer::hops, new_position 
	  integer, allocatable:: copy_X(:,:)



      integer, intent(in):: N 
      double precision, intent(in):: s
      integer, dimension(N,2), intent(inout):: X
	  integer, optional:: max_jumps
	  if(.not. present(max_jumps)) max_jumps=N


	allocate(copy_X(N,2))
	Copy_X=X
	X=0	
	 
	  do i=1, N

		
	    if(copy_X(i,2)>0) then

		    
			do j=1, copy_X(i,2)
			  
			  if(max_jumps==1) then
                hops=1
              else
                hops=zipf(dble(max_jumps), s)
              endif
			 
			    new_position=mod(i-1+hops, N)+1

			 

			  X(new_position,2)=X(new_position,2)+1
			enddo 
		endif


		  if(copy_X(i,1)>0) then

		    
			do j=1, copy_X(i,1)
			  if (max_jumps==1) then
                hops=1
              else
			    hops=zipf(dble(max_jumps), s)
              endif
			  

				new_position=i-hops

				do while(new_position < 1)
					 new_position= N+new_position
				enddo

			    

			  X(new_position,1)=X(new_position,1)+1


			
			enddo 
		 endif
	     
	  enddo

	  deallocate(copy_X)

	


   endsubroutine migration

   
    !Subrutina de migración por bloques
   subroutine block_migration(X, N, max_jumps, s)

   	!variables locales
	  integer:: i
	  double precision:: r
	  integer::hops, new_position 
	  integer, allocatable:: copy_X(:,:)


	!Parámetros de función
      integer, intent(in):: N
      double precision, intent(in):: s
      integer, dimension(N,2), intent(inout):: X
	  integer, optional:: max_jumps
	  if(.not. present(max_jumps)) max_jumps=N

	  
	!Copiamos el vector original para hacer las iteraciones y lo "vaciamos" para empezar a reacomodar cada célula
	
	allocate(copy_X(N,2))
	Copy_X=X
	X=0	
	 	!inicio de ciclo de reacomodo
	  do i=1, N

		!Reacomodo en el canal de velocidad de la derecha
	    if(copy_X(i,2)>0) then

		    
		     if(max_jumps==1) then
               hops=1
              else

			  hops=zipf(dble(max_jumps), s)

			 endif
			    new_position=mod(i-1+hops, N)+1

			 

			  X(new_position,2)=X(new_position,2)+copy_X(i,2)
			
		endif

		  	!Reacomodo en el canal de velocidad de la izquierda
		  if(copy_X(i,1)>0) then

		    
		     if(max_jumps==1) then
              hops=1
              else
	
			  hops=zipf(dble(max_jumps), s)
              endif
			  new_position=i-hops

				do while(new_position < 1)
					 new_position= N+new_position
				enddo

			  X(new_position,1)=X(new_position,1)+copy_X(i,1)

		 endif
	     
	  enddo

	  deallocate(copy_X)

	 


   endsubroutine block_migration


   !_______________________________________________________________________________________________________________
   !Distribución de Zipf; toma el coeficiente 's' y el salto máximo N y arroja un número natural
   !La función utiliza el método iterativo de Newton Raphson

   !Seguimos la sugerencia de la página. Definimos funciones para las potencias
   !______________________________________________________________________________________________________________

   !Cálculo de x^(-s-2)
   double precision function m(x,s)
      implicit none 

     double precision , intent(in):: s,x
     m=x**(-s-2)
   end function m

   !Cálculo de x^(-s-1)
    double precision function mx(x,s)
      implicit none 

     double precision , intent(in):: s,x
     mx=x*m(x,s)
   end function mx

   !Cálculo de x^(-s)
   double precision function mxx(x,s)
      implicit none 

     double precision , intent(in):: s,x
     mxx=x*mx(x,s)
   end function mxx

   !Cálculo de x^(-s+1)
   double precision function mxxx(x,s)
      implicit none 

     double precision , intent(in):: s,x
     mxxx=x*mxx(x,s)
   end function mxxx

   
   !Cálculo del salto con la distribución de Zipf
   integer function Zipf(N,s)
       implicit none 

	  double precision, parameter:: tol=0.01

      double precision, intent(in):: N, s
	  double precision:: a, b
   	  double precision:: Dp
	  real:: p
	  double precision:: xn, xnn
	  
      call random_number(p)
	  Dp=p*(12*(N**(-s+1)-1)/(1-s)+6+6*(N**(-s))+s-s*(N**(-s-1))) 
	  xn=dble(N/2.0)
	  
	  do 
	    a=12*(mxxx(xn,s)-1)/(1-s)+6+6*mxx(xn,s)+s-s*mx(xn,s)-Dp
		b=12*mxx(xn,s)-(6*s*mx(xn,s))+(s*(s+1)*m(xn,s))
		xnn=max(xn-a/b, dble(1))
		xn=xnn
	    if (abs(xnn-xn)< tol) EXIT
	  enddo
      
      Zipf=int(xnn)
	  
	  
   end function Zipf   

end module functions