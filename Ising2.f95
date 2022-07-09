program ising
	implicit none

	integer, parameter :: n = 60 	! Lado de la red, puntos totales = n**2
	integer, parameter :: m = 3 	! Lado de la red, puntos totales = n**2

	integer :: pMC 	! Unidad de tiempo 'paso Monte Carlo' ~ "iteraciones"
	integer :: i, j, k, aux_i, aux_j, l , tt
	integer :: s(0:n+1, 0:n+1), s_init(0:n+1, 0:n+1), delta_E, suma 

	real :: T, p, Temp(1:m), Taux
	real :: r1(1:n, 1:n), r2(1:n, 1:n), eps
 
	character(len = 2) :: path



	call random_seed()

	call init_red(s_init, n)
	call cond_contorno(s_init, n)

	!write(*,*) "Temperatura: "
	!read(*,*) T

	Temp = (/2.2, 2.3, 2.4/)
	! Temperaturas: 0.5, 1.0, 2.0, 2.3, 4.0

	write(*,*) "Pasos MonteCarlo: "
	read(*,*) pMC 

	
	do tt = 1, m
		T = Temp(tt)
		Taux = T*10
		write(path, '(I2.2)') int(Taux)
		!open(10, file="DataRelI2/Data"//trim(path)//".dat", status = "unknown")
		open(11, file="DataExtra/DataRelative"//trim(path)//".dat", status = "unknown")
		s = s_init
		do k = 1, pMC 
			suma = 0 
			call random_number(r1)
			r1 = floor(r1*(n+1))
			call random_number(r2)
			r2 = floor(r2*(n+1))

			do i = 1, n
				do j = 1, n
					aux_i = int(r1(i,j))
					aux_j = int(r2(i,j))
					delta_E = 2*s(aux_i, aux_j)*(s(aux_i+1,aux_j)+s(aux_i-1,aux_j)+s(aux_i,aux_j+1)+s(aux_i,aux_j-1))
					p = min(1.0, exp(-delta_E/T))

					call random_number(eps)

					if (eps < p) then
						s(aux_i, aux_j) = -s(aux_i, aux_j) 
					end if
					if (aux_i==0 .or. aux_i==n .or. aux_j==0 .or. aux_j==n) then
						call cond_contorno(s, n)
					end if
				end do
			end do
			do i = 1,n 
				do j = 1,n
					suma = suma + s(i,j)
				end do 
			end do 

			write(11,*) k, 1.0*suma/n**2

		!	do i = 1, n
		!		do j = 1, n
		!			if ( s(i,j) ==1 ) then
		!				write(10,*) i, j 
		!			end if
		!		end do
		!	end do
		!	write(10,*)
		!	write(10,*)
		end do 
		!close(10)
		close(11)
	end do	
end program ising

subroutine init_red(red, n)
	implicit none
	integer, intent(in) :: n
	integer, intent(inout) ::  red(0:n+1, 0:n+1)

	integer :: i, j
	real :: r 

	do i = 1, n
		do j = 1, n
			call random_number(r)
			if ( r < 0.5 ) then
				red(i,j) = 1
			else 
				red(i,j) = -1
			end if
		end do
	end do

	red(0,0) = 0
	red(n+1,0) = 0
	red(0,n+1) = 0
	red(n+1,n+1) = 0

end subroutine init_red

subroutine cond_contorno(red, n)
	implicit none
	integer,intent(in) :: n
	integer,intent(inout) ::  red(0:n+1, 0:n+1)

	integer :: i, j 

	do i = 1, n
		red(i,0) = red(i,n)
		red(i,n+1) = red(i,1)
		red(0,i) = red(n,i)
		red(n+1,i) = red(1,i)
	end do

end subroutine cond_contorno