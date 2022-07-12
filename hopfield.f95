program ising
	implicit none

	integer, parameter :: n = 60 	! Lado de la red, puntos totales = n**2
	integer, parameter :: m = 3 	! Lado de la red, puntos totales = n**2

	integer :: pMC 	! Unidad de tiempo 'paso Monte Carlo' ~ "iteraciones"
	integer :: i, j, k, aux_i, aux_j, l , tt, paso 
	integer :: s(0:n+1, 0:n+1), s_init(0:n+1, 0:n+1), suma 
	integer :: s_copy(0:n+1, 0:n+1)

	real :: T, p, Temp(1:m), Taux
	real :: r1(1:n, 1:n), r2(1:n, 1:n), eps

	double precision :: w(1:n,1:n,1:n,1:n), theta(1:n, 1:n), H, H_2, delta_H
 
	character(len = 2) :: path

	call patron(s_init, n)
	call cond_contorno(s_init, n)
	call peso_sinap(s_init, n, w)
	call umbral_disparo(w, n, theta)
	call random_red(s, n)

	Temp = (/0.9E-4, 1E-4, 1.1E-4/)

	pMC = 500
	
	! do tt = 1, m
		T = Temp(2)
		Taux = T*10E4
		write(path, '(I2.2)') int(Taux)
		open(10, file="Data/Data1.dat", status = "unknown")

		do paso = 0, pMC
			suma = 0 
			call random_seed()
			call random_number(r1)
			r1 = floor(r1*(n+1))
			call random_number(r2)
			r2 = floor(r2*(n+1))

			do i = 1, n
				do j = 1, n
					aux_i = int(r1(i,j))
					aux_j = int(r2(i,j))
					
					call hamiltoniano(s, w, theta, n, H)

					s_copy = s 
	
					call intercambia(aux_i, aux_j, n, s_copy)
					call hamiltoniano(s_copy, w, theta, n, H_2)
					
					delta_H = H_2 - H 

					! delta_H = 0 
					! do k = 1, n
					! 	do l = 1, n
					! 		delta_H = delta_H + 2.0*w(aux_i,aux_j,k,l)*s(aux_i,aux_j)*s(k,l) + theta(aux_i,aux_j)*s(aux_i,aux_j)
					! 	end do
					! end do
					p = min(1.0, exp(-delta_H/T))

					call random_number(eps)

					if (eps < p) then
						call intercambia(aux_i, aux_j, n, s)
					end if
					if (aux_i==1 .or. aux_i==n .or. aux_j==1 .or. aux_j==n) then
						call cond_contorno(s, n)
					end if
				end do
			end do

			if (mod(paso,10)==0) then
				do i = 1, n
					do j = 1, n
						if (s(i,j) == 1) then
							write(10,*) i, j 
						end if
					end do
				end do
				write(10,*)
				write(10,*)
			end if 
		end do 
		close(10)
	! end do	
end program ising

subroutine random_red(red, n)
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

end subroutine random_red

! Aplico condiciones de contorno de tal manera que la red es equivalente a un toroide.
! (Salir por la derecha = entrar por la izquierda, salir por abajo = entrar por arriba)
subroutine cond_contorno(red, n)
	implicit none
	integer,intent(in) :: n
	integer,intent(inout) ::  red(0:n+1, 0:n+1)

	integer :: i

	do i = 1, n
		red(i,0) = red(i,n)
		red(i,n+1) = red(i,1)
		red(0,i) = red(n,i)
		red(n+1,i) = red(1,i)
	end do

end subroutine cond_contorno


subroutine intercambia(i, j, n, red)
	integer, intent(in) :: i, j, n
	integer, intent(inout) :: red(0:n+1, 0:n+1)

	if (red(i,j) == 1) then 
		red(i,j) = 0
	else 
		red(i,j) = 1
	end if 
	
end subroutine intercambia


! Esto es Eta^mu
! El patron depende del nombre del path
subroutine patron(red, n)
	implicit none 
	integer, intent(in) :: n
	integer, intent(inout) :: red(1:n,1:n)

	integer, parameter :: m = 40
	integer :: patr(1:m,1:m)
	integer :: i, j

	open(1, file='patron1', status='unknown')
		read(1,*) patr
	close(1)

	do i = 1, m
		do j = 1, m
			red(i,j) = patr(i,j)
		end do
	end do

	
end subroutine patron

! Calcula el promedio de los valores de mi red
subroutine promedio(red, n, media)
	implicit none 
	integer, intent(in) :: n 
	integer, intent(in) :: red(0:n+1, 0:n+1)
	double precision, intent(out) :: media 

	integer :: i, j 
	double precision :: suma 

	suma = 0 
	do i = 1, n
		do j = 1, n
			suma = suma + red(i,j)
		end do
	end do

	media = suma/n**2 

end subroutine promedio

subroutine peso_sinap(red, n, w)
	implicit none 
	integer, intent(in) :: n, red(0:n+1, 0:n+1)
	double precision, intent(out) :: w(1:n,1:n,1:n,1:n)

	integer :: i,j,k,l 
	double precision :: media 

	call promedio(red, n, media)

	do i = 1, n
		do j = 1, n
			do k = 1, n
				do l = 1, n
					if (i == k .and. j == l) then 
						w(i,j,k,l) = 0
					else 
						w(i,j,k,l) = (red(i,j) - media)*(red(k,l) - media)
					end if 
				end do
			end do
		end do
	end do


end subroutine peso_sinap

subroutine umbral_disparo(w, n, theta)
	integer, intent(in) :: n 
	double precision, intent(in) :: w(1:n,1:n,1:n,1:n)
	double precision, intent(out) :: theta(1:n,1:n)

	integer :: i,j,k,l 

	theta = 0 
	do i = 1, n
		do j = 1, n
			do k = 1, n
				do l = 1, n
					theta(i,j) = theta(i,k) + 0.5 * w(i,j,k,l)
				end do
			end do
		end do
	end do

end subroutine umbral_disparo

subroutine hamiltoniano(red, w, theta, n, H)
	integer, intent(in) :: red(0:n,0:n), n
	double precision, intent(in) :: w(1:n,1:n,1:n,1:n), theta(1:n,1:n)
	double precision, intent(out) :: H

	integer :: i, j, k, l 

	H = 0 
	do i = 1, n
		do j = 1, n
			do k = 1, n
				do l = 1, n
					H = H - 0.5*w(i,j,k,l)*red(i,j)*red(k,l) + theta(i,j)*red(i,j)
				end do
			end do
		end do
	end do
	
end subroutine hamiltoniano