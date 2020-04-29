!file: orthogonal_collocation.f90
module orthogonal_collocation
    implicit none
    public find_collocation_points, weight, compute_matrix_a, &
    compute_matrix_b, polynomial, d_polynomial
    contains
        subroutine find_collocation_points(collocation_points, number_of_points)
            integer, intent(in) :: number_of_points
            real(8), intent(out) :: collocation_points(number_of_points)  
            real(8) :: alpha(number_of_points), beta(number_of_points)
            integer, parameter :: M = 150
            real(8) :: y_num(M+1)
            real(8) :: y_den(M+1)
            integer :: i, j
            real(8) :: dx, x, x_new, x_old, eps
            real(8) :: summa 
            
            dx = (1.0d0 - 0.0d0)/float(M)
            alpha = 0.d0
            beta = 0.d0
            
            do j = 1, number_of_points
                if (j > 1) then
                do i = 1, M + 1
                    x = float(i-1)*dx                                               
                    y_num(i) = weight(x) * x * &
                     polynomial(j-1, x, alpha, beta,number_of_points)*&
                     polynomial(j-2, x, alpha, beta, number_of_points)

                    y_den(i) = weight(x)*&
                            polynomial(j-2, x, alpha, beta, number_of_points)**2
                end do
                beta(j) = - simpson(y_num, 0.d0, 1.d0, M) / &
                simpson(y_den, 0.d0, 1.d0, M)
            end if
                do i = 1, M + 1
                    x = float(i - 1)*dx
                    y_num(i) = (weight(x) * polynomial(j-1, x, alpha, beta, number_of_points)**2 ) * x
                    y_den(i) = weight(x) * polynomial(j-1, x, alpha, beta, number_of_points)**2
                end do
                alpha(j) = - simpson(y_num, 0.d0, 1.d0,M) /&
                simpson(y_den, 0.d0, 1.d0, M)
            end do

            ! find collocation_points
            x_old = 0.d0
            do i = 1, number_of_points
                x_old = x_old + 1.d0/float(2*number_of_points) 
                summa = 0.d0
                summa =  sum(1/(x_old - collocation_points(1:i-1)))
                do
                    x_new = x_old - & 
                    polynomial(number_of_points, x_old, alpha, beta, &
                    number_of_points)/&
                    d_polynomial(number_of_points, x_old,alpha, beta, &
                     number_of_points)/&
                    (1. - summa*polynomial(number_of_points, x_old, alpha, &
                    beta,number_of_points)/&
                    d_polynomial(number_of_points, x_old, alpha, beta,&
                     number_of_points))
                    
                    eps = abs(x_new - x_old)

                    x_old = x_new
                    if (eps < 1.d-10) exit
                end do
                collocation_points(i) = x_new
            end do
            return
        end subroutine find_collocation_points

        real(8) function weight(x)
            real(8) :: x
                x = x
                weight = 1.
        end function weight

        real(8) recursive function polynomial(j, x, alpha, beta, n) result(polynomial_j)
            integer :: i 
            integer, intent(in) :: j
            integer, intent(in) :: n
            real(8), intent(in) :: x
            real(8), intent(in) :: alpha(n)
            real(8), intent(in) :: beta(n)
            if (j == 0) then
                polynomial_j = 1.d0
                return
            end if

            if (j == -1) then
                polynomial_j = 0.d0
                return
            end if

            do i = 1, j
                polynomial_j = (x + alpha(i))* polynomial(j-1, x, alpha, beta, n) &
                + beta(i) * polynomial(j-2, x, alpha, beta, n)
            end do
            return
        end function polynomial

        real(8) recursive function d_polynomial(j, x, alpha, beta, n) result(d_polynomial_j)
            integer :: i 
            integer, intent(in) :: j
            integer, intent(in) :: n
            real(8), intent(in) :: x
            real(8), intent(in) :: alpha(n)
            real(8), intent(in) :: beta(n)
            if (j == 0 .or. j==-1 ) then
                d_polynomial_j = 0.0d0
                return
            end if
            if (j == 1) then
                d_polynomial_j = 1.0d0
                return
            end if

            do i = 1, j
                d_polynomial_j = polynomial(j-1, x, alpha, beta, n) + (x + alpha(i))* &
                d_polynomial(j-1, x, alpha, beta, n) + beta(i) &
                * d_polynomial(j-2, x, alpha, beta, n)
            end do
            return
        end function d_polynomial

        real(8) function simpson(y, lower, upper, n)
            integer, intent(in) :: n
            real(8), intent(in) :: lower, upper
            real(8), intent(in) :: y(n+1)
            real(8) :: dx 

            integer :: i
            dx = (upper-lower)/n

            simpson = 0.0d0

            simpson = simpson + y(1)
            do i = 2, n-2, 2
                simpson = simpson + 4 * y(i) + 2 * y(i+1)
            end do
            simpson = simpson + 4 * y(n) + y(n+1)

            simpson = dx/3.0d0 * simpson
        end function simpson

        subroutine compute_matrix_b(b, x, n)
            integer, intent(in) :: n
            real(8), intent(in) :: x(n)
            real(8), intent(out) :: b(n, n)

            real(8) :: summa_1, summa_2, prod
            integer :: i, j, k, l, m

            do i = 1, n
                do j = 1, n
                    summa_2 = 0.d0
                    do l = 1, n
                        summa_1 = 0.d0
                        do k = 1, n
                            prod = 1.d0
                            do m = 1, n
                                if (m/=i .and. m/=l .and. m/=k ) & 
                                prod = prod * (x(j)-x(m))/(x(i)-x(m))
                            end do
                            if (k/=i .and. k/=l) & 
                            summa_1 = summa_1 + prod/(x(i)-x(k))
                        end do
                        if (l/=i) & 
                            summa_2 = summa_2 + summa_1/(x(i)-x(l))
                    end do
                    b(j,i) = summa_2    
                end do
            end do
        end subroutine compute_matrix_b

        subroutine compute_matrix_a(a, x, n)
            integer, intent(in) :: n
            real(8), intent(in) :: x(n)
            real(8), intent(out) :: a(n, n)

            real(8) :: summa, prod
            integer :: i,j,l,k
            do i = 1, n
                do j = 1, n
                    summa = 0.0d0
                    do l = 1, n
                        prod = 1.0d0
                        do k = 1, n
                            if (k/=i .and. k/=l) & 
                            prod = prod * (x(j)-x(k))/(x(i)-x(k))
                        end do
                        if (l/=i) & 
                            summa = summa + prod/(x(i)-x(l))
                    end do
                    a(j,i) = summa    
                end do
            end do
        end subroutine compute_matrix_a
end module orthogonal_collocation
!end file: orthogonal_collocation.f90