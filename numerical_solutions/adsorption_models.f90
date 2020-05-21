module adsorption_models
    use types
    use simulation_parameters
    implicit none
    contains
    subroutine fd_eta_3pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar

        integer i ! counter
        real(dp) :: cb
        real(dp) :: cp(0:n+1)
        real(dp) :: der1, der2 ! auxiliar variables for derivatives computation
           
        cp(1:n) = y(1+cc1:n+cc1) !Internal variables definition
        cb = y(neq)
        
        select case(cc2)
        case(0)
            if (bi < 0.d0)  then 
                cp(n+1) = cb
            else    
                cp(n+1) = ( 2.d0*dx*bi*cb + 4.d0*cp(n) - cp(n-1) ) / ( 3.d0 + 2.d0*dx*bi )
            end if
        case(1)
            write(*,*) 'cc2 = 1 is not available'
            stop
        case default
            write(*,*) 'ERROR: CC2 BAD ARGUMENT'
            stop
        end select
    
        select case(cc1)
        case(0)
            write(*,*) 'cc1 = 0 is not available'
            stop
        case(1)
            cp(0) = y(1) 
            der2 = ( 2.d0*cp(1) - 2.d0*cp(0)) / (dx**2) 
            delta(1) = -yprime(1) + (s + 1.d0)*der2
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
        
        !Differencial equations inside the particle
        do i = 1, n
            der1 = ( cp(i+1) - cp(i-1) ) / ( 2.d0*dx ) 
            der2 = ( cp(i+1) - 2.d0*cp(i) + cp(i-1) ) / ( dx**2 ) 
            delta(i+1) = -yprime(i+1) + ( der2 + float(s)/(dx*float(i)) * der1)
        end do
        !Differntial equation for bulk phase.        
        der1 = ( cp(n-1) - 4.d0*cp(n) + 3.d0*cp(n+1) ) / (2.d0 * dx)
        delta(neq) = -yprime(neq) - (s+1)*alfa* der1
        return
    end subroutine fd_eta_3pts

    subroutine fd_eta_5pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar    
        real(dp) :: cb
        
        integer :: i
        real(dp) :: cp(0:n+1)
        real(dp) :: der1, der2
        !	defini��o das vari�veis internas
        cp(1:n) = y(1+cc1:n+cc1)
        cb = y(neq)
        select case(cc2)
        case(0)
            !	condi��o de contorno na superf�cie externa
            if (bi < 0.d0) then
                cp(n+1) = cb
            else
                cp(n+1) = ( 12.d0*dx*bi*cb + 48.d0*cp(n) - 36.d0*cp(n-1) + 16*cp(n-2) - 3.d0*cp(n-3) ) &
                / ( 25.d0 + 12.d0*dx*bi )
            end if
        case(1)
            write(*,*) 'cc2 = 1 is not available'
            stop
        case default
            write(*,*) 'ERROR: CC2 BAD ARGUMENT'
            stop
        end select
    
        select case(cc1)
        case(0)
            write(*,*) 'cc1 = 0 is not available'
            stop
        case(1)
            cp(0) = y(1)
            ! 5 pontos centrais
            der2= ( - 2.d0*cp(2) + 32.d0*cp(1) - 30.d0*cp(0) ) / ( 12.d0*dx**2 )!1/90*y6*h4
            delta(1) = -yprime(1) + (S + 1.d0)*der2
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
        ! Equa��o diferencial no ponto 1
        der1 = ( - cp(3) + 8.d0*cp(2)  + cp(1) - 8.d0*cp(0) ) / ( 12.d0*dx )!1/30*y5*h4
        der2 = ( - cp(3) + 16.d0*cp(2) - 31.d0*cp(1) + 16.d0*cp(0) ) / ( 12.d0*dx**2 )!1/90*y6*h4
        delta(1+cc1) = -yprime(1+cc1) + ( der2 + S / (dx*float(1)) * der1  ) 
        !
        ! Equa��es diferenciais nos pontos 2 a n-1
        do i = 2, n-1
            ! derivadas centrais com 5 pontos
            der1 = ( -cp(i+2) + 8.d0*cp(i+1) - 8.d0*cp(i-1) + cp(i-2) ) / ( 12.d0*dx )!1/30*y5*h4
            der2 = ( -cp(i+2) + 16.d0*cp(i+1) - 30.d0*cp(i) + 16.d0*cp(i-1) - cp(i-2)  ) / ( 12.d0*dx**2 )!1/90*y6*h4
            delta(i+cc1) = -yprime(i+cc1) + ( der2 + S / (dx*float(i)) * der1  ) 
        end do
        ! diferen�as negativas "quasi-centrais" na posi��o n com 6 pontos
        der1 = ( 12.d0*cp(n+1) + 65.d0*cp(n) - 120.d0*cp(n-1) + 60.d0*cp(n-2) - 20.d0*cp(n-3)  + 3.d0*cp(n-4))&
        / ( 60.d0*dx ) !1/30*y6*h5
        
        der2 = ( 10.d0*cp(n+1) - 15.d0*cp(n) - 4.d0*cp(n-1) + 14.d0*cp(n-2) - 6.d0*cp(n-3) + cp(n-4) ) &
        / ( 12.d0*dx**2 ) !13/180*y6*h4
        
        delta(n+cc1) = -yprime(n+cc1) + ( der2 + S / (dx*float(n)) * der1  )        
    
    
        der1 = (3*cp(n-3) - 16*cp(n-2) + 36*cp(n-1) - 48*cp(n) + 25*cp(n+1)) / (12.d0 * dx)
        delta(neq) = -yprime(neq) - (S+1.d0)*alfa*der1
        return
    end subroutine fd_eta_5pts

    subroutine fd_mu_3pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar
        real(dp) :: cb
        integer i ! counter
        real(dp) cp(0:n+1)
        real(dp) der1, der2 ! auxiliar variables for derivatives computation
        
        
        
        !	defini��o das vari�veis internas
        cp(1:n) = y(1+cc1:n+cc1)
        cb = y(neq)
        select case(cc2)
        case(0)
            !	condi��o de contorno na superf�cie externa
            if (bi < 0.d0) then 
                cp(n+1) = cb
            else 
                cp(n+1) = (dx*bi*cb + 4.d0*cp(n) - cp(n-1)) / (3.d0 + dx*bi) !1/30*y3*h2
            end if
        case(1)
            write(*,*) 'cc2 = 1 is not available'
            stop
        case default
            write(*,*) 'ERROR: CC2 BAD ARGUMENT'
            stop
        end select
    
        select case(cc1)
        case(0)
            write(*,*) 'cc1 = 0 is not available'
            stop
        case(1)
            ! define a concentra��o no ponto zero
            cp(0) = y(1)
            ! 3 pontos
            der1 = (- cp(2) + 4.d0*cp(1) - 3.d0*cp(0) ) / ( 2.d0*dx )!1/30*y3*h2
            delta(1) = -yprime(1) +  2.d0*(s + 1.d0)*der1
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
    
        ! Equa��es diferenciais nos pontos 1 a n
        do i = 1, n
            ! derivadas centrais com 3 pontos	
            der1 = ( cp(i+1) - cp(i-1) ) / ( 2.d0*dx )!1/6*y3*h2
            der2 = ( cp(i+1) - 2.d0*cp(i) + cp(i-1) ) / ( dx**2 )!1/12*y4*h2
            delta(i+cc1) = -yprime(i+cc1) + ( 4.d0*dx*float(i)*der2 + 2.d0*(float(s) + 1.d0)*der1  )
        end do
    
        der1 = (cp(n-1) - 4.d0*cp(n) + 3.d0*cp(n+1)) / (2.d0 * dx)
        delta(neq) = -yprime(neq) - (s+1)*alfa*2.d0*sqrt(1.d0)* der1
    RETURn
    end subroutine fd_mu_3pts

    subroutine fd_mu_5pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar

        integer i ! counter
        real(dp) :: cb
        real(dp) cp(0:n+1)
        real(dp) der1, der2 ! auxiliar variables for derivatives computation
        
        !	defini��o das vari�veis internas
        cp(1:n) = y(1+cc1:n+cc1)
        cb = y(neq)
        ! Condi��o de contorno no centro
        select case(cc2)
        case(0)
            !	condi��o de contorno na superf�cie externa
            if (bi < 0.d0) then
                cp(n+1) = cb
            else
                cp(n+1) = ( 6.d0*dx*bi*cb + 48.d0*cp(n) - 36.d0*cp(n-1) + 16*cp(n-2) - 3.d0*cp(n-3) )&
                /( 25.d0 + 6.d0*dx*bi )!1/5*y5*h4
            end if
            case(1)
                write(*,*) 'cc2 = 1 is not available'
                stop
            case default
                write(*,*) 'ERROR: CC2 BAD ARGUMENT'
                stop
            end select
    
        select case(cc1)
        case(0)
            write(*,*) 'cc1 = 0 is not available'
            stop
        case(1)
            ! define a concentra��o no ponto zero
            cp(0) = y(1)
            der1= ( -3.d0*cp(4) + 16.d0*cp(3) - 36.d0*cp(2) + 48.d0*cp(1) - 25.d0*cp(0) ) / ( 12.d0*dx )!1/5*y5*h4
            delta(1) = -yprime(1) +  2.d0*(s + 1.d0)*der1 
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
    
        
        ! Equa��o diferencial no ponto 1
        !
        ! derivadas "quasi-centrais" com 5 pontos
        der1 = ( cp(4) - 6.d0*cp(3) + 18.d0*cp(2) - 10.d0*cp(1) - 3.d0*cp(0) ) / ( 12.d0*dx )!1/20*y5*h4
        der2 = ( cp(5) - 6.d0*cp(4) + 14.d0*cp(3) - 4.d0*cp(2) - 15.d0*cp(1) + 10.d0*cp(0) ) / ( 12.d0*dx**2 )!13/180*y6*h4
        delta(1+cc1) = -yprime(1+cc1) + ( 4.d0*dx*float(1)*der2 + 2.d0*(s + 1.d0)*der1  ) 
        ! Equa��es diferenciais nos pontos 2 a n-1
        do i = 2, n-1
            ! derivadas centrais com 5 pontos
            der1 = ( -cp(i+2) + 8.d0*cp(i+1) - 8.d0*cp(i-1) + cp(i-2) ) / ( 12.d0*dx )!1/30*y5*h4
            der2 = ( -cp(i+2) + 16.d0*cp(i+1) - 30.d0*cp(i) + 16.d0*cp(i-1) - cp(i-2)  ) / ( 12.d0*dx**2 )!1/90*y6*h4
            delta(i+cc1) = -yprime(i+cc1) + ( 4.d0*dx*float(i)*der2 + 2.d0*(s + 1.d0)*der1  ) 
        end do
        ! Equa��o diferencial no ponto n
        ! diferen�as negativas "quasi-centrais" na posi��o n com 5 pontos
        der1 = ( 3.d0*cp(n+1) + 10.d0*cp(n) - 18.d0*cp(n-1) + 6.d0*cp(n-2) - cp(n-3) ) / ( 12.d0*dx )!1/20*y5*h4
        der2 = ( 10.d0*cp(n+1) - 15.d0*cp(n) - 4.d0*cp(n-1) + 14.d0*cp(n-2) - 6.d0*cp(n-3) + cp(n-4) ) / ( 12.d0*dx**2 )!13/180*y6*h4
        delta(n+cc1) = -yprime(n+cc1) + ( 4.d0*dx*float(n)*der2 + 2.d0*(s + 1.d0)*der1  ) 
        !
        der1 = (3*cp(n-3) - 16*cp(n-2) + 36*cp(n-1) - 48*cp(n) + 25*cp(n+1)) / (12.d0 * dx)
        delta(neq) = -yprime(neq) - (s+1)*alfa*2.d0*sqrt(1.d0)* der1
        return
    end subroutine fd_mu_5pts

    subroutine oc_eta(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar
        
        real(dp) :: cb
        real(dp) :: cp(1:n+2)
        real(dp) :: der1, der2 
        integer :: i
        
        
        cp(2:n+1) = y(1+cc1:n+cc1)
        cb = y(neq)
        ! define a concentra��o no ponto zero
        select case(cc1)
        case(0)
            select case(cc2)
            case(0)
                if (bi < 0.d0) then 
                    cp(n+1) = cb
                else
                    cp(n+2) = (bi*cb + sum((A(n+2,1) * A(1,2:n+1)/A(1,1) - A(n+2,2:n+1)) * cp(2:n+1))) /&
                    (A(n+2,n+2) - A(n+2,1)*A(1,n+2)/A(1,1) + bi)
                end if
            case(1)
                write(*,*) 'cc2 = 1  is not available'
                stop
            case default
                write(*,*) 'ERROR: CC2 BAD ARGUMENT'
                stop
            end select
            
            cp(1) = - sum(A(1,2:n+2)*cp(2:n+2)) / A(1,1)
        case(1)
            cp(1) = y(1)
            ! condi��o alg�brica na superf�cie externa � resolvida a priori
            select case(cc2)
            case(0)
                if (bi < 0.d0) then 
                    cp(n+2) = cb
                else
                    cp(n+2) = (bi*cb - sum(A(n+2,1:n+1)*cp(1:n+1))) / ( A(n+2,n+2) + bi)
                end if
            case(1)
                write(*,*) 'cc2 = 1 is not available'
                stop
            case default
                write(*,*) 'ERROR: CC2 BAD ARGUMENT'
                stop
            end select
            
            ! derivada segunda
            der2 = sum(B(1,1:n+2)*cp(1:n+2))
            
            ! Equa��o diferencial no centro
            delta(1) = -yprime(1) +  float(S+1)*der2
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
        
        
        !   Equa��es diferenciais nos pontos 1 a n
        do i = 2, n+1
            der1 = sum(a(i,1:n+2)*cp(1:n+2)) 
            der2 = sum(b(i,1:n+2)*cp(1:n+2)) 
            delta(i+cc1) = -yprime(i+cc1) + (der2 + float(S)/(x(i)) * der1)
        end do 
        
        der1 = sum(A(n+2,1:n+2)*cp(1:n+2))
        delta(neq) = -yprime(neq) - (s+1.0_dp)*alfa*der1
       
    end subroutine

    subroutine oc_mu(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in), dimension(neq) :: y
        real(dp), intent(in), dimension(neq) :: yprime
        real(dp), intent(in) :: cj
        real(dp), intent(out), dimension(neq) :: delta
        integer, intent(in) :: ires
        real(dp), intent(in), dimension(:) :: rpar
        integer, intent(in), dimension(:) :: ipar
        
        real(dp) :: cb
        real(dp) :: cp(1:n+2)
        real(dp) :: der1, der2 
        integer :: i

        cp(2:n+1) = y(1+cc1:n+cc1)

        
        select case(cc1)
        case(0)
            write(*,*) 'cc1 = 0 is not available'
            stop
        case(1)
            cp(1) = y(1)
            ! condi��o alg�brica na superf�cie externa � resolvida a priori
            select case(cc2)
            case(0)
                if (bi < 0.d0) then 
                    cp(n+2) = cb
                else
                    cp(n+2) = (bi*cb - sum(A(n+2,1:n+1)*cp(1:n+1))) / ( A(n+2,n+2) + bi)
                end if
            case(1)
                write(*,*) 'cc2 = 1 is not available'
                stop
            case default
                write(*,*) 'ERROR: CC2 BAD ARGUMENT'
                stop
            end select
            ! derivada segunda
            ! Equa��o diferencial no centro
            der1 = sum(A(1,1:n+2)*cp(1:n+2))
            delta(1) = -yprime(1) +  2.d0*float(s+1)*der1
        
        case default
            write(*,*) 'ERROR: CC1 BAD ARGUMENT'
            stop
        end select
    
    
        ! condi��o de contorno no centro
        ! derivada primeira
        ! equa��es diferenciais nos pontos 1 a n
        do i = 2, n+1
            der1 = sum(A(i,1:n+2)*cp(1:n+2)) 
            der2 = sum(B(i,1:n+2)*cp(1:n+2)) 
            delta(i+cc1) = -yprime(i+cc1) + ( 4.d0*x(i)*der2 + 2.d0*float(s+1) * der1  )
        end do
    
        der1 = sum(A(n+2,1:n+2)*cp(1:n+2)) 
        delta(neq) = -yprime(neq) - (s+1.d0)*alfa* der1
    
        return
    end subroutine oc_mu
end module adsorption_models