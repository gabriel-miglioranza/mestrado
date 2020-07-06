! Solutions for the Diffusion differential equation
module batch
implicit none
public batch_particle, batch_bulk, root_calc, f, df
integer, parameter :: dp = kind(1.d0)

contains
subroutine batch_particle_average(t, u, bi, a, s, eps)
    real(dp), intent(in) :: t(:)
    real(dp), intent(out) :: u(size(t))
    real(dp), intent(in) :: bi, a
    integer, intent(in) :: s
    real(dp), intent(in) :: eps

    real(dp), parameter :: PI = 4.0d0*atan(1.d0)
    real(dp) :: g, g_old, bn, nu
    real(dp) :: term(size(t))
    integer :: n

        nu = real(s-1, dp)/2.0_dp
        u = 1.0_dp/(1.0_dp + a)
        g = 0
        do n = 1, 100
            g_old = g
            if (bi>0 .or. a>0) then
                do while(g_old >= g)
                    g_old = g_old + 0.1_dp
                    call root_calc(g, g_old, bi, a, s, eps)
                end do
                if (a == 0) then
                    bn = -2*bi/(g**2 + bi*(bi - 2*nu))
                else if (bi == -1) then
                    bn = 4*(nu+1)*a/(g**2 + 4*a*(1+a)*(1+nu)**2)
                else if (a > 0 .and. bi > 0) then
                    bn = 2*bi*(2*(nu+1)*a*bi - g**2)/(g**4 + g**2*bi*& 
                (bi-4*a*(nu+1) - 2*nu) + 4*a*(1 + a)*bi**2*(nu + 1)**2)
                end if
            else if (bi == -1 .and. a == 0) then
                select case(s)
                case(0)
                    g = (2*n-1)*PI/2
                case(1)
                    do while (g_old>=g)
                        g_old = g_old + 0.1_dp
                    call root_calc(g, g_old, bi, a, s, eps)
                    end do
                case(2)
                    g = n*PI
                case default
                return
                end select
                    bn = real(-2*((-1)**(s+1))**(n+1), dp)/g
            else
                return
            end if

            select case(s)
            case(0)
                term = bn * tan(g)/g * exp(-g**2 * t)        
            case(1) 
                term =  bn * 1.0_dp/g * bessel_j1(g)/bessel_j0(g) * exp(-g**2 * t)
            case(2)
                term = bn * 1.0_dp/g * (1.0_dp/g - 1.0_dp/tan(g)) * exp(-g**2 * t)
            case default
                return
            end select

            u = u + term
            if (sqrt(maxval(abs(term))) < eps) exit
        end do

end subroutine batch_particle_average

subroutine batch_particle(x, t, u, bi, a, s, eps)
    real(dp), intent(in) :: x(:) 
    real(dp), intent(in) :: t (:)
    real(dp), intent(out) :: u(size(x), size(t))
    real(dp), intent(in) :: bi, a 
    integer, intent(in) :: s 
    real(dp), intent(in) :: eps 

    real(dp), parameter :: PI = 4.0d0*atan(1.d0)
    real(dp) :: g, g_old, bn, nu
    real(dp) :: term(size(x))
    integer :: n, i 

        nu = real(s-1, dp)/2.0_dp
        u = 1.0_dp/(1.0_dp + a)
        g = 0

        do n = 1, 100
            g_old = g
            if (bi>0 .or. a>0) then
                do while(g_old >= g)
                    g_old = g_old + 0.1_dp
                    call root_calc(g, g_old, bi, a, s, eps)
                end do

                if (a == 0) then
                    bn = -2*bi/(g**2 + bi*(bi - 2*nu))
                else if (bi == -1) then
                    bn = 4*(nu+1)*a/(g**2 + 4*a*(1+a)*(1+nu)**2)
                else if (a > 0 .and. bi > 0) then
                    bn = 2*bi*(2*(nu+1)*a*bi - g**2)/(g**4 + g**2*bi*& 
                    (bi-4*a*(nu+1) - 2*nu) + 4*a*(1 + a)*bi**2*(nu + 1)**2)
                end if

            else if (bi == -1 .and. a == 0) then
                    select case(s)
                case(0)
                    g = (2*n-1)*PI/2
                case(1)
                    do while (g_old>=g)
                        g_old = g_old + 0.1_dp
                    call root_calc(g, g_old, bi, a, s, eps)
                    end do
                case(2)
                    g = n*PI
                case default
                return
                end select

                bn = real(-2*((-1)**(s+1))**(n+1), dp)/g
            else
                return
            end if

            do i = 1, size(t)
                if (t(i)==0) then
                    u(:,i) = 0
                    term = 0
                else
                    select case(s)
                    case(0)
                        term = bn * cos(g*x)/cos(g) * exp(-g**2 * t(i))        
                    case(1) 
                        term =  bn * bessel_j0(g*x)/bessel_j0(g) * exp(-g**2 * t(i))
                    case(2)
                        where(x/=0)
                            term = bn * sin(g*x)/(x*sin(g)) * exp(-g**2 * t(i))
                        elsewhere (x==0)
                            term = bn * sin(g*(x+1.0e-6_dp))/((1.0e-6_dp)*sin(g)) * exp(-g**2 * t(i))
                        endwhere
                    case default
                        return
                    end select

                    u(:,i) = u(:,i) + term
                end if
            end do

        if (sqrt(maxval(abs(term))) < eps) exit
        end do
    return
end subroutine batch_particle

subroutine batch_bulk(t, v, bi, a, s, eps)
    real(dp), intent(in) :: t(:) 
    real(dp), intent(out) :: v(size(t))
    real(dp), intent(in) :: bi, a 
    integer, intent(in) :: s 
    real(dp), intent(in) :: eps 

    real(dp) :: term(size(t))
    real(dp) :: g, g_old, cn, nu
    integer :: n
       
    if (a==0) then
        v = 1
    else
        nu = float(s-1)/2.0
        g = 0
        v = 1.0_dp/(1.0_dp + a)
        do n = 1, 100
            g_old = g
            do while(g_old >= g)
                g_old = g + 0.1_dp 
                call root_calc(g, g_old, bi, a, s, eps)
            end do
            if (bi == -1) then 
                cn = 4*(nu+1)*a/(g**2 + 4*a*(1+a)*(nu+1)**2)
            else
                cn = 4.*(nu+1.)*a* bi**2 /(g**4 + g**2*bi*(bi-4.*a*(nu+1.) &
            - 2.*nu) + 4.*a*(1. + a)*bi**2*(nu + 1.)**2)
            end if
            
            term = cn * exp(-g**2 * t)
            v = v + term
            if (maxval(abs(term)) < eps) exit
        end do
    end if

    return
end subroutine batch_bulk

subroutine root_calc(g, g_old, bi, a, s, eps)
    real(dp), intent(in) :: bi, a, eps, g_old
    integer, intent(in) :: s
    real(dp), intent(out) :: g

    real(dp) :: del, b, c
    ! Newton-Raphson Method
    b = g_old
    do  
        c = b + 0.01_dp
        if (f(b, bi, a, s)*f(c, bi, a, s) < 0) exit
        b = c
    end do
    g = c
    do
        del = -f(g, bi, a, s)/df(g, bi, a, s)
        g = g + del
        if(abs(del) < eps) exit
    end do
    return
end subroutine root_calc

real(dp) function f(g, bi, a, s)
    real(dp), intent(in) :: bi, a, g
    integer, intent(in) :: s

    if (a>0 .and. bi>0) then
        select case(s)
        case(0)
            f = (1 - a*bi/g**2)*sin(g) - bi/g*cos(g)
        case(1)
            f = (1 - 2.*a*bi/g**2)*bessel_j1(g) - bi/g*bessel_j0(g)
        case(2)
            f = (3*a*bi/g**3 + (bi-1)/g)*sin(g) + (1 - 3*a*bi/g**2)*cos(g)
        case default
            return
        end select
    
    else if(bi==-1) then
        select case(s)
        case(0)
            f = a/g*sin(g) + cos(g)
        case(1)
            f = 2.*a*bessel_j1(g) + g*bessel_j0(g)
        case(2)
            f = 3*a*(cos(g)/g - sin(g)/g**2)
        case default
            return
        end select

    else if (a==0) then        
        select case(s)
        case(0)
            f = g*tan(g) - bi
        case(1)
            f = g*bessel_j1(g) - bi*bessel_j1(g)
        case(2)
            f = g + (bi + 1)*tan(g)
        case default
            return
        end select
 
    else if (a==0 .and. bi==-1 .and. s==1) then
        f = bessel_j0(g)
    else
        return
    end if
end function

real(dp) function df(g, bi, a, s)
    real(dp), intent(in) :: bi, a, g
    integer, intent(in) :: s

    if (a>0 .and. bi>0) then
        select case(s)
        case(0)
            df = (bi/g + 2*a*bi/g**3)*sin(g) + (1+bi*(1-a)/g**2)*cos(g)
        case(1)
            df = (1-2*a*bi/g**2)*(bessel_j0(g)-bessel_j1(g)/g) + &
                (4*a*bi/g**3 + bi/g)*bessel_j1(g) + bi/g**2 *bessel_j0(g)
        case(2)
            df = (-9*a*bi/g**4 + (3*a*bi-bi+1)/g**2 - 1)*sin(g)+(9*a*bi/g**3 + (bi-1)/g)*cos(g)            
        case default
            return
        end select

    else if (bi == -1) then
        select case(s)
        case(0)
            df = a/g*cos(g) - (1+a/g**2)*sin(g)
        case(1)
            df = (2*a + 1)*bessel_j0(g) - (g**2 + 2*a)/g*bessel_j1(g)
        case(2)
            df = 3*a*(-2*cos(g)/g +(2/g**3 - 1/g)*sin(g))
        case default
            return
        end select

    else if (a==0) then        
        select case(s)
        case(0)
            df = tan(g) + g/(cos(g)**2)
        case(1)
            df = g*bessel_j0(g) + bi*bessel_j1(g)
        case(2)
            df = 1 + (bi - 1)/(sin(g)**2)
        case default
            return
        end select

    else if (a==0 .and. bi==-1 .and. s==1) then
        df = -bessel_j1(g)
    else
        return
    end if
end function 

end module batch
