! file: batch.f90

! Solutions for the Diffusion differential equation
module batch
implicit none
public finite_batch_particle, finite_batch_bulk, root_calc, f, df
contains
subroutine finite_batch_particle(u, x, t, bi, a, s, eps)
    implicit none
    real(8), intent(in) :: x, t, bi, a, eps 
    integer, intent(in) :: s 
    real(8), intent(out) :: u

!f2py intent(in) x, t, bi, a, eps, s
!f2py intent(out) u 

    real(8) :: g, g_old, bn, term, nu
    integer :: n, i

    if (t==0) then
        u = 0.
    else
        nu = float(s-1)/2.0
        u = 1.0/(1.0 + a)
        g_old = 0.1
        g = 0
        do n = 1, 100
            i = 1
            do while (g_old >= g) 
                    g_old = g + float(i)*0.1
                    call root_calc(g, g_old, bi, a, s, eps)
            end do

            bn = 2.*bi*(2.*(nu+1.)*a*bi - g**2)/(g**4 + g**2*bi*& 
            (bi-4.*a*(nu+1.) - 2.*nu) + 4.*a*(1. + a)*bi**2*(nu + 1.)**2)

            select case(s)
            case(0)
                term = bn * cos(g*x)/cos(g) * exp(-g**2 * t)        
            case(1) 
                term =  bn * bessel_j0(g*x)/bessel_j0(g) * exp(-g**2 * t)
            case(2)
                term = bn * sin(g*x)/(x*sin(g)) * exp(-g**2 * t)
            case default
                return
            end select

            g_old = g 
            u = u + term
            !if (term < eps) exit
        end do
    end if
    return
end subroutine finite_batch_particle

subroutine finite_batch_bulk(v, t, bi, a, s, eps)
    implicit none
    real(8), intent(in) :: t, bi, a, eps 
    integer, intent(in) :: s 
    real(8), intent(out) :: v

!f2py intent(in) t, bi, eps, s
!f2py intent(out) v

    real(8) :: g, g_old, cn, term, nu
    integer :: n, i  

       
    if (t==0) then
        v = 1.
    else
        nu = float(s-1)/2.0
        g_old = 0.1
        g = 0
        v = 1.0/(1.0 + a)
        do n = 1, 100
            i = 1
            do while (g_old >= g) 
                    g_old = g + float(i)*0.1
                    call root_calc(g, g_old, bi, a, s, eps)
            end do
            cn = 4.*(nu+1.)*a* bi**2 /(g**4 + g**2*bi*(bi-4.*a*(nu+1.) &
            - 2.*nu) + 4.*a*(1. + a)*bi**2*(nu + 1.)**2)
            
            term = cn * exp(-g**2 * t)
            
            g_old = g 
            v = v + term
            if (term < eps) exit
        end do
    end if
    return
end subroutine finite_batch_bulk

subroutine root_calc(g, g_old, bi, a, s, eps)
    implicit none
    real(8), intent(in) :: bi, a, eps, g_old
    integer, intent(in) :: s
    real(8), intent(out) :: g

!f2py intent(in) bi, a, err, g_old, s
!f2py intent(out) g

    real(8) :: del, b, c
    ! Newton-Raphson Method
    b = g_old
    do  
        c = b + 0.1
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

real(8) function f(g, bi, a, s)
    real(8), intent(in) :: bi, a, g
    integer, intent(in) :: s
!f2py intent(in) bi, a, g, s
!f2py intent(out) f
    select case(s)
    case(0)
        f = (g**2 - a*bi)*sin(g) - g*bi*cos(g)
    case(1)
        f = (g**2 - 2.*a*bi)*bessel_j1(g) - g*bi*bessel_j0(g)
    case(2)
        f = (3.*a*bi/g + g*(bi-1.))*sin(g)+(g**2-3*a*bi)*cos(g)
    case default
        return
    end select
end function

real(8) function df(g, bi, a, s)
    real(8), intent(in) :: bi, a, g
    integer, intent(in) :: s
!f2py intent(in) bi, a, g, s
!f2py intent(out) df
    select case(s)
    case(0)
        df = g*(bi + a)*sin(g) + (g**2 - a + bi)*cos(g)
    case(1)
        df = (g**2 - bi*(1 + 2.*a))*bessel_j0(g) + (g - 2.*a*bi/g &
        + g*bi)*bessel_j1(g)
    case(2)
        df = (bi-1.-3.*a*bi/g**2-g**2 + 3.*a*bi)*sin(g)+ & 
        (3.*a*bi/g + g*bi + g)*cos(g)
    case default
        return
    end select
end function 

end module batch
!end file batch.f90