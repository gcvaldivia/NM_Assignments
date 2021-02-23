!*******************************************************************************
!                  6-dimensional Integral (Gauss-Legendre)           
!                      Gustavo Valdivia, ICTP-EAIFR
!******************************************************dd*************************


!This is the function that will be integrated **********************************
!*******************************************************************************
REAL(8) FUNCTION func(t,x1,y1,z1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: t,x1,y1,z1,x2,y2,z2
    REAL(8) :: part1,part2
    part1=-4d0*DSQRT(x1**2+y1**2+z1**2)-4d0*DSQRT(x2**2+y2**2+z2**2)
    part2=-(t**2)*(x1-x2)**2-(t**2)*(y1-y2)**2-(t**2)*(z1-z2)**2
    res=2d0*DEXP(part1+part2)/DSQRT(3.1415926535d0)
END FUNCTION func
!*******************************************************************************
!*******************************************************************************


!Limits of integration *********************************************************
!*******************************************************************************
REAL(8) FUNCTION tdn(x1,y1,z1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,y1,z1,x2,y2,z2
    res=0
END FUNCTION tdn
REAL(8) FUNCTION tup(x1,y1,z1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,y1,z1,x2,y2,z2
    res=112d0
END FUNCTION tup
REAL(8) FUNCTION z1dn(x1,y1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,y1,x2,y2,z2
    res=-10d0
END FUNCTION z1dn
REAL(8) FUNCTION z1up(x1,y1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,y1,x2,y2,z2
    res=10d0
END FUNCTION z1up
REAL(8) FUNCTION y1dn(x1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,x2,y2,z2
    res=-10d0
END FUNCTION y1dn
REAL(8) FUNCTION y1up(x1,x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x1,x2,y2,z2
    res=10d0
END FUNCTION y1up
REAL(8) FUNCTION x1dn(x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x2,y2,z2
    res=-10d0
END FUNCTION x1dn
REAL(8) FUNCTION x1up(x2,y2,z2) RESULT(res)
    REAL(8), INTENT(IN) :: x2,y2,z2
    res=10d0
END FUNCTION x1up
REAL(8) FUNCTION z2dn(x2,y2) RESULT(res)
    REAL(8), INTENT(IN) :: x2,y2
    res=-10d0
END FUNCTION z2dn
REAL(8) FUNCTION z2up(x2,y2) RESULT(res)
    REAL(8), INTENT(IN) :: x2,y2
    res=10d0
END FUNCTION z2up
REAL(8) FUNCTION y2dn(x2) RESULT(res)
    REAL(8), INTENT(IN) :: x2
    res=-10d0
END FUNCTION y2dn
REAL(8) FUNCTION y2up(x2) RESULT(res)
    REAL(8), INTENT(IN) :: x2
    res=10d0
END FUNCTION y2up
!*******************************************************************************
!*******************************************************************************


!Abscissas and weights of the Gauss-Legendre n-point quadrature formula ********
!*******************************************************************************
SUBROUTINE gauleg(x1,x2,x,w,n)
        INTEGER(8) :: n
        REAL(8), DIMENSION(n) :: x
        REAL(8), DIMENSION(n) :: w
        REAL(8) :: x1,x2
        REAL(8) :: EPS
        INTEGER(8) :: i,j,m
        REAL(8) :: p1,p2,p3,pp,xl,xm,z,z1
        
        EPS=3.d-14
        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)

        DO i=1,m
                z=COS(3.141592654d0*(i-0.25d0)/(n+0.5d0))
1               p1=1.d0
                p2=0.d0
        
                DO j=1,n
                        p3=p2
                        p2=p1
                        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
                END DO

                pp=n*(z*p1-p2)/(z*z-1.d0)
                z1=z
                z=z1-p1/pp
               
                IF (abs(z-z1).gt.EPS) THEN
                        go to 1
                END IF
                
                x(i)=xm-xl*z
                x(n+1-i)=xm+xl*z
                w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
                w(n+1-i)=w(i)
        END DO
END SUBROUTINE gauleg
!*******************************************************************************
!*******************************************************************************


!Subroutines for the integration over each variable ****************************
!*******************************************************************************
SUBROUTINE qgaust(fun,tdn,tup,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: tdn,tup,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(tdn,tup,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgaust
SUBROUTINE qgausx1(fun,x1dn,x1up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: x1dn,x1up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(x1dn,x1up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausx1
SUBROUTINE qgausy1(fun,y1dn,y1up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: y1dn,y1up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(y1dn,y1up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausy1
SUBROUTINE qgausz1(fun,z1dn,z1up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: z1dn,z1up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(z1dn,z1up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausz1
SUBROUTINE qgausx2(fun,x2dn,x2up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: x2dn,x2up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(x2dn,x2up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausx2
SUBROUTINE qgausy2(fun,y2dn,y2up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: y2dn,y2up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(y2dn,y2up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausy2
SUBROUTINE qgausz2(fun,z2dn,z2up,ss)
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: z2dn,z2up,integral,ss
        INTEGER(8) :: i,j,n=18     !******** Here we set the value of n ********
        REAL(8), EXTERNAL :: fun
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gauleg(z2dn,z2up,x,w,n)
        integral=0
        DO i=1,n
        integral=integral+w(i)*fun(x(i))
        END DO
        ss=integral
END SUBROUTINE qgausz2
!*******************************************************************************
!*******************************************************************************


!Partial results ***************************************************************
!*******************************************************************************
REAL(8) FUNCTION d1(tt)
    REAL(8) :: tt,t,x1,y1,z1,x2,y2,z2
    REAL(8), EXTERNAL :: func
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    t=tt
    d1=func(t,x1,y1,z1,x2,y2,z2)
    return
END FUNCTION d1

REAL(8) FUNCTION f1(zz)
    REAL(8) :: zz,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: d1,tdn,tup
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    z1=zz
    call qgaust(d1,tdn(x1,y1,z1,x2,y2,z2),tup(x1,y1,z1,x2,y2,z2),ss)
    f1=ss
    return
END FUNCTION f1

REAL(8) FUNCTION g1(yy)
    REAL(8) :: yy,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: f1,z1dn,z1up
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    y1=yy
    call qgausz1(f1,z1dn(x1,y1,x2,y2,z2),z1up(x1,y1,x2,y2,z2),ss)
    g1=ss
    return
END FUNCTION g1
REAL(8) FUNCTION h1(xx) 
    REAL(8) :: xx,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: g1,y1dn,y1up
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    x1=xx
    call qgausy1(g1,y1dn(x1,x2,y2,z2),y1up(x1,x2,y2,z2),ss)
    h1=ss
    return
END FUNCTION h1
REAL(8) FUNCTION f2(zz1)
    REAL(8) :: zz1,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: h1,x1dn,x1up
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    z2=zz1
    call qgausx1(h1,x1dn(x2,y2,z2),x1up(x2,y2,z2),ss)
    f2=ss
    return
END FUNCTION f2
REAL(8) FUNCTION g2(yy1)
    REAL(8) :: yy1,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: f2,z2dn,z2up
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    y2=yy1
    call qgausz2(f2,z2dn(x2,y2),z2up(x2,y2),ss)
    g2=ss
    return
END FUNCTION g2
REAL(8) FUNCTION h2(xx1) 
    REAL(8) :: xx1,t,x1,y1,z1,x2,y2,z2,ss
    REAL(8), EXTERNAL :: g2,y2dn,y2up
    COMMON /tx1y1z1x2y2z2/ t,x1,y1,z1,x2,y2,z2
    x2=xx1
    call qgausy2(g2,y2dn(x2),y2up(x2),ss)
    h2=ss
    return
END FUNCTION h2
!*******************************************************************************
!*******************************************************************************


!Program ***********************************************************************
!*******************************************************************************
PROGRAM integral6d
    IMPLICIT NONE
    REAL(8) :: ss,x2dn,x2up
    REAL(8), EXTERNAL :: h2
    x2dn=-10d0
    x2up=10d0
    call qgausx2(h2,x2dn,x2up,ss)
    PRINT*,ss
END PROGRAM integral6d
!*******************************************************************************
!*******************************************************************************