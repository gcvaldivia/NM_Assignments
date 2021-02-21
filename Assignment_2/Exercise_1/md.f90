!**********************************************************************
!                     Non-linear equations/Module          
!                   Gustavo Valdivia, ICTP - EAIFR
!**********************************************************************

!**********************************************************************
MODULE md
    IMPLICIT NONE
    CONTAINS
    REAL(8) FUNCTION fun(E) RESULT(res)
        REAL(8), INTENT(IN) :: E
        REAL(8) :: m=938d0
        REAL(8) :: v=60d0
        REAL(8) :: a=1.45d0/197.3d0
        res=DSQRT(m*(v-E))*DCOTAN(a*DSQRT(m*(v-E)))+DSQRT(m*E)
    END FUNCTION fun
    REAL(8) FUNCTION dfun(E) RESULT(res)
        REAL(8), INTENT(IN) :: E
        REAL(8) :: part1,part2
        part1=DSQRT(469d0/2d0)*DCOTAN(0.225083d0*DSQRT(60d0-E))
        part2=3.44678d0/((DSIN(0.225083d0*DSQRT(60d0-E)))**2)
        res=DSQRT(469d0/2d0)/DSQRT(E)-part1/(DSQRT(60d0-E))+part2
    END FUNCTION dfun
    SUBROUTINE bisection(fun,a,b,root,fr,nit)
        IMPLICIT NONE
        INTEGER(8) :: i,nit,N=50
        REAL(8), EXTERNAL :: fun
        REAL(8) :: a,b,p,fa,fp
        REAL(8) :: tol=1d-10
        REAL(8) :: lb=1d-12
        REAL(8) :: root,fr
     
        fa=fun(a)
        DO i=1,N
                p=a+(b-a)/2
                fp=fun(p)
                IF ((abs(fp)<lb).or.((b-a)/2<tol)) THEN
                        root=p
                        fr=fp
                        nit=i
                        exit
                END IF
                IF (fa*fp>0) THEN
                        a=p
                        fa=fp
                ELSE
                        b=p
                END IF                
        END DO
    END SUBROUTINE bisection
    SUBROUTINE newtonraphson(fun,dfun,a,root,fr,nit)
        IMPLICIT NONE
        INTEGER(8) :: i,nit,N=50
        REAL(8), EXTERNAL :: fun,dfun
        REAL(8) :: p,a
        REAL(8) :: tol=1d-10
        REAL(8) :: root,fr
        DO i=1,N
                p=a-fun(a)/dfun(a)
                IF (abs(p-a)<tol)  THEN
                        root=p
                        fr=fun(p)
                        nit=i
                        exit
                END IF
                a=p
        END DO
    END SUBROUTINE newtonraphson
    SUBROUTINE secant(fun,p0,p1,root,fr,nit)
        IMPLICIT NONE
        INTEGER(8) :: i,nit,N=50
        REAL(8), EXTERNAL :: fun
        REAL(8) :: p
        REAL(8) :: P0,p1
        REAL(8) :: q0,q1
        REAL(8) :: tol=1d-10
        REAL(8) :: root,fr
        q0=fun(p0)
        q1=fun(p1)
        DO i=2,N
                p=p1-q1*(p1-p0)/(q1-q0)
                IF (abs(p-p1)<tol)  THEN
                        root=p
                        fr=fun(p)
                        nit=i
                        exit
                END IF
                p0=p1
                q0=q1
                p1=p
                q1=fun(p)
        END DO
    END SUBROUTINE secant
END MODULE md
!**********************************************************************
