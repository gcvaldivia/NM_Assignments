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

REAL(8) FUNCTION f(x) RESULT(res) 
        REAL(8),INTENT(IN) :: x
        res=7/(x**2+1)
END FUNCTION f

PROGRAM ex3
        IMPLICIT NONE
        INTEGER(8), DIMENSION(4) :: m
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: w
        REAL(8) :: x1,x2,integral
        INTEGER(8) :: i,j,n
        REAL(8), EXTERNAL :: f
        
        PRINT*,""
        PRINT*,"Let us compute the integral of the function f(x)=7/(x^2+1)"
        PRINT*,"from 0 to 5 applying the Gauss-Legendre Quadrature method:"
        PRINT*,""
        
        m=(/12,24,36,48/)
        PRINT'(A37)',"Approximate Result"
        PRINT*,""
        
        DO j=1,4
                n=m(j)
                ALLOCATE(x(n))
                ALLOCATE(w(n))
                x1=0
                x2=5
                CALL gauleg(x1,x2,x,w,n)
                integral=0
                DO i=1,n
                        integral=integral+w(i)*f(x(i))
                END DO
                PRINT '(A3,I2,A1,F31.20,F31.20)',"n=",n,":",integral
                DEALLOCATE(x)
                DEALLOCATE(w)
        END DO
        PRINT*,"Exact result: 9.61380536861511102603"
        PRINT*,"" 
END PROGRAM ex3
