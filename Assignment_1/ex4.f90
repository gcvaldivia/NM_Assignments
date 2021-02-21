REAL(8) FUNCTION gammln(xx) RESULT(res)
        REAL(8),INTENT(IN) :: xx
        real(8), DIMENSION(6) :: cof
        REAL(8) :: ser,stp,tmp,x,y
        INTEGER(8) :: j
        
        cof(1)=76.18009172947146d0
        cof(2)=-86.50532032941677d0
        cof(3)=24.01409824083091d0
        cof(4)=-1.231739572450155d0
        cof(5)=0.1208650973866179d-2
        cof(6)=-0.5395239384953d-5
        stp=2.5066282746310005d0
        x=xx
        y=x
        tmp=x+5.5d0
        tmp=(x+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0

        DO j=1,6
                y=y+1.d0
                ser=ser+cof(j)/y
        END DO
        res=tmp+log(stp*ser/x)
END FUNCTION gammln

SUBROUTINE gaulag(x,w,n,alf)
        INTEGER(8) :: i,j,its,n,MAXIT
        REAL(8) :: ai,alf,EPS,p1,p2,p3,pp,z,z1
        REAL(8), DIMENSION(n) :: x
        REAL(8), DIMENSION(n) :: w
        REAL(8),EXTERNAL :: gammln
        EPS=3.d-14
        MAXIT=10
        
        DO i=1,n
                IF(i==1) THEN
                        z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
                ELSE IF (i==2) THEN
                        z=z+(15.+6.25*alf)/(1.+9*alf+2.5*n)
                ELSE
                        ai=i-2
                        z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*(z-x(i-2))/(1.+.3*alf)
                END IF

                DO its=1,MAXIT
                        p1=1.d0
                        p2=0.d0

                        DO j=1,n
                                p3=p2
                                p2=p1
                                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
                        END DO

                        PP=(n*p1-(n+alf)*p2)/z
                        z1=z
                        z=z1-p1/pp
                        IF (abs(z-z1).le.EPS) THEN
                                go to 1
                        END IF
                END DO
1               x(i)=z
                w(i)=-EXP(gammln(alf+n)-gammln(1.d0*n))/(pp*n*p2)
        END DO
END SUBROUTINE gaulag

REAL(8) FUNCTION f(x) RESULT(res) 
        REAL(8),INTENT(IN) :: x
        res=SIN(x)
END FUNCTION f

PROGRAM ex4
        IMPLICIT NONE
        INTEGER(8), DIMENSION(4) :: m
        REAL(8), DIMENSION(:), ALLOCATABLE :: x
        REAL(8), DIMENSION(:), ALLOCATABLE :: W
        REAL(8) :: alf,integral
        REAL(8),EXTERNAL :: f
        INTEGER(8) :: i,j,n

        PRINT*,""
        PRINT*,"Let us compute the integral of the function f(x)=exp(-x)*sin(x)"
        PRINT*,"from 0 to +oo  applying the Gauss-Laguerree Quadrature method:"
        PRINT*,""
        
        m=(/12,24,36,48/)
        PRINT'(A37)',"Approximate Result"
        PRINT*,""

        DO j=1,4
                alf=0
                n=m(j)
                ALLOCATE(x(n))
                ALLOCATE(w(n))
                CALL gaulag(x,w,n,alf)
                integral=0
                DO i=1,n
                        integral=integral+w(i)*f(x(i))
                END DO
                PRINT '(A3,I2,A1,F31.20)',"n=",n,":",integral
                DEALLOCATE(x)
                DEALLOCATE(w)
        END DO
        PRINT*,"Exact result: 0.5"
        PRINT*,""
END PROGRAM ex4
