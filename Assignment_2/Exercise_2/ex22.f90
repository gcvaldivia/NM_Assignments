!********************************************************************************
!                   6-dimensional Integral (Gaus-Laguerre)           
!                       Gustavo Valdivia, ICTP-EAIFR
!********************************************************************************


!Abscissas and weights of the n-point Gauss-Laguerre quadrature formula *********
!for the integration over each variable *****************************************
REAL(8) FUNCTION gammln(xx) RESULT(res)
        REAL(8),INTENT(IN) :: xx
        real(8), DIMENSION(6) :: cof
        REAL(8) :: ser,stp,tmp,x,y
        INTEGER(8) :: j
        
        cof=(/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
        -1.231739572450155d0,0.1208650973866179d-2,-0.5395239384953d-5/)
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
SUBROUTINE gaulag1(x,w,n,alf)
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
END SUBROUTINE gaulag1
SUBROUTINE gaulag2(x,w,n,alf)
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
END SUBROUTINE gaulag2
!********************************************************************************
!********************************************************************************


!This is the function that will be integrated ***********************************
!********************************************************************************
REAL(8) FUNCTION func(r1,r2) RESULT(res) 
        REAL(8),INTENT(IN) :: r1,r2
        res=(3.1415926535d0**2)*(r1+r2-ABS(r1-r2))/(r1*128d0)
END FUNCTION func
!********************************************************************************
!********************************************************************************


!Partial results ****************************************************************
!********************************************************************************
REAL(8) FUNCTION f(rr2)
        REAL(8) :: rr2,r1,r2
        REAL(8), EXTERNAL :: func
        COMMON /r1r2/ r1,r2
        r2=rr2
        f=func(r1,r2)
        return
END FUNCTION f
REAL(8) FUNCTION g(rr1)
        REAL(8) :: rr1,r1,r2,ss,alf
        INTEGER(8) :: n
        REAL(8), EXTERNAL :: f
        REAL(8), DIMENSION(:), ALLOCATABLE :: x,w
        COMMON /r1r2/ r1,r2
        r1=rr1
        n=200            !************** Here we set the value of n *************
        alf=1d0
        ALLOCATE(w(n))
        ALLOCATE(x(n))
        call gaulag1(x,w,n,alf)
        ss=0
        DO i=1,n
            ss=ss+w(i)*f(x(i))
        END DO
        g=ss
        return
END FUNCTION g
!********************************************************************************
!********************************************************************************


!Program ************************************************************************
!********************************************************************************
PROGRAM integral6d
        IMPLICIT NONE
        REAL(8), DIMENSION(:), ALLOCATABLE :: x,w
        REAL(8) :: ss,alf
        INTEGER(8) :: n,i
        REAL(8),EXTERNAL :: g
        n=200             !************** Here we set the value of n *************
        alf=2d0
        ALLOCATE(x(n))
        ALLOCATE(w(n))
        CALL gaulag2(x,w,n,alf)
        ss=0
        DO i=1,n
            ss=ss+w(i)*g(x(i))
        END DO
        PRINT*,ss
END PROGRAM integral6d
!********************************************************************************
!********************************************************************************
