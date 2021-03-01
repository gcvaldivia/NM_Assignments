!************************************************************
!            Gustavo Valdivia, ICTP-EAIFR
!                    Assignment 3
!************************************************************



! RK4 Subroutine ********************************************
subroutine rk4(t0,tf,n,yinits,ts,ys)
    integer(8) :: n
    real(8), dimension(2) :: yinits,yvals,k1,k2,k3,k4
    real(8), dimension(n) :: ts
    real(8), dimension(n,2) :: ys
    real(8) :: h,t0,tf,t1
    
    h=(tf-t0)/(n-1)
    
    ts(1)=t0
    ys(1,:)=yinits
    yvals=yinits
   
    do i=2,n
        t1=t0+h
        
        y0=yvals(1)
        y1=yvals(2)
        k1=(/2*y0-y0*y1,y0*y1-y1/)
        
        y0=yvals(1)+(h/2)*k1(1)
        y1=yvals(2)+(h/2)*k1(2)
        k2=(/2*y0-y0*y1,y0*y1-y1/)
        
        y0=yvals(1)+(h/2)*k2(1)
        y1=yvals(2)+(h/2)*k2(2)
        k3=(/2*y0-y0*y1,y0*y1-y1/)
        
        y0=yvals(1)+h*k3(1)
        y1=yvals(2)+h*k3(2)
        k4=(/2*y0-y0*y1,y0*y1-y1/)

        yvals=yvals+(h/6)*(k1+2*k2+2*k3+k4)
        
        ts(i)=t1
        ys(i,:)=yvals
        
        t0=t1
    end do
end subroutine rk4
!************************************************************


! Program ***************************************************
program a3
integer(8) :: n
real(8), dimension(2) :: yinits
real(8), dimension(:,:),allocatable :: ys
real(8), dimension(:),allocatable :: ts
real(8) :: t0,tf

n=500
t0=0d0
tf=49.9d0

allocate(ts(n))
allocate(ys(n,2))

yinits=(/1d0,2.7d0/)

call rk4(t0,tf,n,yinits,ts,ys)

open(1, FILE="data.dat",status='old')         
do i=1,n
    write(1,*) ts(i),ys(i,:)
end do
close(1)   

end program a3
!************************************************************
