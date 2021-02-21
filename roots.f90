!**********************************************************************
!                         Non-linear equations           
!                   Gustavo Valdivia, ICTP - EAIFR
!**********************************************************************

!**********************************************************************
PROGRAM test
    USE md
    IMPLICIT NONE
    REAL(8) :: a,b,root,fr
    !REAL(8), EXTERNAL :: fun,dfun
    REAL(8) :: rootn,frn,roots,frs
    INTEGER(8) :: nit,nitn,nits
    a=1.3d0
    b=1.5d0
    call newtonraphson(fun,dfun,a,rootn,frn,nitn)
    call bisection(fun,a,b,root,fr,nit)
    call secant(fun,a,b,roots,frs,nits)

    PRINT*,''
    PRINT '(A19)',"1. Bisection Method"
    PRINT '(A10,F23.19)',"Root |E| =",root
    PRINT '(A8,F23.19)',"f(|E|) =",fr
    PRINT '(A31,I3)',"Number of iterations required =",nit
    PRINT*,''
    
    PRINT '(A24)',"2. Newton-Raphson Method"
    PRINT '(A10,F23.19)',"Root |E| =",rootn
    PRINT '(A8,F23.19)',"f(|E|) =",frn
    PRINT '(A31,I3)',"Number of iterations required =",nitn
    PRINT*,''
    
    PRINT '(A16)',"3. Secant Method"
    PRINT '(A10,F23.19)',"Root |E| =",roots
    PRINT '(A8,F23.19)',"f(|E|) =",frs
    PRINT '(A31,I3)',"Number of iterations required =",nits
    PRINT*,''
END PROGRAM test 
!**********************************************************************