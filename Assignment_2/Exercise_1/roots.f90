!**********************************************************************
!                    Non-linear equations/Program          
!                   Gustavo Valdivia, ICTP - EAIFR
!**********************************************************************

!**********************************************************************
!**********************************************************************
PROGRAM test
    USE md
    IMPLICIT NONE
    REAL(8) :: a,b,root,fr
    REAL(8) :: rootn,frn,roots,frs
    INTEGER(8) :: nit,nitn,nits
    a=1.25d0
    b=1.5d0
    call bisection(fun,a,b,root,fr,nit)
    call newtonraphson(fun,dfun,a,rootn,frn,nitn)
    call secant(fun,a,b,roots,frs,nits)

    PRINT*,''
    PRINT '(A19)',"1. Bisection Method"
    PRINT '(A10,F23.19)',"Root |E| =",root
    PRINT '(A6,F23.19)',"f(E) =",fr
    PRINT '(A31,I3)',"Number of iterations required =",nit
    PRINT*,''
    
    PRINT '(A24)',"2. Newton-Raphson Method"
    PRINT '(A10,F23.19)',"Root |E| =",rootn
    PRINT '(A6,F23.19)',"f(E) =",frn
    PRINT '(A31,I3)',"Number of iterations required =",nitn
    PRINT*,''
    
    PRINT '(A16)',"3. Secant Method"
    PRINT '(A10,F23.19)',"Root |E| =",roots
    PRINT '(A6,F23.19)',"f(E) =",frs
    PRINT '(A31,I3)',"Number of iterations required =",nits
    PRINT*,''
END PROGRAM test 
!**********************************************************************
!**********************************************************************
