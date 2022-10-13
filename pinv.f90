PROGRAM PINV 
   INTEGER :: K,ITERS,M,N
   REAL :: T1,T2
   DOUBLE PRECISION :: I(45,45),A(45,30),S(30),X(30,45),A_C(30,30),alpha,NORM
   M = 45
   N = 30
   K = 0

   DO I_=1,M
    DO J_=1,N
       A(I_,J_)= I_**2 + J_**2
     END DO  
   END DO

   DO I_ = 1,M
        DO J_ = 1,M
        IF(I_ .EQ. J_) THEN
                I(I_,J_) = 1
        ELSE
                I(I_,J_) = 0
        END IF
        END DO
   END DO
   A_C = MATMUL(TRANSPOSE(A),A)
   CALL SVD(A_C,N,N,S) 
   alpha=  MAXVAL(S)
   X = (1/alpha**2)*TRANSPOSE(A)
   
   PRINT *, 'Newton-Schulz Scheme'
   CALL CPU_TIME(T1)
   CALL NEWTON_SCHULZ(A,X,M,N,K,100,NORM)
   CALL CPU_TIME(T2)
   print *, 'K = '
   print *,K
   print *, 'ERROR = '
   print *, NORM
   PRINT *, 'Time = '
   PRINT *, T2-T1

   X = (1/alpha**2)*TRANSPOSE(A)
   K = 0
   NORM = 0

   PRINT *, 'Chebyshev Scheme'
   CALL CPU_TIME(T1)
   CALL CHEBYSHEV(A,I,X,M,N,K,100,NORM)
   CALL CPU_TIME(T2)
   PRINT *, 'K = '
   PRINT *, K
   PRINT *, 'ERROR = '
   PRINT *, NORM
   PRINT *, 'TIME = '
   PRINT *, T2-T1
   
   PRINT *,'Homeier Scheme'
   K = 0
   NORM = 0
   X = (1/alpha**2)*TRANSPOSE(A)
   CALL CPU_TIME(T1)
   CALL HOMEIER(A,I,X,M,N,K,100,NORM)
   CALL CPU_TIME(T2)
   PRINT *, 'K = '
   PRINT *, K
   PRINT *, 'ERROR = '
   PRINT *, NORM
   PRINT *, 'TIME = '
   PRINT *, T2-T1

   PRINT *, 'Mid Point scheme'
   X = (1/alpha**2)*TRANSPOSE(A)
   K = 0
   CALL CPU_TIME(T1)
   CALL MIDPOINT(A,I,X,M,N,K,100,NORM)
   CALL CPU_TIME(T2)

   PRINT *, 'K = '
   PRINT *, K
   PRINT *, 'ERROR = '
   PRINT *, NORM
   PRINT *, 'TIME = '
   PRINT *, T2-T1
   CONTAINS
   SUBROUTINE SVD(A,M,N,S)

      DOUBLE PRECISION A(M,N), U(M,M),VT(N,N),S(N),V(N,N)
      DOUBLE PRECISION,ALLOCATABLE :: WORK(:)
      INTEGER LDA,M,N,LWORK,LDVT,INFO
      CHARACTER  JOBU, JOBVT

      JOBU='A'
      JOBVT='A'
      LDA=M
      LDU=M
      LDVT=N

      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))

      ALLOCATE(work(lwork))

      CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO )

      END 
    
    SUBROUTINE NEWTON_SCHULZ(A,X,M,N,K,ITERS,NORM)
            DOUBLE PRECISION :: A(M,N),X(N,M),NORM
            DO I_=1,ITERS
                X = 2*X - MATMUL(MATMUL(X,A),X) 
                NORM  =NORM2(MATMUL(MATMUL(A,X),A)-A)
                IF(NORM<0.00001) THEN

                        EXIT
                END IF   
                K = K + 1
            END DO    
    END

    SUBROUTINE CHEBYSHEV(A,I,X,M,N,K,ITERS,NORM)
        DOUBLE PRECISION :: I(M,M),A(M,N),X(N,M),Y(M,M),NORM
        DO I_ = 1,ITERS
              Y = MATMUL(A,X)  
              X = MATMUL(X,3*I-MATMUL(Y,3*I-Y))
              NORM  =NORM2( MATMUL(MATMUL(A,X),A)-A)
              IF(NORM < 0.00001) THEN
                      EXIT
              END IF

              K = K + 1
              
        END DO
    END

    SUBROUTINE HOMEIER(A,I,X,M,N,K,ITERS,NORM)
        DOUBLE PRECISION :: I(M,M),A(M,N),X(N,M),Y(M,M),NORM
        DO I_ = 1,ITERS
            Y = MATMUL(A,X)
            X = MATMUL(X,I+MATMUL((1/2)*(I-Y),(I+(2*I-Y)**2)))
            NORM = NORM2( MATMUL(MATMUL(A,X),A)-A)
            IF(NORM < 0.00001) THEN
                    EXIT
            END IF
            K = K+1
        END DO
    END 
    SUBROUTINE MIDPOINT(A,I,X,M,N,K,ITERS,NORM)
        DOUBLE PRECISION :: I(N,N),A(M,N),X(N,M),Y(N,N),NORM
        DO I_ = 1,ITERS
                Y = MATMUL(X,A)
                X = MATMUL(I+MATMUL(1/4*(I-Y),(3*I-Y)**2),X)
                NORM = NORM2(MATMUL(MATMUL(A,X),A)-A)
                IF(NORM < 0.00001) THEN
                        EXIT
                END IF
                K = K + 1
        END DO
    
    END
END PROGRAM PINV
