program main
    implicit none
    
    integer N, N2, i, j
    parameter(N = 50, N2 = N**2)
    
    double precision  X(N), Y(N),  H1(N2), H2(N2), B(N2), A(N2, N2), Z(N2)
    
    open(30, FILE="f_1sq.dat")
    do i = 1,N
        do j = 1,N
            read(30, *) X(j), Y(i), B((i-1)*N + j)
        enddo
    enddo

    open(50, FILE="hh1.dat")
    do i = 1,N
        do j = 1,N
            read(50, *) X(i), Y(j), H1((j-1)*N + i)
        enddo
    enddo

    open(60, FILE="hh2.dat")
    do i = 1,N
        do j = 1,N
            read(60, *) X(i), Y(j), H2((j-1)*N + i)
        enddo
    enddo
    CALL myLog('[OK] Data readed')
        
    A = evalA(X, Y, H1, H2, N)
    CALL myLog("[OK] A")
        
    Z = withRegular(A, B, N2)
    CALL myLog("[OK] Z")
     
    open(40, FILE='Z.txt')
    do i = 1,N 
        do j = 1,N
            write(40,*) X(j), Y(i), Z((i-1)*N+j)
        enddo
    enddo 
    CALL myLog("[END] Z.txt written")
    
    contains
    
        function evalA(X, Y, H1, H2, N)
            double precision :: X(N), Y(N), H1(N**2), H2(N**2)
            integer N
            double precision, dimension(N**2, N**2):: evalA
            
            integer i, j, k, l, ai, aj
            double precision  A(N**2,N**2), G, dx, dy
            parameter (G = 0.00667408)
            
            dx = X(2) - X(1)
            dy = Y(2) - Y(1)
            
            do k = 1,N; do l = 1,N
                do i = 1,N; do j = 1,N
                    ai = (l-1)*N + k
                    aj = (i-1)*N + j
                    A(ai, aj) = &
                        1/sqrt( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H1(aj) ** 2 ) - &
                        1/sqrt( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H2(aj) ** 2 )
                        
                    A(ai, aj) = A(ai, aj) * dx * dy * G
                enddo; enddo
            enddo; enddo
            
            evalA = A(:, :)
        end function evalA
    
    
        double precision  function maxL(A, N)
            double precision:: A(:, :)
            integer N, i
            
            double precision XPrev(N), XNext(N), eps, norm
            parameter(eps = 0.001)
            
            XPrev(1) = 1
            XPrev(2:) = 0
            
            norm = 1
            do while (norm > eps)
                XNext = matmul(A, XPrev)
                norm = sqrt(dot_product(XNext, XNext))
                
                do i = 1, N
                    XNext(i) = XNext(i) / norm
                enddo
                
                norm = 0
                do i = 1, N
                    norm = norm + (XNext(i) - XPrev(i)) ** 2
                enddo
                norm = sqrt(norm)
             
                do i = 1, N
                    XPrev(i) = XNext(i)
                enddo
            enddo
            
            XNext = matmul(A, XPrev)
            maxL = sqrt(dot_product(XNext, XNext))
        end function maxL
    
    
        function simpleIteration(A, B, N)
            integer N
            double precision:: A(:, :), B(:)
            double precision, dimension(N):: simpleIteration
            
            integer i, j, k
            double precision l, normB, eps, alpha, ZPrev(N), ZNext(N), diff
            parameter(eps = 0.001, alpha = 0.0001)
            character(len=50) message
            
            l = maxL(A, N)
            write(message, *) "[OK] maxL =", l 
            CALL myLog(message)
            

            normB = sqrt(dot_product(B, B))
            write(message, *) "[OK] normB =", normB 
            CALL myLog(message)
            
            ZPrev(:) = 0
            diff = evalDiff(A, ZPrev, B, N) / normB
            Print*, 0, diff
            
            k = 1
            do while(diff > eps)
                ZNext = matmul(A, ZPrev)
                ZNext = ZNext(:) + alpha * ZPrev(:)
                ZNext = ZPrev(:) - (ZNext(:) - B(:)) / (l + 0.5)
                
                diff = evalDiff(A, ZNext, B, N) / normB
                
                Print*, k, diff
                ZPrev = ZNext(:)
                k = k + 1
            enddo
            
            simpleIteration = ZNext(:)
        end function simpleIteration

        
        double precision function evalDiff(A, Z, B, N)
            integer N
            double precision:: A(:, :), B(:), Z(:)
            
            double precision X(N)
            X = matmul(A, Z) - B(:)
            evalDiff = sqrt(dot_product(X,X))
        end function evalDiff
        
        
        function withRegular(A, B, N)
            integer N
            double precision:: A(:, :), B(:)
            double precision, dimension(N):: withRegular
            
            double precision:: AT(N,N), ATA(N,N), ATB(N)
            
            AT = transpose(A)
            ATA = matmul(AT, A)
            ATB = matmul(AT, B)
            call myLog("[OK] ATA, ATB")
            
            withRegular = simpleIteration(ATA, ATB, N)
        end function withRegular

        SUBROUTINE myLog(message)
            character(LEN=*):: message

            REAL TIME
            call cpu_time(TIME)
            Print*, message, time 
        end SUBROUTINE myLog
end program main