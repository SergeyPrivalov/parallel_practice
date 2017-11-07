program main

real cf
parameter (cf = 0.00667)
parameter (N=50)
parameter (dx = 204, dy = 204)
real A(N*N, N*N), H1(50, 50), H2(50, 50), F(50, 50), X(50, 50), Y(50, 50)


open(4, FILE="data.txt")

do i = 1,N
    do j = 1,N
        READ(4,*) X(i,j), Y(i,j), H1(i,j), H2(i,j), F(i,j)
    enddo
enddo

open(5, FILE='A.txt')
do k = 1,N
    do l = 1,N

        do i = 1,N
            do j = 1,N
                A((l-1)*N+k, (i-1)*N+j) = &
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H1(i,j) ** 2 ) ** 0.5 - &
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H2(i,j) ** 2 ) ** 0.5
                
                A((l-1)*N+k, (i-1)*N+j) = A((l-1)*N+k, (i-1)*N+j) * dx * dy * cf
             enddo
        enddo
        
    enddo
enddo

do i = 1,10
    write(5, 10) (A(i, j), j = 1,10)
    10     FORMAT(10F15.12, ' ')
enddo

end