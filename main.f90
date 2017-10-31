program main

parameter (N=50)
parameter (cf = 0.00667)
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
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H1(i,j) ) ** 0.5 + &
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H2(i,j) ) ** 0.5
            enddo
        enddo
        write(5, 10) (A((l-1)*N+k, i), i = 1,20)
 10     FORMAT(20F15.12, ' ')
    enddo
enddo

end