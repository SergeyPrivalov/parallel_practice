! константы
real cf
parameter (cf = 0.00667)
parameter (N=50)
parameter (dx = 204, dy = 204)

! входные данные
real A(N*N, N*N), H1(N, N), H2(N, N), X(N, N), Y(N, N)

! input
open(4, FILE="data.txt")
do i = 1,N
    do j = 1,N
        READ(4,*) X(i,j), Y(i,j), H1(i,j), H2(i,j)
    enddo
enddo

! вычисление матрицы A
do k = 1,N
    do l = 1,N

        do i = 1,N
            do j = 1,N
                A((k-1)*N+l, (i-1)*N+j) = &
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H1(i,j) ** 2 ) ** 0.5 - &
                1/( (X(i, j) - X(k,l)) ** 2 + (Y(i,j) - Y(k,l)) ** 2 + H2(i,j) ** 2 ) ** 0.5
                
                A((l-1)*N+k, (i-1)*N+j) = A((l-1)*N+k, (i-1)*N+j) * dx * dy * cf
             enddo
        enddo
        
    enddo
enddo

! вывод матрицы A.txt
20 open(5, FILE='A.txt')
do i = 1, N*N
    write(5, 10) (A(i, j), j = 1,N*N)
    10     FORMAT(10F10.7, ' ')
enddo

end