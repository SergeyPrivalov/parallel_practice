program primal

parameter (N=100)
double precision H1(N*N), H2(N*N), DEN(N*N), X(N), Y(N), a, b, F, G, dx, dy
parameter (G = 6.67408)

minX = -200
maxX = -2
step = 2

do i = 1, N
    X(i) = minX + step * (i-1)
    Y(i) = minX + step * (i-1)
enddo

dx = ABS(X(N) - X(1)) / (N - 1)
dy = ABS(Y(N) - Y(1)) / (N - 1)
write (*, *) dx, dy

open(10, FILE="data/H1_100x100.dat")
do i = 1,N*N
    READ(10, *) a, b, H1(i)
enddo

open(20, FILE="data/H2_100x100.dat")
do i = 1,N*N
    READ(20, *) a, b, H2(i)
enddo

open(30, FILE="data/density2_100x100.dat")
do i = 1,N*N
    READ(30, *) a, b, DEN(i)
enddo

open(40, FILE="f_100x100.dat")
! X(l), Y(k)
do k = 1,N
    do l = 1,N
    F = 0
! ite
    ! X(j), Y(i)
    do i = 1,N
        do j = 1,N
        F = F + ( &
            1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H1((i-1) * N + j) ** 2 ) ** 0.5 - &
            1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H2((i-1) * N + j) ** 2 ) ** 0.5 &
            ) * DEN((i-1) * N + j)
        enddo
    enddo
    F = F * G * dx * dy
    WRITE(40, 8) X(l), Y(k), F
    enddo
enddo
8 FORMAT(F25.20, ' ', F25.20, ' ', F25.20)

write(*,*) "Hello, world"

end
