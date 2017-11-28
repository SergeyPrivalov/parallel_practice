program primal
implicit none

integer i,j,k,l,N,maxX,minX,step

parameter (N=100)
double precision H1(N*N), H2(N*N), F(N*N), X(N), Y(N), aa, bb, G, dx, dy, A(N*N, N*N)
double precision XPrev(N*N), XNext(N*N), norm, eps, YCur(N*N), maxL, ZNext(N*N), ZPrev(N*N), normF
double precision nevyaz, az, alpha
parameter (G = 6.67408, eps = 0.001, alpha = 0.1)

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
    READ(10, *) aa, bb, H1(i)
enddo

open(20, FILE="data/H2_100x100.dat")
do i = 1,N*N
    READ(20, *) aa, bb, H2(i)
enddo

open(30, FILE="data/field2_100x100.dat")
do i = 1,N*N
    READ(30, *) aa, bb, F(i)
enddo
WRITE(*,*) "Data read"

!  norm b
normF = 0
do k = 1,N*N
    normF = normF + F(k)**2
enddo
normF = sqrt(normF)
write(*,*) "calculated normF", normF


! X(l), Y(k)
do k = 1,N
    do l = 1,N
        ! X(j), Y(i)
        do i = 1,N
            do j = 1,N
                A((l-1)*N+k, (i-1) * N + j) = &
                1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H1((i-1) * N + j) ** 2 ) ** 0.5 - &
                1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H2((i-1) * N + j) ** 2 ) ** 0.5
                
                A((l-1)*N+k, (i-1) * N + j) = A((l-1)*N+k, (i-1) * N + j) * dx * dy * G
             enddo
        enddo
    enddo
enddo
write(*,*) "calculated A"

! maxL иницилизация
norm = 1
do i = 1,N*N
    XPrev(i) = 0
enddo 
XPrev(1) = 1

! maxL счет
open(8, FILE='maxL.txt')
! READ(8,*) maxL
! close(8)

! WRITE (*,*) "read maxL"

k = 0
do while (norm > eps)
    norm = 0
    do i = 1, N*N
        YCur(i) = 0
        do j = 1, N*N
            YCur(i) = YCur(i) + A(i,j) * XPrev(j)
        enddo
        norm = norm + YCur(i)**2
    enddo
    norm = sqrt(norm)

    do i = 1, N*N
        XNext(i) = YCur(i) / norm
    enddo
    
    norm = 0
    do i = 1, N*N
        norm = norm + (XNext(i) - XPrev(i)) ** 2
    enddo
    norm = sqrt(norm)
    
    do i = 1, N*N
        XPrev(i) = XNext(i)
    enddo
    k = k + 1
enddo

do i = 1, N*N
    YCur(i) = 0
    do j = 1, N*N
        YCur(i) = YCur(i) + A(i,j) * XPrev(j)
    enddo
    norm = norm + YCur(i)**2
enddo
maxL = sqrt(norm)

write(*,*) "calculated maxL", maxL

! write(6,*) maxL


! решение, инициализация
do i = 1,N*N
    ZPrev(i) = 0
enddo

! nevyaz = 0
! do i = 1,N*N
!     az = 0
!     do j = 1,N*N
!         az = az + A(i,j) * ZPrev(j)
!     enddo
!     az = az - F(i)
!     nevyaz = nevyaz + az*az
! enddo
! nevyaz = sqrt(nevyaz) / normF
! write(*,*) nevyaz
!
!      считаем z
!
open(99, FILE='nevyaz.txt')
nevyaz = 1
k = 0
do while (nevyaz > eps)
    do i = 1,N*N
        ZNext(i) = 0
        do j = 1,N*N
            if (i .eq. j) then
                ZNext(i) = ZNext(i) + (A(i,j) + alpha)*ZPrev(j)
            else
                ZNext(i) = ZNext(i) + A(i,j) * ZPrev(j)
            endif
        enddo
        ZNext(i) = ZPrev(i) - (ZNext(i) - F(i)) / (maxL + 0.5)
    enddo

    nevyaz = 0
    do i = 1,N*N
        az = 0
        do j = 1,N*N
            az = az + A(i,j) * ZPrev(j)
        enddo
        az = az - F(i)
        nevyaz = nevyaz + az*az
    enddo

    nevyaz = sqrt(nevyaz) / normF
    
    write(*,*) k, nevyaz
    
    do i = 1, N*N
        ZPrev(i) = ZNext(i)
    enddo
    k = k + 1
enddo
write(*,*) "calculated z"

open(11, FILE='z.txt')
do k = 1,N*N
    write(11,*) X(k), Y(k), ZNext(k)
enddo



! open(40, FILE="f_100x100.dat")
! X(l), Y(k)
! do k = 1,N
!     do l = 1,N
!     F = 0
! ! ite
!     ! X(j), Y(i)
!     do i = 1,N
!         do j = 1,N
!         F = F + ( &
!             1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H1((i-1) * N + j) ** 2 ) ** 0.5 - &
!             1/( (X(j) - X(l)) ** 2 + (Y(i) - Y(k)) ** 2 + H2((i-1) * N + j) ** 2 ) ** 0.5 &
!             ) * DEN((i-1) * N + j)
!         enddo
!     enddo
!     F = F * G * dx * dy
!     WRITE(40, 8) X(l), Y(k), F
!     enddo
! enddo
! 8 FORMAT(F25.20, ' ', F25.20, ' ', F25.20)

write(*,*) "Hello, world"

end
