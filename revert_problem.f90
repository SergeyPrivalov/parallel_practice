program primal
implicit none

integer i,j,k,l,N,maxX,minX

parameter (N=50)
double precision H1(N*N), H2(N*N), F(N*N), X(N), Y(N), aa, bb, G, dx, dy, A(N*N, N*N), ATA(N*N, N*N), ATB(N*N)
double precision XPrev(N*N), XNext(N*N), norm, eps, YCur(N*N), maxL, ZNext(N*N), ZPrev(N*N), normATB
double precision nevyaz, az, alpha, step
parameter (G = 0.00667408, eps = 0.001, alpha = 0.0001)

minX = 0
maxX = 10000
step = 204

do i = 1, N
    X(i) = minX + step * (i-1)
    Y(i) = minX + step * (i-1)
enddo

dx = ABS(X(N) - X(1)) / (N - 1)
dy = ABS(Y(N) - Y(1)) / (N - 1)
write (*, *) dx, dy

! open(10, FILE="data/H1_100x100.dat")
do i = 1,N*N
    ! READ(10, *) aa, bb, H1(i)
    H1(i) = 1000
enddo

! open(20, FILE="data/H2_100x100.dat")
do i = 1,N*N
    ! READ(20, *) aa, bb, H2(i)
    H2(i) = 1500
enddo

! open(30, FILE="data/field2_100x100.dat")
open(30, FILE="f_1sq.dat")
do i = 1,N*N
    READ(30, *) aa, bb, F(i)
enddo
WRITE(*,*) "Data read"


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
! nevyaz = sqrt(nevyaz) / normATB
! write(*,*) nevyaz
!
!      считаем z
!


write(*,*) "calculating ATA and ATB"
do i = 1,N*N
	ATB(i) = 0
	do j = 1,N*N
		ATA(i,j) = 0
		ATB(i) = ATB(i) + A(j,i)*F(j)
		do k = 1,N*N
			ATA(i,j) = ATA(i,j) + A(i,k)*A(k,j)
		enddo
	enddo
enddo
write(*,*) "calculated ATA and ATB"


!  norm b
normATB = 0
do k = 1,N*N
    normATB = normATB + ATB(k)**2
enddo
normATB = sqrt(normATB)
write(*,*) "calculated normATB", normATB

write(*,*) "end calculated ATA"
open(99, FILE='nevyaz.txt')
nevyaz = 1
k = 0
do while (nevyaz > eps)
    do i = 1,N*N
        ZNext(i) = 0
        do j = 1,N*N
            if (i .eq. j) then
                ZNext(i) = ZNext(i) + (ATA(i,j) + alpha)*ZPrev(j)
            else
                ZNext(i) = ZNext(i) + ATA(i,j) * ZPrev(j)
            endif
        enddo
        ZNext(i) = ZPrev(i) - (ZNext(i) - ATB(i)) / (maxL + 0.5)
    enddo

    nevyaz = 0
    do i = 1,N*N
        az = 0
        do j = 1,N*N
            az = az + ATA(i,j) * ZPrev(j)
        enddo
        az = az - ATB(i)
        nevyaz = nevyaz + az*az
    enddo

    nevyaz = sqrt(nevyaz) / normATB
    
    write(*,*) k, nevyaz
    
    do i = 1, N*N
        ZPrev(i) = ZNext(i)
    enddo
    k = k + 1
enddo
write(*,*) "calculated z"

open(11, FILE='z.txt')
do k = 1,N 
    do l = 1,N
        write(11,*) X(l), Y(k), ZNext((l-1)*N+k)
    enddo
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