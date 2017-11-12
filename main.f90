program main

real cf, eps
parameter (cf = 0.00667, eps = 0.001)
parameter (N=50)
parameter (dx = 204, dy = 204)
real A(N*N, N*N), ATA(N*N, N*N), ATAB(N*N), H1(N, N), H2(N, N), F(N, N), X(N, N), Y(N, N), ZNext(N*N), ZPrev(N*N)
real XPrev(N*N), XNext(N*N), YCur(N*N)
real norm, maxL

real alpha
alpha = 0.001

open(4, FILE="data.txt")

do i = 1,N
    do j = 1,N
        READ(4,*) X(i,j), Y(i,j), H1(i,j), H2(i,j), F(i,j)
    enddo
enddo

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


!
!  A^T * A + alpha * E
!

do i = 1,N*N
    do j = 1, N*N
        ATA(i,j) = 0
        do k = 1, N*N
            ATA(i,j) = ATA(i,j) + A(k,i)*A(k,j)
        enddo
        if (i .eq. j) then
            ATA(i,j) = ATA(i,j) + alpha
        endif
    enddo
enddo

!
!  A^T*b
!

do i = 1,N*N
    ATAB(i) = 0
    do j = 1,N*N
        ATAB(i) = ATAB(i) + A(j,i)*F(j / 50, j % 50)
    enddo
enddo

! собств значения

norm = 1
do i = 1,N*N
    if (i .eq. 1) then
        XPrev(i) = 1
    else
        XPrev(i) = 0
    endif
    ZPrev(i) = 0
enddo 

open(6, FILE='maxL.txt')
k = 0
do while (norm > eps)
    norm = 0
    do i = 1, N*N
        YCur(i) = 0
        do j =1, N*N
            YCur(i) = YCur(i) + ATA(i,j) * XPrev(j)
        enddo
        norm = norm + YCur(i)**2
    enddo
    norm = sqrt(norm)
    maxL = norm
    
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
    write(*,*) k, maxL, norm
    k = k + 1
enddo

WRITE(*, *) maxL


!
!      считаем z
!

do while ()
    do i = 1, N*N
        ZNext(i) = 0
        do j = 1, N*N
            ZNext(i) = ZNext(i) + ATA(i,j)*ZPrev(j) - ATAB(j)
        enddo
        ZNext(i) = ZPrev(i) - ZNext(i) / maxL
        ZPrev(i) = ZNext(i)
    enddo
enddo


open(5, FILE='A.txt')
do i = 1,10
    write(5, 10) (A(i, j), j = 1,10)
    10     FORMAT(10F15.12, ' ')
enddo

end

