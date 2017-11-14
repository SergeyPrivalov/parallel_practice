program main

! константы
real cf, eps
parameter (cf = 0.00667, eps = 0.001)
parameter (N=50)
parameter (dx = 204, dy = 204)

! входные данные
real A(N*N, N*N), H1(N, N), H2(N, N), B(N, N), X(N, N), Y(N, N)

! вспомогательные матрицы
real ATA(N*N, N*N), ATB(N*N)

! для счета maxL 
real XPrev(N*N), XNext(N*N), YCur(N*N), norm, maxL

! для метода
real ZNext(N*N), ZPrev(N*N), alpha, az, nevyaz, normB
! коэфициент регуляризации
alpha = 0.001 

! input
open(4, FILE="data.txt")
do i = 1,N
    do j = 1,N
        READ(4,*) X(i,j), Y(i,j), H1(i,j), H2(i,j), B(i,j)
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

!  A^T * A 
do i = 1,N*N
    do j = 1, N*N
        ATA(i,j) = 0
        do k = 1, N*N
            ATA(i,j) = ATA(i,j) + A(k,i)*A(k,j)
        enddo
    enddo
enddo

!  A^T*b
do k = 1,N
    do l = 1,N
        ATB((l-1)*N+k) = 0
        do i = 1,N
            do j = 1,N
                ATB((l-1)*N+k) = ATB((l-1)*N+k) + A((i-1)*N+j,(l-1)*N+k)*B(i, j)
            enddo
        enddo
    enddo
enddo


!  norm b
normB = 0
do k = 1,N
    do l = 1,N
        normB = normB + B(k,l)**2
    enddo
enddo
normB = sqrt(normB)

! maxL иницилизация
norm = 1
do i = 1,N*N
    if (i .eq. 1) then
        XPrev(i) = 1
    else
        XPrev(i) = 0
    endif
enddo 

! maxL счет
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
    write(6,*) k, maxL, norm
    k = k + 1
enddo


! решение, инициализация
do i = 1,N*N
    ZPrev(i) = 0
enddo 

!
!      считаем z
!
nevyaz = 1
do while (nevyaz > eps)
    do i = 1, N*N
        ZNext(i) = 0
        do j = 1, N*N
            if (i .eq. j) then
                ZNext(i) = ZNext(i) + (ATA(i,j) + alpha)*ZPrev(j)
            else
                ZNext(i) = ZNext(i) + ATA(i,j)*ZPrev(j)
            endif
        enddo
        ZNext(i) = ZPrev(i) - (ZNext(i) - ATB(i)) / maxL
    enddo

    nevaz = 0
    do i = N*N
        az = 0
        do j = N*N
            az = az + A(i,j) * ZNext(j)
        enddo
        az = az - B(i / N , MOD(i, N))
        nevyaz = nevyaz + az*az
    enddo

    nevyaz = sqrt(nevyaz) / normB



    do i = 1, N*N
        ZPrev(i) = ZNext(i)
    enddo
enddo


! вывод матрицы A.txt
!20 open(5, FILE='A.txt')
!do i = 1,10
!    write(5, 10) (A(i, j), j = 1,10)
!    10     FORMAT(10F10.7, ' ')
!enddo

end

