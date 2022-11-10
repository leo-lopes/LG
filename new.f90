program lg_m
    implicit none
    integer, parameter :: L=10,npassos=100*1e5
    real(8):: Tempi,se2,varE
    integer,dimension(L+2,L+2) :: k
    real(8) ::E,dE,Emed,tmed,drandom
    integer :: i,j,np,n,ci,cont,cont2,cont3,i3
    character(len=5) :: pos
    character(len=3) :: temp
    integer :: ii,jj,iv,jv
    real(8) :: p,n1,n2,z,beta,zif
    
    call random_seed()
    do i3=1,22
    	write(*,*)i3
    	write(temp,'(I3.3)') i3
    	open(2,file='results.dat',status='unknown')
    	open(1,file='energy'//temp//'.dat',status='unknown')
    	Tempi = 0.40 + float(i3)*0.025
        E = 0.0d0
        tmed=0.8
        se2=0.0
        np = int(L*L/2)
        cont2=0
        n = 0
        cont3=0
        k(:,:) = 0

        ! Definindo a condição inicial aleatoriamente
        
        !do while (n <= np)
        !    call random_number(drandom)
        !    rx = int(L*drandom) + 2
        !    call random_number(drandom)
        !    ry = int(L*drandom) + 2

         !   if (k(rx,ry) == int(0)) then
         !       k(rx,ry) = int(1)
         !       n = n + 1
         !   endif                 
        !enddo
        do i=2,L+1
           do j=2,L+1
           if (mod(i+j,int(2))==int(0)) then 
           k(i,j)=int(1)
           else
           k(i,j)=int(0)
           end if 
           end do 
        end do 
        ! Criando condições de contorno
        k(1,:) = k(L+1,:)
        k(L+2,:) = k(2,:)
        k(:,1) = k(:,L+1)
        k(:,L+2) = k(:,2)  

        ! Calculando a energia do estado inicial
        do i = 2,L+1
            do j = 2,L+1
                if (k(i,j) == 1) then
                    if (k(i,j) == k(i+1,j)) then
                        E = E - 1
                    endif
                    if (k(i,j) == k(i-1,j)) then
                        E = E - 1
                    endif
                    if (k(i,j) == k(i,j+1)) then
                        E = E - 1
                    endif
                    if (k(i,j) == k(i,j-1)) then
                        E = E - 1
                    endif
                endif
            enddo
        enddo
        E = E/2.0d0
        Emed=0.0
        beta = 1/(Tempi)
        cont=99
        do ci=1,npassos
        !! Dada a condição inicial, vamos escolher um sítio aleatóriamente
        call random_number(drandom)
        ii = int(L*drandom) + 2 !! ii são os sítios em x e iv seu vizinho
        call random_number(drandom)
        jj = int(L*drandom) + 2 !! jj são os sítios em y e jv seu vizinho
	    k(1,:) = k(L+1,:)
        k(L+2,:) = k(2,:)
        k(:,1) = k(:,L+1)
        k(:,L+2) = k(:,2) 
        ! Agora, vamos escolher os vizinhos aleatoriamente
        call random_number(drandom)
        p = 4*drandom

        ! Definindo os vizinhos
        if (0 <= p .and. p < 1) then
            if (ii == 2) then ! Condição de fronteira
                iv = L + 1
            else
                iv = ii - 1
            endif
            jv = jj
        elseif (1 <= p .and. p < 2) then
            if (ii == L+1) then ! Condição de fronteira
                iv = 2
            else
                iv = ii + 1
            endif
            jv = jj
        elseif (2 <= p .and. p < 3) then
            if (jj == 2) then ! Condição de fronteira
                jv = L + 1
            else
                jv = jj - 1
            endif
            iv = ii
        else 
            if (jj == L+1) then ! Condição de fronteira
                jv = 2
            else 
                jv = jj + 1
            
            endif
            iv = ii
        endif
        if (k(ii,jj) .ne. k(iv,jv)) then
            if (k(ii,jj) == int(1)) then
                n1 = float(k(ii+1,jj))
                n1 = n1 + float(k(ii-1,jj))
                n1 = n1 + float(k(ii,jj+1))
                n1 = n1 + float(k(ii,jj-1))

                n2 = float(k(iv+1,jv))
                n2 = n2 + float(k(iv-1,jv))
                n2 = n2 + float(k(iv,jv+1))
                n2 = n2 + float(k(iv,jv-1))
		dE = 1. + n1 - n2
		call random_number(drandom)
                z = drandom
                zif = exp(-1.0*beta*dE)
                if (z .le. zif) then
                k(ii,jj) = int(0)
                k(iv,jv) = int(1)
                E = E + dE
                end if 
            else 
                n1 = float(k(iv+1,jv))
                n1 = n1 + float(k(iv-1,jv))
                n1 = n1 + float(k(iv,jv+1))
                n1 = n1 + float(k(iv,jv-1))

                n2 = float(k(ii+1,jj))
                n2 = n2 + float(k(ii-1,jj))
                n2 = n2 + float(k(ii,jj+1))
                n2 = n2 + float(k(ii,jj-1))
                dE = 1. + n1 - n2
                call random_number(drandom)
                z = drandom
                zif = exp(-1.0*beta*dE)
                if (z .le. zif) then
                k(ii,jj) = int(1)
                k(iv,jv) = int(0)
                E = E + dE
                end if 
            endif
            
            ! write(*,*) dE 

            ! Determinando as condições de troca com os vizinhos
            !if (dE > 0) then
            !    call random_number(drandom)
            !    z = drandom
            !    zif = exp(-1.0*beta*dE)
            !    if (z <= zif) then
            !        ! write(*,*) 'Mudou'
            !        if (k(ii,jj) == int(1)) then
            !            k(ii,jj) = int(0)
            !            k(iv,jv) = int(1)
            !        else
            !            k(ii,jj) = int(1)
            !            k(iv,jv) = int(0)
            !        E = E + dE
            !        endif
             !   endif
            !else
                ! write(*,*) 'Mudou'
             !   if (k(ii,jj) == int(1)) then
             !       k(ii,jj) = int(0)
             !       k(iv,jv) = int(1)
             !   else
             !       k(ii,jj) = int(1)
             !       k(iv,jv) = int(0)
             !   E = E + dE
             !   endif 
            !endif
       endif
        cont=cont+1
        cont2 = cont2 + 1
        if (cont2 == int(1e9)) then 
            cont3=cont3+1
            write(pos,'(I5.5)') cont3
            open(15,file='pos'//temp//'_'//pos//'_''.dat',status='unknown')
            do i=2,L+1
                do j=2,L+1   
                    write(15,*) i,j,k(i,j)
                end do 
            end do
            close(15)
            cont2=0
        end if 
        if (cont>int(tmed*float(npassos))) then
        !if (cont==100) then
        	!write(1,*)ci,E
           Emed=Emed+E/((1.0-tmed)*npassos)
           se2 = se2 + (E**2.)/((1.0-tmed)*npassos)
        end if 
    end do 
    varE = se2-Emed**2.
    write(*,*)Tempi,Emed,varE
    write(2,*)Tempi,Emed,varE
    close(1)
    end do 
    close(2)
    end program
