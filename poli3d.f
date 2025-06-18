      program poli
        implicit none
        integer, parameter :: Nmax = 1000000
        real(8), dimension(Nmax) :: obs
        real(8) :: obsmean, tau,O2,sigma
        integer :: i,j,ib,i1,j1,iO,nO

        do i = 440,500,2
          obsmean = 0.0d0
          nO = 0
          do j = 1, Nmax
            read(i,*,end=2) obs(j)
            obsmean = obsmean + obs(j)
            nO = nO + 1
          enddo
2         continue   
          obsmean = obsmean/nO

          ib = 2
          do while(ib .le. nO/2)
            sigma = 0.0d0
            iO = 1
            do i1 = 1, nO-ib+1,ib
              O2 = 0.0d0
              do j1 = i1, i1+ib-1
                O2 = O2 + obs(j1)
              enddo
              O2 = O2/ib
              sigma = sigma + (O2 - obsmean)**2
              iO = iO + 1
            enddo
            sigma = sigma/(iO-1)
            write(1000+i,*)ib, dsqrt(sigma*ib/nO)
            ib=ib+1
          enddo
        enddo

      end program
