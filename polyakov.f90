program wilsonloop
	implicit none
	integer, parameter :: tau = 5000, Ntot = 100000, Npo = Ntot/tau
	integer :: i,j,k1,k2,k3
	real(8) :: ploopmean, sigma, a, b
	real(8), dimension(Npo) :: ploop
	
  open(10,file='data.dat')

	do k1 = 101,172
!     print*, k1
	    k3 = 0
	    ploop = 0.0d0
      k2 = 0
	    ploopmean = 0.0d0
	    do j=1,Ntot-tau+1,tau
	      k3 = k3 + 1
	      do i=j,j+tau-1
	        read(k1,*)b
          k2 = k2 + 1
!         print*, k2
	        ploop(k3) = ploop(k3) + abs(b)
	      enddo
	      ploop(k3) = (ploop(k3))/tau
	      ploopmean = ploopmean + ploop(k3)
	    enddo
	    ploopmean = ploopmean/Npo
	    sigma = 0.0d0
	    do j=1,Npo
	      sigma = (ploop(j) - ploopmean)**2
	    enddo
	    sigma = dsqrt(sigma/(Npo-1))
	    write(10,*)(k1-100)/real(10)+2.90d0,ploopmean,sigma
	enddo
	
end program
