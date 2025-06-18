program su2
	implicit none
	integer, parameter :: Ls = 9
  integer, parameter :: Lt = 4
	integer, parameter :: iseed1 = 10, Nth = 10000, estat = 10000, Nc = 2, tau = 200, dat = Nth/tau
	real(8), parameter :: pi = 4.0d0*atan(1.0d0)
	
	integer :: j,k,k1,k2,mu,n1,n2,n3,n4,nu,cont,norm,ac,total,bin,datum,ibetas
	real(8), dimension(Ls,Ls) :: loop, loopmean, loopli
	real(8), dimension(4) :: x, identidade, ax1,ax1d
	real(8), dimension(4,4,Ls,Ls,Lt) :: U, Uaux
	real(8) :: t1,t2,modr1,det1,t3,T,betas, beta, delta
	real(8) :: eoq,dummy,mod2,Vol,obs,obsmean,obssigma,omf
  real(8), dimension(dat) :: obsval
  real(8), dimension(0:Ls-1) :: pk
	integer, dimension(4,Ls) :: nn
	real(8), dimension(Ls) :: par
	real(8), dimension(5) :: r1
  integer, dimension(3) :: pos
  complex(8) :: z, id

  delta = 0.02d0
  do while (betas .lt. 3.0d0)

  beta = 1.98d0
  betas = beta + delta
  ibetas = 200*betas
    print*, "O valor de betas é", betas
    Vol = Lt*Ls**3
    norm = 12*Vol
    open(11,file='wilsontherm.dat')
    
    
    T = 1.0d0/Lt

    call srand(iseed1)

    do j = 1,Ls
      nn(1,j) = 1 + mod(j,Ls)
      nn(2,j) = Ls - mod(Ls-j+1,Ls)
    enddo
    
    do j = 1,Lt
      nn(3,j) = 1 + mod(j,Lt)
      nn(4,j) = Lt - mod(Lt-j+1,Lt)
    enddo
    
    do n1 = 1, Ls !Configuração inicial com identidade nas matrizes
      do n2 = 1, Ls
          do n4 = 1, Lt
            do mu = 1,4
              U(:,mu,n1,n2,n4) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
            enddo
          enddo
      enddo
    enddo
    
    dummy = 0.0d0
    cont = 0
    ac = 0
    total = 0
    
    call cpu_time(t1)
    do j = 1, Nth
      dummy = 0.0d0
      call calclooppolyakov(U,pos,nn,dummy)
!     call calcloopntnl(U,1,1,nn,dummy)
      call banhotermico(U,betas,nn)
 !    write(11,*)j,dummy/norm
    enddo
    call cpu_time(t2)
    
  !	print*, ac, total, real(ac)/total
    print*, "O tempo gasto pra termalizar com", Nth, "Iterações é", t2-t1
    
!    do j = 1, estat-tau+1, tau
!      loopmean = 0.0d0
!   	  par = 0.0d0
!   	  do bin = j, j + tau - 1
!   	    loop = 0.0d0
!   	    call banhotermico(U,betas,nn)
!   	    call metropolis(U,L,nnp,nnm,beta,ac,total)	  
!   	    do k1 = 1,Ls-2
!   	      do k2 = 1,Ls-2
!   	        call calcloopntnl(U,k1,k2,nn,loop(k1,k2))
!   	      enddo
!   	    enddo
!   	    loopmean = loopmean + loop
!      enddo
!   	  loopmean = loopmean/(tau*norm)
!   	  par = par + real(loopmean(1,:))
!   	  do k1 = 1,1
!   	    do k2 = 1,1
!   	     write(10+k1+(k2-1)*10,*)j,real(loopmean(k1,k2))
!   	   enddo
!   	  enddo
!   	  loopmean = 0.0d0
!   	enddo
!  
!    do j = 1, estat
!      call banhotermico(U, betas, nn)
!      loop = 0.0d0
!      do k1 = 1,Ls-2
!        do k2 = 1, Ls-2
!          call calcloopntnl(U,k1,k2,nn,loop(k1,k2))
!        enddo
!      enddo
!      do k1 = 1, Ls-2
!        do k2 = 1,Ls-2
!          write(10 + k1*10 + k2*1000,*)loop(k1,k2)/norm
!        enddo
!      enddo
!    enddo
!   	call cpu_time(t3)
!   	
!   	print*, "O tempo dos dados estatísticos é", t3-t2
   	
    loop(1,1) = 0.0d0
    call calclooppolyakov(U,pos,nn,loop(1,1))
!   call calcloopntnl(U,1,1,nn,loop(1,1))
    loop(1,1) = loop(1,1)/norm
    
    print*, loop(1,1)
    
    call transform_gauge(U,nn)
    
    loop(1,1) = 0.0d0
    call calclooppolyakov(U,pos,nn,loop(1,1))
!   call calcloopntnl(U,1,1,nn,loop(1,1))
    loop(1,1) = loop(1,1)/norm
    
    print*, loop(1,1)

  ! call losalamos(U,Uaux,nn)
  !    do j = 1, Nth-tau+1, tau
  !      do bin = j, j+tau-1
  !        call banhotermico(U,betas,nn)
  !      enddo
  !      call losalamos(U,Uaux,nn)
  !      call gluonprop(U,nn,pk)
  !      do k = 0, Ls-1
  !        write(k+100,*)k, pk(k)
  !      enddo
  !    enddo
    !
    
    do j = 1, estat
      call banhotermico(U,betas,nn)
      call calclooppolyakov(U,pos,nn,obs)
      write(ibetas,*)obs
    enddo
    delta = delta + 0.1d0
  enddo

contains

subroutine linkmul(a,b,c)
	implicit none
	
	real(8), dimension(4) :: a, b, c
	real(8) :: a1,a2,a3,a4
	
	a4 = -a(1)*b(1) - a(2)*b(2) - a(3)*b(3) + a(4)*b(4)
	a1 = a(1)*b(4) + a(4)*b(1) - a(2)*b(3) + a(3)*b(2)
	a2 = a(2)*b(4) + a(4)*b(2) - a(3)*b(1) + a(1)*b(3)
	a3 = a(3)*b(4) + a(4)*b(3) - a(1)*b(2) + a(2)*b(1)
	
	c = (/a1,a2,a3,a4/)
end subroutine

subroutine dagger(a,c)
	implicit none
	
	real(8), dimension(4) :: a, c
	real(8) :: a1,a2,a3,a4
	
	a1 = -a(1)
	a2 = -a(2)
	a3 = -a(3)
	a4 = a(4)
	
	c = (/a1,a2,a3,a4/)
end subroutine

function det(a)
	implicit none
	
	real(8), dimension(4) :: a, c
	real(8) :: a1,a2,a3,a4
	real(8) :: det
	
	det = a(1)**2 + a(2)**2 + a(3)**2 + a(4)**2
	
	return
end function

subroutine  metropolis(U,nn,beta,ac,total)
	implicit none
	
	real(8), parameter :: eps = 0.1d0
	
	real(8), dimension(4,4,Ls,Ls,Lt) :: U	
	real(8), dimension(4) :: Um, dUm, Am
	integer, dimension(4,Ls) :: nn
	integer, dimension(3) :: pos, n	
	integer :: ac,total,mu,L
	real(8) :: deltaS
	real(8) :: beta
	
	do n1 = 1,L
	  do n2 = 1,L
	      do n4 = 1,L
		      pos = (/n1,n2,n4/)			
	        do mu = 1,4
            r1(4) = 1.0d0 - 2.0d0*rand()	
            r1(1) = 1.0d0 - 2.0d0*rand()
            r1(2) = 1.0d0 - 2.0d0*rand()
            r1(3) = 1.0d0 - 2.0d0*rand()
          
            x(4) = dsign(dsqrt(1.0d0-eps**2),r1(4))
            modr1 = dsqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)				
            x(1:3) = eps*r1(1:3)/modr1
            Um = (/x(1),x(2),x(3),x(4)/)
          
            call linkmul(Um,U(:,mu,n1,n2,n4),Um)
            call grampos(U,Am,mu,pos,nn)
          
            dUm = Um - U(:,mu,n1,n2,n4)
            call linkmul(dUm,Am,dUm)
            deltaS = (-1.0d0*beta/2.0d0)*(2.0d0*dUm(1))
            r1(4) = rand()
          
            if(r1(5) .le. exp(-1.0d0*deltaS)) then
              U(:,mu,n1,n2,n4) = Um/dsqrt(det(U(:,mu,n1,n2,n4)))
              ac = ac + 1
            endif
            total = total + 1
	        enddo
	      enddo
	  enddo
	enddo
end subroutine

subroutine calclooppolyakov(U,pos,nn,obs)
	implicit none
	
	real(8), dimension(4,4,Ls,Ls,Lt) :: U	
	real(8), dimension(4) :: loop, Ud
	integer, dimension(3) :: pos, n
	integer, dimension(4,Ls) :: nn
	real(8) :: obs
	integer :: k
   
  obs = 0.0d0
	do n1 = 1,Ls
	  do n2 = 1,Ls
		      pos = (/n1,n2,1/)
		      loop = (/0.0d0,0.0d0,0.0d0,1.0d0/)
		      n  = pos
		      do k = 1, Lt
		        call linkmul(loop,U(:,4,n(1),n(2),n(3)),loop)
		        n(3) = nn(3,n(3))
		      enddo
		      obs = obs + loop(4)
	  enddo
	enddo
  obs = obs/Lt**3
end subroutine

subroutine calcloopntnl(U,nl,nt,nn,obs)
	implicit none
	
	real(8), dimension(4,4,Ls,Ls,Lt) :: U	
	real(8), dimension(4) :: loop, Ud
	integer, dimension(3) :: pos, n	
	integer, dimension(4,Ls) :: nn
  integer :: dmu, dnu, mu, nu
	integer :: cont,nl,nt,k
	real(8) :: obs
	
	do n1 = 1,Ls
	  do n2 = 1,Ls
	      do n4 = 1,Lt
		      pos = (/n1,n2,n4/)
	        do mu = 1, 3	
      dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
            do nu = 3, mu+1,-1
      dnu = (nu-1)*(nu-2)*3/2-(nu-1)*(nu-3)+(nu-3)*(nu-2)/2
              loop = (/0.0d0,0.0d0,0.0d0,1.0d0/)
              n  = pos
              cont = cont + 1
              do k = 1, nl
                call linkmul(loop,U(:,mu,n(1),n(2),n(3)),loop)
                n(mu) = nn(dmu,n(mu))
              enddo
              do k = 1, nt
                call linkmul(loop,U(:,nu,n(1),n(2),n(3)),loop)
                n(nu) = nn(dnu,n(nu))
              enddo
              do k = 1, nl
                n(mu) = nn(dmu+1,n(mu))
                call dagger(U(:,mu,n(1),n(2),n(3)),Ud)
                call linkmul(loop,Ud,loop)
              enddo
              do k = 1, nt
                n(nu) = nn(dnu+1, n(nu))
                call dagger(U(:,nu,n(1),n(2),n(3)),Ud)
                call linkmul(loop,Ud,loop)
              enddo
              obs = obs + 2.0d0*loop(4)
            enddo
	        enddo
	      enddo
	  enddo
	enddo
end subroutine

subroutine grampos(U,Am,mu,pos,nn)
	implicit none
	
	real(8), dimension(4) :: Um, dUm, Am, Ummu, Umnu, Unu, Am1
	real(8), dimension(4,4,Ls,Ls,Lt) :: U	
	integer, dimension(3) :: pos, n, np
	integer, dimension(4,Ls) :: nn
	integer :: mu, nu, dmu, dnu

  dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
	
	Am = 0.0d0
	np = pos
	np(mu) = nn(dmu,np(mu))
	do nu = 1, 3
	  if (nu .ne. mu) then
  dnu = (nu-1)*(nu-2)*3/2-(nu-1)*(nu-3)+(nu-3)*(nu-2)/2
	    n = np
	    Unu = U(:,nu,n(1),n(2),n(3))
	    n(nu) = nn(dnu,n(nu))
	    n(mu) = nn(dmu+1,n(mu))
	    call dagger(U(:,mu,n(1),n(2),n(3)),Ummu)
	    n(nu) = nn(dnu+1,n(nu))
	    call dagger(U(:,nu,n(1),n(2),n(3)),Umnu)
	    call linkmul(Ummu,Umnu,Am1)
	    call linkmul(Unu,Am1,Am1)
	    Am = Am + Am1	    	  
	    	    
	    n(mu) = nn(dmu,n(mu))
	    n(nu) = nn(dnu+1,n(nu))
	    call dagger(U(:,nu,n(1),n(2),n(3)),Umnu)
	    n(mu) = nn(dmu+1,n(mu))
	    call dagger(U(:,mu,n(1),n(2),n(3)),Ummu)
	    Unu = U(:,nu,n(1),n(2),n(3))
	    call linkmul(Ummu,Unu,Am1)
	    call linkmul(Umnu,Am1,Am1)
	    Am  = Am + Am1
	  endif
	enddo
end subroutine

subroutine transform_gauge(U,nn)
	implicit none
	
	real(8), dimension(4) :: Um, dUm, ident,gd
	real(8), dimension(4,4,Ls,Ls,Lt) :: U	
	real(8), dimension(4,Ls,Ls,Lt) :: g
	integer, dimension(3) :: pos, n
	integer, dimension(4,Ls) :: nn
  integer :: dmu
	
	do n1 = 1,Ls
	  do n2 = 1,Ls
	      do n4 = 1,Lt
	        call random_number(x)
	        x = 1.0d0 - 2.0d0*x
	        x = x/dsqrt(x(4)**2 + x(1)**2 + x(2)**2 + x(3)**2)
	        g(:,n1,n2,n4) = (/x(1),x(2),x(3),x(4)/)
	      enddo
	  enddo
	enddo
	
	do n1 = 1,Ls
	  do n2 = 1,Ls
	      do n4 = 1,Lt
	        pos = (/n1,n2,n4/)
	        do mu = 1,3
  dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
	          n = pos
	          call linkmul(g(:,n1,n2,n4),U(:,mu,n1,n2,n4),U(:,mu,n1,n2,n4))
	          n(mu) = nn(dmu,n(mu))
	          call dagger(g(:,n(1),n(2),n(3)),gd)
		        call linkmul(U(:,mu,n1,n2,n4),gd,U(:,mu,n1,n2,n4))
	        enddo
	      enddo
	  enddo
	enddo
end subroutine

subroutine banhotermico(U,beta,nn)
	implicit none
	
	real(8) :: beta,alfa,lambda2,mod2,a1,a2,a3,a4
	real(8), dimension(4,4,Ls,Ls,Lt) :: U
	real(8), dimension(4) :: Am, ger
	integer, dimension(3) :: pos, n
	integer, dimension(4,Ls) :: nn
	real(8), dimension(1:4) :: val
	integer :: dmu
	
	do n1 = 1,Ls
	  do n2 = 1,Ls
	      do n4 = 1,Lt
	        pos = (/n1,n2,n4/)
	        do mu = 1,3
	          call grampos(U,Am,mu,pos,nn)
	          alfa = sqrt(det(Am))
	          do k = 1, 100
	            call random_number(x)
	            x(4) = 1.0d0 - x(4)
	            a1 = dlog(x(1))
	            a2 = (dcos( 2.0d0*pi*x(2) ))**2
	            a3 = dlog(x(3))
	            a4 = -1.0d0/(2.0d0*alfa*beta)
	            lambda2 = (a1 + a2*a3)*a4
	            if(x(4)**2 .le. (1.0d0 - lambda2)) then
	              val(4) = 1.0d0-2.0d0*lambda2
	              exit
	            endif
	          enddo
	          
	          do k = 1, 100
	            call random_number(x(1:3))
	            x(1:3) = 1.0d0 - 2.0d0*x(1:3)
	            mod2 = dsqrt(x(1)**2 + x(2)**2 + x(3)**2)
	            if(mod2**2 .le. 1.0d0) then
	              val(1:3) = x(1:3)*dsqrt(1.0d0-val(4)**2)/mod2
	              exit
	            endif
	          enddo
		        ger = (/val(1),val(2),val(3),val(4)/)
		        call dagger(Am,Am)
	          call linkmul(ger,Am,U(:,mu,n1,n2,n4))
	          U(:,mu,n1,n2,n4) = U(:,mu,n1,n2,n4)/alfa
	        enddo
	      enddo
	  enddo
	enddo
end subroutine

subroutine losalamos(U,Uaux,nn)
  implicit none

  real(8), dimension(4,4,Ls,Ls,Lt) :: U, Uaux
	real(8), dimension(4) :: g,gd,aux,aux1
	real(8), dimension(3) :: dUx, dUxm, Ufi
	integer, dimension(3) :: pos, n
	integer, dimension(4,Ls) :: nn
	real(8) :: gauge,func,lambda
	integer(8) :: contagem,dmu
	
	Uaux = U
	contagem = 1
	gauge = 1.0d0
		
	do while(gauge .ge. 1.d-12 .and. contagem .le. 1000)
     gauge = 0.0d0
     func = 0.0d0
     do n1 = 1,Ls !Cálculo do gauge usando a Uaux
	    do n2 = 1,Ls
	        do n4 = 1,Lt
	          pos = (/n1,n2,n4/)
	          Ufi = 0.0d0
	          do mu = 1,3
  dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
	            n = pos
	            dUx = Uaux(1:3,mu,n1,n2,n4)
	            n(mu) = nn(dmu+1,n(mu))
	            dUxm = Uaux(1:3,mu,n(1),n(2),n(3))
	            Ufi = Ufi + (dUx - dUxm)
	            func = func - 2.0d0*Uaux(4,mu,n1,n2,n4)
	          enddo
	          gauge = gauge + (Ufi(1)**2 + Ufi(2)**2 + Ufi(3)**2)
	        enddo
	    enddo
	  enddo
	  func = func/(8*Vol)
	  func = func + 1.0d0
	  
	  gauge = gauge/(8*Vol)
	  write(2,*)contagem, gauge
	  write(3,*)contagem, func
	  contagem = contagem + 1
	  
	  do n1 = 1,Ls !Definição da matriz h, g e transformação da matriz auxiliar
	    do n2 = 1,Ls
	        do n4 = 1,Lt
	          pos = (/n1,n2,n4/)
	          g = 0.0d0
	          do mu = 1,3
	            n = pos
  dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
	            call dagger(Uaux(:,mu,n1,n2,n4),aux)
	            n(mu) = nn(dmu+1,n(mu))
	            aux1 = Uaux(:,mu,n(1),n(2),n(3))
	            g = g + aux + aux1
	          enddo
	          g = g/dsqrt(det(g))
	          do mu = 1,3
	            n = pos
  dmu = (mu-1)*(mu-2)*3/2-(mu-1)*(mu-3)+(mu-3)*(mu-2)/2
	            call dagger(g,gd)
	            call linkmul(g,Uaux(:,mu,n1,n2,n4),Uaux(:,mu,n1,n2,n4))
	            n(mu) = nn(dmu+1,n(mu))
	  	      call linkmul(Uaux(:,mu,n(1),n(2),n(3)),gd,Uaux(:,mu,n(1),n(2),n(3)))	            
	          enddo
	        enddo
	    enddo
	  enddo
	enddo
	U = Uaux
end subroutine


subroutine gluonprop(U,nn,pk)
	implicit none
	
	real(8), dimension(3) :: Ux,aux1,aux2,aux3
	real(8), dimension(4,4,Ls,Ls,Lt) :: U
	real(8), dimension(4) :: g,gd,aux
	integer, dimension(3) :: pos, n
  real(8), dimension(0:Ls-1) :: pk
	integer, dimension(4,Ls) :: nn
	real(8) :: gauge,ang
  integer :: kl
  
  do kl = 0,Ls/2
    do mu = 1,3
      aux1 = 0.0d0
      aux2 = 0.0d0
      aux3 = 0.0d0
      do n1 = 1,Ls
        do n2 = 1,Ls
            do n4 = 1,Lt
              ang = 2*pi*kl*n4/Ls
              Ux = U(1:3,mu,n1,n2,n4)
              aux2 = aux2 + Ux*sin(ang)
              aux3 = aux3 + Ux*cos(ang)
            enddo
          enddo
      enddo
      pk(kl) = pk(kl) + aux2(1)**2 + aux2(2)**2 + aux2(3)**2
      pk(kl) = pk(kl) + aux3(1)**2 + aux3(2)**2 + aux3(3)**2
    enddo
  enddo
  pk(1:Ls-1) = pk(1:Ls-1)/(6*Vol)
  pk(0) = pk(0)/(9*Vol)
end subroutine

end program su2
