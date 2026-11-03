C Output from Public domain Ratfor, version 1.05
      subroutine crqdiag(n,p,a1,c1,n1,x,y,c,b,g,m,r,s,d,u,wn,wp,
     *info,tr)
      integer n,p,n1,m,info,nit(3)
      double precision tr(33,m),dsmall,dsum,dabsum
      double precision a1(p,n1),c1(n),x(n,p),y(n),c(n),b(p,m),g(m)
      double precision wn(n1,9),wp(p,p+3),r(p),s(n),d(n),u(n)
      double precision zero,half,one,beta,eps,dh
      parameter( zero = 0.0d0)
      parameter( half = 0.5d0)
      parameter( one = 1.0d0)
      parameter( beta = 0.99995d0)
      parameter( eps = 1.0d-8)
      do23000 k = 2,m 
      tr(25,k-1)=zero
      do 991 i=1,n
      if(d(i).ge.zero)tr(25,k-1)=tr(25,k-1)+one
  991 continue
      dh = -log(one - g(k)) + log(one - g(k-1))
      do23002 i = 1,n 
      u(i) = one
      wn(i,1) = half
      if(d(i) .ge. zero)then
      s(i) = s(i) + dh
      endif
      d(i) = c(i) - s(i)
23002 continue
23003 continue
      call dgemv('T',n,p,one,x,n,d,1,zero,r,1)
      call rqdiag(n1,p,a1,c1,r,d,u,beta,eps,wn,wp,nit,info,
     *tr(1,k-1))
      if(info .ne. 0)then
      goto 23001
      endif
      call dcopy(p,wp,1,b(1,k-1),1)
      call dcopy(n,y,1,d,1)
      call dgemv('N',n,p,one,x,n,b(1,k-1),1,one,d,1)
      do 992 j=26,32
      tr(j,k-1)=zero
  992 continue
      tr(29,k-1)=1.0d300
      do 993 i=1,n
      if(d(i).ge.zero)tr(26,k-1)=tr(26,k-1)+one
      if(d(i).eq.zero)tr(27,k-1)=tr(27,k-1)+one
      dsmall=dabs(d(i))/(one+dabs(y(i)))
      if(dsmall.le.1.0d-10)tr(28,k-1)=tr(28,k-1)+one
      tr(29,k-1)=dmin1(tr(29,k-1),dsmall)
      tr(30,k-1)=tr(30,k-1)+dble(i)*d(i)
      tr(31,k-1)=tr(31,k-1)+dabs(d(i))
      if(dsmall.le.1.0d-8)tr(32,k-1)=tr(32,k-1)+one
  993 continue
23000 continue
23001 continue
      m = k-1
      return
      end
C Output from Public domain Ratfor, version 1.05
      subroutine rqdiag(n,p,a,y,rhs,d,u,beta,eps,wn,wp,nit,info,
     *tr)
      integer n,p,info,nit(3)
      double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,9),wp(p,p+3)
      double precision one,beta,eps,tr(33)
      parameter( one = 1.0d0)
      call lpdiag(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2),wp(1,1),wn(1,
     *3),wn(1,4),wn(1,5),wn(1,6),wp(1,2),wn(1,7),wn(1,8),wn(1,9),wp(1,3
     *),wp(1,4),nit,info,tr)
      return
      end
      subroutine lpdiag(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w,dx,ds,dy,dz,dw,
     *dr,rhs,ada,nit,info,tr)
      integer n,p,pp,i,info,nit(3),maxit
      double precision tr(33)
      double precision a(p,n),c(n),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz,dsdw
      double precision deltap,deltad,beta,eps,mu,gap,g
      double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p
     *)
      double precision dx(n),ds(n),dy(p),dz(n),dw(n),dr(n)
      parameter( zero = 0.0d0)
      parameter( one = 1.0d0)
      parameter( mone = -1.0d0)
      parameter( big = 1.0d+20)
      parameter( maxit = 500)
      tr(20)=0.0d0
      tr(21)=0.0d0
      tr(22)=eps
      tr(23)=dble(maxit)
      nit(1)=0
      nit(2)=0
      nit(3)=n
      pp=p*p
      call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
      do23000 i=1,n
      d(i)=one
23000 continue
23001 continue
      call stepdiag(n,p,a,d,y,ada,info)
      if(info .ne. 0)then
      tr(1)=dble(nit(1))
      tr(2)=dble(nit(2))
      tr(3)=dble(nit(3))
      tr(4)=dble(info)
      tr(5)=0.0d0
      return
      endif
      call dcopy(n,c,1,s,1)
      call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
      do23004 i=1,n
      if(dabs(s(i)).lt.eps)then
      z(i)=dmax1(s(i), zero) + eps
      w(i)=dmax1(-s(i),zero) + eps
      else
      z(i)=dmax1(s(i), zero)
      w(i)=dmax1(-s(i),zero)
      endif
      s(i)=u(i)-x(i)
23004 continue
23005 continue
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
23008 if(gap .gt. eps .and. nit(1).lt.maxit)then
      nit(1)=nit(1)+1
      do23010 i = 1,n
      d(i) = one/(z(i)/x(i) + w(i)/s(i))
      ds(i)=z(i)-w(i)
      dz(i)=d(i)*ds(i)
23010 continue
23011 continue
      call dcopy(p,b,1,dy,1)
      call dgemv('N',p,n,mone,a,p,x,1,one,dy,1)
      call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepdiag(n,p,a,d,dy,ada,info)
      if(info .ne. 0)then
      call lpstate(n,p,a,c,b,x,s,y,z,w,gap,nit,info,1,tr)
      return
      endif
      call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do23014 i=1,n
      dx(i)=d(i)*ds(i)
      ds(i)=-dx(i)
      dz(i)=-z(i)*(dx(i)/x(i) + one)
      dw(i)=-w(i)*(ds(i)/s(i) + one)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23014 continue
23015 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(min(deltap,deltad) .lt. one)then
      nit(2)=nit(2)+1
      mu = ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
      g = mu + deltap*ddot(n,dx,1,z,1)+ deltad*ddot(n,dz,1,x,1) + deltap
     **deltad*ddot(n,dz,1,dx,1)+ deltap*ddot(n,ds,1,w,1)+ deltad*ddot(n,
     *dw,1,s,1) + deltap*deltad*ddot(n,ds,1,dw,1)
      mu = mu * ((g/mu)**3) /dble(2*n)
      do23026 i=1,n
      dr(i)=d(i)*(mu*(1/s(i)-1/x(i))+ dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
23026 continue
23027 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n,one,a,p,dr,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      tr(20)=dble(info)
      if(info.ne.0)tr(21)=tr(21)+1.0d0
      call dgemv('T',p,n,one,a,p,dy,1,zero,u,1)
      deltap=big
      deltad=big
      do23028 i=1,n
      dxdz = dx(i)*dz(i)
      dsdw = ds(i)*dw(i)
      dx(i)= d(i)*(u(i)-z(i)+w(i))-dr(i)
      ds(i)= -dx(i)
      dz(i)= -z(i)+(mu - z(i)*dx(i) - dxdz)/x(i)
      dw(i)= -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23028 continue
23029 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      endif
      call daxpy(n,deltap,dx,1,x,1)
      call daxpy(n,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n,deltad,dz,1,z,1)
      call daxpy(n,deltad,dw,1,w,1)
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
      goto 23008
      endif
23009 continue
      call lpstate(n,p,a,c,b,x,s,y,z,w,gap,nit,info,2,tr)
      call daxpy(n,mone,w,1,z,1)
      call dswap(n,z,1,x,1)
      return
      end
      subroutine stepdiag(n,p,a,d,b,ada,info)
      integer n,p,pp,i,info
      double precision a(p,n),b(p),d(n),ada(p,p),zero
      parameter( zero = 0.0d0)
      pp=p*p
      do23038 j=1,p
      do23040 k=1,p
      ada(j,k)=zero
23040 continue
23041 continue
23038 continue
23039 continue
      do23042 i=1,n
      call dsyr('U',p,d(i),a(1,i),1,ada,p)
23042 continue
23043 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end

C Read-only state inspection before original LP output overwrites x and z.
C Bound sum is one because crqfnb supplies u=1 for every LP.
      subroutine lpstate(n,p,a,c,b,x,s,y,z,w,gap,nit,info,phase,tr)
      use, intrinsic :: ieee_arithmetic
      integer n,p,nit(3),info,phase,i,j
      double precision a(p,n),c(n),b(p),x(n),s(n),y(p),z(n),w(n)
      double precision gap,tr(33),v,den,pr,du,cx,by,sw,bmax,cmax
      tr(1)=dble(nit(1))
      tr(2)=dble(nit(2))
      tr(3)=dble(nit(3))
      tr(4)=dble(info)
      tr(5)=dble(phase)
      tr(6)=gap
      tr(7)=0.0d0
      if(nit(1).ge.500.and.gap.gt.1.0d-8)tr(7)=1.0d0
      pr=0.0d0
      du=0.0d0
      bmax=0.0d0
      cmax=0.0d0
      cx=0.0d0
      by=0.0d0
      sw=0.0d0
      tr(12)=1.0d300
      tr(13)=1.0d300
      tr(14)=0.0d0
      tr(15)=0.0d0
      tr(33)=0.0d0
      do 901 j=1,p
      if(.not.ieee_is_finite(y(j)))tr(33)=tr(33)+1.0d0
      v=-b(j)
      do 902 i=1,n
      v=v+a(j,i)*x(i)
  902 continue
      pr=dmax1(pr,dabs(v))
      bmax=dmax1(bmax,dabs(b(j)))
      by=by+b(j)*y(j)
  901 continue
      do 903 i=1,n
      if(.not.ieee_is_finite(x(i)))tr(33)=tr(33)+1.0d0
      if(.not.ieee_is_finite(s(i)))tr(33)=tr(33)+1.0d0
      if(.not.ieee_is_finite(z(i)))tr(33)=tr(33)+1.0d0
      if(.not.ieee_is_finite(w(i)))tr(33)=tr(33)+1.0d0
      v=c(i)-z(i)+w(i)
      do 904 j=1,p
      v=v-a(j,i)*y(j)
  904 continue
      du=dmax1(du,dabs(v))
      cmax=dmax1(cmax,dabs(c(i)))
      cx=cx+c(i)*x(i)
      sw=sw+w(i)
      tr(12)=dmin1(tr(12),x(i),s(i))
      tr(13)=dmin1(tr(13),z(i),w(i))
      tr(14)=dmax1(tr(14),dabs(x(i)+s(i)-1.0d0))
      tr(15)=dmax1(tr(15),dabs(x(i)*z(i)),dabs(s(i)*w(i)))
  903 continue
      tr(8)=pr
      tr(9)=pr/(1.0d0+bmax)
      tr(10)=du
      tr(11)=du/(1.0d0+cmax)
      tr(16)=cx
      tr(17)=by-sw
      tr(18)=cx-by+sw
      tr(19)=dabs(cx-by+sw)/(1.0d0+dabs(cx)+dabs(by-sw))
      tr(24)=gap/(1.0d0+dabs(cx)+dabs(by-sw))
      return
      end
