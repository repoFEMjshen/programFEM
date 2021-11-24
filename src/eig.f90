

      subroutine rudec(f,r,u)
      implicit real*8(a-h,o-z)
      dimension f(3,3),r(3,3),u(3,3),xs(3),xsi(3),        &
          c(3,3),c2(3,3),uinv(3,3)
      do i=1,3
         do j=1,3
            c(i,j)=0.d0
            do k=1,3
               c(i,j)=c(i,j)+f(k,i)*f(k,j)
            enddo
         enddo
      enddo

      call eig3(c,xs,xsi,ireal)
      do i=1,3
         xs(i)=dsqrt(xs(i))
      enddo
      xi1=xs(1)+xs(2)+xs(3)
      xi2=xs(1)*xs(2)+xs(2)*xs(3)+xs(3)*xs(1)
      xi3=xs(1)*xs(2)*xs(3)
      d=xi1*xi2-xi3
     
      do i=1,3
         do j=1,3
            c2(i,j)=0.d0
            do k=1,3
               c2(i,j)=c2(i,j)+c(i,k)*c(k,j)
            enddo
         enddo
      enddo

      do i=1,3
         do j=1,3
            u(i,j)=-c2(i,j)+(xi1**2.d0-xi2)*c(i,j)
            u(i,j)=u(i,j)/d
         enddo
         u(i,i)=u(i,i)+xi1*xi3/d
      enddo

      do i=1,3
         do j=1,3
            uinv(i,j)=c(i,j)-xi1*u(i,j)
            uinv(i,j)=uinv(i,j)/xi3
         enddo
         uinv(i,i)=uinv(i,i)+xi2/xi3
      enddo

      do i=1,3
         do j=1,3
            r(i,j)=0.d0
            do k=1,3
               r(i,j)=r(i,j)+f(i,k)*uinv(k,j)
            enddo
         enddo
      enddo
      return
      end
         

      
      


      subroutine eig3(a,xlambda,xlambdai,ireal)
      implicit real*8(a-h,o-z)
      dimension a(3,3),ai(3,3),xlambda(3),xlambdai(3)
      xi1=a(1,1)+a(2,2)+a(3,3)
      xi2=0.d0
      do i=1,3
         do k=1,3
            xi2=xi2+a(i,k)*a(k,i)
         enddo         
      enddo
      xi2=0.5d0*(xi1**2.d0-xi2)

      xi3=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)     &
          *a(3,3)+a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)     &
          *a(1,3)*a(2,2))
      call cubsolve(1.d0,-xi1,xi2,-xi3,xlambda,xlambdai,ireal)


      return
      end
      


      subroutine cubsolve(a1,b1,c1,d1,xsol,xsoli,ireal)
      implicit real*8(a-h,o-z)
      dimension xsol(3),xsoli(3)
      a=b1/a1
      b=c1/a1
      c=d1/a1
      q=(3.d0*b-a**2.d0)/9.d0
      r=(9.d0*a*b-27.d0*c-2.d0*a**3.d0)/54.d0
      p=q**3.d0+r**2.d0
      if(p.ge.0.d0)then         
         s=(r+dsqrt(p))**(1.d0/3.d0)
         t=(r-dsqrt(p))**(1.d0/3.d0)
         r1=s+t-a/3.d0
      else
         p=dabs(p)
         p=dsqrt(p)
         xm=r**2.d0+p**2.d0
         xm=dsqrt(xm)
         theta=datan2(p,r)
         r1=2.d0*xm**(1.d0/3.d0)*dcos(theta/3.d0)-a/3.d0
      endif
      temp1=-(a+r1)
      temp2=-3.d0*r1**2.d0-2.d0*a*r1+a**2.d0-4.d0*b
      if(temp2.ge.0.d0)then
         ireal=1
         temp2=dsqrt(temp2)
         xsol(1)=r1
         xsol(2)=(temp1+temp2)/2.d0
         xsol(3)=(temp1-temp2)/2.d0
      else
         ireal=0
         xsol(1)=r1
         xsol(2)=temp1/2.d0
         xsol(3)=temp1/2.d0
         temp2=dabs(temp2)
         temp2=dsqrt(temp2)
         xsoli(2)=temp2/2.d0
         xsoli(3)=-temp2/2.d0
      endif
      return
      end



      subroutine sprinc(s,ps,lstr,ndi,nshr)
      implicit real*8(a-h,o-z)
      dimension s(6),ps(3),sm(3,3),v(3,3),xli(3)
      do i=1,ndi
         sm(i,i)=s(i)
      enddo
      fac=1.d0/(1.d0*lstr)
      sm(1,2)=fac*s(4)
      if(nshr.gt.1)then
         sm(2,3)=fac*s(5)
         sm(1,3)=fac*s(6)
      endif
      do i=1,3
         do j=i+1,3
            sm(j,i)=sm(i,j)
         enddo
      enddo
      call eig3(sm,ps,xli,ireal)
!      call jacobi(sm,ndi,3,ps,v,nrot)
      return
      end

      subroutine jacobi(a,n,np,d,v,nrot)
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      dimension a(np,np),d(np),v(np,np),b(nmax),z(nmax)
      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.d0
         enddo
         v(ip,ip)=1.d0
      enddo
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.d0
      enddo
      nrot=0
      do i=1,50
         sm=0.d0
         do ip=1,n-1
            do iq=i+1,n
               sm=sm+dabs(a(ip,iq))
            enddo
         enddo
         if(dabs(sm).lt.1d-16)return
         tresh=1d-16
         do ip=1,n-1
            do iq=ip+1,n
               if((i.gt.400).and.(dabs(a(ip,iq)/d(ip))).lt.1.d-10 &               
                      .and.(dabs(a(ip,iq)/d(iq)).lt.1.d-10))then
                  a(ip,iq)=0.d0
               elseif(dabs(a(ip,iq)).gt.tresh)then
                  h=d(iq)-d(ip)
                  if(dabs(a(ip,iq)/h).lt.1.d-12)then
                     t=a(ip,iq)/h
                  else
                     theta=0.5d0*h/a(ip,iq)
                     t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
                     if(theta.lt.0)t=-t  
                  endif
                  c=1.d0/sqrt(1.d0+t**2.d0)
                  s=t*c
                  tau=s/(1.d0+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.d0
                  do j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
                  enddo
                  do j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
                  enddo
                  nrot=nrot+1
               endif
            enddo
         enddo
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.d0
         enddo
      enddo
      write(*,*)'Too many iterations for eigenvalues'
      return
      end
