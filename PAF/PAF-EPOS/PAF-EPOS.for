c**********************************************************************
c                                                                     *
c                Inver P A F    (25/03/2017)                          * 
c                                                                     *
c     A software tool for inversion of terrain deformations           *
c        as free-geometry extended bodies for anomalous pressure.     *
c                                                                     *
c     A.G. Camacho (1), J. Fernández (1), F. Cannavó (2)              *
c     (1) Institute of Geosciences (CSIC-UCM), Madrid, Spain          *
c     (2) Osservatorio Etneo, INGV, Catania, Italy.                   *
c                                                                     *
c     See PAFmanual.txt for description of the input and output files *
c      and operation approach                                         *
c                                                                     *
c**********************************************************************
      use ifport
c      use ifqwin

      implicit real*8(a-h,o-z)
      parameter (ms=5000,mc=100000,mb=2000,mm=5000)

      dimension fs(2),x(ms),y(ms),z(ms),rab(ms),
     -dg(ms),wg(ms),rg(ms),
     -dz(ms),wz(ms),zcp(ms),zcd(ms),rz(ms),zo(ms),za(ms),
     -dx(ms),wx(ms),xcp(ms),xcd(ms),rx(ms),xo(ms),xa(ms),
     -dy(ms),wy(ms),ycp(ms),ycd(ms),ry(ms),yo(ms),ya(ms),
     -qg(mc),qz(mc),xc(mc),yc(mc),zc(mc),dxc(mc),dyc(mc),dzc(mc),vc(mc),
     -zas(ms),di(ms),sz(ms),sen(mc),
     -cdz(ms),cdx(ms),cdy(ms),cpz(ms),cpx(ms),cpy(ms),cdg(ms),cpg(ms)

      integer*4 m(mc),p(mc),nor(mc),js(2),ipa(mc),ifa(mc)
      integer*2 nada,idz(ms),idh(ms),ks(2),ld(50),lp(50)
      character*50 obs,mod,gri,fil,flos,tex(12),raya
      character*9 hoy,hora1,hora2, fic(mm)
      character*130 texto

      data raya/'-'/,pi/3.1141592/
      call seed(2016)
      idum=1000
      blun=9.0
      nuf=1
      sig=0.25
      dmu=10
      depli=3000.         

c---------------------------Input data and covariance matrix-------------!
      j=time()
      idum=j-int(j/1.e3)*1.e3
      iu=1
      iu0=0
      ih0=0
      los=0
      dpr=0
      is3=0
      lx=560
      ly=190
      lz=460
      tx=360
      ty=360
      lx1=lx-tx/2
      lx2=lx+tx/2
      lx3=lx-tx/2+380
      lx4=lx+tx/2+380
      ly1=ly-ty/2
      ly2=ly+ty/2

      open(1,file='PAF-Parameters.txt')
      read(1,'(///a7)') hora2
      read(1,*) pem
      read(1,*) jsmo
      read(1,*) jseli
      read(1,*) nepo
      read(1,*) iep1
      if(nc.gt.mc) call dim(mc)
      close(1)                
      
      write(*,200)  pem,jsmo,jseli,nepo,iep1
  200 format('      PAF-modeling of pressure sources'
     -//'   Parameters'
     -/' -----------------------------------------' 
     -/ f7.1,3x,'....Pressure change (MPa)        '
     -/   i7,3x,'....Smoothing coeff. (0<smo<1000)'
     -/   i7,3x,'....Significance limit (0<sig<10)'
     -/   i7,3x,'....Number of data epochs        ' 
     -/   i7,3x,'....First epoch identification   ' 
     -/' -----------------------------------------') 


      dlap=jsmo
      seli=jseli
      iz=1
      ih=1

      open(2,file='DeforData.txt')
      i=0
      ax=9.d9
      ay=ax
      bx=-ax
      by=-ay
      xm=0
      ym=0
      zm=0
    1 read(2,*,err=2,end=2) xp,yp,zp
      if(xp.eq.0.and.yp.eq.0) go to 2
      if(xp.lt.ax) ax=xp
      if(xp.gt.bx) bx=xp
      if(yp.lt.ay) ay=yp
      if(yp.gt.by) by=yp
      xm=xm+xp
      ym=ym+yp
      zm=zm+zp
      i=i+1
      if(i.gt.ms) write(*,*) ' ***** Num of data site > ',i6
      if(i.gt.ms) stop
      x(i)=xp
      y(i)=yp
      z(i)=zp
      go to 1
    2 ns=i
      close(2)
              write(*,'(/10x,a,i5)') 'Num. data points =',ns
      xm=xm/ns
      ym=ym/ns
      zm=zm/ns
      if(ns.eq.0) stop
      na=4
      if(ns.gt.500) na=9
      if(ns.lt.30) na=2
      fe=(bx-ax+by-ay)/2./30000.
      do 3 i=1,ns
      x(i)=x(i)-xm
      y(i)=y(i)-ym
    3 z(i)=z(i)-zm
                                          
      do 6 i=1,nepo
      no=iep1-1+i
      if(no.le.99999) write(fic(i),'(i5)') no
      if(no.le.9999) write(fic(i),'(a1,i4)') '0',no
      if(no.le.999) write(fic(i),'(a2,i3)') '00',no
      if(no.le.99) write(fic(i),'(a3,i2)') '000',no
      if(no.le.9) write(fic(i),'(a4,i1)') '0000',no
    6 continue
                 
            call cells(ms,mc,ns,x,y,z,nc,xc,yc,zc,dxc,dyc,dzc)
            
            write(*,'(10x,a,i7)') 'Num. volume cells=',nc
            write(*,'(//20x,a)') ' Wait please ....'
      
      siz=0
      do 5 j=1,nc
      p(j)=0
      sen(j)=0
      vo=dxc(j)*dyc(j)*dzc(j)/fe/fe/fe
      r=(dxc(j)+dyc(j)+dzc(j))/3./fe
      siz=siz+r
      if(r.lt.500) r=500
      r2=r*r
      at2=0
      dv=0.
      sw=9.d9
         cc=0
      do 4 i=1,ns
      cx=(x(i)-xc(j))/fe                         
      cy=(y(i)-yc(j))/fe
      cz=(z(i)-zc(j))/fe
      d2=cx*cx+cy*cy+cz*cz
      if(d2.lt.sw) dv=cz
      if(d2.lt.sw) sw=d2
         if(d2.lt.r2) d2=r2
      c=cz*cz/d2/d2/d2
      at2=at2+c
         if(c.gt.cc) cc=c
    4 continue
         sen(j)=sqrt((at2-cc)/(ns-1))*1.e9
         if(dv.lt.depli) sen(j)=sen(j) *dv/depli
      at2=at2/ns
      qz(j)=vo*vo*at2                    
      vc(j)=vo*fe*fe*fe                  
    5 continue
      siz=siz/nc

c-------------Inversion process--------------------------------------

      do 88 ii=1,nepo                                      
      write(obs,'(a1,a5,a4)') 'D',fic(ii),'.txt'           
      write(mod,'(a1,a5,a4)') 'M',fic(ii),'.txt'           

      if(nepo.eq.1) write(obs,'(a13)') 'DeforData.txt'    
      if(nepo.eq.1) write(mod,'(a10)') 'ModPAF.txt'       
                           
                      write(*,*) 'Mod.',ii
      
      nz=0
      nx=0
      ny=0
      swz=0
      swx=0
      swy=0
      sw=0
      ez=1
      ex=1
      ey=1
      if(iz.eq.0) ez=0
      if(ih.eq.0) ex=0
      if(ih.eq.0) ey=0
      open(2,file=obs)
      k=0
      do 13 i=1,ns
      read(2,*,err=14,end=14) xp,yp,zp,(tex(j),j=1,6)
      k=k+1
      if(xp.eq.0.and.yp.eq.0.) go to 14
      idz(i)=0
      dz(i)=0.
      if(iz.ne.1.or.tex(1).eq.raya.or.tex(2).eq.raya) go to 12
       read(tex(1),*) dz(i)
       read(tex(2),*) erz
       if(erz.eq.0) go to 12
        idz(i)=1
        wz(i)= erz
        swz=swz+wz(i)
        nz=nz+1
   12 idh(i)=0
      dx(i)=0.
      dy(i)=0.
      if(ih.eq.0) go to 13
      if(tex(3).eq.raya.or.tex(4).eq.raya.or.tex(5).eq.raya) go to 13
       read(tex(3),*) dx(i)
       read(tex(4),*) erx
       read(tex(5),*) dy(i)
       read(tex(6),*) ery
       if(ery.eq.0.and.erx.eq.0) go to 13
       idh(i)=1
       wx(i)= erx   !1
       wy(i)= ery   !1
       swx=swx+wx(i)
       swy=swy+wy(i)
       nx=nx+1
       ny=ny+1
   13 continue
   14 close(2)
      if(k.ne.ns) write(*,*) 'Data error'
      if(k.ne.ns) stop
      if(nz.gt.0) swz=swz/nz
      if(nx.gt.0) swx=swx/nx
      if(ny.gt.0) swy=swy/ny
      cx=swx
      cy=swy
      cz=swz

      sdz=0
      sdx=0
      sdy=0
      emp=0.
      swu=0
      swz=0
      swx=0
      swy=0
      dzm=0
      do 15 i=1,ns
      dzm=dzm+dz(i)
      if(idz(i).eq.1) wz(i)=wz(i)/cz
      if(idz(i).eq.2) wz(i)=wz(i)/cz2
      if(idz(i).eq.3) wz(i)=wz(i)/cz3
      if(cx.ne.0) wx(i)=wx(i)/cx
      if(cy.ne.0) wy(i)=wy(i)/cy
      if(sw.ne.0) wg(i)=wg(i)/sw
      zcp(i)=0.
      zcd(i)=0.
      xcp(i)=0.
      xcd(i)=0.
      ycp(i)=0.
      ycd(i)=0.
      if(idz(i).eq.1) then
       wz(i)=ez*wz(i)
       sdz=sdz+dz(i)*dz(i)
       emp=emp+dz(i)*dz(i)*wz(i)
       swu=swu+wz(i)
       swz=swz+wz(i)
      endif
       dx(i)=dx(i)-dx0
       dy(i)=dy(i)-dy0
       if(idh(i).eq.1) then
       wx(i)=ex*wx(i)
       wy(i)=ey*wy(i)
       sdx=sdx+dx(i)*dx(i)
       sdy=sdy+dy(i)*dy(i)
       emp=emp+dx(i)*dx(i)*wx(i)
       swu=swu+wx(i)
       swx=swx+wx(i)
       emp=emp+dy(i)*dy(i)*wy(i)
       swu=swu+wy(i)
       swy=swy+wy(i)
      endif
   15 continue
      dzm=dzm/ns
      if(nz.gt.0) sdz=sqrt(sdz/nz)
      if(ii.eq.1) tdef=nint(sdz*2.)/1.0
      if(nx.gt.0) sdx=sqrt(sdx/nx)
      if(ny.gt.0) sdy=sqrt(sdy/ny)
      if(swu.gt.0.) emp1=sqrt(emp/swu)
      nrej=0 
      c=5 
      do 16 i=1,ns
      cx=0
      cy=0
      cz=0
      s=0
      do 17 j=1,ns    
      if(i.eq.j) go to 17
      d=1./((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
      s=s+d
      cx=cx+dx(j)*d
      cy=cy+dy(j)*d
      cz=cz+dz(j)*d
   17 continue   
      cx=cx/s    
      cy=cy/s    
      cz=cz/s    
      if(abs(dz(i)-cz).gt.c*sdz) idz(i)=0
      if(abs(dx(i)-cx).gt.c*sdx) idh(i)=0
      if(abs(dy(i)-cy).gt.c*sdy) idh(i)=0
      k=idz(i)*idh(i)
   16 if(k.eq.0) nrej=nrej+1
      

      nu5=2
      nrem=0
      svp=0.
      fps1=0
      ozs=0
      oxs=0
      oys=0
      ol2=0
      ol3=0
      nu=0
      npr=0
      do 28 j=1,nc
      p(j)=-j
   28 nor(j)=0
          nsig=0
      kpf=1
      ls=1
   20 continue
      fitu=999.d9
      ncs=nc
      if(na.gt.1.and.nu.gt.10) ncs=nc/na
      do 21 i=1,ns
      if(idz(i).eq.1) then
        zo(i)=dz(i) -ozs
        zas(i)=zcp(i)
      endif
      if(idz(i).eq.2) then
        zo(i)=dz(i)-ol2
        zas(i)=-xcp(i)*cs1+ycp(i)*ss1+zcp(i)*cl1
      endif
      if(idz(i).eq.3) then
        zo(i)=dz(i) -ol3
        zas(i)=-xcp(i)*cs2+ycp(i)*ss2+zcp(i)*cl2
      endif
      if(idh(i).eq.0) go to 21
        xo(i)=dx(i) -oxs
        yo(i)=dy(i) -oys
   21 continue
              dlp=dlap*1.e-4/siz *(1.-1./(nu+1))     
      do 30 ja=1,ncs
      j=ja
      if(ncs.lt.nc) call RANDOM_NUMBER(ranval)
      if(ncs.lt.nc) j=(nc-1)*ranval+1
      vo=vc(j)
      sv=svp+qz(j)*dlp
      if(p(j).ge.0) go to 30
      call deformod(dmu,sig,fra,xc(j),yc(j),zc(j),vo,ms,ns,x,y,z,1,1,
     -cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
      do 24 k=1,2
      pb=2*k-3
      std=0.
      sdd=0.
      do 23 i=1,ns
      xp=cpx(i)
      yp=cpy(i)
      zp=cpz(i)
      if(idz(i).ne.0) then
       if(idz(i).eq.1) c=zp
       if(idz(i).eq.2) c=-xp*cs1+yp*ss1+zp*cl1
       if(idz(i).eq.3) c=-xp*cs2+yp*ss2+zp*cl2
       za(i)=zas(i)+c*pb
       zp=za(i)
       ww=wz(i)*zp
       std=std+ww*zo(i)
       sdd=sdd+ww*zp
      endif
      if(idh(i).eq.0) go to 23
       xa(i)=xcp(i)+xp*pb
       xp=xa(i)
       ww=wx(i)*xp
       std=std+ww*xo(i)
       sdd=sdd+ww*xp
       ya(i)=ycp(i)+yp*pb
       yp=ya(i)
       ww=wy(i)*yp
       std=std+ww*yo(i)
       sdd=sdd+ww*yp
   23 continue
      fp=std/(sdd+sv)
      if(fp.le.0.) go to 24
      em=0.
      do 22 i=1,ns
      if(idz(i).eq.0) go to 22
      c=zo(i)-za(i)*fp
      em=em+c*c*wz(i)
      if(idh(i).eq.0) go to 22
      c=xo(i)-xa(i)*fp
      em=em+c*c*wx(i)
      c=yo(i)-ya(i)*fp
      em=em+c*c*wy(i)
   22 continue
      em=em/swu
      c=em+fp*fp*sv
      if(c.ge.fitu) go to 24
       fitu=c
       fs(1)=fp
       js(1)=j
       ks(1)=k
       kps=1
       emp=em
   24 continue
   30 continue
      j=js(ls)
      if(j.eq.0) go to 50
      vo=vc(j)

       call deformod(dmu,sig,fra,xc(j),yc(j),zc(j),vo,ms,ns,x,y,z,
     - 1,1,cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
       fitu1=fitu
       npr=npr+1
       emp1=emp
       fps1=fs(1)
       jp=js(1)
       p(jp)=ks(1)
       nor(jp)=fps1*10
       if(npr.le.10) lp(npr)=jp
       pb=2*ks(1)-3
       svp=svp+qz(jp)*dlp
       ipa(jp)=1
       ifa(jp)=1
       sdx=0.
       sdy=0.
       sdz=0.
       sdl2=0
       sdl3=0
       cz2=0.
       cz3=0.
       cz=0
       cx=0
       cy=0

       do 33 i=1,ns
       xp=cpx(i)
       yp=cpy(i)
       zp=cpz(i)
       gp=cpg(i)
       zcp(i)=zcp(i)+zp*pb
       xcp(i)=xcp(i)+xp*pb
       ycp(i)=ycp(i)+yp*pb
       d=zcd(i)
       if(idz(i).eq.2) d=-xcd(i)*cs1+ycd(i)*ss1+zcd(i)*cl1
       if(idz(i).eq.3) d=-xcd(i)*cs2+ycd(i)*ss2+zcd(i)*cl2
       zz=ozs
       if(idz(i).eq.2) zz=ol2
       if(idz(i).eq.3) zz=ol3
       c=dz(i)-zz
       d=zcp(i)
       if(idz(i).eq.2) d=-xcp(i)*cs1+ycp(i)*ss1+zcp(i)*cl1
       if(idz(i).eq.3) d=-xcp(i)*cs2+ycp(i)*ss2+zcp(i)*cl2
       rz(i)=c-d*fps1
       c=c-d
       if(idz(i).eq.1) then
         sdz=sdz+wz(i)*rz(i)*rz(i)
         cz=cz+wz(i)*rz(i)
       endif
       if(idz(i).eq.2) then
         sdl2=sdl2+wz(i)*rz(i)*rz(i)
         cz2=cz2+wz(i)*rz(i)
       endif
       if(idz(i).eq.3) then
         sdl3=sdl3+wz(i)*rz(i)*rz(i)
         cz3=cz3+wz(i)*rz(i)
       endif
       rx(i)=dx(i)-oxs -xcp(i)*fps1
       ry(i)=dy(i)-oys -ycp(i)*fps1
       if(idh(i).eq.0) go to 33
       sdx=sdx+wx(i)*rx(i)*rx(i)
       sdy=sdy+wy(i)*ry(i)*ry(i)
       cx=cx+wx(i)*rx(i)
       cy=cy+wy(i)*ry(i)
   33  continue
       if(swz.gt.0.) sdz=sqrt(sdz/swz)
       if(swl2.gt.0.) sdl2=sqrt(sdl2/swl2)
       if(swl3.gt.0.) sdl3=sqrt(sdl3/swl3)
       if(swx.gt.0.) sdx=sqrt(sdx/swx)
       if(swy.gt.0.) sdy=sqrt(sdy/swy)
       if(swz2.gt.0) cz2=cz2/swl2
       if(swz3.gt.0) cz3=cz3/swl3
       if(swx.gt.0.) cx=cx/swx
       if(swy.gt.0.) cy=cy/swy
       if(swz.gt.0.) cz=cz/swz
       if(iu0.eq.1) ozs=ozs+cz*0.5
       if(ih0.eq.1) oxs=oxs+cx*0.5
       if(ih0.eq.1) oys=oys+cy*0.5
       if(il0.eq.1) then
        ol2=ol2+cz2*0.5
        ol3=ol3+cz3*0.5
       endif
      
      tt=fps1/pem
      call dmedian(ms,ns,rz,em)
      em=em/0.6745*blun *tt*tt*tt
      mout=0
      d=0.
      do 34 i=1,ns
   34 if(wz(i).gt.d) d=wz(i)
      do 35 i=1,ns
      c=rz(i)/em
      if(c.le.1.) go to 35
       c=c*c-1.
       if(c.lt.0.5) mout=mout+1
       c=d/(1.+c*c)
       if(c.lt.wz(i)) wz(i)=c
   35 continue
           if(sen(jp).ge.seli) nsig=nsig+1
      nu=npr

                                write(*,'(i6,f6.0\)') nu,fps1/pem
      
   36 continue
      l=ls
      ls=l
      js(1)=0
      js(2)=0
      if(fps1.gt.0.and.fps1.lt.pem.and.npr.gt.2) kpf=0
      if(kpf.ne.0) go to 20

c----------------------------------------------------------------------+
c--------------  End of process ---------------------------------------+

   50 continue
      call time(hora2)

      open(1,file=mod)
      open(2,file='PAF-Parameters.txt')
      do 59 i=1,10
      read(2,'(a130)') texto
   59 write(1,'(a130)') texto
      write(1,'(/a//a/a)')
     -' Cells: location (UTM, m) , sides (m)  and  press (MPa)',
     -'    X      Y       Z      sx   sy   sz     Press   Signi',
     -'---------------------------------------------------------'
      nn=0
      do 51 j=1,nc
   51 if(p(j).gt.0.or.m(j).gt.0) nn=nn+1

      totmas=0
      totpre=0
      swz=0
        npp=0
        npn=0
      do 52 j=1,nc
      jx=xc(j)+xm
      jy=yc(j)+ym
      jz=zc(j)+zm
      kx=dxc(j)
      ky=dyc(j)
      kz=dzc(j)
      vol=kx*ky*kz
      pb=0.
      i=p(j)
      if(i.le.0) go to 52
      pb=(2*i-3)*fps1    !presion
        if(i.eq.1) npn=npn+1
        if(i.eq.2) npp=npp+1
   58 totpre=totpre+vol*abs(pb)
      swz=swz+jz*vol*abs(pb)
      k=sen(j)
      write(1,'(i7,i8,i7,2x,3i5,f9.2,i5)') jx,jy,jz,kx,ky,kz,pb,k
   52 continue
      swz=swz/totpre
      write(1,'(a,a9,a)') 
     -'--------------------- Date:',hoy,'------------------'
      write(1,202) ns,nc,npn+npp,npn,npp,sig,dmu,na
  202 format(/' Num. data points =',i6
     -/' Num. total cells=',i6
     -/' Num.filled cells=',i5,4x,'Neg, Pos =',2i5
     -/' Medium param.:',6x,'Poisson=',f4.2,4x,'Share Mod.=',f4.0,'GPa'
     -/' Random explor.coeff.=',i4)
      write(1,203) fps1,totpre/1.e9,swz,sdz,sdx,sdy,fitu1
  203 format(' Pressure model:',6x,'Press. contrast (-+)=',f7.2,' MPa'
     -/22x,'Press*vol:',f7.1,' MPA*Km3',4x,' Mean mod.depth =',f7.0
     -/' RMS residuals (cm):  Up=',f5.2,4x,'WE=',f5.2,4x,'SN=',f5.2,
     -5x,'Misfit=',f9.4)
   53 close(2)
      write(1,'(a,2a12)') ' Initial and final exec.times:',hora1,hora2

      write(1,'(///20x,a)') 'Observed, modeled, and residual values'
      write(1,'(/a,10x,a,12x,a,16x,a/2x,a,6x,a,7x,a,8x,a,10x,a,10x,a
     -/a,a)')
     -'Data point loc(UTM, m)','dz (cm)','dx(cm)','dy(cm)','X','Y',
     -'Z','obs  mod  res','obs  mod  res','obs mod res    Weight'
     -,'----------------------------------------------',
     -'-----------------------------------------------'
      k=0
      emg=0
      emu=0
      do 55 i=1,ns
      if(idz(i).eq.0) go to 55
      emu=emu+dz(i)*dz(i)
      k=k+1
   55 continue
      if(k.gt.0) emu=sqrt(emu/k)
      cg=0
      cz=0
      c2=0
      c3=0
      do 56 i=1,ns
      xp=xm+x(i)
      yp=ym+y(i)
      zp=zm+z(i)
      ie=0
      if(emg.ne.0.) ie=abs(rg(i))/emg
      if(d.lt.0.01) d=0.01
      up=0
      if(emu.ne.0.) up=abs(rz(i))/emu
      zcal=dz(i)-rz(i)
      gcal=dg(i)-rg(i)
      xcal=dx(i)-rx(i)
      ycal=dy(i)-ry(i)
      d=0
      if(idz(i).eq.2) d=dl2+dxl2*x(i)/1000.+dyl2*y(i)/1000.
      if(idz(i).eq.3) d=dl3+dxl3*x(i)/1000.+dyl3*y(i)/1000.
      write(1,210) nint(xp),nint(yp),nint(zp),
     -dz(i)+d,zcal+d,rz(i),
     -dx(i)+dx0,xcal+dx0,rx(i),dy(i)+dy0,ycal+dy0,ry(i),wz(i)
  210 format(i6,i8,i5,2x,3(3f7.2,1x),f5.1)
   56 continue
      write(1,'(a,a//60x,a,i4)') 
     -'----------------------------------------------',
     -'-----------------------------------------------',
     -'Number of outlier values=',mout
      close(1)

   88 continue                                  !<<<

      write(*,'(//20x,a)') ' End of process'
      stop
      end

c***********************************************************************
c   Dimension                                                          *
c***********************************************************************
      subroutine dim(m)
      write(*,200) m
  200 format(/6x,'*** Error: data size > max. dimension',i6,' !!!'/
     -10x,'Reduce the data side or change dimension in the code'/)
      stop
      end
C***********************************************************************
c Surface modeled elevation changes (cm)                               *
c  for pressure and mass = 1 MPa, 1 Kg                                 *
c***********************************************************************
      subroutine deformod(dmu,sig,fra,xc,yc,zc,vo,ms,ns,x,y,z,kt,kf,
     -cdz,cdx,cdy,cdg,cpz,cpx,cpy,cpg)
      implicit real*8(a-h,o-z)
      dimension x(ms),y(ms),z(ms),
     -cdz(ms),cdx(ms),cdy(ms),cpz(ms),cpx(ms),cpy(ms),cdg(ms),cpg(ms),
     -cx(14),cy(14),cz(14)
      data cx/0, 1, 0, 0, 0.7071, 0.7071,0.7071,-0.7071, 0.0000, 0.0000,
     - 0.57735,-0.57735, 0.57735, 0.57735/
      data cy/0, 0, 1, 0, 0.7071,-0.7071,0.0000, 0.0000, 0.7071,-0.7071,
     - 0.57735, 0.57735,-0.57735, 0.57735/
      data cz/0, 0, 0, 1, 0.0000, 0.0000,0.7071, 0.7071, 0.7071, 0.7071,
     - 0.57735, 0.57735, 0.57735,-0.57735/
      csp=0.577
      pi=3.14159265359
      dera=pi/180.
      u1=9.8066/4./pi/dmu/1.e8          
      ra3=0.750/pi*vo
      siz=vo**0.333333
      u2=(1.-sig)/dmu*ra3/10.   
      do 10 i=1,ns
      xp=x(i)-xc
      yp=y(i)-yc
      zp=z(i)-zc
      x2=xp*xp
      y2=yp*yp
      z2=zp*zp
      d2=x2+y2+z2 +siz
      d=sqrt(d2)
      d3=d2*d
      cdx(i)=0
      cdy(i)=0
      cdz(i)=0
      cdg(i)=0
      cpx(i)=0
      cpy(i)=0
      cpz(i)=0
      cpg(i)=0
      if(kt.eq.0) then   
        cdz(i)=-u1/d*(2.*(1.-sig)+z2/d2)*vo        
        pp=-u1/d*(zp/d2-(1.-2.*sig)/(d+zp))*vo     
        cdx(i)=pp*xp                               
        cdy(i)=pp*yp                               
        cdg(i)=vo*zp/d3*6.672d-3 - cdz(i)*fra      
      endif
      if(kf.le.1) go to 2
      if(kt.eq.2.and.kf.gt.1) then
        if(kf.eq.2.and.yp.gt.0) go to 10
        if(kf.eq.3.and.yp.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.3.and.kf.gt.1) then
        if(kf.eq.2.and.xp.gt.0) go to 10
        if(kf.eq.3.and.xp.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.4.and.kf.gt.1) then
        if(kf.eq.2.and.xp.gt.0) go to 10
        if(kf.eq.3.and.xp.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.5.and.kf.gt.1) then
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.6.and.kf.gt.1) then
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.zp.gt.0) go to 10
        if(kf.eq.5.and.zp.lt.0) go to 10
      endif
      if(kt.eq.7.and.kf.gt.1) then
        s=xp-zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.8.and.kf.gt.1) then
        s=xp+zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.yp.gt.0) go to 10
        if(kf.eq.5.and.yp.lt.0) go to 10
      endif
      if(kt.eq.9.and.kf.gt.1) then
        s=yp-zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.xp.gt.0) go to 10
        if(kf.eq.5.and.xp.lt.0) go to 10
      endif
      if(kt.eq.10.and.kf.gt.1) then
        s=yp+zp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        if(kf.eq.4.and.xp.gt.0) go to 10
        if(kf.eq.5.and.xp.lt.0) go to 10
      endif
      if(kt.eq.11.and.kf.gt.1) then        
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp-zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.12.and.kf.gt.1) then    
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp-zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.13.and.kf.gt.1) then     
        s=xp+yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp+zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
      if(kt.eq.14.and.kf.gt.1) then     
        s=xp-yp
        if(kf.eq.2.and.s.gt.0) go to 10
        if(kf.eq.3.and.s.lt.0) go to 10
        s=yp+zp
        if(kf.eq.4.and.s.gt.0) go to 10
        if(kf.eq.5.and.s.lt.0) go to 10
      endif
    2 p=u2/d3  
      if(kt.eq.1) then    
       cpx(i)=p*xp *1.28
       cpy(i)=p*yp *1.28
       cpz(i)=p*zp
       go to 3
      endif
      if(kt.gt.1) p=p*1.0   
      if(kf.gt.1) p=p*1.0   
      s=p* d*(1-csp)        
      c=d2-cx(kt)*cx(kt)*x2-cy(kt)*cy(kt)*y2-cz(kt)*cz(kt)*z2
      c=p* (cx(kt)*xp+cy(kt)*yp+cz(kt)*zp)/sqrt(c)*csp
      cpx(i)=s*cx(kt)+c*(1-cx(kt))*xp
      cpy(i)=s*cy(kt)+c*(1-cy(kt))*yp
      cpz(i)=s*cz(kt)+c*(1-cz(kt))*zp
    3 cpg(i)=-cpz(i)*fra      
   10 continue
      return
      end

c***********************************************************************
      subroutine ORDEN(mx,nx,x,ix)
      real*8 x(mx)
      integer ix(mx)
      ix(1)=1
      if(nx.eq.1) return
      ia=1
      ib=nx
      do 1 i=1,nx
      if(x(i).gt.x(ib)) ib=i
    1 if(x(i).lt.x(ia)) ia=i
      ix(1)=ia
      do 2 i=2,nx
      i1=i-1
      im=ix(i1)
      ia=ib
      do 3 j=1,nx
      if(x(j).le.x(im)) go to 3
      if(x(j).lt.x(ia)) ia=j
    3 continue
      ix(i)=ia
    2 continue
      return
      end
c***********************************************************************
c   median of absolute values for n<m different numbers x(i)           c
c***********************************************************************
      subroutine dmedian(m,n,x,xa)
      implicit real*8(a-h,o-z)
      dimension x(m)
      xi=x(1)
      xs=xi
      do 1 i=1,n
      xx=x(i)
      if(xx.lt.xi)  xi=xx
      if(xx.gt.xs)  xs=xx
        do 6 j=i+1,n
    6   if(xx.eq.x(j)) x(j)=x(j)*1.00001
    1 continue
      j=0
    4 xa=xs
      do 2 i=1,n
      xx=x(i)
      if(xx.lt.xa.and.xx.gt.xi) xa=xx
    2 continue
      j=j+1
      if(j.eq.n) go to 9
      xi=xa
      xb=xi
      do 3 i=1,n
      xx=x(i)
      if(xx.gt.xb.and.xx.lt.xs) xb=xx
    3 continue
      j=j+1
      xa=(xa+xb)/2.
      xs=xb
      if(j.lt.n) go to 4
    9 return
      end
      
c***********************************************************************
c                                                                     *
c   Determination of a 3D partition of the subsurface volume          *
c             into nc prismatic cells                                 *
c                                                                     *     
c**********************************************************************
      subroutine cells(ms,mc,ns,x,y,z,nc,xc,yc,zc,dx,dy,dz)
      implicit real*8(a-h,o-z)
      dimension at(mc),di(ms),x(ms),y(ms),z(ms),
     -xc(mc),yc(mc),zc(mc),dx(mc),dy(mc),dz(mc)
      
      xm=0.
      ym=0.
      zm=0.
      ax=9.d9
      bx=-ax
      ay=ax
      by=-ay
      techo=-9.d9 
      do 1 i=1,ns
      if(x(i).lt.ax) ax=x(i)
      if(x(i).gt.bx) bx=x(i)
      if(y(i).lt.ay) ay=y(i)
      if(y(i).gt.by) by=y(i)
      if(z(i).gt.techo) techo=z(i)    
      xm=xm+x(i)
      ym=ym+y(i)
    1 zm=zm+z(i)
      xm=nint(xm/ns/10.)*10.
      ym=nint(ym/ns/10.)*10.
      zm=nint(zm/ns/10.)*10.
      fe=(bx-ax+by-ay)/2. /30000
      pp=800*fe
      base=-15000*fe+zm
      avm=pp*pp*pp*pp*1.d-10   
      nc=0
      kr=0
   10 ncr=nc
      kr=kr+1
      zs=techo
      ax=ax-17000*fe
      ay=ay-17000*fe
      bx=bx+17000*fe
      by=by+17000*fe
   12 k=100-(zs-base)/(techo-base)*100
      tz=pp*0.8
      call especapa(zs,avm,ms,ns,x,y,z,tz,fe) 
      if(tz.gt.1) go to 11
         zs=zs-pp*0.5
         if(zs.gt.base) go to 12
   11 if(kr.eq.2) tz=tz*2 
      if(kr.eq.3) tz=tz*2.7
      call celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,ncr,xc,yc,zc,dx,dy,dz,at)
      zs=zs-tz                      
      if(zs.gt.base) go to 12  
      l=nc
      do 15 i=nc+1,ncr
      if(nc.eq.0) go to 16
      xx=dx(i) /2.       ! 1.9
      yy=dy(i) /2.       ! 1.9    
      zz=dz(i) /1.9
      k=0
      do 14 j=1,nc
      if(abs(xc(j)-xc(i)).gt.xx) go to 14
      if(abs(yc(j)-yc(i)).gt.yy) go to 14
      if(abs(zc(j)-zc(i)).gt.zz) go to 14
      k=1  
   14 continue
      if(k.eq.1) go to 15
   16 l=l+1
      xc(l)=xc(i)
      yc(l)=yc(i)
      zc(l)=zc(i)
      dx(l)=dx(i)
      dy(l)=dy(i) 
      dz(l)=dz(i)
   15 continue
      nc=l
      if(kr.lt.1) go to 10   
      
      return
      end
      
c***********************************************************************
c     Sensitivity calculus
c***********************************************************************
      subroutine att(ms,ns,x,y,z,xx,yy,zz,dx,dy,dz,av,zp ,fe)
      implicit real*8(a-h,o-z)
      dimension x(ms),y(ms),z(ms)
      av=0.
      t=(dx+dy+dz)/3. 
      if(t.lt.500*fe) t=500*fe   
      t2=t*t
      zp=0
      sp=0 
      do 1 i=1,ns
      zr=zz-z(i)
      xr=xx-x(i)
      yr=yy-y(i)
      d=xr*xr+yr*yr+zr*zr
      if(d.lt.t2) d=t2
      c=zr*zr/d/d/d
      av=av+c
      d=1/d/d
      sp=sp+d
      zp=zp+z(i)*d
    1 continue
      av=av/ns
      t2=dx*dy*dz
      av=av*t2*t2  
      zp=zp/sp
      return
      end
c***********************************************************************
c     Step dr for autocorrelation analysis
c***********************************************************************
      subroutine step(m,n,x,y,d,dr)
      implicit real*8(a-h,o-z)
      integer x(m),y(m)
      real*4 d(m)
      dm=0.
      do 1 i=1,n
      d4=9.d9
      d3=9.d9
      d2=9.d9
      d1=9.d9
      do 2 j=1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      dd= xx*xx+yy*yy
      if(dd.gt.d4.or.dd.le.1.) go to 2
      if(dd.lt.d1) d1=dd
      if(dd.lt.d2.and.dd.gt.d1) d2=dd
      if(dd.lt.d3.and.dd.gt.d2) d3=dd
      if(dd.lt.d4.and.dd.gt.d3) d4=dd
    2 continue
      d1=sqrt(d1)
      d2=sqrt(d2)
      d3=sqrt(d3)
      d4=sqrt(d4)
    1 d(i)=(d1+d2+d3+d4)/4.0
      call dmedian(m,n,d,dr)
      return
      end
c***********************************************************************
c          Cell design 1                                               *
c***********************************************************************
      subroutine celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,nc,xb,yb,zb,dx,dy,dz,at)
      implicit real*8(a-h,o-z)
      dimension xb(mc),yb(mc),zb(mc),dx(mc),dy(mc),dz(mc),
     -x(ms),y(ms),z(ms) 
      character*50 texto
      real*4 at(mc)
      zi=zs-tz                       
      zz=zs-tz/2
      tx=tz     
      ty1=tz/3  
      ty2=tz*3  
      tt=tz*tz 
      do 7 i=1,999                
      xx=ax+tx*(i-1)
      if(xx.gt.bx) go to 9
      ty=tz
      yr=ay
      do 6 j=1,999                 
      if(yr.ge.by) go to 7
      c=ty
      l=0
    3 yy=yr+c/2
      call att(ms,ns,x,y,z,xx,yy,zz,tx,c,tz,am,zp, fe)
      if(zp.lt.zs.or.am.le.0) go to 6      
      av=am/avm
      if(abs(av-1.).le.0.1) go to 4        
      c=c/av**0.41
      if(c.lt.ty1.or.c.gt.ty2) go to 6    
      l=l+1
      if(l.le.10) go to 3
      go to 6
    4 if(c.gt.9999) go to 6
      ty=c                           
      nc=nc+1
      if(nc.ge.mc) then 
        write(*,*) ' >> Warning: number of cells >',mc 
        stop   
      endif  
      xb(nc)=xx
      yb(nc)=yy
      zb(nc)=(zi+zs)/2.
      dx(nc)=tx
      dy(nc)=ty
      dz(nc)=(zs-zi)
      at(nc)=av
    6 yr=yr+ty
    7 continue
    9 return
      end
c***********************************************************************
c               Cells design 2                                         * 
c***********************************************************************
      subroutine especapa(zs,avm,ms,ns,x,y,z,tz,fe)
      implicit real*8(a-h,o-z)
      dimension x(ms),y(ms),z(ms)
      l=0
      c=tz                            
    1 zz=zs-c*0.5                     
      k=0
      am=-9.d9
      do 2 i=1,ns
      call att(ms,ns,x,y,z,x(i)*1.5,y(i)*1.5,zz,c,c,c,aa,zp,fe)
      if((zp-zz).le.100.*fe) go to 2
      if(aa.gt.am) am=aa
      k=k+1
    2 continue
      if(k.le.1) go to 9
      am=am/avm
      c=c/am**0.13
      l=l+1
      if(abs(am-1.).gt.0.1.and.l.lt.12) go to 1
      if(c.gt.tz) tz=c
    9 return
      end
c************************************************************************
c*************************End code **************************************
