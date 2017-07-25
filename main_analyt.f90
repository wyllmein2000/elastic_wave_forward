! -------------------------------------------
! Date: 19 May 2014
! -------------------------------------------
  program main
  implicit none
  integer,parameter :: lbyte=4
  real,parameter :: pi=4.0*atan(1.0)
  real :: vp, vs, rho
  real :: mu, area, slip, mo
  real :: mxx, myy, mzz, mxy, myz, mzx, mt(3,3)
  real :: strike, dip, rake
  real :: xs, ys, zs, tr, td, nv(3)
  real :: x, y, z
  real,allocatable :: s(:,:)
  real,allocatable :: u(:,:,:),f(:,:,:)
  real,allocatable :: xr(:),yr(:),zr(:)

  real :: dt
  integer :: nt
  integer :: n,ng
  integer :: flag,istype
  integer :: i,j,k
  character(len=30) :: frec,str

! -----------------------------------
!
! -----------------------------------
  open (1,file='elas3d0.in')
  read (1,*) flag
! flag=1            ! compute complete wavefield
! flag=2            ! far-field approximation

! -----------------------------------
! model parameters
! -----------------------------------
  read (1,*) vp, vs
  read (1,*) rho
  read (1,*) nt, dt

  mu=rho*vs*vs

! -----------------------------------
! receiver information
! -----------------------------------
  read (1,*) ng
  if (ng>0) then

     allocate(xr(ng))
     allocate(yr(ng))
     allocate(zr(ng))

     do i=1, ng
        read (1,*) xr(i),yr(i),zr(i)
     enddo

  else
     read (1,*) frec

     open (2,file=frec)
     read (2,*) ng

     allocate(xr(ng))
     allocate(yr(ng))
     allocate(zr(ng))

     do i=1, ng
        read (2,*) xr(i),yr(i),zr(i)
     enddo
     close(2)
  endif

  read (1,*) nv
! nv=0.0; nv(3)=1.0    ! normal direction

! -----------------------------------
! source function
! -----------------------------------
  read (1,*) area, slip
  read (1,*) tr, td
! slip=0.01
! tr=0.0     
! td=0.005

  allocate(s(0:nt,2))
  call sourcefun (s,nt,dt,tr,td)
  call writedata ('accsour.bin',s(:,2),nt+1,1,lbyte)

! -----------------------------------
! source mechanism
! itype = 1, only shear slip
!       = 2, arbitrary source
! -----------------------------------
  read (1,*) xs,ys,zs
  read (1,*) istype

  if (istype<=2) then
     read (1,*) strike,dip,rake        ! degree

     mo=mu*area*slip

     strike=strike*pi/180.0
        dip=   dip*pi/180.0
       rake=  rake*pi/180.0

  else
     read (1,*) mxx,myy,mzz,mxy,myz,mzx
  endif

  read (1,*) str
  close(1)


! ------------------------------------
! ------------------------------------
  allocate(u(ng,nt,3))
  allocate(f(ng,nt,3))

  
! ------------------------------------------------
! shear + complete (N+I+F)
! ------------------------------------------------
  if (istype==1) then

  do k=  1, ng
     if (mod(k,1000)==1) write(6,*) 'i=',k

     x=xr(k)-xs
     y=yr(k)-ys
     z=zr(k)-zs

     call synthetic (u(k,:,:),s,x,y,z,vp,vs,dt,nt)

  enddo

  u=u*mo/(4.0*pi*rho)


! ------------------------------------------------
! arbitrary + far-field approximation (F)
! ------------------------------------------------
  elseif (istype>=2) then

     if (istype==2) then
        call calmoment (mxx,myy,mzz,mxy,myz,mzx,strike,dip,rake)
     endif

     mt(1,1)=mxx
     mt(1,2)=mxy
     mt(1,3)=mzx
     mt(2,1)=mxy
     mt(2,2)=myy
     mt(2,3)=myz
     mt(3,1)=mzx
     mt(3,2)=myz
     mt(3,3)=mzz
     write(6,*) mt

     xr=xr-xs
     yr=yr-ys
     zr=zr-zs

!    write(6,*) strike,dip,rake
!    write(6,*) mo,myz
!    write(6,*) xr,yr,zr

     call farfield (u,f,s(:,2),mt,tr+td,xr,yr,zr,nv,vp,vs,ng,nt,dt)

     u= u*mo/(4.0*pi*rho)
     f= f*mo/(4.0*pi)

     call writedata (trim(str)//'tx.bin',f(:,:,1),ng,nt,lbyte)
     call writedata (trim(str)//'ty.bin',f(:,:,2),ng,nt,lbyte)
     call writedata (trim(str)//'tz.bin',f(:,:,3),ng,nt,lbyte)

  endif

  call writedata (trim(str)//'vx.bin',u(:,:,1),ng,nt,lbyte)
  call writedata (trim(str)//'vy.bin',u(:,:,2),ng,nt,lbyte)
  call writedata (trim(str)//'vz.bin',u(:,:,3),ng,nt,lbyte)

  deallocate(s)
  deallocate(u,f)
  deallocate(xr,yr,zr)
  stop
  end program main
! ======================================================

! ======================================================
  subroutine synthetic (u,s,xr,yr,zr,vp,vs,dt,nt)
  implicit none
  integer,intent(in) :: nt
  real(kind=8),parameter :: pi=4.0*atan(1.0)
  real,intent(in) :: s(0:nt,2),xr,yr,zr,vp,vs,dt
  real,intent(out) :: u(nt,3)
  real(kind=8) :: r, rh, take, azi
  real(kind=8) :: p1, p2, p3
  real(kind=8) :: conv(3,3),anf(3),aip(3),ais(3),afp(3),afs(3)
  real :: cosazi, sinazi, costak, sintak
  integer :: i,j,k,np,ns,it2

  real,allocatable :: si(:)

  u=0.0
  r=sqrt(xr*xr+yr*yr+zr*zr)
  if (r==0.0) return

  rh=sqrt(xr*xr+yr*yr)
! if (zr==0.0) then
!    take=0.5d0*pi
! elseif (zr>0d0) then
!    take=datan(rh/zr)
! else
!    take=pi+datan(rh/zr)
! endif
! if (rh==0d0) then
!    azi=0d0
! elseif (xr>0d0) then
!    azi=dasin(yr/rh)
! else
!    azi=pi-dasin(yr/rh)
! endif

! write(6,*) take*180.0/pi, azi*180.0/pi
! write(6,*) take*180.0/pi

  if (r==0.0) then
     costak=0.0
     sintak=0.0
  else
     costak=zr/r
     sintak=rh/r
  endif
  if (rh==0.0) then
     cosazi=0.0
     sinazi=0.0
  else
     cosazi=xr/rh
     sinazi=yr/rh
  endif

! p1=dsin(2d0*take)*dcos(azi)
! p2=dcos(2d0*take)*dcos(azi)
! p3=dcos(take)*dsin(azi)
  p1=2.0*sintak*costak*cosazi
  p2=(costak*costak-sintak*sintak)*cosazi
  p3=costak*sinazi

! spherical => xyz
! conv(1,1)=dsin(take)*dcos(azi)   !    r->x
! conv(2,1)=dsin(take)*dsin(azi)   !    r->y
! conv(3,1)=dcos(take)
! conv(1,2)=dcos(take)*dcos(azi)   ! take->x
! conv(2,2)=dcos(take)*dsin(azi)   ! take->y
! conv(3,2)=-dsin(take)
! conv(1,3)=-dsin(azi)            !  azi->x
! conv(2,3)= dcos(azi)            !  azi->y
  conv(1,1)=sintak*cosazi         !    r->x
  conv(2,1)=sintak*sinazi         !    r->y
  conv(3,1)=costak
  conv(1,2)=costak*cosazi         ! take->x
  conv(2,2)=costak*sinazi         ! take->y
  conv(3,2)=-sintak
  conv(1,3)=-sinazi               !  azi->x
  conv(2,3)= cosazi               !  azi->y
  conv(3,3)=0.0



  do k=1, 3
     anf(k)= 9d0*p1*conv(k,1)-6d0*p2*conv(k,2)+6d0*p3*conv(k,3)
     aip(k)= 4d0*p1*conv(k,1)-2d0*p2*conv(k,2)+2d0*p3*conv(k,3)
     ais(k)=-3d0*p1*conv(k,1)+3d0*p2*conv(k,2)-3d0*p3*conv(k,3)
     afp(k)=     p1*conv(k,1)
     afs(k)=                      p2*conv(k,2)-    p3*conv(k,3)
  enddo
  anf=anf/(r**4)
  aip=aip/(r*vp)**2
  ais=ais/(r*vs)**2
  afp=afp/(r*vp**3)
  afs=afs/(r*vs**3)

! write(6,*) anf
! write(6,*) aip
! write(6,*) ais
! write(6,*) afp
! write(6,*) afs
 
  np=int(r/(vp*dt))
  ns=int(r/(vs*dt))

  allocate(si(0:nt))
  if (np<=nt) then
  do i=np, nt
     it2=i
     if (it2>ns) it2=ns
  do j=np, it2
     si(i)=si(i)+(j-1)*s(i-j,1)
  enddo
  enddo
  si=si*dt*dt
  endif

! ------------------------------------------------------
! near field
! ------------------------------------------------------
  do k=1, nt
  do i=1, 3
     u(k,i)=anf(i)*si(k)
  enddo
  enddo

! ------------------------------------------------------
! intermediate field
! ------------------------------------------------------
  do k=np, nt
  do i=1, 3
     u(k,i)=u(k,i)+aip(i)*s(k-np,1)+ais(i)*s(k-ns,1)
  enddo
  enddo

! ------------------------------------------------------
! far field
! ------------------------------------------------------
  do k=np, nt
  do i=1, 3
     u(k,i)=u(k,i)+afp(i)*s(k-np,2)+afs(i)*s(k-ns,2)
  enddo
  enddo

  deallocate(si)
  return
  end subroutine synthetic
! ======================================================
! s(:,1) --- s(t),      slip
! s(:,2) --- s1=ds/dt,  velocity
! s(:,3) --- s2=ds1/dt, acceleration
! ======================================================
  subroutine sourcefun (s,nt,dt,tr,td)
  implicit none
  real,parameter :: pi2=8.0*atan(1.0)
  integer,intent(in) :: nt
  real,intent(in) :: dt,tr,td
  real,intent(out) :: s(0:nt,2)
  real :: t
  integer :: i,j,i1,i2
  s=0.0
  i1=max(tr,0.0)/dt
  i2=(tr+td)/dt
  write(6,*) i1, i2
  if (i2>nt) i2=nt
  do i=i1, i2
     t=(i*dt-tr)/td
     s(i,2)=pi2/(td*td)*sin(pi2*t)
     s(i,1)=(1.0-cos(pi2*t))/td
  enddo
  return
  end subroutine sourcefun
! ======================================================
  subroutine farfield (vv,ff,s,mo,te,xr,yr,zr,nv,vp,vs,nr,nt,dt)
  implicit none
  integer,intent(in) :: nr,nt
  real,intent(in) :: vp,vs,te,dt
  real,intent(in) :: s(0:nt),mo(3,3),nv(3),xr(nr),yr(nr),zr(nr)
  real,intent(out) :: vv(nr,nt,3),ff(nr,nt,3)
  real :: ro(3),delta(3,3),vp2,vs2,vp3,vs3
  real :: pv(3),qv(3),pf(3),qf(3),ftmp,p,r
  integer :: i,j,k,l,m,n
  integer :: itp,its,itbeg,itend,ite

  delta=0.0
  do k=1, 3
     delta(k,k)=1.0
  enddo
  vp2=vp*vp
  vs2=vs*vs
  vp3=vp*vp2
  vs3=vs*vs2
  p=2.0*(vs/vp)**2

  vv=0.0
  ff=0.0

  do n=1, nr
     if (mod(n,1000)==1) &
        write(6,*) 'n=',n
     r=sqrt(xr(n)*xr(n)+yr(n)*yr(n)+zr(n)*zr(n))
     ro(1)=xr(n)/r
     ro(2)=yr(n)/r
     ro(3)=zr(n)/r
     itp=r/vp/dt
     if (itp>nt) cycle
     ftmp=ro(1)*nv(1)+ro(2)*nv(2)+ro(3)*nv(3)
!    write(6,*) itp,itend

     pv=0.0; qv=0.0
     pf=0.0; qf=0.0
     do k=1, 3
     do i=1, 3
     do j=1, 3
        pv(k)=pv(k)+(ro(k)*ro(i)-delta(k,i))*ro(j)*mo(i,j)
        qv(k)=qv(k)+ ro(k)*ro(i)            *ro(j)*mo(i,j)
        pf(k)=pf(k)+((1.0-p)*nv(k)+p*ro(k)*ftmp)*ro(i)*ro(j)*mo(i,j)
        qf(k)=qf(k)+((2.0*ro(k)*ro(i)-delta(k,i))*ro(j)*ftmp-ro(k)*ro(j)*nv(i))*mo(i,j)
     enddo
     enddo
     enddo

!    goto 10
     itbeg=itp
     if (itp<1) itbeg=1
     itend=(r/vp+te)/dt
     if (itend>nt) itend=nt
     do m=itbeg, itend
     do k=1, 3
        vv(n,m,k)= qv(k)/vp3*s(m-itp)
        ff(n,m,k)=-pf(k)/vp2*s(m-itp)
     enddo
     enddo

10   continue
!    goto 5
     its=r/vs/dt
     if (its>nt) goto 5
     itbeg=its
     if (its<1) itbeg=1
     itend=(r/vs+te)/dt
     if (itend>nt) itend=nt
     do m=itbeg, itend
     do k=1, 3
        vv(n,m,k)=vv(n,m,k)-pv(k)/vs3*s(m-its)
        ff(n,m,k)=ff(n,m,k)+qf(k)/vs2*s(m-its)
     enddo
     enddo

5    continue
     vv(n,:,:)=vv(n,:,:)/r
     ff(n,:,:)=ff(n,:,:)/r
  enddo
! call writedata ('s.bin',s,nt+1,1,4)

  return
  end subroutine farfield
! ======================================================
  subroutine calmoment (mxx,myy,mzz,mxy,myz,mzx,strike,dip,rake)
  implicit none
  real,intent(in) :: strike,dip,rake
  real,intent(out) :: mxx,myy,mzz,mxy,myz,mzx
     mxx=-sin(dip)*cos(rake)*sin(2.0*strike)-sin(2.0*dip)*sin(rake)*sin(strike)**2
     myy=+sin(dip)*cos(rake)*sin(2.0*strike)-sin(2.0*dip)*sin(rake)*cos(strike)**2
     mzz=+sin(2.0*dip)*sin(rake)
     mxy=+sin(dip)*cos(rake)*cos(2.0*strike)+sin(2.0*dip)*sin(rake)*sin(2.0*strike)/2.0
     myz=-cos(dip)*cos(rake)*sin(strike)+cos(2.0*dip)*sin(rake)*cos(strike)
     mzx=-cos(dip)*cos(rake)*cos(strike)-cos(2.0*dip)*sin(rake)*sin(strike)
  return
  end subroutine calmoment
! ======================================================
