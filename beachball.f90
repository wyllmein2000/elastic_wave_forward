! Date:        19 Mar 2011
! Last Update: 23 May 2014

  program main
  implicit none
  integer,parameter :: N1=360
  integer,parameter :: N2=181
  real :: strike,dip,rake
  real :: azimul,takeoff,pi
  real :: x(N1,N2)
  real :: y,z,z0
  real :: da,dt
  integer :: na,nt,npar
  integer :: i,j
  real,external :: firstmotion
  character(len=10) :: str

  npar=iargc()
  if (npar<3) stop 'Arguments: strike, dip, rake (degree)'

  call getarg (1, str)
  read (str,*) strike
  call getarg (2, str)
  read (str,*) dip
  call getarg (3, str)
  read (str,*) rake
  pi=4.*atan(1.)/180.0

  na=N1
  nt=N2

  call beachball (x,na,nt,strike,dip,rake)

  da=360.0/na
  dt=180.0/(nt-1)

  open (10,file='focal.dat')
  do i=1, na
     azimul=(i-1)*da
  do j=1, (nt+1)/2
     takeoff=-(j-1)*dt
     y=tan(takeoff*pi/2.0)*sin(azimul*pi)
     z=tan(takeoff*pi/2.0)*cos(azimul*pi)
     write(10,*) y,z,x(i,j)
  enddo
  enddo
  close(10)

! azimul=210.0
! takeoff=90.0+atan(6.0/12.0)/pi
! azimul=azimul+180.0
! takeoff=180.0-takeoff
! z0=firstmotion (azimul,takeoff,strike,dip,rake)
! y=tan(takeoff*pi/2.0)*sin(azimul*pi)
! z=tan(takeoff*pi/2.0)*cos(azimul*pi)
! write(6,*) 'first motion =', y,z,z0

  stop
  end program main

  subroutine beachball (x,na,nt,strike,dip,rake)
  implicit none
  real,parameter :: pi=4.*atan(1.0)/180.0
  integer,intent(in) :: na,nt
  real,intent(in) :: strike,dip,rake
  real,intent(out) :: x(na,nt)
  real :: azimul,takeoff
  real :: s,d,r
  real :: da,dt
  integer :: i,j
  real,external :: firstmotion
  s=strike*pi
  d=dip*pi
  r=rake*pi
  da=360.0/na*pi
  dt=180.0/(nt-1)*pi
  do i=1, na
     azimul=(i-1)*da
     do j=1, nt
        takeoff=(j-1)*dt
        x(i,j)=firstmotion (takeoff,azimul,s,d,r)
     enddo
  enddo
  return
  end subroutine beachball

  real function firstmotion (t,a,s,d,r)
  implicit none
  real :: s,d,r,t,a
  real :: y
  y=cos(r)*sin(d)*sin(t)*sin(t)*sin(2.*(a-s))
  y=y-cos(r)*cos(d)*sin(2.*t)*cos(a-s)
  y=y+sin(r)*sin(2.*d)*(cos(t)*cos(t)-sin(t)*sin(t)*sin(a-s)*sin(a-s))
  y=y+sin(r)*cos(2.*d)*sin(2.*t)*sin(a-s)
  firstmotion=y
  return
  end function firstmotion
