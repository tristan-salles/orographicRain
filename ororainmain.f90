!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Compute the precipitation using the linear orographic precipitation model of
! Smith and Barstad (2004).
program orographic_rain

  use smithmodel
  implicit none

  integer :: i,j
  integer :: nx1, ny1
  integer :: i1, j1, i2, j2
  real :: el

  ! Conversion time -- 200 to 2000 s
  tauc = 1000.
  ! Fallout time -- 200 to 2000 s
  tauf = 1000.
  ! Background precipitation rate -- 0 to 5 m/a
  prbgd = 2.
  ! Rain min/max range
  minRain = 1.
  maxRain = 3.
  ! Horizontal wind components (m/s) -- 1 to 100 m/s
  u = 5.
  v = 0.
  ! Moist stability frequency -- 0 to 0.01 s 1
  nm = 0.005
  ! Uplift sensitivity factor -- 0.001 to 0.02 kg/m3
  cw = 0.005
  ! Water vapor scale height -- 1 to 5 km in m
  hw = 2500.
  ! Grid size
  nx1 = 534
  ny1 = 719
  dx = 250.
  nx1 = 300
  ny1 = 300
  dx = 1000

  ! Define computational grid
  nx = int(nx1/2.) + nx1
  ny = int(ny1/2.) + ny1
  i1 = int(nx1/4.)
  j1 = int(ny1/4.)
  i2 = int(nx1/4.) + nx1 - 1
  j2 = int(ny1/4.) + ny1 - 1

  ! Read elevation file
  allocate(elev(nx,ny))
  open(10,file='data/elev.csv',status='old')
  !open(10,file='nodes.csv',status='old')
  elev(:,:) = 0.
  do j=j1,j2
    do i=i1,i2
      read(10,*) el
      !read(10,*) a,b,el
      elev(i,j)=el
      if(el<0) elev(i,j)=0.
    enddo
  enddo
  close(10)

  ! Extrapolate on borders
  call extrapolate_border(i1,j1,i2,j2)

  ! Compute orographic rain
  call get_rain

  ! Define precipitation using user range
  allocate(rainval(nx1,ny1))
  call clip_rain(i1,j1,i2,j2,nx1,ny1)

  ! Write output
  open(10,file='data/precip.csv',status='unknown')
  write(10,*)'i,j,z,p'
  do j=1,ny1
    do i=1,nx1
      write(10,*) (i)*dx,',',(j)*dx,',',elev(i+i1-1,j+j1-1),',',rainval(i,j)
    enddo
  enddo
  close(10)

  return

end program orographic_rain
