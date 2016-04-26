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
module smithmodel

  use fourrier
  implicit none

  real :: pi
  real :: minRain, maxRain

  ! Physical grid number
  integer :: nx, ny

  ! Physical grid space
  real :: dx

  ! Wind velocity
  real :: u, v

  ! Moist static stability (/s)
  real :: nm

  ! Cloud time scale
  real :: tauc

  ! Precipitation fallout time scale
  real :: tauf

  ! Uplift sensitivity factor (kg.m-3)
  real :: cw

  ! Water vapor scale height (km)
  real :: hw

  ! Background precipitation (mm/hr)
  real :: prbgd

  ! Topography (real + imaginary)
  real,dimension(:,:),allocatable :: hhr, hhi

  ! Precipitation (real + imaginary)
  real,dimension(:,:),allocatable :: prr, pri

  ! Wavenumbers
  real,dimension(:),allocatable :: wk, wl

  ! Elevation value
  real,dimension(:,:),allocatable :: elev

  ! Precipitation value
  real,dimension(:,:),allocatable :: rainval

contains
  ! =====================================================================================
  subroutine get_rain !(nx, ny, dx, prbgd, u, v, nm2, tauc, tauf, cw, hw)

    integer :: i,j,p
    real :: k1,k2,hr,hi,sigma,mi,mr

    real :: funcr, funci, funcr2, funci2
    real :: nm2

    pi = 4.*atan(1.)

    nm2 = nm * nm

    ! Get topography
    if(allocated(hhr)) deallocate(hhr,hhi)
    allocate(hhr(nx,ny),hhi(nx,ny))
    p = 1
    do j=1,ny
      do i=1,nx
        hhr(i,j) = elev(i,j)
        hhi(i,j) = 0.
        p = p+1
      enddo
    enddo

    ! FFT terrain
    call fasts(hhr,hhi,nx,ny,1,2,1,2)

    ! Wavenumbers
    if(allocated(wk)) deallocate(wk,wl)
    allocate(wk(nx),wl(ny))

    do i=1,nx
      if(i<=nx/2+1)then
        wk(i) = 2.*pi*(i-1)/(nx*dx)
      else
        wk(i) = -2.*pi*(nx-(i-1))/(nx*dx)
      endif
    enddo

    do j=1,ny
      if(j<=ny/2+1)then
        wl(j)=2.*pi*(j-1)/(ny*dx)
      else
        wl(j)=-2.*pi*(ny-(j-1))/(ny*dx)
      endif
    enddo

    ! FFT precipitation
    if(allocated(prr)) deallocate(prr,pri)
    allocate(prr(nx,ny),pri(nx,ny))
    do j=1,ny
      do i=1,nx
        k1 = wk(i)
        k2 = wl(j)
        hr = hhr(i,j)
        hi = hhi(i,j)
        sigma = u*k1 + v*k2
        if(sigma /= 0.)then
          if((nm2/(sigma*sigma)-1.)>=0)then
            mr = sqrt((nm2/(sigma*sigma)-1.)*(k1*k1+k2*k2))*sigma/abs(sigma)
            mi = 0.
          else
            mi = sqrt((1.-nm2/(sigma*sigma))*(k1*k1+k2*k2))
            mr = 0.
          endif
        else
          mi = 0.
          mr = 0.
        endif

        funcr = -cw*sigma*((hi*(1+mi*hw)+hr*mr*hw)*(1-sigma*sigma*tauc*tauf)-(hr*(1+mi*hw)-hi*mr*hw)*sigma*(tauc+tauf))/ &
                   (((1+mi*hw)**2+mr*mr*hw*hw)*(1+sigma*sigma*tauc*tauc)*(1+sigma*sigma*tauf*tauf))
        funci = cw*sigma*((hr*(1+mi*hw)-hi*mr*hw)*(1-sigma*sigma*tauc*tauf)+(hi*(1+mi*hw)+hr*mr*hw)*sigma*(tauc+tauf))/ &
                   (((1+mi*hw)**2+mr*mr*hw*hw)*(1+sigma*sigma*tauc*tauc)*(1+sigma*sigma*tauf*tauf))

        funcr2 = funcr-funci
        funci2 = funci+funcr

        ! Conversion in mm/hr
        funcr2=3600*funcr2
        funci2=3600*funci2

        ! Precipiation
        prr(i,j)=funcr2
        pri(i,j)=funci2
      enddo
    enddo

    ! Calculate the precipitation field in physical space by taking the inverse
    ! Fourier transform of the function calculated above.
    call fasts(prr,pri,nx,ny,1,2,-1,2)

    ! Conversion from mm/hr to m/yr
    do j=1,ny
      do i=1,nx
         prr(i,j) = prr(i,j)*24.*365./1000.
       enddo
    enddo

    return

  end subroutine get_rain
  ! =====================================================================================
  subroutine clip_rain(i1,j1,i2,j2,nx1,ny1)
    ! Calculate the precipitaion in each grid point. This is defined as the precipitation
    ! calculated hourly with the present model plus the background precipitation. If the
    ! precipition is below zero (it can happen with the present model), it is set to zero.
    real :: minprr, maxprr, aa, bb, aa2, bb2
    integer :: i, j, i1, j1, i2, j2, nx1, ny1

    minprr = 1000.
    maxprr = -1000.

    do j=j1,j2
      do i=i1,i2
        ! Conversion mm/hr to m/yr
        rainval(i-i1+1,j-j1+1) = prr(i,j) + prbgd
        minprr = min(rainval(i-i1+1,j-j1+1),minprr)
        maxprr = max(rainval(i-i1+1,j-j1+1),maxprr)
      enddo
    enddo

    aa =  ( prbgd - minRain) / (prbgd - minprr)
    bb = minRain - aa * minprr

    aa2 =  (maxRain - prbgd) / (maxprr - prbgd)
    bb2 = prbgd - aa2 * prbgd

    do j=1,ny1
      do i=1,nx1
        if(rainval(i,j)<=prbgd)then
          rainval(i,j) = aa * rainval(i,j) + bb
        else
          rainval(i,j) = aa2 * rainval(i,j) + bb2
        endif
      enddo
    enddo

    return

  end subroutine clip_rain
  ! =====================================================================================
  subroutine extrapolate_border(i1,j1,i2,j2)

    integer :: i,j,i1,j1,i2,j2

    do j=1,j1-1
      do i=i1,i2
        elev(i,j) = elev(i,j1)
      enddo
    enddo

    do j=j2+1,ny
      do i=i1,i2
        elev(i,j) = elev(i,j2)
      enddo
    enddo

    do i=1,i1-1
      do j=j1,j2
        elev(i,j) = elev(i1,j)
      enddo
    enddo

    do i=i2+1,nx
      do j=j1,j2
        elev(i,j) = elev(i2,j)
      enddo
    enddo

    do i=1,i1-1
      do j=1,j1-1
        elev(i,j) = elev(i1,j1)
      enddo
    enddo

    do i=1,i1-1
      do j=j2+1,ny
        elev(i,j) = elev(i1,j2)
      enddo
    enddo

    do i=i2+1,nx
      do j=1,j1-1
        elev(i,j) = elev(i2,j1)
      enddo
    enddo

    do i=i2+1,nx
      do j=j2+1,ny
        elev(i,j) = elev(i2,j2)
      enddo
    enddo

    return

  end subroutine extrapolate_border
  ! =====================================================================================

end module smithmodel
