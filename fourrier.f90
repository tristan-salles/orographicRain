!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!
!!                                                                                   !!
!!  This file forms part of the Badlands surface processes modelling application.    !!
!!                                                                                   !!
!!  For full license and copyright information, please refer to the LICENSE.md file  !!
!!  located at the project root, or contact the authors.                             !!
!!                                                                                   !!
!!~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~!!

! Fast Fourrier Transform.
module fastfourrier

  public :: fasts
  private :: fft, reals, fftmx, bkdt, realt, istkgt, istkrl

contains

  ! =====================================================================================
  subroutine fasts(a, b, n1, n2, n3, ndim, rdi, rrc)
    ! FFT for:
    ! complex serie a(i) = Re[z(i)], b(i) = Im[z(i)] of dimension n1 x n2 x n3
    ! real serie a(i), b(i)
    ! + ndim=1 1D Serie 1D (n2=n3=1)                                                                                               !
    ! + ndim=2 2D Serie 2D (n3=1)                                                                                                  !
    ! + ndim=3 3D Serie 3D
    ! + rdi=1 direct transform
    ! + rdi=-1 inverse transform
    ! + rrc=2 Complex serie
    ! + rcc=1 Real serie

    real :: a(*),b(*)
    common/cstak/dstak(3000)
    double precision dstak
    integer :: n1, n2, n3, ndim
    integer :: istak(6000), lout, lused, rdi, rrc
    equivalence (dstak(1),istak(1))
    equivalence (istak(1),lout)
    equivalence (istak(3),lused)

    ! Initialisation of cstak
    call bkdt

    ! Parameters initialisation
    isn = -1
    if(rdi == -1) isn = 1
    if(rrc == 2) goto 10

    ! FFT real serie
    nn1 = n1/2
    nn = nn1 * n2 * n3
    goto (1,2,3) ndim

    ! 1D FFT
1   if( rdi /= -1) call fft(a,b,1,nn1,1,-1)
    call reals(a, b, nn1, isn)
    if( rdi == -1) call fft(a,b,1,nn1,1,1)
    goto 100

    ! 2D FFT
2   if(rdi == -1) goto 4
    call fft(a, b, n2, nn1, 1, -1)
    call realt(a, b, n2, nn1, 1, -1)
    call fft(a(nn+1), b(nn+1), 1, n2, 1, -1)
    call fft(a, b, 1, n2, nn1, -1)
    goto 100

4   call fft(a, b, 1, n2, nn1, 1)
    call fft(a(nn+1), b(nn+1), 1, n2, 1, 1)
    call realt(a, b, n2, nn1, 1, 1)
    call fft(a, b, n2, nn1, 1, 1)
    goto 100

    ! 3D FFT
3   if(rdi == -1) goto 5
    call fft(a, b, n3*n2, nn1, 1, -1)
    call realt(a, b, n3*n2, nn1, 1, -1)
    call fft(a(nn+1), b(nn+1), n3, n2, 1, -1)
    call fft(a(nn+1), b(nn+1), 1, n3, n2, -1)
    call fft(a, b, n3, n2, nn1, -1)
    call fft(a, b, 1, n3, n2*nn1, -1)
    goto 100

5   call fft(a,b,1,n3,n2*nn1,1)
    call fft(a,b,n3,n2,nn1,1)
    call fft(a(nn+1),b(nn+1),1,n3,n2,1)
    call fft(a(nn+1),b(nn+1),n3,n2,1,1)
    call realt(a,b,n3*n2,nn1,1,1)
    call fft(a,b,n3*n2,nn1,1,1)
    goto 100

    ! FFT complex serie
10  nn = n1 * n2 * n3
    goto (11,12,13) ndim

    ! 1D FFT
11  call fft(a, b, 1, n1, 1, isn)
    goto 100

    ! 2D FFT
12  call fft(a, b, n2, n1, 1, isn)
    call fft(a, b, 1, n2, n1, isn)
    goto 100

    ! 3D FFT
13  call fft(a, b, n2*n3, n1, 1, isn)
    call fft(a, b, n3, n2, n1, isn)
    call fft(a, b, 1, n3, n1*n2, isn)

100 continue

    j = lout
    if(j /= 0) call istkrl(j)

    return

  end subroutine fasts
  ! =====================================================================================
  subroutine fft(a, b, nseg, n, nspn, isn)

    dimension nfac(15)
    real a(*),b(*)
    common/cstak/dstak(3000)
    double precision dstak
    integer istak(6000)
    real rstak(6000)
    equivalence (dstak(1),istak(1))
    equivalence (dstak(1),rstak(1))

    nfac(1:15) = 0
    m = 0
    nf = iabs(n)
    k = nf
    if(nf == 1) return
    nspan = iabs(nf * nspn)
    ntot = iabs(nspan * nseg)
    if(isn * ntot /= 0) goto 20
    write(7,9999) nseg, n, nspn, isn
    return

10  m = m + 1
    nfac(m) = 4
    k = k/16

20  if(k - (k/16) * 16 == 0) goto 10
    j = 3
    jj = 9
    goto 40

30  m = m + 1
    nfac(m) = j
    k = k/jj

40  if(mod(k,jj) == 0) goto 30
    j = j + 2
    jj = j**2
    if( jj <= k) goto 40
    if(k > 4) goto 50
    kt = m
    nfac(m + 1) = k
    if(k /= 1) m = m + 1
    goto 90

50  if(k-(k/4)*4 /= 0) goto 60
    m=m+1
    nfac(m)=2
    k=k/4

    ! All square factores out now, but k.ge.5 still
60  kt=m
    maxp=max0(kt+kt+2,k-1)
    j=2

70  if(mod(k,j) /= 0) goto 80
    m=m+1
    nfac(m)=j
    k=k/j

80  j = ((j+1)/2)*2+1
    if(j<=k) goto 70

90  if(m<=kt+1) maxp=m+kt+1
    if(m+kt>15) goto 120
    if(kt==0) goto 110
    j=kt

100 m=m+1
    nfac(m)=nfac(j)
    j=j-1
    if(j/=0) goto 100

110 maxf=m-kt
    maxf=nfac(maxf)
    if(kt>0) maxf=max0(nfac(kt),maxf)
    j=istkgt(maxf*4,3)
    jj=j+maxf
    j2=jj+maxf
    j3=j2+maxf
    k=istkgt(maxp,2)
    call fftmx(a,b,ntot,nf,nspan,isn,m,kt,rstak(j),rstak(jj),rstak(j2),rstak(j3),istak(k),nfac)
    call istkrl(2)

    return

120 write(0,9998) n

    return
9998    format(' Erro - parametro n(fft) tem mais de 15 factores ',i10)
9999    format(' Erro - zero nos parametros da fft ',4i10)

  endsubroutine fft
  ! =====================================================================================
  subroutine reals(a,b,n,isn)

    real :: a(*),b(*)

    inc=iabs(isn)
    nf=iabs(n)
    if(nf*isn.ne.0) goto 10
    write(0,9999) n,isn
    return

10  nk=nf*inc+2
    nh=nk/2
    rad=atan(1.0)
    dr=-4.0/float(nf)
    cd=2.0*sin(0.5*dr*rad)**2
    sd=sin(dr*rad)
    lim=32
    mm=lim
    ml=0
    sn=0.0
    if(isn.gt.0) goto 40
    cn=1.0
    a(nk-1)=a(1)
    b(nk-1)=b(1)

20  continue
    do 30 j=1,nh,inc
      k=nk-j
      aa=a(j)+a(k)
      ab=a(j)-a(k)
      ba=b(j)+b(k)
      bb=b(j)-b(k)
      re=cn*ba+sn*ab
      em=sn*ba-cn*ab
      b(k)=(em-bb)*0.5
      b(j)=(em+bb)*0.5
      a(k)=(aa-re)*0.5
      a(j)=(aa+re)*0.5
      ml=ml+1
      if(ml.eq.mm) goto 50
      aa=cn-(cd*cn+sd*sn)
      sn=(sd*cn-cd*sn)+sn
      cn=0.5/(aa**2+sn**2)+0.5
      sn=cn*sn
      cn=cn*aa
      goto 30
50    mm=mm+lim
      sn=float(ml)*dr*rad
      cn=cos(sn)
      if(isn.gt.0) cn=-cn
      sn=sin(sn)
30  continue

    return

40  cn=-1.0
    sd=-sd

    goto 20

9999    format(' Error - parameter in reals function is equal to zero',2i10)

  end subroutine reals
  ! =====================================================================================
  subroutine fftmx(a,b,ntot,n,nspan,isn,m,kt,at,ck,bt,sk,np,nfac)

    integer :: np(*),nfac(*)
    real    :: at(*),ck(*),bt(*),sk(*)
    real    :: a(*),b(*)

    inc=iabs(isn)
    nt=inc*ntot
    ks=inc*nspan
    rad=atan(1.0)
    s72=rad/0.625
    c72=cos(s72)
    s72=sin(s72)
    s120=sqrt(0.75)
    if(isn>0) goto 10
    s72=-s72
    s120=-s120
    rad=-rad
    goto 30

    ! Scale by 1/n para isn.gt.0
10  ak=1.0/float(n)
    do j=1,nt,inc
      a(j)=a(j)*ak
      b(j)=b(j)*ak
    enddo

30  kspan=ks
    nn=nt-inc
    jc=ks/n
    lim=32
    klim=lim*jc
    i=0
    jf=0
    maxf=m-kt
    maxf=nfac(maxf)
    if(kt>0) maxf=max0(nfac(kt),maxf)

40  dr=8.0*float(jc)/float(kspan)
    cd=2.0*sin(0.5*dr*rad)**2
    sd=sin(dr*rad)
    kk=1
    i=i+1
    if(nfac(i)/=2) goto 110

    ! Factor 2
    kspan=kspan/2
    k1=kspan+2

50  k2=kk+kspan
    ak=a(k2)
    bk=b(k2)
    a(k2)=a(kk)-ak
    b(k2)=b(kk)-bk
    a(kk)=a(kk)+ak
    b(kk)=b(kk)+bk
    kk=k2+kspan
    if(kk<=nn) goto 50
    kk=kk-nn
    if(kk<=jc) goto 50
    if(kk>kspan) goto 350

60  c1=1.0-cd
    s1=sd
    mm=min0(k1/2,klim)
    goto 80

70  ak=c1-(cd*c1+sd*s1)
    s1=(sd*c1-cd*s1)+s1
    c1=0.5/(ak**2+s1**2)+0.5
    s1=c1*s1
    c1=c1*ak

80  k2=kk+kspan
    ak=a(kk)-a(k2)
    bk=b(kk)-b(k2)
    a(kk)=a(kk)+a(k2)
    b(kk)=b(kk)+b(k2)
    a(k2)=c1*ak-s1*bk
    b(k2)=s1*ak+c1*bk
    kk=k2+kspan
    if(kk.lt.nt) goto 80
    k2=kk-nt
    c1=-c1
    kk=k1-k2
    if(kk.gt.k2) goto 80
    kk=kk+jc
    if(kk.le.mm) goto 70
    if(kk.lt.k2) goto 90
    k1=k1+inc+inc
    kk=(k1-kspan)/2+jc
    if(kk.le.jc+jc) goto 60
    goto 40

90  s1=float((kk-1)/jc)*dr*rad
    c1=cos(s1)
    s1=sin(s1)
    mm=min0(k1/2,mm+klim)
    goto 80

    ! Factor 3
100 k1=kk+kspan
    k2=k1+kspan
    ak=a(kk)
    bk=b(kk)
    aj=a(k1)+a(k2)
    bj=b(k1)+b(k2)
    a(kk)=ak+aj
    b(kk)=bk+bj
    ak=-0.5*aj+ak
    bk=-0.5*bj+bk
    aj=(a(k1)-a(k2))*s120
    bj=(b(k1)-b(k2))*s120
    a(k1)=ak-bj
    b(k1)=bk+aj
    a(k2)=ak+bj
    b(k2)=bk-aj
    kk=k2+kspan
    if(kk.lt.nn) goto 100
    kk=kk-nn
    if(kk.le.kspan) goto 100
    goto 290

    ! Factor 4
110 if(nfac(i).ne.4) goto 230
    kspnn=kspan
    kspan=kspan/4

120 c1=1.0
    s1=0.0
    mm=min0(kspan,klim)
    goto 150

130 c2=c1-(cd*c1+sd*s1)
    s1=(sd*c1-cd*s1)+s1
    c1=0.5/(c2**2+s1**2)+0.5
    s1=c1*s1
    c1=c1*c2

140 c2=c1**2-s1**2
    s2=c1*s1*2.0
    c3=c2*c1-s2*s1
    s3=c2*s1+s2*c1

150 k1=kk+kspan
    k2=k1+kspan
    k3=k2+kspan
    akp=a(kk)+a(k2)
    akm=a(kk)-a(k2)
    ajp=a(k1)+a(k3)
    ajm=a(k1)-a(k3)
    a(kk)=akp+ajp
    ajp=akp-ajp
    bkp=b(kk)+b(k2)
    bkm=b(kk)-b(k2)
    bjp=b(k1)+b(k3)
    bjm=b(k1)-b(k3)
    b(kk)=bkp+bjp
    bjp=bkp-bjp
    if(isn.lt.0) goto 180
    akp=akm-bjm
    akm=akm+bjm
    bkp=bkm+ajm
    bkm=bkm-ajm
    if(s1.eq.0.0) goto 190

160 a(k1)=akp*c1-bkp*s1
    b(k1)=akp*s1+bkp*c1
    a(k2)=ajp*c2-bjp*s2
    b(k2)=ajp*s2+bjp*c2
    a(k3)=akm*c3-bkm*s3
    b(k3)=akm*s3+bkm*c3
    kk=k3+kspan
    if(kk.le.nt) goto 150

170 kk=kk-nt+jc
    if(kk.le.mm) goto 130
    if(kk.lt.kspan) goto 200
    kk=kk-kspan+inc
    if(kk.le.jc) goto 120
    if(kspan.eq.jc) goto 350
    goto 40

180 akp=akm+bjm
    akm=akm-bjm
    bkp=bkm-ajm
    bkm=bkm+ajm
    if(s1.ne. 0.0) goto 160

190 a(k1)=akp
    b(k1)=bkp
    a(k2)=ajp
    b(k2)=bjp
    a(k3)=akm
    b(k3)=bkm
    kk=k3+kspan
    if(kk.le.nt) goto 150
    goto 170

200 s1=float((kk-1)/jc)*dr*rad
    c1=cos(s1)
    s1=sin(s1)
    mm=min0(kspan,mm+klim)
    goto 140

    ! Transform for factor of 5 (opcional code)
210 c2=c72**2-s72**2
    s2=2.0*c72*s72

220 k1=kk+kspan
    k2=k1+kspan
    k3=k2+kspan
    k4=k3+kspan
    akp=a(k1)+a(k4)
    akm=a(k1)-a(k4)
    bkp=b(k1)+b(k4)
    bkm=b(k1)-b(k4)
    ajp=a(k2)+a(k3)
    ajm=a(k2)-a(k3)
    bjp=b(k2)+b(k3)
    bjm=b(k2)-b(k3)
    aa=a(kk)
    bb=b(kk)
    a(kk)=aa+akp+ajp
    b(kk)=bb+bkp+bjp
    ak=akp*c72+ajp*c2+aa
    bk=bkp*c72+bjp*c2+bb
    aj=akm*s72+ajm*s2
    bj=bkm*s72+bjm*s2
    a(k1)=ak-bj
    a(k4)=ak+bj
    b(k1)=bk+aj
    b(k4)=bk-aj
    ak=akp*c2+ajp*c72+aa
    bk=bkp*c2+bjp*c72+bb
    aj=akm*s2-ajm*s72
    bj=bkm*s2-bjm*s72
    a(k2)=ak-bj
    a(k3)=ak+bj
    b(k2)=bk+aj
    b(k3)=bk-aj
    kk=k4+kspan
    if(kk.lt.nn) goto 220
    kk=kk-nn
    if(kk.le.kspan) goto 220
    goto 290

    ! Transform for odd factors
230 k=nfac(i)
    kspnn=kspan
    kspan=kspan/k
    if(k.eq.3) go to 100
    if(k.eq.5) go to 210
    if(k.eq.jf) go to 250
    jf=k
    s1=rad/(float(k)/8.0)
    c1=cos(s1)
    s1=sin(s1)
    ck(jf)=1.0
    sk(jf)=0.0
    j=1

240 ck(j)=ck(k)*c1+sk(k)*s1
    sk(j)=ck(k)*s1-sk(k)*c1
    k=k-1
    ck(k)=ck(j)
    sk(k)=-sk(j)
    j=j+1
    if(j.lt.k) goto 240

250 k1=kk
    k2=kk+kspnn
    aa=a(kk)
    bb=b(kk)
    ak=aa
    bk=bb
    j=1
    k1=k1+kspan

260 k2=k2-kspan
    j=j+1
    at(j)=a(k1)+a(k2)
    ak=at(j)+ak
    bt(j)=b(k1)+b(k2)
    bk=bt(j)+bk
    j=j+1
    at(j)=a(k1)-a(k2)
    bt(j)=b(k1)-b(k2)
    k1=k1+kspan
    if(k1.lt.k2) goto 260
    a(kk)=ak
    b(kk)=bk
    k1=kk
    k2=kk+kspnn
    j=1

270 k1=k1+kspan
    k2=k2-kspan
    jj=j
    ak=aa
    bk=bb
    aj=0.0
    bj=0.0
    k=1

280 k=k+1
    ak=at(k)*ck(jj)+ak
    bk=bt(k)*ck(jj)+bk
    k=k+1
    aj=at(k)*sk(jj)+aj
    bj=bt(k)*sk(jj)+bj
    jj=jj+j
    if(jj.gt.jf) jj=jj-jf
    if(k.lt.jf) goto 280
    k=jf-j
    a(k1)=ak-bj
    b(k1)=bk+aj
    a(k2)=ak+bj
    b(k2)=bk-aj
    j=j+1
    if(j.lt.k) goto 270
    kk=kk+kspnn
    if(kk.le.nn) goto 250
    kk=kk-nn
    if(kk.le.kspan) goto 250

290 if(i.eq.m) goto 350
    kk=jc+1

300 c2=1.0-cd
    s1=sd
    mm=min0(kspan,klim)
    goto 320

310 c2=c1-(cd*c1+sd*s1)
    s1=s1+(sd*c1-cd*s1)
    c1=0.5/(c2**2+s1**2)+0.5
    s1=c1*s1
    c2=c1*c2

320 c1=c2
    s2=s1
    kk=kk+kspan

330 ak=a(kk)
    a(kk)=c2*ak-s2*b(kk)
    b(kk)=s2*ak+c2*b(kk)
    kk=kk+kspnn
    if(kk.le.nt) goto 330
    ak=s1*s2
    s2=s1*c2+c1*s2
    c2=c1*c2-ak
    kk=kk-nt+kspan
    if(kk.le.kspnn) goto 330
    kk=kk-kspnn+jc
    if(kk.le.mm) goto 310
    if(kk.lt.kspan) goto 340
    kk=kk-kspan+jc+inc
    if(kk.le.jc+jc) goto 300
    goto 40

340 s1=float((kk-1)/jc)*dr*rad
    c2=cos(s1)
    s1=sin(s1)
    mm=min0(kspan,mm+klim)
    goto 320

350 np(1)=ks
    if(kt.eq.0) goto 440
    k=kt+kt+1
    if(m.lt.k) k=k-1
    j=1
    np(k+1)=jc

360 np(j+1)=np(j)/nfac(j)
    np(k)=np(k+1)*nfac(j)
    j=j+1
    k=k-1
    if(j.lt.k) goto 360
    k3=np(k+1)
    kspan=np(2)
    kk=jc+1
    k2=kspan+1
    j=1
    if(n.ne.ntot) goto 400

    ! Permutation for single variate transform (opcional code)
370 ak=a(kk)
    a(kk)=a(k2)
    a(k2)=ak
    bk=b(kk)
    b(kk)=b(k2)
    b(k2)=bk
    kk=kk+inc
    k2=kspan+k2
    if(k2.lt.ks) goto 370

380 k2=k2-np(j)
    j=j+1
    k2=np(j+1)+k2
    if(k2.gt.np(j)) goto 380
    j=1

390 if(kk.lt.k2) goto 370
    kk=kk+inc
    k2=kspan+k2
    if(k2.lt.ks) goto 390
    if(kk.lt.ks) goto 380
    jc=k3
    goto 440

    ! Permutation for multivariate transform
400 k=kk+jc

410 ak=a(kk)
    a(kk)=a(k2)
    a(k2)=ak
    bk=b(kk)
    b(kk)=b(k2)
    b(k2)=bk
    kk=kk+inc
    k2=k2+inc
    if(kk.lt.k) goto 410
    kk=kk+ks-jc
    k2=k2+ks-jc
    if(kk.lt.nt) goto 400
    k2=k2-nt+kspan
    kk=kk-nt+jc
    if(k2.lt.ks) goto 400

420 k2=k2-np(j)
    j=j+1
    k2=np(j+1)+k2
    if(k2.gt.np(j)) goto 420
    j=1

430 if(kk.lt.k2) goto 400
    kk=kk+jc
    k2=kspan+k2
    if(k2.lt.ks) goto 430
    if(kk.lt.ks) goto 420
    jc=k3

440 if(2*kt+1.ge.m) return
    kspnn=np(kt+1)

    j=m-kt
    nfac(j+1)=1

450 nfac(j)=nfac(j)*nfac(j+1)
    j=j-1
    if(j.ne.kt) goto 450
    kt=kt+1
    nn=nfac(kt)-1
    jj=0
    j=0
    goto 480

460 jj=jj-k2
    k2=kk
    k=k+1
    kk=nfac(k)

470 jj=kk+jj
    if(jj.ge.k2) goto 460
    np(j)=jj

480 k2=nfac(kt)
    k=kt+1
    kk=nfac(k)
    j=j+1
    if(j.le.nn) goto 470

    j=0
    goto 500

490 k=kk

    kk=np(k)
    np(k)=-kk
    if(kk.ne.j) goto 490
    k3=kk

500 j=j+1
    kk=np(j)
    if(kk.lt.0) goto 500
    if(kk.ne.j) goto 490
    np(j)=-j
    if(j.ne.nn) goto 500
    maxf=inc*maxf

    goto 570

510 j=j-1
    if(np(j).lt.0) goto 510
    jj=jc

520 kspan=jj
    if(jj.gt.maxf) kspan=maxf
    jj=jj-kspan
    k=np(j)
    kk=jc*k+i+jj
    k1=kk+kspan
    k2=0

530 k2=k2+1
    at(k2)=a(k1)
    bt(k2)=b(k1)
    k1=k1-inc
    if(k1.ne.kk) goto 530

540 k1=kk+kspan
    k2=k1-jc*(k+np(k))
    k=-np(k)

550 a(k1)=a(k2)
    b(k1)=b(k2)
    k1=k1-inc
    k2=k2-inc
    if(k1.ne.kk) goto 550
    kk=k2
    if(k.ne.j) goto 540
    k1=kk+kspan
    k2=0

560 k2=k2+1
    a(k1)=at(k2)
    b(k1)=bt(k2)
    k1=k1-inc
    if(k1.ne.kk) goto 560
    if(jj.ne.0) goto 520
    if(j.ne.1) goto 510

570 j=k3+1
    nt=nt-kspnn
    i=nt-inc+1
    if(nt.ge.0) goto 510

    return

  end subroutine fftmx
  ! =====================================================================================
  subroutine bkdt

    double precision dstak
    integer istak(6000),isize(5),lout,lnow,lused,lmax,lbook
    common/cstak/dstak(3000)
    equivalence (dstak(1),istak(1))
    equivalence (istak(1),lout)
    equivalence (istak(2),lnow)
    equivalence (istak(3),lused)
    equivalence (istak(4),lmax)
    equivalence (istak(5),lbook)
    equivalence (istak(6),isize(1))

    isize(1)=1
    isize(2)=1
    isize(3)=1
    isize(4)=2
    isize(5)=2
    lout=0
    lnow=10
    lused=10
    lmax=6000
    lbook=10

    return

  end subroutine bkdt
  ! =====================================================================================
  subroutine realt(a,b,nseg,n,nspn,isn)

    real a(*),b(*)

    inc=iabs(isn)
    ks=iabs(nspn)*inc
    nf=iabs(n)
    ns=ks*nf
    nt=iabs(ns*nseg)
    if(isn*nt.ne.0) goto 10
    write(0,9999) nseg,n,nspn,isn

    return

10  jc=ks
    k2=iabs(ks*nseg)-inc
    kd=ns
    nh=ns/2+1
    nn=nt-inc
    nt=nt+1
    kk=1
    rad=atan(1.0)
    dr=-4.0/float(nf)
    cd=2.0*sin(0.5*dr*rad)**2
    sd=sin(dr*rad)

    lim=32
    klim=lim*ks
    mm=min0(nh,klim)
    sn=0.0
    if(isn.gt.0) goto 70

20  aa=a(kk)
    ba=b(kk)
    b(kk)=0
    a(kk)=aa+ba
    a(nt)=aa-ba
    b(nt)=0
    nt=nt+jc
    kk=kk+ns
    if(kk.le.nn) goto 20
    nt=nt-k2
    kk=kk-nn
    if(kk.le.jc) goto 20
    cn=1.0

30  if(nf.eq.1) return

40  aa=cn-(cd*cn+sd*sn)
    sn=(sd*cn-cd*sn)+sn

    cn=0.5/(aa**2+sn**2)+0.5
    sn=cn*sn
    cn=cn*aa
50  jc=jc+ks
    kd=kd-ks-ks

60  k2=kk+kd
    aa=a(kk)+a(k2)
    ab=a(kk)-a(k2)
    ba=b(kk)+b(k2)
    bb=b(kk)-b(k2)
    re=cn*ba+sn*ab
    em=sn*ba-cn*ab
    b(k2)=(em-bb)*0.5
    b(kk)=(em+bb)*0.5
    a(k2)=(aa-re)*0.5
    a(kk)=(aa+re)*0.5
    kk=kk+ns
    if(kk.le.nn) goto 60
    kk=kk-nn
    if(kk.le.jc) goto 60
    if(kk.le.mm) goto 40
    if(kk.gt.nh) return
    sn=float(jc/ks)*dr*rad
    cn=cos(sn)
    if(isn.gt.0) cn=-cn
    sn=sin(sn)
    mm=min0(nh,mm+klim)
    goto 50

70  aa=a(kk)
    ba=a(nt)
    a(kk)=(aa+ba)*0.5
    b(kk)=(aa-ba)*0.5
    nt=nt+jc
    kk=kk+ns
    if(kk.le.nn) goto 70
    nt=nt-k2
    kk=kk-nn
    if(kk.le.jc) goto 70
    cn=-1.0
    sd=-sd

    goto 30

9999 format(' Error - one zero parameter in function realt ',3i10,i9)

  end subroutine realt
  ! =====================================================================================
  integer function istkgt(nitems,itype)

    common/cstak/dstak(3000)
    double precision dstak
    integer istak(6000)
    integer isize(5),lout,lnow,lused,lmax,lbook
    equivalence (dstak(1),istak(1))
    equivalence (istak(1),lout)
    equivalence (istak(2),lnow)
    equivalence (istak(3),lused)
    equivalence (istak(4),lmax)
    equivalence (istak(5),lbook)
    equivalence (istak(6),isize(1))

    istkgt=(lnow*isize(2)-1)/isize(itype)+2
    i=((istkgt-1+nitems)*isize(itype)-1)/isize(2)+3
    if(i.gt.lmax) go to 10
    istak(i-1)=itype
    istak(i)=lnow
    lout=lout+1
    lnow=i
    lused=max0(lused,lnow)

    return

10  write(0,9999) i
    write(0,9998) (istak(j),j=1,10),istak(lnow-1),istak(lnow)

9998    format(12i6)
9999    format(' Overflow of common array istak - need ',i10)

    stop

  end function istkgt
  ! =====================================================================================
  subroutine istkrl(k)

    common/cstak/dstak(3000)
    double precision dstak
    integer istak(6000),lout,lnow,lmax,lbook,lused
    equivalence(dstak(1),istak(1))
    equivalence(istak(1),lout)
    equivalence(istak(2),lnow)
    equivalence(istak(3),lused)
    equivalence(istak(4),lmax)
    equivalence(istak(5),lbook)

    in=k

    if(lbook.le.lnow .and. lnow.le.lused .and. lused.le.lmax) goto 10

    write(0,9999)
    write(0,9997) (istak(j),j=1,10),istak(lnow-1),istak(lnow)

10  if(in.le.0) return

    if(lbook.gt.istak(lnow) .or. istak(lnow).ge.lnow-1) go to 20
    lout=lout-1
    lnow=istak(lnow)
    in=in-1

    goto 10

20  write(0,9998)

9997    format(1x,12i6)
9998    format(' Warning: pointer at istak(lnow) overwritten',11x,' allocation not completed')
9999    format(' Warning: istak(2),istak(3),istak(4) or istak(5) hit')

    return
  end subroutine istkrl
  ! =====================================================================================

end module fastfourrier
