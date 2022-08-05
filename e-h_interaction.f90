! code that reads gaussian cube file for HOMO and LUMO wavefunctions (WF)
! integrates the WFs radially and calculates their interaction energy

  implicit none
  double precision, allocatable :: hwf(:,:,:) ! values for the density
  double precision, allocatable :: lwf(:,:,:) ! values for the density
  double precision, allocatable :: hamp(:)    ! sum of homo intensities (squared amplitudes) within the same shell
  double precision, allocatable :: lamp(:)    ! sum of homo intensities (squared amplitudes) within the same shell
  double precision r0_br, r0, r_c(3), r_p(3), d, ch, dx, dy, dz, incr, pi, bhr, norm_f(2)
  double precision hpsi, lpsi, theta, r
  integer n, n1, n2, n3, i, j, k, nat, m, at, limit, pos, np, norb
  
  pi=3.141592653589d0
  bhr=0.52917721092d0   ! bohr radius = 0.529... angstroms
!  print*,'normalization factors'
!  read*, norm_f(1), norm_f(2)
  read(210,*)
  read(210,*)
  read(210,*) n, r0_br
  read(210,*) n1, dx, dy, dz
  read(210,*) n2, dx, dy, dz
  read(210,*) n3, dx, dy, dz
  ! converting bohr radius into angstroms
  r0=r0_br*bhr
  j=nint(r0)
  r0=dble(j)
  dx=dz*bhr
  dy=dx
  dz=dx
  i=nint(dx*1000)
  incr=dble(i)/1000.
  read(220,*)
  read(220,*)
  read(220,*)
  read(220,*)
  read(220,*)
  read(220,*)

  limit=nint(abs(r0)/incr)
  read(210,*) ch, at, (r_c(j),j=1,3)  ! reading the center of the NC
  read(220,*)
  norb=15
  do i=2,n
    read(210,*) ch
    if (ch==48 .or. ch==30) then
      norb=norb+15
    elseif(ch==16) then
      norb=norb+13
    else
      norb=norb+5
    endif
    read(220,*) ch
  enddo

  allocate(hwf(n3,n2,n1),lwf(n3,n2,n1),hamp(0:limit),lamp(0:limit))
  hamp(:)=0    ! the intensities variable
  lamp(:)=0    ! the intensities variable
  theta=0
  ! reading the WF amplitudes
  do i=1,n1
    do j=1,n2
      read(210,'(6f13.5)') (hwf(k,j,i),k=1,n3) ! reads through the whole set of points along the third axes
      read(220,'(6f13.5)') (lwf(k,j,i),k=1,n3) ! reads through the whole set of points along the third axes
      do m=1,n3
        r_p(1)=dble(i-1)*incr+r0
        r_p(2)=dble(j-1)*incr+r0
        r_p(3)=dble(m-1)*incr+r0
        d=sqrt(dot_product((r_p(:)-r_c(:)),(r_p(:)-r_c(:))))   ! calculating the distance from the center of the NC and the WF point
        pos=nint(d/incr)
        if (pos<=limit .and. pos>=0.) then
!          hamp(pos)=hamp(pos)+hwf(m,j,i)*hwf(m,j,i)*incr**3  ! summing the intensities at a given distance
          lamp(pos)=lamp(pos)+lwf(m,j,i)*lwf(m,j,i)*incr**3  ! summing the intensities at a given distance
        endif
      enddo
    enddo
  enddo
  np=1
!  hpsi=hamp(0)
!  lpsi=lamp(0)
  do i=1,limit
    np=np+1
!    lpsi=lpsi+lamp(i)
!    hpsi=hpsi+hamp(i)
    if (lamp(i)==0) goto 10
  enddo
10  r=dble(np-1)*incr
  hpsi=0.
  lpsi=0
  do i=1,n1
    do j=1,n2
      do m=1,n3
        r_p(1)=dble(i-1)*incr+r0
        r_p(2)=dble(j-1)*incr+r0
        r_p(3)=dble(m-1)*incr+r0
        d=sqrt(dot_product((r_p(:)-r_c(:)),(r_p(:)-r_c(:))))   ! calculating the distance from the center of the NC and the WF point
        if (d<(r+0.01) .and. d>(-0.1)) then
          hpsi=hpsi+hwf(m,j,i)*hwf(m,j,i)*incr**3  ! calculating the integral over HOMO states
          lpsi=lpsi+lwf(m,j,i)*lwf(m,j,i)*incr**3  ! calculating the integral over HOMO states
        endif
      enddo
    enddo
  enddo
  do i=1,n1
    do j=1,n2
      do m=1,n3
        theta=theta+(hwf(m,j,i)*hwf(m,j,i)/hpsi)*(lwf(m,j,i)*lwf(m,j,i)/lpsi)  ! calculating the overlap integral theta
      enddo
    enddo
  enddo
  print('(3f13.8)'), hpsi, lpsi, theta

end