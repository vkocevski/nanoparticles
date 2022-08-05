! Program that surrounds Si unit cell with its copies and cuts away a sphere with given cut-off
! Removes Si atoms with 1 neighbor from the nanoparticle, and hydrogenates the surface

  implicit none

  double precision, allocatable :: r_np(:,:)    ! coordinates of cut nanoparticle
  double precision, allocatable :: r(:,:)       ! coordinates of trimed nanoparticle
  double precision, allocatable :: rr(:,:)      ! repeated r
  double precision                 rnn(3,100)   ! positions of nearest neighbors
  integer                          n            ! number of atoms
  integer                          nnn          ! number of nearest neighbors
  character, allocatable        :: elem*2(:)    ! element in nanoparticle

  double precision c_o, r_s(3), r0(3), a, si_uc(3,8), dis, cutoff
  integer h, k, l, i, j, nat
  character el*2

  external add1, add2, add3, add4

  si_uc = reshape(  (/ 0.50, 0.00, 0.00, &
                       0.00, 0.50, 0.00, &
                       0.75, 0.25, 0.25, &
                       0.25, 0.75, 0.25, &
                       0.00, 0.00, 0.50, &
                       0.50, 0.50, 0.50, &
                       0.25, 0.25, 0.75, &
                       0.75, 0.75, 0.75 /), shape(si_uc) )  ! Si unit cell in fractional coordinates
					   

  print*, 'Cutt-off for NP: '
  read*, c_o

  n=8    ! read number of atoms
  a=5.47 ! read lattice parameter
  cutoff = 2.6  ! cutoff radius for Si nearest neighbors
  el='Si'

  allocate(r(3,n*27),elem=(n*27))
  nat=0
  do i=1,n
    read(1,*) el, (r0(j),j=1,3)  ! read element type and coordiates
  ! making a 7x7x7 supercell
    do h=-3,3                              
      r_s(1)=(r0(1)+h)*a   
      do k=-3,3                              
        r_s(2)=(r0(2)+k)*a
        do l=-3,3                              
          r_s(3)=(r0(3)+l)*a
          ! cutting out a nanoparticle with radius of c_o
          if (r_s(3) < c_o .and. r_s(3) > dble(-1)*c_o) then
            if (r_s(3) < c_o .and. r_s(3) > dble(-1)*c_o) then
              if (r_s(3) < c_o .and. r_s(3) > dble(-1)*c_o) then
                nat=nat+1
                r_np(1,nat)=r_s(1)
                r_np(2,nat)=r_s(2)
                r_np(3,nat)=r_s(3)
                elem(nat)=el
              end if   
            end if 
          end if 
        end do                              
      end do                          
    end do                                      
  end do

! removing atoms with only 1 neareast neighbor
  n=0
  do i=1,nat
    rr(:,i)=r_np(:,i)
  enddo
  do i=1,nat
    nnn = 0
    do k=1,nat
      dis = sqrt(dot_product(r(:,i)-rr(:,k),r(:,i)-rr(:,k)))  ! distance between atoms
      if (dis <= cutoff .and. dis > 0.01) nnn=nnn + 1
    enddo
    if (nnn>1.5) then
      n=n+1
      r(:,i)=r_np(:,i)
    endif
  enddo

! write nanopartile coordinates
  write(10,'(i5)') n
  write(10,*)
  do i=1,n
    rr(:,i)=r(:,i)
    write(10,'(a4,3f12.6)') el(i), r(:,i)
  enddo

  do i=1,nat
    nnn = 0
    do k=1,nat
      dis = sqrt(dot_product(r(:,i)-rr(:,k),r(:,i)-rr(:,k)))  ! distance between atoms
      if (dis <= cutoff .and. dis > 0.01) then
        nnn=nnn + 1
        rnn(:,nnn)=rr(:,k)     ! store atom coordinate
      endif
    enddo
    write(11,'(i4,i4)') i, nnn   ! print number of nearest neighbor for debuging purpose
! write main atom
    write(20,'(a4,3f12.6)') 'Si  ', r(:,i)
! write hydrogen if nearest neighbor is less than 4
    if(nnn==3) call add1(r(1,i),t(i),rnn)
    if(nnn==2) call add2(r(1,i),t(i),rnn)
    if(nnn==1) call add3(r(1,i),t(i),rnn)
    if(nnn==0) call add4(r(1,i),t(i),rnn)
  end do

  end

! ---------------------------------------
! three neighbors around, one to be added
subroutine add1(r0,t0,r)

  implicit none
  
  double precision r0(3),r(3,3)
  integer          t0, t(3)
  
! local variables
  double precision u(3), unorm, h(3)

! calculate the a+b+c vector and normalize it
  u(:) = (r(:,1)-r0(:)) + (r(:,2)-r0(:)) + (r(:,3)-r0(:))
  unorm = sqrt(dot_product(u,u))
  u(:) = u(:) / unorm

! print out the H position in xyz-style format
  h = 1.5  ! Si-H distance
  pointh(:) = r0(:) - h*u(:)
  write(20,'(a4,3f12.6)') 'H ', pointh(:)

end subroutine


! -------------------------------------
! two neighbors around, two to be added
subroutine add2(r0,t0,r)

  implicit none
  
  double precision r0(3),r(3,2)
  integer          t0, t(3)
  real h
  double precision a(3), b(3), c(3), d(3), cnorm, dnorm
  double precision h1(3), h2(3), pointh1(3), pointh2(3)

  a(:) = (r(:,1)-r0(:))
  b(:) = (r(:,2)-r0(:))
  d(:) = (a(:)+b(:))/2
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
  cnorm  = sqrt(dot_product(c,c))
  dnorm = sqrt(dot_product(d,d))
  c(:) = c(:) / cnorm
  d(:) = d(:) / dnorm
  h = 1.5  ! Si-H distance
  h1(:) = -h*cos(0.95547)*d(:) + h*sin(0.95547)*c(:)
  h2(:) = -h*cos(0.95547)*d(:) - h*sin(0.95547)*c(:)
  pointh1(:) = r0(:) + h1(:)
  pointh2(:) = r0(:) + h2(:)

  write(20,'(a4,3f12.6)') 'H ', pointh1(:)
  write(20,'(a4,3f12.6)') 'H ', pointh2(:)

end subroutine


! ---------------------------------------
! one neighbors around, three to be added
subroutine add3(r0,t0,r)

  implicit none
  
  double precision r0(3),r(3,1)
  integer          t0, t(3)
  real h
  double precision a(3), b(3), c(3), d(3), u(3)
  double precision anorm, bnorm, dnorm 
  double precision pointh1(3), pointh2(3), pointh3(3), pointd(3)

  h = 1.5  ! Si-H distance
  dnorm = h*cos(1.2322)
  a(:) = (r0(:)-r(:,1))
  anorm  = sqrt(dot_product(a,a))
  a(:) = a(:) / anorm
  d(:) = (anorm + dnorm)*a(:)
  if (a(1) == 0) then
  if (a(2) == 0) then
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  else if (a(3) == 0) then
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  else   
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  end if
  else if (a(2) == 0) then
  if (a(1) == 0) then
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  else if (a(3) == 0) then  
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  else   
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  end if
  else
  u(1) = 0.0
  u(2) = 0.0
  u(3) = 1.0
  end if
  b(1) = a(2)*u(3) - a(3)*u(2)
  b(2) = a(3)*u(1) - a(1)*u(3)
  b(3) = a(1)*u(2) - a(2)*u(1)
  bnorm  = sqrt(dot_product(b,b))
  b(:) = b(:) / bnorm
  pointd(:) = (d(:) + r(:,1))
  pointh1(:) = pointd(:) + h*sin(1.2322)*b(:)
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
  pointh2(:) =  pointd(:) - h*sin(1.2322)*sqrt(3.0)*c(:)/2 - h*sin(1.2322)*b(:)/2.0
  pointh3(:) =  pointd(:) + h*sin(1.2322)*sqrt(3.0)*c(:)/2 - h*sin(1.2322)*b(:)/2.0

  write(20,'(a4,3f12.6)') 'H ', pointh1(:)
  write(20,'(a4,3f12.6)') 'H ', pointh2(:)
  write(20,'(a4,3f12.6)') 'H ', pointh3(:)

end subroutine


! ---------------------------------------
! totally sad orphan atom
subroutine add4(r0,t0,r)

  implicit none
  
  double precision r0(3),r(3,3)
  integer          t0, t(3)
  real h
  double precision a(3), b(3), c(3), d(3), u(3)
  double precision anorm, bnorm, dnorm, r0norm
  double precision pointh1(3), pointh2(3), pointh3(3), pointh4(3), pointd(3)

  r0norm = sqrt(dot_product(r0,r0))
  r0(:) = r0(:)/r0norm
  h = 1.5  ! Si-H distance
  pointh4(:) = (h + r0norm) * r0(:)
  dnorm = h*cos(1.2322)
  a(:) = (r0(:)-r(:,1))
  anorm = sqrt(dot_product(a,a))
  a(:) = a(:) / anorm
  d(:) = (anorm + dnorm)*a(:)
  if (a(1) == 0) then
  if (a(2) == 0) then
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  else if (a(3) == 0) then
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  else   
     u(1) = 1.0
     u(2) = 0.0
     u(3) = 0.0
  end if
  else if (a(2) == 0) then
  if (a(1) == 0) then
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  else if (a(3) == 0) then  
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  else   
     u(1) = 0.0
     u(2) = 1.0
     u(3) = 0.0
  end if
  else
  u(1) = 0.0
  u(2) = 0.0
  u(3) = 1.0
  end if
  b(1) = a(2)*u(3) - a(3)*u(2)
  b(2) = a(3)*u(1) - a(1)*u(3)
  b(3) = a(1)*u(2) - a(2)*u(1)
  bnorm  = sqrt(dot_product(b,b))
  b(:) = b(:) / bnorm
  pointd(:) = (d(:) + r(:,1))
  pointh1(:) = pointd(:) + h*sin(1.2322)*b(:)
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
  pointh2(:) =  pointd(:) - h*sin(1.2322)*sqrt(3.0)*c(:)/2 - h*sin(1.2322)*b(:)/2.0
  pointh3(:) =  pointd(:) + h*sin(1.2322)*sqrt(3.0)*c(:)/2 - h*sin(1.2322)*b(:)/2.0

  write(20,'(a4,3f12.6)') 'H ', pointh1(:)
  write(20,'(a4,3f12.6)') 'H ', pointh2(:)
  write(20,'(a4,3f12.6)') 'H ', pointh3(:)
  write(20,'(a4,3f12.6)') 'H ', pointh4(:)

end subroutine