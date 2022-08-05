! Program that replaces H atoms with C2H4NH2 molecules
! asking if only give H atoms should be replaced, all H atoms
! or partial replacement of H atoms (by removing molecules that are close to each other)

  implicit none

  double precision, allocatable :: rsi(:,:)      ! coordinates of the Si atoms
  double precision, allocatable :: rh(:,:)       ! coordinates of the H atoms
  double precision, allocatable :: rhnotc(:,:)   ! coordinates of the H atoms that are not substituted
  double precision, allocatable :: rhyesc(:,:)   ! coordinates of the H atoms that are substituted
  double precision, allocatable :: point2(:,:,:) ! position of the Si and H atom
  double precision point1(3,2), vec1(3), vec2(3), rot_m(3,3)
  double precision mol(3,9), r0(3), d, rc(3), vec_mol(3,8), rot_mol(3), temp(3)
  integer, allocatable :: phtos(:)
  integer n, nh, nsi, i, j, k, nnh, pos(4), same, h2c, hnotc, nhtos, found
  character el*4, elm(9)*4, sub_type*4

  print*, 'should all H atoms be substituted (a), only partially (p) or H atom at givem position (g)'
  read*, sub_type
  if (sub_type=='g') then
    print*, 'give the number and the position of H atoms which are being substituted'
    read*, nhtos
    allocate(phtos(nhtos))
    read*, (phtos(i),i=1,nhtos)
  endif
! coordinates of the relaxed molecule
  mol = reshape(  (/ 0.000000,  0.000000,  0.000000, &
                     0.189898,  0.984368, -0.559833, &
                     1.018924, -0.504845,  0.158412, &
                    -0.669389,  0.275489,  1.358042, &
                    -0.904974, -0.728513,  1.864309, &
                    -1.697426,  0.757024,  1.184163, &
                     0.098646,  1.106285,  2.297656, &
                     1.067541,  0.641616,  2.342362, &
                     0.315015,  2.010073,  1.757590  /), shape(mol) )
  elm = reshape( (/ 'C   ', 'H   ', 'H   ', 'C   ', 'H   ', 'H   ', 'N   ', 'H   ', 'H   ' /), shape(elm) )
! calculating the vector from the first C atom (0,0,0) to all other atoms
  do i=1,8
    vec_mol(:,i)=mol(:,i+1)-mol(:,1)
  enddo
! setting up the vector which will be used to calculate the rotational matrix
! this is the direction of the molecule
  vec1 = reshape( (/ 1.286161, 1.296682, 1.226543 /), shape(vec1) )

! reading the xyz file with the NC structure
  read(20,*) n
  read(20,*)
  allocate(rsi(3,n),rh(3,n),point2(3,2,n),rhnotc(3,n),rhyesc(3,n))
  nsi=0
  nh=0
  hnotc=0
  h2c=0
  do i=1,n
    read(20,*) el, (r0(j),j=1,3)
! saving the coordinates of Si and H atoms in separate arrays
    if (el=='Si') then
      nsi=nsi+1
      rsi(:,nsi)=r0(:)
    elseif (el=='H') then
      found=0
      if (sub_type=='g') then
        do j=1,nhtos
          if (phtos(j)==i) then
            h2c=h2c+1
            found=1
            rhyesc(:,h2c)=r0(:)   ! saving the coordinates of the H atom which is being substituted
            point2(:,1,h2c)=rsi(:,nsi)  ! the first point (Si) of the vector from Si atom to the H atom
            goto 30
          endif
        enddo
30      if (found==0) then
          hnotc=hnotc+1 
          rhnotc(:,hnotc)=r0(:)   ! saving the coordinates of not substituted H atoms
        endif
      else
        nh=nh+1
        rh(:,nh)=r0(:)
      endif
    endif
  enddo

  if (sub_type=='g') then
    do i=1,h2c
      rc(:)=0.   ! coordinates of the C atom from C2H4NH2
      call calc_rc(rc(:),point2(:,1,i),rhyesc(:,i))
      point2(:,2,i)=rc(:)  ! the second point (H) of the vector from Si atom to the H atom
    enddo
  elseif (sub_type=='a') then
    h2c=nh
    do i=1,nh
      rc(:)=0.   ! coordinates of the C atom from C2H4NH2
      do j=1,nsi
      ! calculating the distance between the H and Si atom
        d=sqrt(dot_product((rsi(:,j)-rh(:,i)),(rsi(:,j)-rh(:,i))))
      ! considering only the Si atom bonded to H
        if (d<1.8) then
          point2(:,1,i)=rsi(:,j)  ! the first point (Si) of the vector from Si atom to the H atom
          call calc_rc(rc(:),rsi(:,j),rh(:,i))
          point2(:,2,i)=rc(:)  ! the second point (H) of the vector from Si atom to the H atom
          goto 10
        endif
      enddo
10    hnotc=0
    enddo
  elseif (sub_type=='p') then
    hnotc=0
    h2c=0
    do i=1,nsi
      rc(:)=0.   ! coordinates of the C atom from C2H4NH2
      nnh=0
      do j=1,nh
      ! calculating the distance between the H and Si atom
        d=sqrt(dot_product((rh(:,j)-rsi(:,i)),(rh(:,j)-rsi(:,i))))
      ! considering only the H atom bonded to Si
        if (d<1.8) then
          nnh=nnh+1    ! counting the number of H atoms bonded to the Si atom
          pos(nnh)=j   ! saving the position of atoms sharing the same Si atom
        endif
      enddo
    ! if the Si atom has 2 or more NN H atoms
      if (nnh>1.5) then
        do k=1,nnh
          same=0
          do j=1,nh
          ! calculating the distance between the H atoms, NN to Si, with the other H atoms
            d=sqrt(dot_product((rh(:,j)-rh(:,pos(k))),(rh(:,j)-rh(:,pos(k)))))
          ! checking if there are H atoms in the vicinity
            if (d<1.9 .and. d>0.5) then
              same=1
            endif
          enddo
        ! if there is NOT a H atom in vicinity, the particular H atom is substituted
          if (same==0) then
            h2c=h2c+1   ! counting the number of H atoms to be substituted
            point2(:,1,h2c)=rsi(:,i)  ! the first point (Si) of the vector from Si atom to the H atom
            rhyesc(:,h2c)=rh(:,pos(k))
            call calc_rc(rc(:),rsi(:,i),rh(:,pos(k)))
            point2(:,2,h2c)=rc(:)  ! the second point (H) of the vector from Si atom to the H atom      
        ! if there is a H atom in vicinity, the particular H atom is kept (not substituted)
          elseif (same==1) then
            hnotc=hnotc+1  ! counting the number of H atom not substituted
            rhnotc(:,hnotc)=rh(:,pos(k))
          endif
        enddo
    ! if the Si atom has only 1 NN H atom, then this atom is substituted
      elseif (nnh>0.5) then
        h2c=h2c+1
        point2(:,1,h2c)=rsi(:,i)  ! the first point (Si) of the vector from Si atom to the H atom
        rhyesc(:,h2c)=rh(:,pos(1))
        call calc_rc(rc(:),rsi(:,i),rh(:,pos(1)))
        point2(:,2,h2c)=rc(:)  ! the second point (H) of the vector from Si atom to the H atom
      endif
    enddo
  endif

  print'(a9,f10.4,a2)', 'coverage:', dble(h2c)*100./dble(h2c+hnotc), ' %'
  write(30,'(i4)') nsi+h2c*9+hnotc
  write(30,*)
! writing the coordinates of the Si atoms
  do i=1,nsi
    write(30,'(a4,3f12.6)') 'Si  ', rsi(:,i)
  enddo
! writing the coordinates of the H atoms that are being substituted
  write(200,'(i5)') nh
  write(200,*)
  do i=1,h2c
    write(200,'(a4,3f12.6)') 'H   ', rhyesc(:,i)
  enddo
! writing the coordinates of the unsubstituted H atoms
  if (hnotc>0.5) then
    do i=1,hnotc
      write(30,'(a4,3f12.6)') 'H   ', rhnotc(:,i)
      write(200,'(a4,3f12.6)') 'He  ', rhnotc(:,i)
    enddo
  endif

! substituting the H atoms with C2H4NH2
  do i=1,h2c
    rot_m(:,:)=0.
    vec2(:)=point2(:,2,i)-point2(:,1,i)  ! vector between the Si and the H atom i.e. direction that the molecule should aquire
    call calc_rot_matrix(vec1,vec2,rot_m)
    mol(:,1)=point2(:,2,i)  ! placing the C atom at the position of the H atom
    do j=1,8
    ! rotating the vectors of the molecule to be in direction of the vec1 (from C atom to the middle of the H atoms)
      do k=1,3
        rot_mol(k)=dot_product(rot_m(:,k),vec_mol(:,j))
      enddo
      mol(:,j+1)=rot_mol(:)+mol(:,1) ! coordinates of the rotated atoms
    enddo
    do j=1,9
      write(30,'(a4,3f12.6)') elm(j), mol(:,j)
    enddo
  enddo

end

subroutine calc_rc(rc,r1,r2)
implicit none
  double precision, intent (in) :: r1(3),r2(3)
  double precision u(3), mod_u
  double precision, intent (inout) :: rc(3)
  integer i

  u=r2-r1   ! vector between the Si and H atom
  mod_u=sqrt(dot_product(u,u))  ! module of the vector
! calculating the coordinates of a point 1.905 A from the Si atom
! later the C atom will be placed on this point
  rc=r1+1.905*u/mod_u

end subroutine

! calculating the rotational matrix using the Rodrigues' rotation formula:
! R=e^{A*theta}=I+sin(theta)*A+(1-cos(theta))*A^2
subroutine calc_rot_matrix(vec1,vec2,rot_m)
implicit none
  double precision, intent (in) :: vec1(3), vec2(3)
  double precision, intent (inout) :: rot_m(3,3)
  double precision theta, x(3), A(3,3), mvec1, mvec2, mx, crossx(3), I(3,3)

  mvec2=sqrt(dot_product(vec2,vec2))  ! module of vec2
  mvec1=sqrt(dot_product(vec1,vec1))  ! module of vec1
! calculating a vector, x, perpendicular to vec1 and vec2
  ! cross product of vec1 and vec2
  crossx(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
  crossx(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
  crossx(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
  mx=sqrt(dot_product(crossx,crossx)) ! module of crossx
  x=crossx/mx  ! the vector perpendicular to vec1 and vec2

! calculating a angle, theta, between vec1 and vec2
  theta=acos(dot_product(vec1,vec2)/(mvec1*mvec2))

! assigning the values for the skew-symetric matrix corresponding to x
  A(:,:)=0.
  A(2,1)=(-1)*x(3)
  A(3,1)=x(2)
  A(1,2)=x(3)
  A(3,2)=(-1)*x(1)
  A(1,3)=(-1)*x(2)
  A(2,3)=x(1)

! identity matrix I
  I(:,:)=0.
  I(1,1)=1.
  I(2,2)=1.
  I(3,3)=1.

! calculating the rotational matrix
  rot_m=I+sin(theta)*A+(1-cos(theta))*matmul(A,A)

end subroutine
