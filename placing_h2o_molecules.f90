! Program that places H2O molecules at specific O positions
! Direction of the H2O molecule is given by the O and Na atoms

  implicit none
  double precision, allocatable :: rna(:,:)      ! coordinates of the Na atoms
  double precision, allocatable :: ro(:,:)       ! coordinates of the O atoms
  double precision, allocatable :: point2(:,:,:) ! position of the Si and H atom
  double precision point1(3,2), vec1(3), vec2(3), rot_m(3,3)
  double precision mol(3,3), vec_mol(3,2), rot_mol(3)
  integer n, no, nna, i, j, k, o_to_h2o
  character el*4, elm(3)*4

! coordinates of an H2O molecule
  mol = reshape(  (/ 0.000000, 0.000000, 0.000000, &
                    -0.756491, 0.591035, 0.000000, &
                     0.756491, 0.591035, 0.000000  /), shape(mol) )
  elm = reshape( (/ 'O   ', 'H   ', 'H   '/), shape(elm) )
! setting the vectors from O atom (0,0,0) to H atoms
  do i=1,2
    vec_mol(:,i)=mol(:,i+1)
  enddo
! setting up the vector which will be used to calculate the rotational matrix
! this is the direction of the molecule
  vec1 = reshape( (/ 0.000000, 1.000000, 0.000000  /), shape(vec1) )

! reading the xyz file with the Na and O atoms
  read(20,*) n
  read(20,*)
  allocate(rna(3,n),ro(3,n),point2(3,2,n))
  nna=n/2
  no=n/2
  o_to_h2o=n/2
  do i=1,no
    read(20,*) el, (ro(j,i),j=1,3)
    point2(:,2,i)=ro(:,i)
  enddo
  do i=1,nna
    read(20,*) el, (rna(j,i),j=1,3)
    point2(:,1,i)=rna(:,i)
  enddo

! substituting the O atoms with H2O
  do i=1,o_to_h2o
    rot_m(:,:)=0.
    vec2(:)=point2(:,2,i)-point2(:,1,i)  ! vector between the O and the Na atom i.e. direction that the molecule should aquire
    call calc_rot_matrix(vec1,vec2,rot_m)
    mol(:,1)=point2(:,2,i)  ! placing the H2O atom at the position of the O atom
    do j=1,2
    ! rotating the vectors of the molecule to be in direction of the vec1 (from O atom to the middle of the H atoms)
      do k=1,3
        rot_mol(k)=dot_product(rot_m(:,k),vec_mol(:,j))
      enddo
      mol(:,j+1)=rot_mol(:)+mol(:,1) ! coordinates of the rotated atoms
    enddo
    do j=1,3
      write(30,'(a4,f14.9,2f20.9)') elm(j), mol(:,j)
    enddo
  enddo

end

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
  A(2,1)=dble(-1)*x(3)
  A(3,1)=x(2)
  A(1,2)=x(3)
  A(3,2)=dble(-1)*x(1)
  A(1,3)=dble(-1)*x(2)
  A(2,3)=x(1)

! identity matrix I
  I(:,:)=0.
  I(1,1)=1.
  I(2,2)=1.
  I(3,3)=1.

! calculating the rotational matrix
  rot_m=I+sin(theta)*A+(1-cos(theta))*matmul(A,A)

end subroutine
