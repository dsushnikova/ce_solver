module fast_direct

!##########################################################################



type, public :: node
    integer :: i
    type(node), pointer :: next => null()
end type
        
type, public :: row
  integer :: sz, nclose, nfar=0, ind, start
  real(8), pointer :: block_row(:)
  integer :: diag_block, close_block, far_block
  integer, pointer :: stencil(:, :)
  integer, pointer :: col_index(:)
  !integer, pointer :: block_col(:)
  integer, pointer :: clos(:)
  integer :: close_col, far_col=0
  integer, pointer :: far(:)
end type
!The code is very-simple: gnew = g - LL^{-1} fcur
!first we extract diagonal (triangular) block of the matrix L and then solve the linear system with it. 
!Then we multiply the remaining block by the resulting vector and then (?) redistribute it all around the matrix.
!How do we store such matrix? The size of the block, the rank, the blocks themselves (long array or list of blocks?) 
type, public :: point2d
    real(8), pointer, contiguous :: p(:, :)
end type
type, public :: point1d
    real(8), allocatable :: p(:)
end type
type, public :: block_diag
    integer :: M
    integer, allocatable :: sz(:)
    type(point2d), pointer :: u(:)
end type

type, public :: block_triangular !This matrix is organized as follows: it is block-square 
    integer :: M
    integer, pointer :: sz_red(:), sz(:)
    real(8), allocatable :: blocks(:) !The blocks
    integer, pointer :: i_block(:), j_block(:), s_block(:)
    integer, pointer :: i_block_csc(:), j_block_csc(:), s_block_csc(:)
end type

type, public :: block_sparse_matrix
    integer :: M
    type(row), pointer :: rows(:)
    integer, pointer :: i_block(:), i_block_csc(:)
    integer, pointer  :: sz(:), nc_col(:), nf_col(:)
    integer, pointer :: j_block(:), s_block(:),j_block_csc(:), s_block_csc(:)
    real(8), allocatable :: blocks(:)
end type
    !real(8), allocatable :: blocks(:) !We have to have an allocatable array for contigious memory storage. 
    !Thus, it goes to the module data
integer :: level_global=0
type(block_sparse_matrix), allocatable :: matrix_global(:)
type(block_triangular), allocatable :: Lmat_global(:)
type(block_diag), allocatable :: umat_global(:)

contains

subroutine matrix_to_python(level,matrix, M,i_block, i_block_csc,sz, nc_col, nf_col, j_block,&
 s_block, j_block_csc, s_block_csc, blocks)
    type(block_sparse_matrix) :: matrix(level)
    integer :: level
    integer :: M(level)
    integer, pointer :: i_block(:,:), i_block_csc(:,:)
    integer, pointer  :: sz(:,:), nc_col(:,:), nf_col(:,:)
    integer, pointer :: j_block(:,:), s_block(:,:),j_block_csc(:,:), s_block_csc(:,:)
    real(8), allocatable :: blocks(:,:)
    !M(1) = matrix(1)% M
    !print *, matrix(1)% i_block
    do i = 1,level
        M(i) = matrix(i)% M
    end do
    
end subroutine matrix_to_python
    subroutine init_block_matrix(matrix, M, rows, i_block, i_block_csc, sz, nc_col, nf_col, j_block, s_block, & 
             j_block_csc, s_block_csc)

     type(block_sparse_matrix) :: matrix
     integer :: M
     type(row), pointer :: rows(:)
     integer, pointer :: i_block(:), i_block_csc(:)
     integer, pointer  :: sz(:), nc_col(:), nf_col(:)
     integer, pointer :: j_block(:), s_block(:),j_block_csc(:), s_block_csc(:)
     real(8), allocatable :: blocks(:)
     matrix % M = M
     matrix % sz => sz
     matrix % rows => rows
     matrix % i_block  => i_block
     matrix % j_block => j_block
     matrix % s_block => s_block
     matrix % i_block_csc => i_block_csc
     matrix % j_block_csc => j_block_csc
     matrix % s_block_csc => s_block_csc
     return 
     end subroutine 
    
    subroutine dtransp(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(n,m)
    real(8), intent(inout):: B(m,n)
    integer i,j
    do i=1,n
       call dcopy(m, A(i,1), n, B(1,i),1)
    end do
  end subroutine dtransp

  subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end
  
  subroutine isparse_transpose(n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      implicit none
      integer :: n, n2
      integer :: job, ipos
      integer :: ia(n+1),iao(n2+1),ja(*),jao(*)
      integer ::  a(*),ao(*)
      integer :: i, j, k, next
!  Taken from SPARSEKIT subroutine csrcsc2 by Y. Saad
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n	= number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job	= integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!	  for any other normal usage, enter ipos=1.
! a	= real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of!	  the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao	= real array of size nzz containing the "a" part of the transpose
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the "ia" index array of
!	  the transpose. 
!
!----------------------------------------------------------------------- 
!----------------- compute lengths of rows of transp(A) ----------------
      do i=1,n2+1
         iao(i) = 0
      end do
      do i=1, n
         do k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
         end do 
      end do
!---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
      end do
!--------------- now do the actual copying ----------------------------- 
      do i=1,n
         do k=ia(i), ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
         end do
      end do
         !-------------------------- reshift iao and leave ---------------------- 
      do i=n2,1,-1
         iao(i+1) = iao(i)
      end do
      iao(1) = ipos
!--------------- end of csrcsc2 ---------------------------------------- 
      end subroutine

  integer function chop(s,tol,rmax,err) result (r)
  implicit none
  double precision,intent(in) :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  double precision :: nrm,er,er2,bound
  double precision,external :: dnrm2
  r=size(s); er2=0.d0
  if(present(rmax))then
      if(rmax.lt.r)then
          er2=dot_product(s(rmax+1:r),s(rmax+1:r))
          r=rmax
      end if 
  end if 
  if(present(tol))then
      bound=tol*tol
      er=er2+s(r)*s(r)
      do while(er.lt.bound .and. r .ge. 1)
         er2=er; r=r-1; er=er+s(r)*s(r)
      end do
  end if
  if(present(err))err=dsqrt(er2)
  return
end function
!##########################################################################

subroutine compute_u_mult(m, n, far_block,n_c, close_block, blocks, f_pos, c_pos, u, r, r_fix, c_l, mbc,eps)
! Here we compute, u, apply u, compute r, return r 
use dispmodule
implicit none
real(8), pointer :: u(:,:)
real(8) :: blocks(:)
integer :: m, n, i, j, r,k, n_c, f_pos, c_pos, mbc, c_l(2,mbc), r_fix
integer :: lwork, info
integer, allocatable :: jpvt(:)
real(8), allocatable :: work(:), tau(:),cur_bl(:,:)!, tmp_array(:)
real(8) :: far_block(m, n), close_block(m, n_c)!, tmp(m, n_c)
real(8) :: origin_far_block(m, n), tmp_mat(n, m)
real(8) :: r_row_norm(n), eps, s(min(m, n))
integer :: loc_sz(2), pos, nfar
real(8), allocatable  :: tmp(:, :)
real(8), external :: dnrm2
real(8) :: nrm, cf, orcond, sval(4)
integer :: rk, rmax

allocate(tmp(m, max(maxval(c_l(1, :)), n_c)))
nrm = dnrm2(m * n_c, close_block, 1)
rmax = 200
allocate(jpvt(n))
allocate(tau(min(m, n)))
lwork = 256 * max(m, n)
allocate(work(lwork))
cf = dnrm2(m * n, far_block, 1)
!print *, 'Far row norm:', cf
allocate(u(m, m))
u(:, :) = 0d0
do i = 1, m
   u(i, i) = 1d0
end do
!print *, "mn", m,n
if ( cf < 1.0E-006 ) then
    !print *, 'ZERO ROW', cf, eps, nrm
    r = 0 !Heh
    return
end if
!print *, 'cf', cf
    !cf = (nrm/cf) * eps 
    cf = eps !* nrm
    !call dcopy(n_c * m, close_block, 1, u(1, m + 1), 1)
    !Do RRQR
    !call dgeqpx(3, m, n, m, origin_far_block, m, u, m, jpvt, cf, orcond, r, sval, work, lwork, info)  
    !if ( info .ne. 0 ) then
    !   !print *, 'RRQR failed with info=', info
    !   stop
    !end if 
    !U *  [ C F ] 
    !Let us compute column norms of the far block
    !Compute the norms and skip them if it is small.
    nfar = 0
    origin_far_block(:, :) = 0d0
    do i = 1, n
       if (dnrm2(m, far_block(:, i), 1) >= 1e-10) then
           nfar = nfar + 1
           call dcopy(m, far_block(1, i), 1, origin_far_block(1, nfar), 1)
       end if
    end do
    nfar = max(m, nfar) !Otherwise it works not very well
    !print *, m, n, nfar
    !call dcopy(m * n, far_block, 1, origin_far_block, 1)
    !nfar = n
    !print *, 'far matrix of size:', m, nfar
    !call disp(origin_far_block(1:m, 1:nfar))
    !print *, 'and diag block:'
    !call disp(close_block(1:m, 1:m))
    !pause
    call dgesvd('o', 'n', m, nfar, origin_far_block, m, s, origin_far_block, m, origin_far_block, m, work, lwork, info)
    !r = chop(s, eps * nrm)
    !print *,r,r_fix
    if (r_fix > 0) then
        if (r_fix > m) then
            r = m
        else
            r = r_fix
        endif
    else
        r = chop(s, eps * nrm)
    endif
    r = min(r, rmax)
    u = origin_far_block(1:m, 1:m)
    !call disp(u)
    if ( info .ne. 0 ) then
        print *, 'DGESVD failed in compute_mult_u with info=', info
        stop
    end if
!call dgeqpb(3, m, n, m, far_block, m, u, m, jpvt, cf, orcond, r, sval, work, lwork, info)  
if (r == 0) then
   return 
end if
!Do the block scaling
!origin_far_block(1:m, jpvt) = far_block(1:m, 1:n)
!do i = 1, m
!   do j = 1, n
!       !origin_far_block(i, j) = far_block(i, jpvt(j))
!       origin_far_block(i, jpvt(j)) = far_block(i, j)
!   end do 
!end do
!origin_far_block(1:m, 1:n) = far_block(1:m, jpvt)
!call dcopy(m * n, origin_far_block, 1, blocks(f_pos), 1)
!call dcopy(m * n_c, u(1, m + 1), 1, blocks(c_pos), 1)
!print *, 'Old far zone'
!do i = 1, m
!   !print *, dnrm2(n, origin_far_block(i, 1), m)
!end do
call dcopy(n * m, far_block, 1, origin_far_block, 1)
call dgemm ('t', 'n', m, n, m, 1d0, u, m, origin_far_block, m, 0d0, blocks(f_pos), m) 
call dgemm('t', 'n', m, n_c, m, 1d0, u, m, close_block, m, 0d0, tmp, m)
call dcopy(m * n_c, tmp, 1, blocks(c_pos), 1)
!call dcopy(m * n_c, u(1, m + 1), 1, blocks(c_pos), 1)
do i = 1, mbc
    !print *, m, c_l(1, i)
    if (c_l(1, i) .NE. 0 ) then
        call dgemm('n', 'n', c_l(1, i), m, m, 1d0, blocks(c_l(2, i)), c_l(1, i), u, m, 0d0, tmp, c_l(1, i)) 
        call dcopy(c_l(1, i) * m, tmp, 1, blocks(c_l(2, i)), 1)
    end if
end do
!print *, 'diag_block after transform:'
!call disp(close_block(1:m, 1:m))
nfar = 0
origin_far_block(:, :) = 0d0
do i = 1, n
   if (dnrm2(m, far_block(:, i), 1) >= 1e-10) then
       nfar = nfar + 1
       call dcopy(m, far_block(1, i), 1, origin_far_block(1, nfar), 1)
   end if
end do
!print *, 'far zone after transform:'
!call disp(origin_far_block(1:m, 1:nfar))
!print *, 'singular values (unscaled)', s
!print *, 'Determined rank:', r
!pause
end subroutine

!##########################################################################
subroutine copy_L(sz, r, close_block, ccs, L)
integer :: ccs
integer :: sz, r
real(8) :: close_block(sz, ccs),L(sz-r,ccs)
close_block(r+1:sz, :) = L
end subroutine 
!##########################################################################
subroutine copy_Lblock(sz, r, close_block, ccs, L)
integer :: ccs
integer :: sz, r
real(8) :: close_block(sz, ccs),L(sz-r,ccs)
L(1:sz-r, 1:ccs) = close_block(r+1:sz, :)
end subroutine 
!##########################################################################
subroutine local_eliminate(sz, r, close_block, ccs, L)
use dispmodule
implicit none
!Cholesky-factorization of the diagonal block, apply to close row; -> new L columns (save)
integer :: info,i, ccs
integer :: sz, r
real(8) :: close_block(sz, ccs),L(sz-r,ccs)
real(8) :: diag_block(sz - r, sz - r)
real(8) :: c_b(sz, ccs)
c_b(:, :) = close_block(1:sz, 1:ccs)
diag_block(:, :) = close_block(r+1:sz, r+1:sz)
call dpotrf( 'l', (sz-r) , diag_block, sz-r, info) !
if (info /= 0 ) then
   print *, 'Cholesky failed, info=', info
   !print *,"bl", diag_block, 'sz', sz, 'r', r
   stop
end if
call dtrtrs('l','n','n',(sz-r),ccs,diag_block,sz-r,c_b(r+1,1),sz,info)
if (info /= 0 ) then
   print *, 'Triangular solver failed, info=', info
end if
L =  c_b(r+1:sz,:)
!Put the L block inside the close block at the corresponding rows
end subroutine
!##########################################################################
!The prepare_data gets the "raw" data, and then computes the matrix. Probably, we need a structure for that.
!We only use the "join" parameter to create the close list. 
subroutine compute_block_structure(matrix, mem)
    use dispmodule
    implicit none
    type(block_sparse_matrix) :: matrix 
    type(row), pointer :: cur_row
    integer :: N, M, i, k, kr, c_col, jr, j, nnz_blocks, Nblocks, f_col, mem
    integer :: b_pos, pos_m, max_close, pos, pos1
    integer, allocatable :: msk(:), stencil_tmp(:, :)
    integer, allocatable :: all_far(:)
    M = matrix%M
    do i = 1, M
        allocate(matrix%rows(i)%col_index(matrix%rows(i)%nclose))
        !Reorder to make the diagonal block the first
        do j = 1, matrix%rows(i)%nclose
           if (matrix%rows(i)%clos(j) == i) then
               exit
           end if
        end do
        k = matrix%rows(i)%clos(1)
        matrix%rows(i)%clos(1) = i
        matrix%rows(i)%clos(j) = k
        !Allocate col_index
        allocate(matrix%rows(i)%col_index(matrix%rows(i)%nclose+1))
        matrix%rows(i)%col_index(1) = 1
        do j = 1, matrix%rows(i)%nclose
           matrix%rows(i)%col_index(j+1) = matrix%rows(i)%col_index(j) + matrix%sz(matrix%rows(i)%clos(j))
        end do
    end do 
    ! [La Lb]' [La Lb] = [La La' La' Lb]; Lb * Lb corresponds to the new fill-in, where as La' Lb does not? i.close x i.close2 ->
    !
    !Allocate the block structure 
    allocate(msk(M))
    allocate(all_far(M))
    msk(1:M) = 0
    !Find out the block structure
    Nblocks = 0 !Total number of blocks
    nnz_blocks = 0 !Elements in the blocks array
    do i = 1, M
        do k = 1, matrix%rows(i)%nclose
            kr = matrix%rows(i)%clos(k)
            msk(kr) = 1 !Close block
        end do
        matrix%rows(i)%nfar = 0 
        !Go through all the close blocks and create, if necessary, the (i, j) block
        do k = 1, matrix%rows(i)%nclose
            kr = matrix%rows(i)%clos(k)
            do j = 1, matrix%rows(kr)%nclose
                jr = matrix%rows(kr)%clos(j)
                if (msk(jr) == 0) then
                    msk(jr) = 1
                    matrix%rows(i)%nfar = matrix%rows(i)%nfar + 1
                    all_far(matrix%rows(i)%nfar) = jr
                end if      
            end do
        end do
        allocate(matrix%rows(i)%far(matrix%rows(i)%nfar))
        matrix%rows(i)%far(1:matrix%rows(i)%nfar) = all_far(1:matrix%rows(i)%nfar)
        c_col = 0
        do k = 1, matrix%rows(i)%nclose
            kr = matrix%rows(i)%clos(k)
            msk(kr) = 0 
            c_col = c_col + matrix%sz(kr) 
        end do
        nnz_blocks = nnz_blocks + matrix%sz(i) * c_col
        matrix%rows(i)%close_col = c_col
        f_col = 0
        do k = 1, matrix%rows(i)%nfar
            kr = matrix%rows(i)%far(k)
            msk(kr) = 0
            f_col = f_col + matrix%sz(kr)
        end do 
        nnz_blocks = nnz_blocks + matrix%sz(i) * f_col
        matrix%rows(i)%far_col = f_col
        Nblocks = Nblocks + matrix%rows(i)%nclose + matrix%rows(i)%nfar   !To reference the k-th block, we look
    end do
    !Allocate the blocks
    
    allocate(matrix%blocks(nnz_blocks))
    mem = mem + nnz_blocks 
    !print *, "Non-zero elements in the blocks structure", nnz_blocks
    matrix%blocks(1:nnz_blocks) = 0d0
    allocate(matrix%i_block(M+1))
    allocate(matrix%j_block(Nblocks))
    allocate(matrix%s_block(Nblocks))
    
    pos = 1
    pos_m = 1
    b_pos = 1
    matrix%i_block(1) = 1 !The start of the block-CSR format
    pos1 = 1
!############# Fill block csr and blocks ##################
!First we should fill the i_block structure, and then copy the blocks
max_close = 0
do i = 1, M
      cur_row => matrix%rows(i)
      max_close = max(max_close, cur_row%nclose) 
      !First the diagonal block, then the close blocks, then the far blocks
      cur_row%diag_block = pos 
      cur_row%close_block = pos 
      pos = pos + matrix%sz(i) * cur_row%close_col
      cur_row%far_block = pos 
      pos = pos + matrix%sz(i) * cur_row%far_col
      do k = 1, cur_row%nclose !Cycle through close blocks
          matrix%j_block(b_pos) = cur_row%clos(k)
          matrix%s_block(b_pos) = pos_m !Where it starts
          b_pos = b_pos + 1
          kr = cur_row%col_index(k+1)-cur_row%col_index(k)    
          pos_m = pos_m + matrix%sz(i) * kr
      end do
       do k = 1, cur_row%nfar
          matrix%j_block(b_pos) = cur_row%far(k)
          matrix%s_block(b_pos) = pos_m
          b_pos = b_pos + 1
          kr = matrix%sz(cur_row%far(k))
          pos_m = pos_m + matrix%sz(i) * kr
       end do
      matrix%i_block(i+1) = matrix%i_block(i) + cur_row%nclose + cur_row%nfar
end do
      
    Nblocks = matrix%i_block(M+1) - 1
allocate(matrix%i_block_csc(M+1), matrix%j_block_csc(Nblocks), matrix%s_block_csc(Nblocks))
!############# Fill block csc ##################
call isparse_transpose(M, M, 1, 1, matrix%s_block, matrix%j_block, matrix%i_block, matrix%s_block_csc, &
    matrix%j_block_csc, matrix%i_block_csc)
!Now the problem is to evaluate the stencil
!Probably, the simplest idea is to store the block matrix in the CSR format, and then extract the stencil
!##################Now find the stencil###########################
    allocate(stencil_tmp(max_close, M))
    stencil_tmp(:, :) = 0
    do i = 1, M
        allocate(matrix%rows(i) % stencil(matrix%rows(i)%nclose, matrix%rows(i)%nclose))
        cur_row => matrix%rows(i)
        do j = 1, cur_row%nclose
           !Cycle through all the block row
           do k = matrix%i_block(cur_row%clos(j)), matrix%i_block(cur_row%clos(j)+1) - 1
                stencil_tmp(j, matrix%j_block(k)) = matrix%s_block(k)
           end do
        end do
        !Read the stencil
        do j = 1, cur_row%nclose
           do k = 1, cur_row%nclose
               cur_row % stencil(j, k) = stencil_tmp(j, cur_row%clos(k))
           end do
        end do
        !Clean the stencil
        do j = 1, cur_row%nclose
           do k = matrix%i_block(cur_row%clos(j)), matrix%i_block(cur_row%clos(j)+1) - 1
                stencil_tmp(j, matrix%j_block(k)) = 0
           end do      
        end do
    end do
!      if ( matrix%M < 256 ) then
!          do i = 1, matrix%M
!             !print *, 'filling block:', matrix%rows(i)%clos(1), i, matrix%i_block(i), b_pos, matrix%j_block(matrix%i_block(i))
!          end do
!      end if 
!    return 
end subroutine

subroutine prepare_data_csr(matrix, n, M, ia, ja, sa, bl_s, bl_ar,mem)
    use dispmodule
    implicit none
    type(block_sparse_matrix) :: matrix 
    integer :: M, n, nnz, nclose, ia(n+1), ja(*), bl_ar(M)!,r_row
    real(8) :: sa(*)
    type(row), pointer :: cur_row
    integer :: s1, s2, mem
    real(8), allocatable :: block_row(:, :)
    integer :: i, j, k, s, Nblocks, nnz_blocks, pos, c_col, f_col
    integer :: pos_m, b_pos, mbc
    integer :: ir, jr, kr, k_cl, max_close, pos1, nc
    integer, allocatable :: bc_num(:), iw(:), all_close(:)
    integer :: bl_s
    integer :: q
    !M = n / bl_s !Number of block rows
    !print *, 'Block rows:', M
    !print *, 'N', n
    matrix%M = M
    allocate(matrix%sz(M), matrix%rows(M)) 
    allocate(bc_num(n)) !Mapping of block coumns
    !Big loop, a block row each type
    pos = 1
    !Fill the block numbers of each block, compute the sizes of each larger block
    !print *, bl_ar(1)
    do i = 1, M
       matrix%sz(i) = bl_ar(i)
       matrix%rows(i)%start = pos-1 !!!!!
       matrix%rows(i)%sz = bl_ar(i)
       matrix%rows(i)%ind = i
       do k = 1, bl_ar(i)
           bc_num(pos) = i
           pos = pos + 1  
       end do
    end do
    !print *, "It work!!"
    !Compute the list of close block
    allocate(all_close(M)) 
    allocate(iw(M))
    iw(:) = 0
    kr = 1 
    do i = 1, M
       nc = 0 !Number of close blocks
       do k = 1, bl_ar(i)
          do j = ia(kr), ia(kr+1)-1
              q = bc_num(ja(j))
              if ( iw(q) .eq. 0 ) then
                  nc = nc + 1
                  iw(q) = 1
                  all_close(nc) = q
              end if
          end do
          kr = kr + 1
       end do
       matrix%rows(i)%nclose = nc
       allocate(matrix%rows(i)%clos(nc))
       matrix%rows(i)%clos(1:nc) = all_close(1:nc)
       !And clean iw
       do j = 1, nc
          iw(matrix%rows(i)%clos(j)) = 0
       end do
    end do
    !print *, 'Max block:', maxval(matrix%sz) 
    call compute_block_structure(matrix, mem)
  !Actually, fill the blocks
    allocate(block_row(maxval(matrix%sz), n))
    block_row(:, :) = 0d0
    pos = 1
    pos_m = 1
    pos1 = 1
   do i = 1, M
       cur_row => matrix%rows(i)
      !Now fill the block row
      do j = 1, matrix%sz(i)
          jr = cur_row%start + j
          do k = ia(jr), ia(jr+1)-1
              block_row(j, ja(k)) = sa(k)
          end do
      end do
      !Copy the block row
      do k = 1, cur_row%nclose !Cycle through close blocks
          kr = cur_row%col_index(k+1)-cur_row%col_index(k)    
          k_cl = cur_row%clos(k)
          do s1 = 1, matrix%sz(i)
             do s2 = 1, kr
                 matrix%blocks(pos_m + s1 - 1 + (s2 - 1) * matrix%sz(i)) = &
                 block_row(s1, matrix%rows(k_cl)%start + s2)
             end do
          end do
          !call dcopy(matrix%sz(i)*kr,block_row(1:matrix%sz(i), & 
          !matrix%rows(k_cl)%start+1:matrix%rows(k_cl)%start+matrix%sz(k_cl)),1, &
          !matrix%blocks(pos_m),1)      
          pos_m = pos_m + matrix%sz(i) * kr
      end do
      !Zero out the block row
      do j = 1, matrix%sz(i)
          jr = cur_row%start + j
          do k = ia(jr), ia(jr+1)-1
              block_row(j, ja(k)) = 0
          end do
       end do
       !Zero out the far blocks (Hmm, we can just allocate memory for them, and fill them with zeros, no need to cycle)
       do k = 1, cur_row%nfar
          kr = matrix%sz(cur_row%far(k))
          pos_m = pos_m + matrix%sz(i) * kr
          matrix%blocks(cur_row%far_block:cur_row%far_block+matrix%sz(i)*kr-1) = 0d0
       end do
    end do
  return  
 end subroutine
subroutine  make_rhs(n, rhs, sol0, sa, ja, ia)
    use dispmodule
    implicit none
    integer :: n, ia(n+1), ja(*)
    real(8) :: sa(*)
    real(8) :: sol0(n)
    real(8) :: rhs(n)
    !call 
    rhs(2) = 2d0
    call amux(n, rhs, sol0, sa,ja,ia)
    !print *,"TEST!"
    return 
end subroutine
subroutine prepare_data(matrix, M, n, nnz, nclose, ia, ja, sa, nodes, clos, mem)
    use dispmodule
    implicit none
    type(block_sparse_matrix) :: matrix 
    integer :: M, n, nnz, nclose, ia(n+1), ja(nnz), mem!,r_row
    real(8) :: sa(nnz)
    integer :: nodes(M, 4) !number of particles, start of the block, number of close blocks, position in the clos array
    integer :: clos(nclose) 
    integer, pointer :: nc_col(:), nf_col(:)
    type(row), pointer :: cur_row
    integer :: s1, s2
    real(8), allocatable :: block_row(:, :)
    integer :: i, j, k, s, Nblocks, nnz_blocks, pos, c_col, f_col
    integer :: msk(M), pos_m, b_pos, mbc
    integer :: ir, jr, kr, k_cl, max_close, pos1
    integer, allocatable :: stencil_tmp(:, :)
    matrix%M = M !Importan/
    allocate(matrix%rows(M)) !Allocate the rows
    allocate(matrix%sz(M))    
    !Copy the raw data into the structure: the minimal info are the blocks sizes and the close lists
    !We can actually get rid of it in future, since we can just look at the "block size", to get the information about the block
    !structure.
    do i = 1, M
        matrix%rows(i)%sz = nodes(i, 1)
        matrix%sz(i) = nodes(i, 1)
        matrix%rows(i)%start = nodes(i, 2)
        matrix%rows(i)%ind = i
        matrix%rows(i)%nclose = nodes(i, 3)
    end do 
    !Then we can allocate the memory for the 
    do i = 1, M
        allocate(matrix%rows(i) % clos(matrix%rows(i) % nclose))
        matrix%rows(i)%clos(1:matrix%rows(i)%nclose) = clos(nodes(i, 4):nodes(i, 4) + nodes(i, 3))
        !Reorder to make the diagonal block the first
    end do
    !print *, 'Max block:', maxval(matrix%sz) 
    !Crucial: compute the final close-far structure (and allocate the arrays)
    call compute_block_structure(matrix, mem)
    
  !Actually, fill the blocks
    allocate(block_row(maxval(matrix%sz), n))
    block_row(:, :) = 0d0
    pos = 1
    pos_m = 1
    pos1 = 1
   do i = 1, M
      cur_row => matrix%rows(i)
      !Now fill the block row
      do j = 1, matrix%sz(i)
          jr = cur_row%start + j
          do k = ia(jr), ia(jr+1)-1
              block_row(j, ja(k)) = sa(k)
          end do
      end do
      !Copy the block row
      do k = 1, cur_row%nclose !Cycle through close blocks
          kr = cur_row%col_index(k+1)-cur_row%col_index(k)    
          k_cl = cur_row%clos(k)
          do s1 = 1, matrix%sz(i)
             do s2 = 1, kr
                 matrix%blocks(pos_m + s1 - 1 + (s2 - 1) * matrix%sz(i)) = &
                 block_row(s1, matrix%rows(k_cl)%start + s2)
             end do
          end do
          !call dcopy(matrix%sz(i)*kr,block_row(1:matrix%sz(i), & 
          !matrix%rows(k_cl)%start+1:matrix%rows(k_cl)%start+matrix%sz(k_cl)),1, &
          !matrix%blocks(pos_m),1)      
          pos_m = pos_m + matrix%sz(i) * kr
      end do
      !Zero out the block row
      do j = 1, matrix%sz(i)
          jr = cur_row%start + j
          do k = ia(jr), ia(jr+1)-1
              block_row(j, ja(k)) = 0
          end do
       end do
       !Zero out the far blocks (Hmm, we can just allocate memory for them, and fill them with zeros, no need to cycle)
       do k = 1, cur_row%nfar
          kr = matrix%sz(cur_row%far(k))
          pos_m = pos_m + matrix%sz(i) * kr
          matrix%blocks(cur_row%far_block:cur_row%far_block+matrix%sz(i)*kr-1) = 0d0
       end do
    end do
    
       !print *, 'Total number of blocks written:', pos_m-1
  return  
 end subroutine

subroutine sol_ll(n,x,sol)

    implicit none
    integer :: n!,level
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: sol(n)
    real(8) :: rhs(n)
    !type(block_sparse_matrix), allocatable :: matrix(:)
    !type(block_triangular), allocatable :: Lmat(:)
    !type(block_diag), allocatable :: umat(:)

    rhs = x
    !level = level_global
    !allocate(matrix(level))
    !allocate(umat(level))
    !allocate(Lmat(level))
    !matrix = matrix_global
    !umat = umat_global
    !Lmat = Lmat_global
    !sol(1:n) = 34. 
    !print *, level
    call solve_system(n, level_global, matrix_global, Lmat_global, umat_global, rhs, sol)
    !deallocate(matrix)
    !deallocate(umat)
    !deallocate(Lmat)
    !print *, "done"
end subroutine sol_ll

subroutine fact_ll(n,M,nnz,ia,ja,sa,level,eps,r_fix,bl_s,bl_ar,join,tot_hyp,hyper_join,m_j,mem,f_time)
    use dispmodule
    implicit none
    integer :: n, nnz,join
    integer :: ia(n+1), ja(nnz), bl_ar(M)
    real(8) :: sa(nnz), eps
    integer :: level, r_fix, bl_s, tot_hyp
    integer :: m_j(level), hyper_join(tot_hyp)
    real(8) :: rhs(n)
    real(8) :: sol(n)
    integer :: M_py(level),mem
    integer, pointer :: i_block(:,:), i_block_csc(:,:)
    integer, pointer  :: sz(:,:), nc_col(:,:), nf_col(:,:)
    integer, pointer :: j_block(:,:), s_block(:,:),j_block_csc(:,:), s_block_csc(:,:)
    real(8), allocatable :: blocks(:,:)
    integer :: j_p1, j_p2, len_hj
    type(block_sparse_matrix), allocatable :: matrix(:)
    !integer :: n, ia(n+1), ja(*)!,r_row
    !real(8) :: sa(*) 
    integer :: M
    integer, allocatable :: sz_new(:)
    integer :: i, j
    real(8) :: dnrm2, nrm
    type(block_triangular), allocatable :: Lmat(:)
    type(block_diag), allocatable :: umat(:)
    real(8), allocatable :: sol0(:), rhs1(:)
    real(8) :: t0, t1, fact_time,f_time
    !level = 4
    !eps = 1e-2
    allocate(matrix(level))
    allocate(umat(level))
    allocate(Lmat(level))

    !print *, "len hyper_join",tot_hyp
    !print *, "m_j",m_j
    !print *, " hyper_join", hyper_join
    mem = 0
    !call prepare_data(matrix(1), M, n, nnz, nclose, ia, ja, sa, nodes, clos)
    call prepare_data_csr(matrix(1), n, M, ia, ja, sa, bl_s, bl_ar,mem)
    M = matrix(1)%M
    allocate(sz_new(M))
    j_p1 = 1
    j_p2 = m_j(2)
    call cpu_time(t0)
    do i = 1, level-1
       call eliminate_level(matrix(i), sz_new, umat(i),eps,r_fix,mem)  !Here it contains both the new matrix and the L matrix
       call extract_L(matrix(i), sz_new, Lmat(i), mem)
       !print *, "TEST MEM", mem
       len_hj = j_p2 - j_p1 + 1
       call next_level(matrix(i), matrix(i+1), sz_new, mem, join,len_hj, hyper_join(j_p1:j_p2)) !Here it should contain only the new matrix
       j_p1 = j_p2 + 1
       j_p2 = j_p2 + m_j(i+2)
    end do
    call cpu_time(t1)
    f_time = t1 - t0
    !print *, "Fact time", f_time
    !fact_time = fact_time + t0
    level_global = level
    mem = mem/131072
    !print *, 'Memory ll:', mem, 'Mb'
    allocate(matrix_global(level_global))
    allocate(umat_global(level_global))
    allocate(Lmat_global(level_global))
    matrix_global = matrix
    umat_global = umat
    Lmat_global = Lmat
   !call matrix_to_python(level,matrix, M_py,i_block, i_block_csc,sz, nc_col, nf_col, j_block,&
! s_block, j_block_csc, s_block_csc, blocks)
    !print *, 'Factorization time:', fact_time
    !call write_ll_to_file(level,matrix,lmat,umat)
    
end subroutine fact_ll




!########################
subroutine solve_problem(n,M,nnz,ia,ja,sa,rhs,level,eps,r_fix,sol,bl_s,bl_ar,join,tot_hyp,hyper_join,m_j,mem, &
                        f_time,sol_time,res,error)
    use dispmodule
    implicit none
    integer :: n, nnz, M, join
    integer :: ia(n+1),ja(nnz), bl_ar(M)
    real(8) :: sa(nnz),eps
    integer :: level, r_fix, bl_s, mem, tot_hyp
    integer :: m_j(level), hyper_join(tot_hyp)
    real(8) :: rhs(n)
    real(8) :: sol(n)
    real(8) :: f_time, sol_time
    real(8) :: res, error
    type(block_sparse_matrix), allocatable :: matrix(:)
    !integer :: n, ia(n+1), ja(*)!,r_row
    !real(8) :: sa(*) 
    !integer :: M
    !integer :: mem
    integer, allocatable :: sz_new(:)
    integer :: i, j, i1,jr, pos_m
    integer :: j_p1, j_p2, len_hj
    real(8) :: dnrm2, nrm
    type(block_triangular), allocatable :: Lmat(:)
    type(block_diag), allocatable :: umat(:)
    real(8), allocatable :: sol0(:), rhs1(:)
    real(8) :: t0, t1, fact_time
    !level = 4
    !eps = 1e-2
    !print *, tot_hyp
    !print *, m_j
    !print *, hyper_join
    allocate(matrix(level))
    allocate(umat(level))
    allocate(Lmat(level))
    mem = 0
    !print *, "Befor prepare_data_csr"
    !call prepare_data(matrix(1), M, n, nnz, nclose, ia, ja, sa, nodes, clos)
    call prepare_data_csr(matrix(1), n, M, ia, ja, sa, bl_s, bl_ar, mem)
    !M = matrix(1)%M
    allocate(sz_new(M))
    !print *, "Befor factorization loop"
    call cpu_time(t0)
    j_p1 = 1
    j_p2 = m_j(2)
    do i = 1, level-1
       !print *, "l ", i
       call eliminate_level(matrix(i), sz_new, umat(i),eps,r_fix, mem)  !Here it contains both the new matrix and the L matrix
       !print *, sz_new(1), sz_new(2), sz_new(3)
       !print *, "eliminate_level is done"
       call extract_L(matrix(i), sz_new, Lmat(i), mem)
       !print *, "extract_L is done"
       !print *, "SZ_N", sz_new
       len_hj = j_p2 - j_p1 + 1
       call next_level(matrix(i), matrix(i+1), sz_new, mem, join, len_hj, hyper_join(j_p1:j_p2)) !Here it should contain only the new matrix
       j_p1 = j_p2 + 1
       j_p2 = j_p2 + m_j(i+2)
       !print *, " matrix(i)",  matrix(i)%M,  matrix(i)%sz(1:matrix(i)%M)
       !print *, " matrix(i+1)",  matrix(i+1)%M,  matrix(i+1)%sz(1:matrix(i+1)%M)
       !print *, "next_level is done"
    end do
    call cpu_time(t1)
    !print *, matrix(1)%M
    !print *, matrix(1)%sz(1:matrix(1)%M)
    !print *, matrix(2)%M
    !print *, matrix(2)%sz(1:matrix(2)%M)
    !do j = matrix(2)%i_block(12), matrix(2)%i_block(13) - 1
    !   jr = matrix(2)%j_block(j)
    !   pos_m  = matrix(2)%s_block(j)
    !   if (jr .eq. 12) then
    !      !print *, matrix(2)%blocks(pos_m)
    !   end if
    !end do
    f_time = t1 - t0
    !print *, "Fact time", f_time
    !fact_time = fact_time + t0
    !print *, 'Factorization time:', fact_time
    !pause
    !print *, "Solution"
    if ( 1 < 0 ) then
        allocate( rhs1(n),sol0(n))
        !rhs(1:n) = 1d0
        sol0(1:n) = 1d0
        !call make_rhs(n, rhs, sol0, sa, ja, ia)
        call amux(n, sol0, rhs, sa,ja,ia)
        !print *, 'rhs', rhs(1:10)
        !Here we compute memory requirements:
        
        call solve_system(n, level, matrix, Lmat, umat, rhs, sol) 
        call amux(n, sol, rhs1, sa, ja, ia)
        !print *, 'residual:', dnrm2(n, rhs1(1:n) - rhs(1:n), 1)/dnrm2(n, rhs(1:n), 1)
        !print *, 'error:', dnrm2(n, sol(1:n) - sol0(1:n), 1)/dnrm2(n, sol0(1:n), 1)
    else
        allocate( rhs1(n),sol0(n))
        rhs(1:n) = 1d0
        sol0(1:n) = 1d0
        call make_rhs(n, rhs, sol0, sa, ja, ia)
        call amux(n, sol0, rhs, sa,ja,ia)
        !print *, 'rhs', rhs(1:10)
        call solve_system(n, level, matrix, Lmat, umat, rhs, sol) 
        mem = mem/131072
        !print *, 'Memory:', mem, 'Mb'
        call amux(n, sol, rhs1, sa, ja, ia)
        res = dnrm2(n, rhs1(1:n) - rhs(1:n), 1)/dnrm2(n, rhs(1:n), 1)
        error = dnrm2(n, sol(1:n) - sol0(1:n), 1)/dnrm2(n, sol0(1:n), 1)
        !print *, 'residual:', dnrm2(n, rhs1(1:n) - rhs(1:n), 1)/dnrm2(n, rhs(1:n), 1), res
        !print *, 'error:', dnrm2(n, sol(1:n) - sol0(1:n), 1)/dnrm2(n, sol0(1:n), 1), error

    end if
    call cpu_time(t0)
    sol_time = t0 - t1
    !print *, "Sol time", sol_time
    !sol_time = 0.1
    return
    end subroutine


!########################
 subroutine test_factorization(n, level, matrix0, matrix, Lmat, umat)
      use dispmodule
      implicit none
      integer :: level, n
      integer, allocatable :: nall(:)
      type(block_triangular), target :: Lmat(level)
      type(block_diag) :: umat(level)
      type(block_sparse_matrix) :: matrix0, matrix(level)
      real(8), allocatable :: tmp(:), rhs_loc(:, :)
      integer :: i, q, k
      type(point1d) :: rhs_all(level)
      real(8), pointer, contiguous :: tmp_mat(:, :), full_mat(:, :), full_delta(:, :), full_u(:, :), full_L(:, :)
      real(8), allocatable ::  tmp2(:, :), Ldiag(:, :), tmp_mat_mod(:, :), sol_true(:)
      real(8), pointer, contiguous :: tmp1(:, :)
      real(8) :: dnrm2
      integer, allocatable :: prm(:), tmp_prm(:), pos, pos1, nf(:)
      real(8), allocatable :: rhs(:), rhs1(:), sol(:), rhs2(:)
      real(8), allocatable :: x(:), y(:), f(:), g(:), L0(:, :), L1(:, :), f1(:), g1(:)
      integer :: info, N1, N2
      call full_matrix(matrix0, full_mat) 
      call full_umat(n, umat(1), full_u)
      call full_Lmat(matrix(1), Lmat(1)%sz_red, full_L, full_delta) 
      
      !Permutation part
      allocate(prm(n))
      call get_permutation(Lmat(1), prm)
      full_L(:, :) = full_L(:, prm)
      full_u(:, :) = full_u(:, prm)
      !full_mat = full_mat(prm, prm)
      full_delta(:, :) = full_delta(prm, prm)
      !call disp(full_delta(1:20, 1:20))
      !pause
      !So we have U' * A * U - L' * L - delta = 0, and we permute column of L, i.e. L -> L P, L' * L -> P' L' L * P
      !P' * U' * A * U * P - P' * L' * L * P - P' * Delta * P = 0
      allocate(tmp1(n, n), tmp2(n, n))
      !call dgemm('t', 'n', n, n, n, 1d0, full_L, n, full_L, n, 0d0, tmp1)
      call dgemm('t', 'n', n, n, n, 1d0, full_u, n, full_mat, n, 0d0, tmp2, n)
      call dgemm('n', 'n', n, n, n, 1d0, tmp2, n, full_u, n, 0d0, tmp1, n)
      allocate(tmp_mat_mod(n, n)) 
      call dcopy(n * n, tmp1, 1, tmp_mat_mod, 1)
      call dgemm('t', 'n', n, n, n, -1d0, full_L, n, full_L, n, 1d0, tmp1, n)
      !print *, dnrm2(n * n, tmp1 - transpose(tmp1), 1)
      !call disp(tmp1(n-10:n, n-10:n)-full_delta(n-10:n, n-10:n))
      !print *, 'Approximation error:', dnrm2(n * n, tmp1 - full_delta, 1)/dnrm2(n * n, tmp1, 1)
      allocate(rhs(n), rhs1(n), sol(n), rhs2(n))
      rhs(1:n) = 1d0
     ! call dcopy(n, rhs, 1, sol, 1)
      
      !call dposv('u', n, 1, full_mat, n, sol, n, info)
      !if ( info .ne. 0 ) then 
      !   !print *, 'dposv in test_factorization failed with info=', info
      !   stop
      !end if
      deallocate(tmp1)
      allocate(tmp1(n, n))
      tmp1(1:N, 1:N) = tmp_mat_mod(1:n, 1:n)
      call dcopy(n, rhs, 1, sol, 1)
      call dposv('u', n, 1, tmp1, n, sol, n, info)
      if ( info .ne. 0 ) then 
         print *, 'dposv in test_factorization failed with info=', info
         stop
      end if

      !call dgemv('n', n, n, 1d0, tmp_mat_mod, n, sol, 1, 0d0, rhs1, 1)
      !print *, rhs1
      !pause
      !Now we can try to do it approximately
      !call dgemv('t', n, n, 1d0, full_u, n, rhs, 1, 0d0, rhs1, 1)  
      call dcopy(n, rhs, 1, rhs1, 1)
      !Then we have to do the "forward step"
       N1 = sum(Lmat(1)%sz_red)
       N2 = N - N1
       !Just do a mean Schur complement computation
       deallocate(tmp1)
       allocate(tmp1(N2, N2))
       tmp1(1:N2, 1:N2) = tmp_mat_mod(1:N2, 1:N2)
       call dcopy(n, rhs1, 1, rhs2, 1)
       call dposv('u', N2, 1, tmp1, N2, rhs2, N2, info)
       if ( info .ne. 0 ) then 
           print *, 'dposv in test_factorization failed with info=', info
           stop
       end if
       deallocate(tmp1)
       allocate(tmp1(N1, N2))
       tmp1(1:N1, 1:N2) = tmp_mat_mod(N2+1:N, 1:N2)
       call dgemv('n', N1, N2, -1d0, tmp1, N1, rhs2, 1, 1d0, rhs1(N2+1), 1)
       deallocate(tmp1)
       call forward_step(Lmat(1), rhs) !the upper part of rhs should coincide with rhs1 starting from N2+1
       !call disp(rhs(1:25))
       !print *, 'and: (true one)'
       !call disp(rhs1(N2+1:N2+25))
       !call disp(rhs1(N2+1:N) - rhs(1:N1))
       !stop
       !print *, 'Forward step rhs error: (from Schur complement solver)', dnrm2(N1, rhs1(N2+1:N)-rhs(1:N1), 1)
       !print *, 'Forward step rhs error: (from the part that should not have changed)', dnrm2(N2, rhs1(1:N2)-rhs(N1+1:N), 1)
       !call disp(rhs(N1+1:N1+20))
       !print *, 'and: (true one)'
       !call disp(rhs1(1:20))
       !allocate(tmp1(N1, N1))
       
       ! L L'^-* B 

       allocate(tmp1(N2, N2)) !Check if computed rhs satisfies L' * x = z
       tmp1(1:N2, 1:N2) = full_L(1:N2, 1:N2)
        
       !call disp(rhs2(1:N2))
       !stop
       !rhs2(1:N2) = rhs1(N1+1:N) !The L^-* or L^- f 
       !Now we need to check what happens if we take the right rhs
       ![L0, L1]' * [L0, L1] = [ L0' L0 | L0' L1] 
      
       !#Doing MEGA-TEST
       deallocate(tmp1)
       call full_matrix(matrix(2), tmp1)
       call solve_full(N1, tmp1, rhs)
       deallocate(tmp1)
       allocate(x(N2), y(N1), f(N2), g(N1), L0(N2, N2), L1(N2, N1), f1(N2))
       L0(1:N2, 1:N2) = full_L(1:N2, 1:N2)
       L1(1:N2, 1:N1) = full_L(1:N2, N2+1:N)
       x(1:N2) = sol(1:N2)
       y(1:N1) = sol(N2+1:N)
       f(1:N2) = 1d0 !The true rhs 
       g(1:N1) = 1d0 
       f1(1:N2) = rhs(N1+1:N)
       !print *, dnrm2(N1, rhs(1:N1) - y, 1)
       !print *, dnrm2(N2, matmul(transpose(L0), f1) - f, 1)
       !print *, dnrm2(N2, matmul(matmul(transpose(L0), L0), x) + matmul(matmul(transpose(L0), L1), y) - f, 1)
       !print *, dnrm2(N2, matmul(L0, x) + matmul(L1, y) - f1, 1)
       f1 = f1 - matmul(L1, y)
       !print *, dnrm2(N2, matmul(L0, x) - f1, 1) 
       call backward_step(Lmat(1), rhs) 
       !print *, dnrm2(N1, rhs(1:N1) - y, 1) 
       !print *, dnrm2(N2, f1 - rhs(N1+1:N), 1)
       !print *, dnrm2(N2, x - rhs(N1+1:N), 1)


       !call dgemv('n', N2, N1, -1d0, full_L(1:N2, N2+1:N), N2, rhs(1:N1), 1, 1d0, rhs2, 1) 
       !call dgemv('t', N2, N2, -1d0, tmp1, N2, sol(1:N2), 1, 1d0, rhs2, 1)
       !call disp(rhs2(1:5))
       stop
      
        !print *, 'And the upper part:', dnrm2(N2, rhs1(1:N2) - rhs(N1+1:N), 1)
       allocate(tmp1(N1, N1))
       tmp1(1:N1, 1:N1) = full_delta(N2+1:N, N2+1:N)
       call dposv('u', N1, 1, tmp1, N1, rhs1(N2+1), N1, info)
       !print *, 'Low-order solution error:', dnrm2(N1, rhs1(N2+1:N)-sol(N2+1:N), 1)
       if ( info .ne. 0 ) then 
           print *, 'dposv in test_factorization failed with info=', info
           stop
       end if
       !print *, 'Forward step rhs error: (from Schur complement solver)', dnrm2(N1, rhs1(N2+1:N)-rhs(1:N1), 1)
       !print *, 'Forward step rhs error: (from the part that should not have changed)', dnrm2(N2, rhs1(1:N2)-rhs(N1+1:N), 1)
       !call disp(rhs(1:10))
       !pause
       call disp(sol(N2-4:N2))
       !print *, 'AND:'
       !print *, dnrm2(N, rhs, 1), dnrm2(N, sol, 1)
       !stop
       !print *, 'Backward step rhs error: (from Schur complement solver)', dnrm2(N1, sol(N2+1:N)-rhs(1:N1), 1)
       !print *, 'Backward step rhs error: (from the remaining part)', dnrm2(N2, sol(1:N2)-rhs(N1+1:N), 1)
       
       !Check if the small part of the rhs solves L * sol = rhs
       allocate(tmp1(N2, N2))
       tmp1(1:N2, 1:N2) = full_L
       call dgemv('n', N2, N2, -1d0, full_L, N2, rhs(N1+1:N), 1, 1d0, sol, 1)
       call disp(sol(1:5))
       pause
       deallocate(tmp1)

       !call disp(sol(1:5))
       !print *, 'AND:'
       !call disp(rhs(N1+1:N1+5))
       call disp(rhs(N-4:N))
       stop
       !Then do the backward step, just (?) add the result back 
       deallocate(tmp1)
       allocate(tmp1(N2, N1))
       tmp1(1:N2, 1:N1) = tmp_mat_mod(1:N2, N2+1:N)
       call dgemv('n', N2, N1, -1d0, tmp1, N2, rhs1(N2+1), 1, 1d0, rhs1, 1)
       !N2, N1 
       !And finally solve it
       deallocate(tmp1)
       allocate(tmp1(N2, N2))
       tmp1(1:N2, 1:N2) = tmp_mat_mod(1:N2, 1:N2)
       call dposv('u', N2, 1, tmp1, N2, rhs1, N2, info)
       if ( info .ne. 0 ) then 
           print *, 'dposv in test_factorization failed with info=', info
           stop
       end if
       !print *, 'full error:', dnrm2(n, sol(1:n) - rhs1(1:n), 1)  
       pause


       !call disp(rhs1(1:40)) 
      !print *, dnrm2(n * n, tmp1, 1), dnrm2(n * n, full_mat, 1)
      !Now do the funny permutation column permutation for L
      !In the permutation, we take the first to the 
      !do i = 1, matrix0%M !
      !    tmp_prm(1:(Lmat(1)%sz(i)-Lmat(1)%sz_red(i))) = prm(pos+Lmat(1)%sz_red(i):pos+Lmat(1)%sz(i)-1)
      !    tmp_prm(Lmat(1)%sz(i)-Lmat(1)%sz_red(i)+1:Lmat(1)%sz(i)) = prm(pos:pos+Lmat(1)%sz_red(i)-1)
      !    prm(pos:pos+Lmat(1)%sz(i)-1) = tmp_prm(1:Lmat(1)%sz(i))
      !    pos = pos + Lmat(1)%sz(i)
      !end do
      !print *, prm
      !pause
      !call disp(full_L(1:8, 1:10))
      !pause
 end subroutine 
 

 subroutine get_permutation(L, prm)
     use dispmodule
     implicit none
     type(block_triangular) :: L
     integer prm(*)
     integer, allocatable :: nf(:), tmp_prm(:)
     integer :: i, pos, pos1, N, N1
      N = sum(L%sz(1:L%M))
      do i = 1, N
          prm(i) = i
      end do
      allocate(nf(L%M+1), tmp_prm(N))
      nf(1) = 1
      do i = 1, L%M
         nf(i+1) = nf(i) + L%sz(i)
      end do
      
      N1 = sum(L%sz_red(1:L%M))
      pos = 1
      pos1 = N1+1
      do i = 1, L%M
         tmp_prm(pos:pos+L%sz_red(i)-1) = prm(nf(i):nf(i)+L%sz_red(i)-1)
         tmp_prm(pos1:pos1+L%sz(i) - L%sz_red(i)-1) = prm(nf(i)+L%sz_red(i):nf(i+1)-1)
         pos = pos + L%sz_red(i)
         pos1 = pos1 + L%sz(i) - L%sz_red(i)
      end do
      prm(1:N-N1) = tmp_prm(N1+1:N)
      prm(N-N1+1:N) = tmp_prm(1:N1)
      return 
     end subroutine 
 subroutine full_Lmat(matrix, sz_red, a, b)
      use dispmodule
      implicit none
      type(block_sparse_matrix) :: matrix
      real(8), pointer, contiguous  :: a(:, :), b(:, :)
      real(8), allocatable :: tmp(:, :)
      integer :: i, j, jr, pos_m, N
      integer :: sz_red(*)
      integer, allocatable :: nf(:), nf_red(:)
      N = sum(matrix%sz(1:matrix%M))
      allocate(a(N, N), b(N, N), nf(matrix%M), nf_red(matrix%M)) !Then size of the small system
      a(1:N, 1:N) = 0d0
      b(1:N, 1:N) = 0d0
      nf(1) = 1
      nf_red(1) = 1
      do i = 1, matrix%M-1
         nf(i+1) = nf(i) + matrix%sz(i)
         nf_red(i+1) = nf_red(i) + matrix%sz(i) - sz_red(i)
      end do
      do i = 1, matrix%M
         do j = matrix%i_block(i), matrix%i_block(i+1) - 1
             jr = matrix%j_block(j)
             pos_m  = matrix%s_block(j)
             allocate(tmp(matrix%sz(i), matrix%sz(jr)))
             call dcopy(matrix%sz(i) * matrix%sz(jr), matrix%blocks(pos_m), 1, tmp, 1)

             !call copy_submatrix(matrix%sz(i)-sz_red(i), matrix%sz(jr), N, N, nf(i)+sz_red(i), nf(jr), &
             !    tmp( sz_red(i)+1:matrix%sz(i),1:matrix%sz(jr)), a)
             !call copy_submatrix(matrix%sz(i)-sz_red(i), matrix%sz(jr), N, N, nf(i)+sz_red(i), nf(jr), &
             !    tmp( sz_red(i)+1:matrix%sz(i),1:matrix%sz(jr)), a)
             call copy_submatrix(matrix%sz(i)-sz_red(i), matrix%sz(jr), N, N, nf_red(i), nf(jr), &
                 tmp( sz_red(i)+1:matrix%sz(i),1:matrix%sz(jr)), a)
             
             call copy_submatrix(sz_red(i), sz_red(jr), N, N, nf(i), nf(jr), &
                 tmp( 1:sz_red(i),1:sz_red(jr)), b)
             deallocate(tmp)
         end do
      end do
  end subroutine


 subroutine full_umat(n, umat, full_u)
 use dispmodule
 implicit none
 type(block_diag) :: umat
 real(8), pointer, contiguous :: full_u(:, :)
 integer :: i, n, pos
 allocate(full_u(n, n))
 full_u(1:n, 1:n) = 0d0
 pos = 1
 do i = 1, umat%M
     full_u(pos:pos+umat%sz(i)-1, pos:pos+umat%sz(i)-1) = umat%u(i)%p(1:umat%sz(i), 1:umat%sz(i))
     pos = pos + umat%sz(i)
 end do
 end subroutine 
!########################
! This subroutine, given the factorization, computes the (approximate) solution of a system of linear equations.
!########################
!We can check, if a small right-hand side for 16x16 problem coincides with the Python code.  
 subroutine solve_system(n, level, matrix, Lmat, umat, rhs, sol)
      use dispmodule
      implicit none
      integer :: level, n
      integer, allocatable :: nall(:)
      type(block_triangular) :: Lmat(level)
      type(block_diag) :: umat(level)
      type(block_sparse_matrix) :: matrix(level)
      real(8) :: rhs(*), sol(*)
      real(8), allocatable :: tmp(:), rhs_loc(:, :)
      integer :: i, q, k
      type(point1d) :: rhs_all(level)
      real(8), pointer, contiguous :: full_mat(:, :)
      real(8) :: dnrm2, st, ed, solver_time
      solver_time = 0
      allocate(tmp(n))
      allocate(nall(level))
      call cpu_time(st)
      do i = 1, level
         nall(i) = sum(matrix(i)%sz(1:matrix(i)%M))
         allocate(rhs_all(i)%p(nall(1)))
      end do
      call dcopy(n, rhs, 1, rhs_all(1)%p, 1) 
      !call disp(rhs_all(k)%p(1:nall(k)))
      !pause
      do i = 1, level-1
         call apply_umat_transpose(umat(i), rhs_all(i)%p, rhs_all(i+1)%p) 
         call forward_step(Lmat(i), rhs_all(i+1)%p) !Then rhs_loc will contain (as a part) the "right" right-hand side for a system 
      end do
      !print *, "Level ok"
      !print *, 'Level', level
      call sparse_matrix(matrix(level), full_mat)
      !print *, "make_full ok"
      call solve_full(nall(level), full_mat, rhs_all(level)%p)
      !print *, "Full ok"
      !return
      !Now given sol we have to do "the backwards" step
      do i = level-1, 1, -1
          !The "g" component of the solution is stored where it should be stored (in the first components) 
          call backward_step(Lmat(i), rhs_all(i+1)%p)
          call apply_umat(umat(i), rhs_all(i+1)%p, rhs_all(i)%p)
      end do
      !print *, "backvard ok"
      !print *, dnrm2(nall(1), rhs_all(1)%p, 1)
      call dcopy(nall(1), rhs_all(1)%p, 1, sol, 1)
      call cpu_time(ed)
      st = ed - st
      solver_time = solver_time + st
      !print *, 'Solver time:', solver_time
      deallocate(full_mat)
  end subroutine
!########################
! Forward step is just a cycle through the lower-triangular blocks
!########################
  subroutine forward_step(L, rhs)
      use dispmodule
      implicit none
      type(block_triangular) :: L
      integer :: pos_m, pos, pos1, i, j, jr, k, sz_max, N, N1
      integer, allocatable :: nf(:)
      integer :: info
      real(8) :: rhs(*)
      real(8), allocatable :: Ldiag(:, :), rhs_loc(:), tmp(:)
      real(8), allocatable :: rhs_old(:) 
      sz_max = maxval(L%sz) 
      allocate(nf(L%M+1))
      nf(1) = 1
      do i = 1, L%M
         nf(i+1) = nf(i) + L%sz(i)
      end do
      N = sum(L%sz(1:L%M))
      allocate(rhs_old(N))
      call dcopy(N, rhs, 1, rhs_old, 1)
      !What we have to do, is having L^T B^T = L^T L B^T L compute A^{-1} B^T = (L' B)' = (L/B') (L' B) = (LL' L B)  (LL')^{-1} L B
      != L'^{-1} B -> B' * L^{-1} -> L is lower triangular (that means we start from the first equation and **should** get the same
      !rhs!!!!, however, indeed, since we store the transpose of L, we should start from the bottom, i.e. first solve the last
      !equation in our setting, when [L/B']*L^{-1} = 
      N1 = sum(L%sz_red(1:L%M))
      allocate(tmp(N))
!      pos1 = N1+1
!      do i = 1, L%M
!         tmp(pos1:pos1+L%sz(i) - L%sz_red(i)-1) = rhs(nf(i)+L%sz_red(i):nf(i+1)-1)
!         pos1 = pos1 + L%sz(i) - L%sz_red(i)
!      end do
      do i = 1, L%M
         !Extract diag block   The row is L^{top}, i.e. it is "upper" triangular in this form
         ! B^t A^{-1} i.e. multiply the row A B from the left by A^{-1}, i.e. the db is factorized as L^{\top} L i
         ! and L^{\top} LL i.e. we multiply by L^{\top}^{-1} if we transpose we get L LL^{\top} and we have to multiply by L^{-1}
         ! from the left, i.e. we have transpose solve!
         pos_m = L%s_block(L%i_block(i))
         k = L%sz(i) - L%sz_red(i)
         !Go through all the close
         !call dcopy(k, rhs(nf(i)+L%sz_red(i)), 1, rhs_loc, 1)  !Copy the current rhs to a local storage
         !Solve
         if ( k > 0 ) then
             allocate(Ldiag(k, k))
             call get_diag_block(L%sz(i), L%sz_red(i), L%blocks(pos_m), Ldiag)
             !call disp(Ldiag)
             !pause
             !rhs_loc(1:k) = rhs(nf(i)+L%sz_red(i):nf(i+1)-1)
             call dtrtrs('u','t','n',k,1,Ldiag,k,rhs(nf(i)+L%sz_red(i)),k,info)
             !call dtrtrs('u','t','n',k,1,Ldiag,k,rhs_loc,k,info)
             deallocate(Ldiag)
             if (info .ne. 0) then
                 print *, 'dtrtrs in forward_step failed with info=', info, k
                 stop
             end if
             !Modify the rhs as required (and it seems it will exactly produce the right rhs for a smaller system)
             do j = L%i_block(i), L%i_block(i+1)-1
             pos_m = L%s_block(j)
             jr = L%j_block(j)
             !print *, i, jr, L%sz_red(i)
             !call dgemv('t', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs(nf(i)+L%sz_red(i)), 1, 1d0, rhs(nf(jr)), 1) 
             !call dcopy(L%sz(jr), rhs(nf(jr)), 1, rhs_loc, 1)  
             !call dgemv('t', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs(nf(i)+L%sz_red(i)), 1, 1d0, rhs_loc, 1)

             !The DGEMV here should be done only for the "outer part", i.e. we eliminate B - r rows, thus r rows have to be
             !However, at each step we eliminate only $B - r$ variables, we can also try to "reorder" the variable to avoid rhs
             !pollution
             if ( jr > i) then
                 !call dgemv('t', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs_loc, 1, 1d0, rhs(nf(jr)), 1)
                 call dgemv('t', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs(nf(i)+L%sz_red(i)), 1, 1d0, rhs(nf(jr)), 1)
             else
                 !call dgemv('t', k, L%sz_red(jr), -1d0, L%blocks(pos_m), k, rhs_loc, 1, 1d0, rhs(nf(jr)), 1)
                 call dgemv('t', k, L%sz_red(jr), -1d0, L%blocks(pos_m), k, rhs(nf(i)+L%sz_red(i)), 1, 1d0, rhs(nf(jr)), 1)
             end if
             end do
     end if
      end do
      pos1 = N1+1
      do i = 1, L%M
         tmp(pos1:pos1+L%sz(i) - L%sz_red(i)-1) = rhs(nf(i)+L%sz_red(i):nf(i+1)-1)
         pos1 = pos1 + L%sz(i) - L%sz_red(i)
      end do
      pos = 1
      do i = 1, L%M
         tmp(pos:pos+L%sz_red(i)-1) = rhs(nf(i):nf(i)+L%sz_red(i)-1)
         !tmp(pos1:pos1+L%sz(i) - L%sz_red(i)-1) = rhs(nf(i)+L%sz_red(i):nf(i+1)-1)
         pos = pos + L%sz_red(i)
         !pos1 = pos1 + L%sz(i) - L%sz_red(i)
      end do
      call dcopy(N, tmp, 1, rhs, 1)
  end subroutine 
!########################
  !The backward step is obtained by reversing the operations in the forward step, and also replacing "-" by "+"
!########################
  subroutine backward_step(L, rhs)
      use dispmodule
      implicit none
      type(block_triangular) :: L
      integer :: pos_m, pos, pos1, i, j, jr, k, sz_max, N, N1
      integer, allocatable :: nf(:)
      integer :: info, ind_diag
      real(8) :: rhs(*)
      real(8), allocatable :: Ldiag(:, :), rhs_loc(:), tmp(:)
      sz_max = maxval(L%sz) 
      allocate(nf(L%M+1))
      nf(1) = 1
      do i = 1, L%M
         nf(i+1) = nf(i) + L%sz(i)
      end do
      N = sum(L%sz(1:L%M))
      N1 = sum(L%sz_red(1:L%M))
      !Copy "top" rhs back to the right plance
      allocate(tmp(N))
      !We have to take the y vector and do f - B y and then (f - A^{-1} B y)
      call dcopy(N, rhs, 1, tmp, 1)
      pos = 1
      pos1 = N1+1
      do i = 1, L%M
         rhs(nf(i):nf(i)+L%sz_red(i)-1) = tmp(pos:pos+L%sz_red(i)-1) !Put the solution *into* the rhs
         rhs(nf(i)+L%sz_red(i):nf(i+1)-1) = tmp(pos1:pos1+L%sz(i)-L%sz_red(i)-1) !This is L^-* f here
         pos = pos + L%sz_red(i)
         pos1 = pos1 + L%sz(i) - L%sz_red(i)
      end do 
     !We have found $y$ (correctly), thus we have Lx = L^-* f - L_B y - this is what we have to form in the backward_step
     ! Probably, if we do not solve, we have correctly formed the rhs for Lx = rhs -- and that is what we have to think about - for
     ! the diag part it is (B - sz(i)) x (B - sz(i))
     ! matrix is stored as:
     ! L11 L12
     !     L22 ; we have to solve
     ! L11' 0     
     ! L12' L22'
     ! So we first form the rhs.
     ! Then solve for 
     ! At first we multiply y and modify the rhs for f. That should give the right rhs for the last $x$ 
     do i = L%M, 1, -1
         k = L%sz(i) - L%sz_red(i)
         if ( k > 0 ) then
             do j = L%i_block(i), L%i_block(i+1)-1
             pos_m = L%s_block(j)
             jr = L%j_block(j)
             allocate(rhs_loc(k))
             call dcopy(k, rhs(nf(i) + L%sz_red(i)), 1, rhs_loc, 1)    
             if ( jr > i ) then !Off diagonal terms 
                 !call dgemv('n', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs(nf(jr)), 1, 1d0, rhs(nf(i) + L%sz_red(i)), 1)
                 call dgemv('n', k, L%sz(jr), -1d0, L%blocks(pos_m), k, rhs(nf(jr)), 1, 1d0, rhs_loc, 1)
             else
                 !call dgemv('n', k, L%sz_red(jr), -1d0, L%blocks(pos_m), k, rhs(nf(jr)), 1, 1d0, rhs(nf(i) + L%sz_red(i)), 1)
                 call dgemv('n', k, L%sz_red(jr), -1d0, L%blocks(pos_m), k, rhs(nf(jr)), 1, 1d0, rhs_loc, 1)
             end if
             !call dgemv('n', k, L%sz_red(jr), -1d0, L%blocks(pos_m), k, rhs(nf(jr)), 1, 1d0, rhs_loc, 1)
             call dcopy(k, rhs_loc, 1, rhs(nf(i)+L%sz_red(i)), 1)
             deallocate(rhs_loc)
             !end if
             end do
             pos_m = L%s_block(L%i_block(i))
             allocate(Ldiag(k, k))
             call get_diag_block(L%sz(i), L%sz_red(i), L%blocks(pos_m), Ldiag)
             call dtrtrs('u','n','n',k,1,Ldiag,k,rhs(nf(i)+L%sz_red(i)),k,info)
             deallocate(Ldiag)
         end if
      end do
      !This step was only for debugging purposes
      !pos1 = N1+1
      !pos = 1
      !do i = 1, L%M
      !   tmp(pos:pos+L%sz_red(i)-1) = rhs(nf(i):nf(i)+L%sz_red(i)-1)
      !   tmp(pos1:pos1+L%sz(i) - L%sz_red(i)-1) = rhs(nf(i)+L%sz_red(i):nf(i+1)-1)
      !   pos = pos + L%sz_red(i)
      !   pos1 = pos1 + L%sz(i) - L%sz_red(i)
      !end do
      !call dcopy(N, tmp, 1, rhs, 1)
      return
  end subroutine 
!########################
!Solve the full system matrix (can remove later, since it is just a copy + dposv)  
  subroutine solve_full(n, a, rhs)
      use dispmodule 
      implicit none
      integer :: n, i
      real(8) :: a(n, n), b(n, n), rhs(n)
      integer :: info
      real(8) :: dnrm2
      !print *,"N_sm", n
      !print *, "A", a(n-5:n, n-5:n)
      !do i = 1, n
      !  !print *, i, a(i,i)
      !end do
      call dcopy(n * n, a, 1, b, 1)
      !print *, a
      call dposv('u', n, 1, b, n, rhs, n, info)
      if ( info .ne. 0 ) then
         print *, 'Solve_full, dposv failed with info=', info
         stop
      end if
      !call dgemv('n', n, n, -1d0, a, n, x, 1, 1d0, rhs, 1)
      !print *, dnrm2(n, rhs, 1)
      !call disp(a)
      end subroutine 
!########################
  subroutine full_matrix(matrix, a)
      use dispmodule
      implicit none
      type(block_sparse_matrix) :: matrix
      real(8), pointer, contiguous  :: a(:, :)
      integer :: i, j, jr, pos_m, N
      integer, allocatable :: nf(:)
      N = sum(matrix%sz(1:matrix%M))
      allocate(a(N, N), nf(matrix%M)) !Then size of the small system
      a(1:N, 1:N) = 0d0
      nf(1) = 1
      do i = 1, matrix%M-1
         nf(i+1) = nf(i) + matrix%sz(i)
      end do
      
  end subroutine

  subroutine sparse_matrix(matrix, a)
      use dispmodule
      implicit none
      type(block_sparse_matrix) :: matrix
      real(8), pointer, contiguous  :: a(:, :), b(:,:)
      integer :: i, j, jr, pos_m, N, p, k, bn_i = 0, bn_j = 0
      integer, allocatable :: nf(:)
      !print *, matrix%sz(1:matrix%M)
      N = sum(matrix%sz(1:matrix%M))
      allocate(a(N, N), nf(matrix%M)) !Then size of the small system
      a(1:N, 1:N) = 0d0
      nf(1) = 1
      do i = 1, matrix%M-1
         nf(i+1) = nf(i) + matrix%sz(i)
      end do
      !print *, matrix%M
      do i = 1, matrix%M
         do j = matrix%i_block(i), matrix%i_block(i+1) - 1
             jr = matrix%j_block(j)
             pos_m  = matrix%s_block(j)
             !if (i .eq. jr) then
             !   !print *,i,jr, matrix%sz(i),  matrix%sz(jr)
             !   !print *, matrix%blocks(pos_m)
              
             !end if
             call copy_submatrix(matrix%sz(i), matrix%sz(jr), N, N, nf(i), nf(jr), matrix%blocks(pos_m), a)
         end do
      end do
  end subroutine


!########################
  subroutine print_diag_block(n, matrix)
      use dispmodule
      implicit none
      integer :: n 
      real(8) :: matrix(n, n)
      call disp(matrix)
      end subroutine 

!########################
  subroutine copy_submatrix(n1, m1, n2, m2, s1, s2, a, b)
      !Copy matrix a into the submatrix of the matrix b starting at s1, s2
      use dispmodule
      implicit none
      integer :: n1, m1, n2, m2, s1, s2
      real(8) :: a(n1, m1), b(n2, m2)
      b(s1:s1+n1-1, s2:s2+m1-1) = a(1:n1, 1:m1)
  end subroutine 
! We actually have (B - r) x ncols matrix, where each block corresponds to the "remaining" part. 
! The 
! A B    x
! B^t C  y says C y - B^T A^{-1} B y = g - B^t A^{-1} f = g - L^{\top} L^{-1} f so we solve for f (which we eliminate) and modify g
! at each step L^{-1} f is a vector of length (B-r) so the components of g are all-modified B
!########################
  subroutine get_diag_block(sz, r, Lblock, Ldiag)
      implicit none
      integer :: sz, r
      real(8) :: Ldiag(sz-r, sz-r), Lblock(sz-r, sz)
      Ldiag(:, :) = Lblock(:, r+1:sz)
  end subroutine

!########################
  subroutine apply_umat(umat, rhs, res)
      use dispmodule
      implicit none
      type(block_diag) :: umat
      real(8) :: rhs(*), res(*)
      integer :: i, pos
      pos = 1
      do i = 1, umat%M
          call dgemv('n', umat%sz(i), umat%sz(i), 1d0, umat%u(i)%p, umat%sz(i), rhs(pos), 1, 0d0, res(pos), 1)
          !call disp(rhs(pos:pos+umat%sz(i)-1))
          !call disp(res(pos:pos+umat%sz(i)-1))
          !call disp(umat%u(i)%p)
          !pause
          pos = pos + umat%sz(i)
      end do
  end subroutine 
!########################
  subroutine apply_umat_transpose(umat, rhs, res)
      implicit none
      type(block_diag) :: umat
      real(8) :: rhs(*), res(*)
      integer :: i, pos
      pos = 1
      do i = 1, umat%M
          call dgemv('t', umat%sz(i), umat%sz(i), 1d0, umat%u(i)%p, umat%sz(i), rhs(pos), 1, 0d0, res(pos), 1)
          pos = pos + umat%sz(i)
      end do
  end subroutine 
!########################
subroutine extract_L(matrix, sz_new, L, mem)
    use dispmodule
    implicit none
    type(block_sparse_matrix) :: matrix
    type(block_triangular) :: L
    integer :: sz_new(*), mem
    integer :: i, j, jr, M, pos, pos_m, ncols, nblocks, nnz
    !integer :: M
    !integer, pointer :: sz_red(:), sz(:)
    !real(8), allocatable :: blocks(:) !The blocks
    !integer, pointer :: i_block(:), j_block(:), s_block(:)
    M = matrix%M
    L%M = M
    allocate(L%sz_red(M), L%sz(M), L%i_block(M+1))
    !Count the number of close blocks
    nblocks = 0
    nnz = 0
    do i = 1, M
       L%sz(i) = matrix%sz(i)
       L%sz_red(i) = sz_new(i)
       nblocks = nblocks + matrix%rows(i)%nclose
       do j = 1, matrix%rows(i)%nclose
           jr = matrix%rows(i)%clos(j)
           nnz = nnz + (L%sz(i) - L%sz_red(i)) * (matrix%sz(jr)) !The Lmatrix takes into account only el vars     
       end do
    end do
    allocate(L%j_block(nblocks))
    allocate(L%s_block(nblocks))
    allocate(L%blocks(nnz))
    mem = mem + nnz 
    !Fill the sparsity pattern
    L%i_block(1) = 1
    pos = 1
    pos_m = 1
    do i = 1, M
        L%i_block(i+1) = L%i_block(i) + matrix%rows(i)%nclose
        ncols = 0
        do j = 1, matrix%rows(i)%nclose
           jr = matrix%rows(i)%clos(j)
           L%j_block(pos) = jr
           L%s_block(pos) = pos_m
           pos_m = pos_m + (L%sz(i) - L%sz_red(i)) * (L%sz(jr))
           ncols = ncols + (L%sz(jr))
           pos = pos + 1
        end do
        jr = L%s_block(L%i_block(i))
        call copy_Lblock(L%sz(i),L%sz_red(i),matrix%blocks(matrix%rows(i)%close_block),ncols,L%blocks(jr)) !The matrix L is saved in the
     end do
    allocate(L%i_block_csc(M+1), L%j_block_csc(nblocks), L%s_block_csc(nblocks))
    !############# Fill block csc ##################
    call isparse_transpose(M, M, 1, 1, L%s_block, L%j_block, L%i_block, L%s_block_csc, &
    L%j_block_csc, L%i_block_csc)
    return 
    end subroutine 
!########################

subroutine next_level(matrix, matrix_small, sz_new, mem, join, len_hj, hyper_join) !Modify the matrix according to the level compression
    !What we have to do, is to recompute all the stuff. 
    use dispmodule
    implicit none
    type(block_sparse_matrix) :: matrix
    type(block_sparse_matrix) :: matrix_small
    integer :: sz_new(*), mem, join
    integer :: i, j, k, pos_m, pos_k,  pos, irow, kr, jr, nc, pos_i, pos_j
    integer :: N, M, q, s1, s2, loc_sz(2), sm_join, len_hj
    !integer :: M_1
    integer :: hyper_join(len_hj)
    real(8) :: dnrm2, nrm, nrm1, nrm2
    real(8), allocatable :: block_row(:, :)
    integer, allocatable :: nf(:), nf_small(:), iw(:), bc_num(:), all_close(:)
    !integer, parameter :: join  !How many block rows we need to merge; 
    !print *, "len_hj", len_hj
    !print *, "hyper_join", hyper_join

    ! Before:
    !M = (matrix%M) / join !Number of block rows
    !sm_join = MOD(matrix%M,join)
    !if (sm_join .ne. 0) then
    !    M = (matrix%M) / join + 1
    !else
    !    sm_join = join
    !end if

    ! After:
    M = len_hj

    !print *, M
    allocate(bc_num(matrix%M)) !Mapping of block coumns
    !Big loop, a block row each type
    matrix_small%M = M
    allocate(matrix_small%sz(M), matrix_small%i_block(M+1), matrix_small%rows(M))
    pos = 1
    !Fill the block numbers of each block, compute the sizes of each larger block
    
    ! Before:
    !do i = 1, M-1

    !After:
    do i = 1, M
       !print *, i
       matrix_small%sz(i) = 0

       ! Before:
       !do k = 1, join
        
       ! After:
       do k = 1, hyper_join(i)
           bc_num(pos) = i
           matrix_small%sz(i) = matrix_small%sz(i) + sz_new(pos)
           !print *, matrix_small%sz(i), sz_new(pos)
           !print *, pos, sz_new(pos)
           pos = pos + 1  
       end do
       !print *, matrix_small%sz(i), join, pos
    end do
    
    !Before:
    !i = M
    !    matrix_small%sz(i) = 0
    !    do k = 1, sm_join
    !       bc_num(pos) = i
    !       matrix_small%sz(i) = matrix_small%sz(i) + sz_new(pos)
    !       pos = pos + 1  
    !   end do
    !Compute the list of close block
    !print *, "Compute the list of close block"
    allocate(all_close(M)) 
    allocate(iw(M))
    iw(:) = 0d0
    kr = 1 
    ! Before:
    !do i = 1, M-1

    !After:
    do i = 1, M

       !print *, "I", i
       nc = 0 !Number of close blocks
       do k = 1, hyper_join(i)
          !print *, "k",k
          do j = 1, matrix%rows(kr)%nclose
              q = bc_num(matrix%rows(kr)%clos(j))
              if ( iw(q) .eq. 0 ) then
                  nc = nc + 1
                  iw(q) = 1
                  all_close(nc) = q
              end if
          end do
          kr = kr + 1
       end do
       matrix_small%rows(i)%nclose = nc
       allocate(matrix_small%rows(i)%clos(nc))
       matrix_small%rows(i)%clos(1:nc) = all_close(1:nc)
       !And clean iw
       do j = 1, nc
          iw(matrix_small%rows(i)%clos(j)) = 0
       end do
    end do
    
    !Before:    
    !i = M
    !   !print *, "I", i
    !   nc = 0 !Number of close blocks
    !   do k = 1, sm_join
    !      !print *, "k",k
    !      do j = 1, matrix%rows(kr)%nclose
    !          q = bc_num(matrix%rows(kr)%clos(j))
    !          if ( iw(q) .eq. 0 ) then
    !              nc = nc + 1
    !              iw(q) = 1
    !              all_close(nc) = q
    !          end if
    !      end do
    !      kr = kr + 1
    !   end do
    !   matrix_small%rows(i)%nclose = nc
    !   allocate(matrix_small%rows(i)%clos(nc))
    !   matrix_small%rows(i)%clos(1:nc) = all_close(1:nc)
    !   !And clean iw
    !   do j = 1, nc
    !      iw(matrix_small%rows(i)%clos(j)) = 0
    !   end do

    !Now given the matrix_small_rows structure we can do all
    !print *, "Now given the matrix_small_rows structure we can do all"
    call compute_block_structure(matrix_small, mem)
    !print *, "Done"
    !Now we have to fill the blocks by copying them in a right way. We just copy the parts of the bl
    !kr = 1
    !We have to merge several block rows into a block.
    !We also (!) have to work with the full sparsity pattern, not only with close blocks
    !we can actually 
    !The size of block row should be (max_sz) x total number of unknowns in the reduced matrix (i.e., sum of sz_new)
    N = sum(sz_new(1:matrix%M))
    !print *, 'Total number of unknowns in the reduced system is', N, 'maximal rank:', maxval(sz_new(1:matrix%M)), & 
    !         'Mean rank:', N*1.0/matrix%M
    allocate(block_row(maxval(matrix_small%sz), sum(matrix%sz(1:matrix%M))))
    block_row(:, :) = 0d0
    allocate(nf(matrix%M))
    allocate(nf_small(matrix_small%M))
    nf(1) = 1
    do i = 2, matrix%M
       nf(i) = nf(i-1) + sz_new(i-1) 
    end do
    kr = 1

    ! Before:
    !do i  = 1, matrix_small%M - 1

    !After:
    do i  = 1, matrix_small%M 

       nf_small(i) = nf(kr)
       
       ! Before:
       !kr = kr + join
       
       !After:
       kr = kr + hyper_join(i)
    end do
    
    ! Before:
    !i = matrix_small%M
    !   nf_small(i) = nf(kr)
    !   kr = kr + sm_join
    
    !We just have to write the blocks in right positions taking into the account their reduced size.
    !And then just reread the block row
    kr = 1

    ! Before:    
    !do i = 1, matrix_small%M - 1 !Go through all the new block rows

    !After:
    do i = 1, matrix_small%M
        pos_k = 0

        ! Before:
        !do k = 1, join

        !After:
        do k = 1, hyper_join(i)
           do j = matrix%i_block(kr), matrix%i_block(kr+1)-1 !Through all the blocks in this block row
                jr = matrix%j_block(j)
                pos_m = matrix%s_block(j) 
                do s1 = 1, sz_new(kr)
                   do s2 = 1, sz_new(jr)
                      block_row(pos_k + s1, nf(jr) + s2 - 1) = matrix%blocks(pos_m + s1 + (s2 - 1) * matrix%sz(kr) - 1)  
                   end do
                end do
            end do
            pos_k = pos_k + sz_new(kr)
            kr = kr + 1
        end do
        !Read the block row to the matrix 
        do j = matrix_small%i_block(i), matrix_small%i_block(i+1)-1
            pos_m = matrix_small%s_block(j)
            jr = matrix_small%j_block(j)
            do s1 = 1, matrix_small%sz(i)
               do s2 = 1, matrix_small%sz(jr) 
                   matrix_small%blocks(pos_m + s1 + (s2 - 1) * matrix_small%sz(i) - 1) = block_row(s1, nf_small(jr) + s2 - 1)
               end do
            end do
        end do
        pos_k = 0
        
        ! Before:
        !kr = kr - join
        !do k = 1, join
        
        !After:
        kr = kr - hyper_join(i)
        do k = 1, hyper_join(i)
            do j = matrix%i_block(kr), matrix%i_block(kr+1)-1 !Through all the blocks in this block row
                jr = matrix%j_block(j)
                do s1 = 1, sz_new(kr)
                   do s2 = 1, sz_new(jr)
                      block_row(pos_k + s1, nf(jr) + s2 - 1) = 0d0!matrix%blocks(pos_m + s1 + (s2 - 1) * matrix%sz(jr))  
                   end do
                end do
            end do
            pos_k = pos_k + sz_new(kr)
            kr = kr + 1
        end do
    end do

    ! Before:
    !i = matrix_small%M
    !    pos_k = 0
    !    do k = 1, sm_join
    !       do j = matrix%i_block(kr), matrix%i_block(kr+1)-1 !Through all the blocks in this block row
    !            jr = matrix%j_block(j)
    !            pos_m = matrix%s_block(j) 
    !            do s1 = 1, sz_new(kr)
    !               do s2 = 1, sz_new(jr)
    !                  block_row(pos_k + s1, nf(jr) + s2 - 1) = matrix%blocks(pos_m + s1 + (s2 - 1) * matrix%sz(kr) - 1)  
    !               end do
    !            end do
    !        end do
    !        pos_k = pos_k + sz_new(kr)
    !        kr = kr + 1
    !    end do
    !    !Read the block row to the matrix 
    !    do j = matrix_small%i_block(i), matrix_small%i_block(i+1)-1
    !        pos_m = matrix_small%s_block(j)
    !        jr = matrix_small%j_block(j)
    !        do s1 = 1, matrix_small%sz(i)
    !           do s2 = 1, matrix_small%sz(jr) 
    !               matrix_small%blocks(pos_m + s1 + (s2 - 1) * matrix_small%sz(i) - 1) = block_row(s1, nf_small(jr) + s2 - 1)
    !           end do
    !        end do
    !    end do
    !    pos_k = 0
    !    kr = kr - sm_join
    !    do k = 1, sm_join
    !        do j = matrix%i_block(kr), matrix%i_block(kr+1)-1 !Through all the blocks in this block row
    !            jr = matrix%j_block(j)
    !            do s1 = 1, sz_new(kr)
    !               do s2 = 1, sz_new(jr)
    !                  block_row(pos_k + s1, nf(jr) + s2 - 1) = 0d0!matrix%blocks(pos_m + s1 + (s2 - 1) * matrix%sz(jr))  
    !               end do
    !            end do
    !        end do
    !        pos_k = pos_k + sz_new(kr)
    !        kr = kr + 1
    !    end do
    !print *, "UPD version"
end subroutine


subroutine eliminate_level(matrix, sz_new, umat,eps, r_fix, mem)
    use dispmodule
    implicit none
    integer :: sz_new(*), mem
    type(block_sparse_matrix) :: matrix
    type(block_diag) :: umat
    !Now we will have to do all the dirty work here; first, allocate the blocks array (symbolic part)
    !real(8),  allocatable:: u(:)  
    integer, allocatable:: col_list(:,:)
    integer :: i, j, k, s1, s2
    integer :: mbc
    integer :: r, rmax, r_fix
    real(8) :: eps     
    real(8), allocatable :: L(:,:)
    integer :: f_pos, c_pos, pos_j, pos_j_n, pos_i, pos_i_n, col_close_size 
    real(8) :: ed, st, local_time, prepare_time, subtract_time
    real(8) :: mult_time, dnrm2
    local_time = 0
    subtract_time = 0
    prepare_time = 0
    mult_time = 0
    rmax = 999
!############################################################################!
    allocate(umat%u(matrix%M)) !Allocate pointers
    allocate(umat%sz(matrix%M)) 
    umat%M = matrix%M 
    umat%sz(1:matrix%M)  = matrix%sz(1:matrix%M)
    do i = 1, matrix%M
        mem = mem + umat%sz(i)*umat%sz(i)
        mbc = matrix%i_block_csc(i+1) - matrix%i_block_csc(i)
        allocate(col_list(2,mbc))
        k=1
        do j = matrix%i_block_csc(i),matrix%i_block_csc(i+1) -1
            col_list(1,k) = matrix%sz(matrix%j_block_csc(j))
            col_list(2,k) = matrix%s_block_csc(j)
            k=k+1
        end do
        f_pos = matrix%rows(i)%far_block
        c_pos = matrix%rows(i)%close_block
        col_close_size = 0
        do j = 1,matrix%rows(i)%nclose
            col_close_size = col_close_size + matrix%sz(matrix%rows(i)%clos(j))
        end do
        call cpu_time(st)
        !print *, "I", i,matrix%sz(i),matrix%M, col_list(1,1)
        call compute_u_mult(matrix%sz(i), matrix%rows(i)%far_col,matrix%blocks(f_pos),matrix%rows(i)%close_col, &
        matrix%blocks(c_pos), & 
        matrix%blocks,f_pos,c_pos,umat%u(i)%p,r,r_fix,col_list,mbc,eps)
        !print *, "I", i, "r",r
        call cpu_time(ed)
        mult_time = mult_time + (ed - st)
        
        !Local eliminate (Cholesky + subtract from the blocks) 
        sz_new(i) = r !The new size
        !print *, "IN ELIM", i, r
        if (r < matrix%sz(i)) then
            allocate (L((matrix%sz(i)-r) , col_close_size)) 
            call cpu_time(st)
            !if (i .eq. 353) then
            !   !print *, i, matrix%sz(i), r, matrix%blocks(matrix%rows(i)%close_block),col_close_size    
            !end if
            !print *,i, matrix%sz(i),r, matrix%blocks(matrix%rows(i)%close_block),col_close_size 
            call local_eliminate(matrix%sz(i),r, matrix%blocks(matrix%rows(i)%close_block),col_close_size,L)
            call cpu_time(ed)
            st = ed - st
            local_time = local_time + st
            !######################## subtract #######################
            call cpu_time(st)
            pos_i = 1
            pos_j = 1
            do j = 1,matrix%rows(i)%nclose
                pos_i_n = pos_i + matrix%sz(matrix%rows(i)%clos(j))
                do k = 1,matrix%rows(i)%nclose
                    pos_j_n = pos_j + matrix%sz(matrix%rows(i)%clos(k))
                        !print *, i, j, k, 'sz0',matrix%sz(matrix%rows(i)%clos(j)), 'sz1', matrix%sz(matrix%rows(i)%clos(k))
                        if  (matrix%sz(matrix%rows(i)%clos(j)) .NE. 0) then
                        call dgemm('t','n', matrix%sz(matrix%rows(i)%clos(j)),matrix%sz(matrix%rows(i)%clos(k)),&
                                    matrix%sz(i)-r,-1d0,L(1, pos_i), matrix%sz(i)-r, &
                                    L(1, pos_j), matrix%sz(i)-r,1d0, &
                                    matrix%blocks(matrix%rows(i)%stencil(j,k)), matrix%sz(matrix%rows(i)%clos(j)))
                         pos_j = pos_j_n
                         end if
                     !end if
                end do
                pos_i = pos_i_n
                pos_j = 1
            end do
            call copy_L(matrix%sz(i),r,matrix%blocks(matrix%rows(i)%close_block),col_close_size,L) !The matrix L is saved in the
            !full matrix
            call cpu_time(ed)
            st = ed - st
            subtract_time = subtract_time + st
            deallocate(L)
       end if
       deallocate(col_list)
    end do
    !print *, 'prepare time:', prepare_time
    !print *, 'Mult time:', mult_time
    !print *, 'subtract time:', subtract_time
    !print *, 'local_eliminate time:', local_time
return
end subroutine
end module
