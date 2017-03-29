subroutine write_to_file(a)
    implicit none
    real(8) :: a(5)
    open(unit=1,file="data_fd.bin",form="Unformatted",status="OLD",action="write")
    a(1:5) = 2.0
    write(1) a(1:5)
    print *, a(1:5)
    close(1)
end subroutine write_to_file

subroutine test(arr, brr)
    integer :: arr
    integer :: brr
    brr = 3 + arr
    !print *, 'testsub executed'
 end subroutine test


module prec
contains 
subroutine free_mem()
    use fast_direct
    if (level_global == 0 ) then
        print *, "First in"
    else
        deallocate(matrix_global)
        deallocate(umat_global)
        deallocate(Lmat_global)
    end if
end subroutine free_mem
subroutine factor_sol(n,M,nnz,ia,ja,sa,x,level,eps,r,sol,bl_s,bl_ar,join, &
                      tot_hyp,hyper_join,m_j,mem,f_time,sol_time,res,error)
    use fast_direct
    implicit none
    integer, intent(in) :: n, nnz, M, join
    integer, intent(in) :: ia(n+1), ja(nnz), bl_ar(M)
    real(8), intent(in) :: sa(nnz)
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: sol(n)
    integer, intent(out) :: mem
    real(8), intent(out) :: res,error
    real(8), intent(out) :: f_time,sol_time
    integer, intent(in) :: level, r, bl_s, tot_hyp
    integer, intent(in) :: m_j(level), hyper_join(tot_hyp)
    real(8), intent(in) :: eps
    real(8) :: rhs(n)
    !print *,"In prec", m_j
    !print *,"In prec", hyper_join
    !print *, bl_ar
    !print *, join
    !integer :: mem
    rhs = x
    !print *, level
    !print *,eps
    !level = 4
    !eps = 1e-2
    !call fact_sol(n,nnz,ia,ja,sa,x,level,eps,sol)
    call solve_problem(n, M, nnz, ia, ja, sa,rhs,level,eps,r,sol,bl_s,bl_ar, &
                       join,tot_hyp, hyper_join,m_j,mem,f_time,sol_time,res,error)
    !print *, 'Memory:', mem, 'MB'
    !print *, "done"
end subroutine factor_sol
subroutine factor_ll(n,M,nnz,ia,ja,sa,level,eps,r,bl_s,bl_ar,join,tot_hyp,hyper_join,m_j,mem,f_time)
    use fast_direct
    implicit none
    integer, intent(in) :: n, nnz, M, join
    integer, intent(in) :: ia(n+1),ja(nnz),bl_ar(M)
    real(8), intent(in) :: sa(nnz)
    integer, intent(in) :: level, r, bl_s, tot_hyp
    integer, intent(in) :: m_j(level), hyper_join(tot_hyp)
    real(8), intent(in) :: eps
    integer , intent(out):: mem
    real(8), intent(out) :: f_time
    !print *, "level", level
    !print *,"In ll", tot_hyp
    !print *,"In ll", m_j
    !print *,"In ll", hyper_join
    call fact_ll(n,M,nnz,ia,ja,sa,level,eps,r,bl_s,bl_ar,join,tot_hyp,hyper_join,m_j,mem,f_time)
    !print *, 'Memory:', mem, 'Mb'
!    !print *,M_global
!    !call precond(n, ia, ja, sa)
!    !call test_test(n,nnz, ia, ja, sa)
end subroutine factor_ll
subroutine solve_ll(n,x,sol)
    use fast_direct
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8), intent(out) :: sol(n)
    real(8) :: rhs(n)
    rhs = x
    call solve_system(n, level_global, matrix_global, Lmat_global, umat_global, rhs, sol)
!    !print *, matrix_global(1)%M
!    call sol_ll(n,x,sol)
!    !print *,M_global(1)
!    !call precond(n, ia, ja, sa)
!    !call test_test(n,nnz, ia, ja, sa)
end subroutine solve_ll

end module prec
