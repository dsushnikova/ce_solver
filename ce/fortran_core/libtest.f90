module testmod
implicit none
contains

  subroutine testsub(arr, brr)
    integer, intent(in) :: arr
    integer, intent(out) :: brr
    brr = 3 + arr
    !print *, 'testsub executed'
  end subroutine testsub

end module testmod
