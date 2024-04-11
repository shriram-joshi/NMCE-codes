! p=15 : double, p=17 : extended, p=19 : quadruple
Module kind
    integer*4, parameter :: rk =selected_real_kind(p=15), rk2=selected_real_kind(p=15), ik=4
    !integer*4, parameter :: rk =selected_real_kind(p=17), rk2=selected_real_kind(p=15), ik=4
    integer*4, parameter :: ths=4
end Module kind
! p = 7 Single precision: 8 significant digits
! p = 15 Double precision: 16 significant digits
! p = 17 Extended precision: 19 significant digits
! p = 19 Quadratic precision: 32 significant digits
