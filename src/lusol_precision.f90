!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lusol_precision.f90
!
! LUSOL module for defining integer(ip), real(rp).
!
! ip  huge
!  4  2147483647
!  8  9223372036854775807
!
! rp  huge
!  4  3.40282347E+38
!  8  1.79769313486231571E+308
! 16  1.18973149535723176508575932662800702E+4932
!
! rp  eps
!  4  1.19209290E-07
!  8  2.22044604925031308E-016
! 16  1.92592994438723585305597794258492732E-0034
!
! We don't need selected_int_kind or selected_real_kind now.
! Previously we used these values:
! ip = integer precision    int_kind( 7) = integer(4)
!                           int_kind(15) = integer(8)
! rp = real precision      real_kind( 6) = real(4)   (not used in SNOPT)
!                          real_kind(15) = real(8)
!                          real_kind(30) = real(16)
!
! 11 Mar 2008: First version.
! 20 Apr 2012: First quad version.
! 22 Apr 2012: Made three versions of snPrecision.f90.
!              See README.QUAD for use with configure and make.
! 17 Jul 2013: Changed name to lusol_precision.f90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lusol_precision

  implicit none
  private
  public                :: ip, rp

  integer(4), parameter :: ip = 8, rp = 8

end module lusol_precision
