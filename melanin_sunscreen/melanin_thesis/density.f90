subroutine density(x,y,z,zmax,layer)
  !***subroutine to setup density grid***!
  !***for a given x,y,z, returns layer***!
  !***all distances in CM


  implicit none

  real*8 x,y,z

  real*8 sine,pi
  real*8 a,w
  real*8 b_epi_top, b_epi_mid, b_mel,b_basal,b_sc
  real*8 zmax
  real*8 um_to_cm
  integer layer

  !******calculate depths of layers from the surface of the skin
  !all distances in CM
  !**'wavelength' for sinosoidal undulation (dermal paillae)
  a=0.003  !amplitude
  w=0.0015!

  um_to_cm=0.0001

  sine=a*(SIN(x/w)*COS(y/w))

!Thicknesses TYPE I
! b_sc=zmax-(0.00148)
!
! b_epi_top=(2.d0/3.d0)*sine+(zmax-0.00565+a/2.d0)
! b_epi_mid=sine+(zmax-0.00765)
! b_mel=sine+(zmax-0.00865)
! b_basal=sine+(zmax-0.00965)



!Thicknesses TYPE VI
  b_sc=zmax-(0.00148)

  b_epi_top=(2.d0/3.d0)*sine+(zmax-0.00565+a/2.d0)
  b_epi_mid=sine+(zmax-0.00665)
  b_mel=sine+(zmax-0.00865)
  b_basal=sine+(zmax-0.00965)



    !for given x,y,z; return layer
  if (z.ge.b_sc)then
     layer=1
  elseif (z.ge.b_epi_top) then
     layer=2

  elseif (z.ge.b_epi_mid) then
     layer=3

  elseif (z.ge.b_mel) then
    layer=4



  elseif (z.ge.b_basal)then
     layer=5
  else
     layer=6
  endif


  return
end subroutine density
