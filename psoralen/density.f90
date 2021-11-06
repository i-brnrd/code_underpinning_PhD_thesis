subroutine density(x,y,z,zmax,layer)
  !***subroutine to setup density grid***!
  !***for a given x,y,z, returns layer***!
  !***all distances in CM
  implicit none

  real*8 x,y,z

  real*8 sine,pi
  real*8 a,w,b_epi,b_meli,b_basal,b_sc
  real*8 zmax
  real*8 um_to_cm
  integer layer


  !******calculate depths of layers from the surface of the skin
  a=0.003  !amplitude
  w=0.0015!

  um_to_cm=0.0001

  sine=a*(SIN(x/w)*COS(y/w))


  ! Epiderman thinkness Jane Sandy Moller ! Sc = 14.8 cellular epidermis= 83.7

  b_sc=zmax-(0.00296)
  b_epi=sine+(zmax-0.00637)
  b_meli=sine+(zmax-0.00737)
  b_basal=sine+(zmax-0.00837)



  !for given x,y,z; return layer
  if (z.ge.b_sc)then
     layer=1
  elseif (z.ge.b_epi) then
     layer=2
  elseif (z.ge.b_meli) then
     layer=3
  elseif (z.ge.b_basal)then
     layer=4
  else
     layer=5
  endif


  return
end subroutine density
