module grid_mod
  use optical_properties_mod, only:nwl,nlayer
  implicit none
  save
  !contains
  integer, parameter :: nxg=50,nyg=50,nzg=400
  !integer, parameter :: nxg=100,nyg=100,nzg=100
  real*8 :: xface(-3:nxg+4),yface(-3:nyg+4),zface(-3:nzg+4)
  real*8 :: MASK(nxg,nyg,nzg,nlayer)
  real*8 :: rhokap(nxg,nyg,nzg)

  REAL*8 :: PL_SUM(nwl,NXG,NYG,NZG)

  real*8,parameter :: xmax=0.025,ymax=0.025,zmax=0.02
  !real*8,parameter :: xmax=0.05,ymax=0.05,zmax=0.05
  real*8 :: grid_max(3)
  real*8 :: delta

end module grid_mod