program mcsphere
implicit none
integer iseed,npk
integer i,j,ipk,nDepth
real r(3),u(3)
real pi
real ran
real rmag,rmax,l,taumax,a,nScatt,sCount

open(8, FILE='mc_iso.dat', status='replace')
! initalise variables
pi=4.*ATAN(1.)

iseed=-17284!keep this throughout-CHECK this is that simple!
! using cartesian coords r_x,r_y,r_z; with |r|=rmag
rmag=0.! |r(x,y,z)|
rmax=1.
a=1.!albedo: uniform density=>albedo not a function of r(x,y,z)

! npk: no. of photon packets per optical depth tau
!....... counted using ipk- the ith photon packet
! scount: total no. of scatterings per optical depth
!........counted using scount: scatter count
nDepth=20 !number of tau_maxes to sample
npk=20000
nScatt=0
call random_seed()
do j=1,nDepth
  print*,j
  ! set a max optical depth; and set scatter count sCount=0
  taumax=real(j)
  sCount=0
!       for this depth, simulate npk photons to be emitted
  do ipk=1,npk

    call emit_photon(iseed,pi,r,u)
    do !exit cond. in loop after |r| is updated
!            ...exit when |r|>rmax
!           get a new l (from a random tau), to find new position r
    call random_number(ran)
      l=-log(ran)/taumax*rmax
      do i=1,3
        r(i)=r(i)+l*u(i)
      enddo
!          use r to update rmag, and exit loop if rmag>rmax
      rmag=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
      if (rmag.gt.rmax) exit
!          if pkt is still in sphere, scatter isotropically
        call scatter(iseed,pi,u)
        sCount=sCount+1.
    enddo

enddo

  nScatt=sCount/(real(npk))
  write(8,*)taumax,nScatt
enddo
close(8)

end program

!-----SUBROUTINES------------------------------------------------------

subroutine emit_photon(iseed,pi,r,u)
implicit none
integer iseed, i
real r(3),u(3),pi
!     give packets initial position and direction
!     corresponding to emission of isotropic point source at origin
!     position is easy (x,y,z)=(0,0,0)
!     assuming velocity is same mag, any direction

do i=1,3
  r(i)=0.
enddo
call isotropic(iseed,pi,u)

return
end
!-------------
subroutine scatter(iseed,pi,u)
implicit none
integer iseed
real pi, u(3)
!     albedo not function of position so just return new direction
!     and increment scatter count

call isotropic(iseed,pi,u)

return
end
!------------
subroutine isotropic(iseed,pi,u)
implicit none
integer iseed
real pi, u(3),theta,phi,ran
call random_number(ran)
theta=ACOS((2.*ran-1.))

call random_number(ran)
phi=2.*pi*ran

u(1)=SIN(theta)*COS(phi)
u(2)=SIN(theta)*SIN(phi)
u(3)=COS(theta)

return
end
