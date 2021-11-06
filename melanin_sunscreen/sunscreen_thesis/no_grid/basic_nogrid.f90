! MCRT code to swiftly simulate transmission through a sunscreen film
program mcsphere
  !use get_dim_mod
!  use load_mod
!  use load_spec_mod
!  use search_bisec_mod
!
implicit none
integer iseed,npk,sz,ipk,nDepth
integer i,j
real r(3),u(3)
real pi
real ran
real rmag,zmax,l,taumax,a,nScatt,sCount

character*30 :: fname
real*8, allocatable :: input_spec(:,:),Lambda(:),cdf(:),transmitted(:) !

pi=4.*ATAN(1.)

zmax=1.
a=0.!albedo: pure absorber


fname='e_eff_sim.txt'

!call load_spec(fname,sz,input_spec)

allocate(lambda(sz),cdf(sz),transmitted(sz))
lambda(:)=input_spec(1,:)
!call get_cdf(cdf,input_spec,sz)


do i=1,sz
   print*, lambda(i),input_spec(1,i),cdf(i)
enddo

stop
      npk=10.
      nScatt=0

      do j=1,nDepth
        print*,j
! set a max optical depth; and set scatter count sCount=0
        taumax=real(j)
        sCount=0
!       for this depth, simulate npk photons to be emitted
        do ipk=1,npk

          call emit_photon(iseed,pi,r,u)
          do !exit cond. in loop after |r| is updated
!             ...exit when |r|>rmax
             !           get a new l (from a random tau), to find new position r
             call random_number(ran)
            l=-log(ran)!/taumax*rmax
            do i=1,3
              r(i)=r(i)+l*u(i)
            enddo
!           use r to update rmag, and exit loop if rmag>rmax
            if (r(3).gt.zmax) exit
!           if pkt is still in film..
          enddo

        enddo

        !nScatt=sCount/(real(npk))
        !write(8,*)taumax,nScatt
      enddo
      close(8)
      stop
      end
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
!    albedo not function of position so just return new direction
!    and increment scatter count

      call isotropic(iseed,pi,u)

      return
      end
!------------
      subroutine isotropic(iseed,pi,u)
      implicit none
      integer iseed
      real pi, u(3),theta,phi,ran

      call random_number(ran)
      theta=ACOS((2.*ran)-1.)
      call random_number(ran)
      phi=2.*pi*ran

      u(1)=SIN(theta)*COS(phi)
      u(2)=SIN(theta)*SIN(phi)
      u(3)=COS(theta)

      return
      end
