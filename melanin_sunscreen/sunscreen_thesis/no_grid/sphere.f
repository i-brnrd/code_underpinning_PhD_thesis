





c Monte Carlo code to simulate emission from an isotropic point source
c at origin of a uniform density sphere where the scattering albedo=1
      program mcsphere
      implicit none
      integer iseed,npk
      integer i,j,ipk,nDepth
      real r(3),u(3)
      real pi
      real ran2
      real rmag,rmax,l,taumax,a,nScatt,sCount

      open(8, FILE='mc_iso.dat', status='replace')
c initalise variables
      pi=4.*ATAN(1.)

      iseed=-17284!keep this throughout-CHECK this is that simple!
c using cartesian coords r_x,r_y,r_z; with |r|=rmag
      rmag=0.! |r(x,y,z)|
      rmax=1.
      a=1.!albedo: uniform density=>albedo not a function of r(x,y,z)

c npk: no. of photon packets per optical depth tau
c....... counted using ipk- the ith photon packet
c scount: total no. of scatterings per optical depth
c........counted using scount: scatter count

      nDepth=20 !number of tau_maxes to sample
      npk=20000
      nScatt=0

      do j=1,nDepth
        print*,j
c set a max optical depth; and set scatter count sCount=0
        taumax=real(j)
        sCount=0
c       for this depth, simulate npk photons to be emitted
        do ipk=1,npk

          call emit_photon(iseed,pi,r,u)
          do !exit cond. in loop after |r| is updated
c             ...exit when |r|>rmax
c           get a new l (from a random tau), to find new position r
            l=-log(ran2(iseed))/taumax*rmax
            do i=1,3
              r(i)=r(i)+l*u(i)
            enddo
c           use r to update rmag, and exit loop if rmag>rmax
            rmag=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            if (rmag.gt.rmax) exit
c           if pkt is still in sphere, scatter isotropically
            call scatter(iseed,pi,u)
            sCount=sCount+1.
          enddo

        enddo

        nScatt=sCount/(real(npk))
        write(8,*)taumax,nScatt
      enddo
      close(8)
      stop
      end
c-----SUBROUTINES------------------------------------------------------

      subroutine emit_photon(iseed,pi,r,u)
      implicit none
      integer iseed, i
      real r(3),u(3),pi
c     give packets initial position and direction
c     corresponding to emission of isotropic point source at origin
c     position is easy (x,y,z)=(0,0,0)
c     assuming velocity is same mag, any direction

      do i=1,3
         r(i)=0.
      enddo
      call isotropic(iseed,pi,u)

      return
      end
c-------------
      subroutine scatter(iseed,pi,u)
      implicit none
      integer iseed
      real pi, u(3)
c     albedo not function of position so just return new direction
c     and increment scatter count

      call isotropic(iseed,pi,u)

      return
      end
c------------
      subroutine isotropic(iseed,pi,u)
      implicit none
      integer iseed
      real pi, u(3),theta,phi,ran2

      theta=ACOS((2.*ran2(iseed))-1.)
      phi=2.*pi*ran2(iseed)

      u(1)=SIN(theta)*COS(phi)
      u(2)=SIN(theta)*SIN(phi)
      u(3)=COS(theta)

      return
      end
