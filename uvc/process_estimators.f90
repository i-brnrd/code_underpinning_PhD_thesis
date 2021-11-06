!ok trying to replicate paper outputs:::
program process_estimators
use optical_properties_mod
use optical_properties_init_mod
use load_spec2_mod
use search_bisec_mod
use interpolate_mod

use grid_mod, only: nxg,nyg,nzg,xface,yface,zface,xmax,ymax,zmax, PL_SUM,MASK,zmax
use constants_mod

implicit none
real*8, dimension(nxg,nyg,nzg) :: j_mean
integer :: i,j,k,p,m
real*8 :: depth_val(nzg)

integer :: mask_snapshot(nyg,nzg)

integer :: CASE

real*8 :: dna_abs(nwl)
real*8 spec(nwl)

character*70 :: filename
integer :: basal, n_basal
real*8 :: wl
real*8 :: energy_norm
real*8 :: energy_basal(nwl),fluence_basal(nwl)

real*8 :: top_epi(nwl), mid_epi(nwl)

real*8 :: energy_uva(2)
real*8 :: energy_uvb(2)
real*8 :: energy_total(2)

real*8 :: picture_slice(nxg,nzg)
real*8 :: abs

real*8 :: x,y,z
real*8 :: sine,a,w

PL_SUM=0.d0
MASK=0.d0


energy_uva=0.d0
energy_uvb=0.d0
energy_total=0.d0


!---MASK EXTRACTION

!extracting some mask slices

    !call optical_properties_init()

    !open(8,file='optical_properties_check.txt',status='unknown')
    !do i=1,nwl
    !  read(8,*) l(i),u_s_all(1:5,i), u_a_all(1:5,i),g_all(i)
    !enddo
  !  close(8)
  a=0.003  !amplitude
  w=0.0015!



! !  b_meli=sine+(zmax-0.00737)
! !  b_basal=sine+(zmax-0.00837)
!   do i=1,nxg
!      do j=1,nyg
!         do k=1,nzg
!
!            x=xface(i)-xmax+xmax/nxg
!            y=yface(j)-ymax+ymax/nyg
!            z=zface(k)-zmax+zmax/nzg
!
!            sine=a*(SIN(x/w)*COS(y/w))
!            !check this condition....
!            if ((z.lt.(sine+(zmax-0.00737))).and. (z.ge.(sine+(zmax-0.00747))))then
!              print*, i,j,k
!           endif
!
!
!
!         end do
!      end do
!   end do
!
! print*, 'went through the xfacey thing'
!
!
!
!
!
! STOP
    open(10, file='./data_out/densitygrid.dat',status='unknown',form='unformatted')
     read(10) MASK
    close(10)


i=50
do j=1,NYG
  do k=1,NZG
    do p=1,nlayer
      if (MASK(i,j,k,p)==1) then
          mask_snapshot(j,k)=p
      endif
    end do
  end do
enddo

open(10, file='./data_out/pic.dat')
do k=1,nzg
  write(10,*) mask_snapshot(:,k)
enddo
  close(10)



print*, '...ABOUT TO calculateTO THE basal incident fluence....'
open(10, file='./data_out/pathlengths.dat',status='unknown',form='unformatted')
 read(10) PL_SUM
close(10)

PL_SUM=PL_SUM/10000.d0  !luminosity was in m^2; need to put it in cm^2

basal=4

fluence_basal=0.d0
n_basal=0
 do m=1,nwl
   !wl= wl_start+real(m-1)
   do i=1,nxg
     do j=1,nyg

       do k=1,nzg
         if (MASK(i,j,k,basal)==1) then
           n_basal=n_basal+1
           fluence_basal(m)=fluence_basal(m)+ PL_SUM(m,i,j,k)
         endif
       enddo
     enddo
   enddo
 enddo
print*, n_basal
print*, '...ABOUT TO WRITE TO THE basal incident FILE'
open(10,file='./data_out/basal.dat',status='replace')

do m=1,nwl
  write(10,*) wl_start+real(m-1),fluence_basal(m)

enddo
close(10)
print*, 'written it out!!!!'








STOP
!Oh lordy


! does my code actually work lols

  open(10, file='./data_out/direct_unity.dat',status='unknown',form='unformatted')
   read(10) PL_SUM
  close(10)


  basal=4
  top_epi=0.d0
  mid_epi=0.d0
  fluence_basal=0.d0
  n_basal=0
   do m=1,nwl
     wl= wl_start+real(m-1)
     do i=1,nxg
       do j=1,nyg
         do k=1,nzg
           if (MASK(i,j,k,basal)==1) then
             !n_basal=n_basal+1
             fluence_basal(m)=fluence_basal(m)+ PL_SUM(m,i,j,k)
           endif
           if (k==100) then !very top of epidemis!
             top_epi(m)= top_epi(m) + PL_SUM(m,i,j,k)
           endif

           if (k==97) then !very top of epidemis!

             mid_epi(m)= mid_epi(m) + PL_SUM(m,i,j,k)
           endif

            enddo
         enddo
       enddo
     enddo

  print*, '...ABOUT TO WRITE TO THE final files....'
  open(10,file='./data_out/basal_datasets/direct_basal.dat',status='replace')
  do m=1,nwl
    write(10,*) wl_start+real(m-1),fluence_basal(m)
  enddo

  open(10,file='./data_out/basal_datasets/direct_top.dat',status='replace')
  do m=1,nwl
    write(10,*) wl_start+real(m-1),top_epi(m)
  enddo


  open(10,file='./data_out/basal_datasets/direct_mid.dat',status='replace')
  do m=1,nwl
    write(10,*) wl_start+real(m-1),mid_epi(m)
  enddo


  print*, 'written it out!!!!'




















STOP
!   print*,l_sum
!   ! print out values at each layer
! !Total jmean over all wavelengths.....
    do i=1,nxg
      do j=1,nyg
        do k=1,nzg
          j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))
        enddo
      enddo
    enddo

    open(10,file='depths_from_postprocess.dat',status='replace')

    depth_val=0.0d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
        enddo
      enddo
      depth_val(k)=depth_val(k)/(nxg*nyg)
      write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)

    enddo

    close(10)
end program
