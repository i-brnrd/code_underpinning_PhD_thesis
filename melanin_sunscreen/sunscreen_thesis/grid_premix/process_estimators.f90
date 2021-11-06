!ok trying to replicate paper outputs:::
program process_estimators
use optical_properties_mod
use optical_properties_init_mod
use load_spec2_mod
use search_bisec_mod
use interpolate_mod

use grid_mod, only: nxg,nyg,nzg, PL_SUM,MASK,zmax
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

real*8 :: picture_depths(nwl,nzg)

PL_SUM=0.d0
MASK=0.d0


energy_uva=0.d0
energy_uvb=0.d0
energy_total=0.d0


!---MASK EXTRACTION

! !extracting some mask slices
!
!     !call optical_properties_init()
!
!     open(8,file='optical_properties_check.txt',status='unknown')
!     do i=1,nwl
!       read(8,*) l(i),u_s_all(1:5,i), u_a_all(1:5,i),g_all(i)
!     enddo
!     close(8)
!
!
     open(10, file='./data_out/densitygrid.dat',status='unknown',form='unformatted')
      read(10) MASK
     close(10)
!
!
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

print*, mask_snapshot(1,1)
print*, '...writing snapshot of density...'
open(10,file='./data_out/density.dat',status='replace')
do k=1,NZG
  write(10,*)  mask_snapshot(:,k)
enddo
close(10)
print*, 'written it out!!!!'

!printing out the slice pictures
print*, 'nlayer', nlayer
!

! open(10, file='./data_out/.dat',status='unknown',form='unformatted')
!  read(10) PL_SUM
! close(10)
!
! picture_slice=0.d0
! abs=0.d0
! j=50
! do k=1,nzg
!     do i=1,nxg
!         do m=1,40
!           do p=1,nlayer
!             if (MASK(i,j,k,p)==1) then
!
!               picture_slice(i,k) = picture_slice(i,k) + u_a_all(p,m)*PL_SUM(m,i,j,k)
!             endif
!           enddo
!         enddo
!     enddo
! enddo
!
! print*, '...ABOUT TO WRITE TO THE SLICE FILE'
! open(10,file='./data_out/slice/sb_I_uvb.dat',status='replace')
! do k=1,NZG
!     write(10,*)  picture_slice(:,k)
! enddo
! close(10)
! print*, 'written it out!!!!'
!
!
!
!
!   do CASE=1,2
!
!      if (case==1) then !SOLAR
!        open(10, file='./data_out/sol_II.dat',status='unknown',form='unformatted')
!         read(10) PL_SUM
!        close(10)
!      else !SUNBED
!        open(10, file='./data_out/sb_II.dat',status='unknown',form='unformatted')
!         read(10) PL_SUM
!        close(10)
!      end if
!
! !reinititalise op_props so we can extract the DNA abs
!     print*, h,c, h*c
!
! !------BASAL EXTRACTION---------------------------------
! basal=4
!
! energy_basal=0.d0
! n_basal=0
!  do m=1,nwl
!    wl= wl_start+real(m-1)
!    do i=1,nxg
!      do j=1,nyg
!        do k=1,nzg
!          if (MASK(i,j,k,basal)==1) then
!            n_basal=n_basal+1
!            !energy_norm= ((10.**(-9))/(h*c*wl))*dna_abs(m)
!            energy_basal(m)=energy_basal(m)+ dna_abs(m)*PL_SUM(m,i,j,k)
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!
! !split up into uva & energy_uvb
! energy_uva(case)=0.d0
! energy_uvb(case)=0.d0
! energy_total(case)=0.d0
!
! do m=1,nwl
!    wl= wl_start+real(m-1)
!
!    if (wl.le.315) then
!      energy_uvb(case)=energy_uvb(case)+energy_basal(m)
!    else
!      energy_uva(case)=energy_uva(case)+energy_basal(m)
!    endif
! end do
!
! energy_total(case)=energy_uvb(case)+energy_uva(case)
! print*, case
! enddo
!
!
!
! print*, energy_uva, energy_uvb
! print*,energy_uvb(2)/energy_uvb(1),energy_uva(2)/energy_uva(1)!,energy_total
!
!
filename='./data_out/MEL_RES/type6.dat'


print*, '...ABOUT TO calculate basal incident fluence....'
open(10, file=filename,status='unknown',form='unformatted')
 read(10) PL_SUM
close(10)

basal=5

fluence_basal=0.d0
n_basal=0
 do m=1,nwl
   wl= wl_start+real(m-1)
   do i=1,nxg
     do j=1,nyg

       do k=1,nzg
         if (MASK(i,j,k,basal)==1) then
          ! print*, k
           !n_basal=n_basal+1
           fluence_basal(m)=fluence_basal(m)+ PL_SUM(m,i,j,k)
         endif
       enddo
     enddo
   enddo
 enddo

print*, '...ABOUT TO WRITE TO THE basal incident FILE'
open(10,file='./data_out/MEL_RES/basal_t6.dat',status='replace')

do m=1,nwl
  write(10,*) wl_start+real(m-1),fluence_basal(m)

enddo
do m=1,nwl
  print*, wl_start+real(m-1), fluence_basal(m)
  enddo
close(10)
print*, 'written it out!!!!'


STOP
    picture_depths=0.d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          do m=1,nwl
          picture_depths(m,k)= picture_depths(m,k) + PL_SUM(m,i,j,k)
          enddo
        enddo
      enddo
      picture_depths(:,k)=picture_depths(:,k)/(nxg*nyg)
      !write(11,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    enddo

    print*, '...ABOUT TO WRITE TO THE PICTURE FILE'
    open(10,file='./data_out/mel_results/mel_6.dat',status='replace')
    do k=1,NZG
      write(10,*)  picture_depths(:,k)
    enddo
    close(10)
    print*, 'written it out!!!!'
!Oh lordy

STOP

  basal=5
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
  open(10,file='./data_out/mel_results/basal_type6_50_v2.dat',status='replace')
  do m=1,nwl
    write(10,*) wl_start+real(m-1),fluence_basal(m)
  enddo

  open(10,file='./data_out/mel_results/top_type6_50.dat',status='replace')
  do m=1,nwl
    write(10,*) wl_start+real(m-1),top_epi(m)
  enddo


  open(10,file='./data_out/mel_results/mid_type6_50.dat',status='replace')
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
