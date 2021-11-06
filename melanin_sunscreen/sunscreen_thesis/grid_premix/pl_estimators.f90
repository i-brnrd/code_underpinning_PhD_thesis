  module pl_estimators_mod
  implicit none
  save
contains


!--------------------------------------------------------------------------------
  SUBROUTINE PL_ESTIMATORS(pkt_count,incident_spec_irr)
    use optical_properties_mod, only: n_pkt_wl, l, nwl, wl_start
    use grid_mod, only: nxg,nyg,nzg,xmax,ymax,zmax, PL_SUM,MASK
    use constants_mod
    implicit none

    integer, intent(in) :: pkt_count
    real*8, intent(in) :: incident_spec_irr(nwl)

    real*8 :: Area


    real*8 :: wls(nwl)
    real*8 :: lumin(nwl)
    real*8 :: V
    real*8, dimension(nxg,nyg,nzg) :: j_mean
    real*8 :: l_norm, l_sum
    integer :: i,j,k,m
    integer:: n_basal, basal
    real*8 :: depth_val(nzg)
    real*8 :: picture_depths(nwl,nzg)
    real*8 :: fluence_basal(nwl)


    print*, nxg,nyg,nzg

    wls=l
    lumin=incident_spec_irr
    l_sum=sum(incident_spec_irr)
    print*, '****'
    !print*, lumin
    Area=(2.d0*(xmax))*(2.d0*(ymax))
    V= (2.d0*(XMAX/NXG))*(2.d0*(YMAX/NYG))*(2.d0*(ZMAX/NZG))

    print*,'Area',Area, ' in cm**2'
    !L=L/(100.d0*100.d0) !W/cm
    !gprint*,'Luminosity',L,'W/m^2'
    print*,'Volume of Voxel',v
    print*, 'cells', nxg,nyg,nzg

    j_mean=0.d0
    do m=1,nwl
      if (n_pkt_wl(m).eq.0) cycle
      l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
      do i=1,nxg
        do j=1,nyg
          do k=1,nzg
            j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))*l_norm
          enddo
        enddo
      enddo
    enddo

    open(10,file='./data_out/depths1.dat',status='replace')
    depth_val=0.0d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
        enddo
      enddo
      depth_val(k)=depth_val(k)/(nxg*nyg)
      write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    !print*, (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    enddo
    close(10)

    do m=1,nwl
      if (n_pkt_wl(m).eq.0) then
        PL_SUM(m,:,:,:)=0.d0
      else
        l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
        PL_SUM(m,:,:,:)=l_norm*PL_SUM(m,:,:,:)
      endif
    enddo

    print*, 'Writing out the MASK'
    open(10,file='./data_out/densitygrid.dat',status='replace',form='unformatted')
    write(10) MASK
    close(10)


    print*, '...ABOUT TO WRITE TO THE PL FILE'
      open(10,file='./data_out/pathlengths.dat',status='replace',form='unformatted')

        write(10)  PL_SUM

    close(10)
  print*, 'written it out!!!!'

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
      write(11,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    enddo

    print*, '...ABOUT TO WRITE TO THE PICTURE FILE'
    open(10,file='./data_out/wavelengths_depths.dat',status='replace')
    do k=1,NZG
      write(10,*)  picture_depths(:,k)
    enddo
    close(10)
    print*, 'written it out!!!!'



    basal=5
    print*, '...ABOUT TO WRITE TO THE BASAL FILE layer::',basal
    fluence_basal=0.d0
    n_basal=0
     do m=1,nwl
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

    fluence_basal=fluence_basal/n_basal
    print*, '...ABOUT TO WRITE TO THE basal incident FILE'
    open(10,file='./data_out/basal.dat',status='replace')

    do m=1,nwl
      write(10,*) wl_start+real(m-1),fluence_basal(m)
    enddo

    print*, 'written it out!!!!'


!------------------------------------------------------------------

    basal=4
    print*, '...ABOUT TO WRITE TO THE MID FILE layer::',basal
    fluence_basal=0.d0
    n_basal=0
     do m=1,nwl
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

    fluence_basal=fluence_basal/n_basal
    print*, '...ABOUT TO WRITE TO THE basal incident FILE'
    open(10,file='./data_out/mid.dat',status='replace')

    do m=1,nwl
      write(10,*) wl_start+real(m-1),fluence_basal(m)
    enddo

    print*, 'written it out!!!!'





!------------------------------------------------------------------

basal=3
print*, '...ABOUT TO WRITE TO THE MID FILE layer::',basal
fluence_basal=0.d0
n_basal=0
 do m=1,nwl
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

fluence_basal=fluence_basal/n_basal
print*, '...ABOUT TO WRITE TO THE basal incident FILE'
open(10,file='./data_out/upper.dat',status='replace')

do m=1,nwl
  write(10,*) wl_start+real(m-1),fluence_basal(m)
enddo

print*, 'written it out!!!!'







    RETURN
  END SUBROUTINE PL_ESTIMATORS

end module pl_estimators_mod
