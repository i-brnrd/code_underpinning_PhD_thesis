module optical_properties_init_mod
  implicit none
  save
contains

  subroutine optical_properties_init()
    use optical_properties_mod, ONLY: nwl, nlayer,lcount,l,u_a_all,u_s_all,g_all,n_all,wl_start
    use load_spec2_mod
    use search_bisec_mod
    use interpolate_mod
    implicit none

    real*8,dimension (nwl):: u_a_base, u_s_base,u_a_mel, u_a_epidermis,oxy_hgb,deoxy_hgb !base values from omlc
    real*8 spec(nwl),wl

    real*8 nonmel_epi(nwl)
    integer :: i_wl, layer
    real*8 :: jaques_basal

    character*70 :: filename
    real*8 lg10
    real*8 :: V_m, v_m_epi
    real*8 n_mel !'how many' cells
    real*8 v_b,sO2, hgb

    real*8 upper_frac, middle_frac, lower_frac



    do i_wl=1,nwl
      l(i_wl)= wl_start + real(i_wl) -1.
    enddo


        u_s_all=0.
        u_a_all=0.
        n_all=1.38

        g_all=0.
          lg10 = Log(10.)

          v_m_epi=0.04

          v_b=0.02   !blood volume only 2%
          SO2=0.75
          hgb=150./66500.



!MELANIN!!!
!volumetric calculation of melanin!!
!**************************



!values from Dami's FASEB paper
!for type VI

!*****************************************
! !Type 1
! upper_frac=0.0
! middle_frac=0.0
! lower_frac=1.0
! v_m=0.03 !total melanin fraction

!Type VI
! upper_frac=0.05
! middle_frac=0.15
! lower_frac=0.8
! v_m=0.03 !total melanin fraction
!v_m=v_m*8.d0

! !Type I v2
! upper_frac=0.05
! middle_frac=0.15
! lower_frac=0.8
! v_m=0.03 !total melanin fraction

! !Type I v2
upper_frac=0.13
middle_frac=0.27
lower_frac=0.60
v_m=0.03 !total melanin fraction
v_m=v_m*8.d0








n_mel=((lcount(2)+lcount(3)+lcount(4))*v_m) ! number of 'melanin containing voxels' in total epidermis


print*, 'layers', lcount




!print*,'TYPE VI', v_m, upper_frac,middle_frac,lower_frac
print*, n_mel
print*, n_mel*(upper_frac+middle_frac+lower_frac)


      do i_wl=1,nwl
        wl=l(i_wl)
        g_all(i_wl)=0.62+0.29*wl*(10.**(-3))!
        u_a_base(i_wl)=7.84*(10.**8)*(wl**(-3.255))  ! Jaques et al 1998 skin optics OMLC News
        u_s_base(i_wl)=1.752*(10.**8)*(wl**(-2.33)) + 134.67*(wl**(-0.494))  ! Jaques et al 1998 skin optics OMLC News
        u_a_mel(i_wl)=6.6*(10.**11)*(wl**(-3.33))  !omlc

      nonmel_epi=0.d0
      filename='./spectra/absorption_epi.txt'!Van Gemert
      call load_spec2(filename,u_a_epidermis,10.d0)
      nonmel_epi=u_a_epidermis-(v_m_epi*u_a_mel)





    layer=1 !STRATUM CORNEUM
    !Absorption
    filename='./spectra/absorption_sc.txt' !Van Gemert et al, 1989 (taken from Everett....)
    call load_spec2(filename,spectrum=spec,unit_to_cm=10.d0)
    u_a_all(layer,:)=spec
    !scattering
    filename='./spectra/scattering_sc.txt'  !Van Gemert et al, 1989
    call load_spec2(filename,spec,10.d0)
    u_s_all(layer,:)=spec
    !--------------------------------------------------

    layer =2 !EPIDERMIS (upper)
    !Absorption
    u_a_all(layer,:)=nonmel_epi + (n_mel*upper_frac/real(lcount(2)))*u_a_mel

    !Scattering
    filename='./spectra/scattering_epi.txt' !Van Gemert
    call load_spec2(filename,spec,10.d0)
    u_s_all(layer,:)=spec


    !-------------------------------------------------------

    layer=3 !EPIDERMIS (middle)
    !Absorption
    u_a_all(layer,:)=nonmel_epi + (n_mel*middle_frac/real(lcount(3)))*u_a_mel
    u_s_all(layer,:)=u_s_all(2,:)
    !-------------------------------------------------------

    layer=4 !Melanin (bottom)
    !Absorption
    u_a_all(layer,:)=u_a_epidermis + (n_mel*lower_frac/real(lcount(4))) * u_a_mel
    !Scattering
    u_s_all(layer,:)=u_s_all(2,:)


    layer=5 !DNA LAYER DNA: units?? is it out by 1000???
    !Absorption
    !filename='./spectra/dna_molar_extinct.txt' !Mouret et al...
    !call load_spec2(filename,spec,1.d0)
    !u_a_all(layer,:)=lg10*0.0185*spec
    u_a_all(layer,:)=u_a_epidermis
    !Scattering
    u_s_all(layer,:)=u_s_all(2,:)
    !-------------------------------------------------------
    layer = 6 !DERMIS
    !Absorption
    filename='./spectra/ohb_as.txt'
    call load_spec2(filename,spec,1.d0)
    oxy_hgb=spec

    filename='./spectra/dhb_as.txt'
    call load_spec2(filename,spec,1.d0)
    deoxy_hgb=spec
    u_a_all(layer,:)=(v_b*hgb*(SO2*oxy_hgb + (1.d0-sO2)*deoxy_hgb)*lg10)+u_a_base*(1.d0-v_b)

    !scattering

    filename='./spectra/scattering_epi_derm.txt' !Van Gemert
    call load_spec2(filename,spec,10.d0)
    !u_s_all(layer,:)=u_s_all(2,:)
    u_s_all(layer,:)=spec

    return
  END SUBROUTINE optical_properties_init

end module optical_properties_init_mod
