module optical_properties_init_mod
  implicit none
  save
contains
  !!!!!!N
   !!!!ns=1.39-((wl-279)*0.0001)
  subroutine optical_properties_init()
    use optical_properties_mod, ONLY: nwl, nlayer,lcount,l,u_a_all,u_s_all,g_all,n_all,wl_start
    use load_spec2_mod
    use search_bisec_mod
    use interpolate_mod
    implicit none

    real*8,dimension (nwl):: u_a_base, u_s_base,u_a_mel, u_a_epidermis,oxy_hgb,deoxy_hgb !base values from omlc "CHECK THIS CHECK THIS urgent****"
    real*8, dimension(nwl):: u_s_rayleigh
    real*8 spec(nwl),wl
    integer :: i_wl, layer
    real*8 :: jaques_basal

    character*70 :: filename
    real*8 lg10
    real*8 :: V_m, v_m_epi
    real*8 n_mel
    real*8 v_b,sO2, hgb



    do i_wl=1,nwl
      l(i_wl)= wl_start + real(i_wl) -1.
    enddo

        u_s_all=0.
        u_a_all=0.
        n_all=1.38

        g_all=0.
          lg10 = Log(10.)
          V_m=0.02 !p. 135 patricks thesis
          v_m_epi=0.04
          n_mel=(lcount(2)*v_m) ! number of 'melanin voxels' in epidermis
          v_b=0.02   !blood volume only 2%
          SO2=0.75
          hgb=150./66500.


      do i_wl=1,nwl
        wl=l(i_wl)
        g_all(i_wl)=0.62+0.29*wl*(10.**(-3))!
        u_a_base(i_wl)=7.84*(10.**8)*(wl**(-3.255))  !Jaques et al 1998 skin optics OMLC News
        u_s_base(i_wl)=1.752*(10.**8)*(wl**(-2.33)) + 134.67*(wl**(-0.494))  !Jaques et al 1998 skin optics OMLC News
        u_a_mel(i_wl)=6.6*(10.**11)*(wl**(-3.33))  !omlc

      enddo


    layer=1 !STRATUM CORNEUM
    !Absorption
    filename='./spectra/uvc/absorption_sc_uvc.txt' !Van Gemert et al, 1989 (taken from Everett....)
    call load_spec2(filename,spectrum=spec,unit_to_cm=10.d0)
    u_a_all(layer,:)=spec
    !scattering
    filename='./spectra/scattering_sc.txt'  !Van Gemert et al, 1989
    call load_spec2(filename,spec,10.d0)
    u_s_all(layer,:)=spec
    !u_s_all(layer,:)=u_s_rayleigh
    !--------------------------------------------------

    layer =2 !EPIDERMIS
    !Absorption
    filename='./spectra/uvc/absorption_epi_uvc.txt'!Van Gemert
    call load_spec2(filename,u_a_epidermis,10.d0)
    u_a_all(layer,:)=u_a_epidermis-(v_m_epi*u_a_mel)
    !Scattering
    filename='./spectra/scattering_epi_derm.txt' ! !Van Gemert
    call load_spec2(filename,spec,10.d0)
    u_s_all(layer,:)=spec
    !u_s_all(layer,:)=u_s_rayleigh
    !-------------------------------------------------------

    layer=3 !MELANIN LAYER
    !Absorption
    u_a_all(layer,:)= u_a_epidermis + (n_mel/real(lcount(3))) * u_a_mel
    !Scattering
    u_s_all(layer,:)=u_s_all(2,:)
    !-------------------------------------------------------

    layer=4 !DNA LAYER 
    u_a_all(layer,:)=u_a_all(2,:)
    !Scattering
    u_s_all(layer,:)=u_s_all(2,:)
    !-------------------------------------------------------
    layer = 5 !DERMIS
    !Absorption
    filename='./spectra/uvc/absorption_derm_uvc.txt'
    call load_spec2(filename,spec,10.d0)

    u_a_all(layer,:)=spec

    !scattering

    filename='./spectra/scattering_epi_derm.txt' !Van Gemert
    call load_spec2(filename,spec,10.d0)
    u_s_all(layer,:)=u_s_all(2,:)
    !u_s_all(layer,:)=u_s_rayleigh

    return
  END SUBROUTINE optical_properties_init

end module optical_properties_init_mod
