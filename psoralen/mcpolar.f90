!************************************************************
program mcpolar
  include 'modules_full_list.txt'
  implicit none
  include 'photon.txt'

  !*****Parameter declarations ****************************************

  integer nphotons   !number of photon packets sent into simulation
  real*8 xmax,ymax,zmax
  real*8 diff !diffuse fractions of incident light
  character*30 :: fname_incident_irradiation
  real*8, parameter :: pi = 4.*atan(1.)
  real*8, parameter:: twopi=2.*pi
  real*8, parameter:: fourpi=4.*pi

  !initialise input arrays
  real*8 incident_spec_irr(nwl)
  real*8 tot_irr
  real*8 :: cdf(nwl)
  integer pkt_count,scatter_count   !Counts packets that actually travel through medium (not reflected)

  real*8 hgg,g2  !Henyey Greenstein phase function variables
  real*8 delta
  real*8 u_s_old(nlayer),u_a_old(nlayer),g
  real*8 ns ! refractive index of skin

 ! (loop) indices, counters and seeds & FLAGS
  integer iseed
  integer j
  integer xcell,ycell,zcell
  integer tflag, seg_flag, r_flag
  INTEGER I,k,p
  integer d_flag
  integer sz
  integer b_wl

  real*8 kappa,albedo,pl,pc,sc
  real*8 e

  real*8 WL
  real*8 dnaval
  real*8 ran
  integer  A,B
  integer CASE
  real*8 :: tot(2),uva(2)

  print*,'**********************'
  print*,'Simulation of MC-UVRT through upper layers of skin'

  include 'parameters.txt'
  diff=1.d0
  fname_incident_irradiation='./spectra/bb_uva.csv'

  print*, 'model size(cm)', xmax,ymax,zmax
  print*, 'sim in',nlayer, ' layers'
  print*, 'using photons', nphotons

  call iarray !Initialise the arrays shared by the modules to 0

  ! Initialise counter
    pkt_count=0
  !Skin parameters/ Optical Properties parameters/ Tissue parameters
    ns=1.38
    kappa=1.
    pl=0.
    pc=0.
    sc=1.

!INITIALISE GRID
    call gridset(xmax,ymax,zmax,kappa)
!INITIALISE PACKETS OPTICAL PROPERTIES
print*, wl_start
do i=1,nwl
  l(i)= wl_start + real(i) -1.
enddo

 if (nwl.gt.1) then

  !**Initialise Spectra**
    call random_seed()
    !Check sum, luminosity, and CDF all work ok.
    call load_spec2(fname_incident_irradiation,incident_spec_irr,1.d0)
    call get_cdf(cdf,l,incident_spec_irr,nwl)
    tot_irr=sum(incident_spec_irr)
    print*,'Incident spectral irradiance loaded from ', fname_incident_irradiation
    print*,'Total irradiance:',tot_irr
    print*,'Diffuse Fraction',diff


   else
     !JAQUES VERIFICATION
     print*, 'Jaques verification mode on nwls:', nwl
     Diff = 0.d0
     tot_irr=1.d0
     l(1)= 630.d0
     incident_spec_irr(1)=1.d0
     print *, tot_irr, l
     print*, 'diffuse fraction', diff
     call random_seed()
     g=0.9
     wl=630.d0
   do i=1,nlayer
     u_a_old(i)=0.23
     u_s_old(i)=(21.)/(1.-g)
  enddo

 endif

    call op_prop_set(u_a,u_s,g_skin)

    print*,'layer 2 cell count',lcount(2),'layer 3 cell count',lcount(3)
!---------------------------------------------------------------------------------------------------------------
!MCRT
!-------------------------------------------------------------------------------------------------------------------------

     !*****Set small distance for use in optical depth integration routines
     !*****for roundoff effects when crossing cell walls
     delta=1.e-6*(2.*zmax/nzg)
!!$
     !**** Loop over nph photons from each source *************************
     scatter_count=0

     do j=1,nphotons
        if(mod(j,100000).eq.0)then
           print *, j,' scattered photons completed'
        end if

        !*****Release photon from point source *******************************
        call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
             fi,fq,fu,fv,xmax,ymax,zmax,twopi,xcell,ycell,zcell,iseed,L,&
             cdf,DELTA,WL,diff,d_flag)


        if (d_flag.eq.1) then  !if diffuse photon, check fresnel
           !direct photons ~2% misdiagnosed as diffuse ...
           call reflect(1.d0,ns,cost,pi,r_flag) !fresnel reflection
           if (r_flag.eq.1) then
              CYCLE
           endif
        endif

        pkt_count=pkt_count+1

        seg_flag=0
        e=1./wl

        !bin wavelength

        if (nwl.gt.1) then
        b_wl = ceiling( ( wl - l(1) ) * nwl / ( l(nwl) - l(1) ) )
      else
          b_wl = 1
        endif
        !add one to wavelength bin
        n_phot_wl(b_wl)= n_phot_wl(b_wl) + 1

        ! obtain all the optical properties for the photon

        u_a_old=u_a(:,b_wl)
        u_s_old=u_s(:,b_wl)
        g=g_skin(b_wl)


        dnaval=u_a_old(4)/0.0185
        hgg=g
        g2=hgg**2      ! Henyey-Greenstein parameter, hgg^2

        !******Find FIRST scattering location
       ! if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
        call tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
             xcell,ycell,zcell,tflag,iseed,delta,u_s_old,u_a_old,b_wl,e,seg_flag,dnaval)


        if (seg_flag.eq.1) EXIT
        !********Photon scatters in grid until it exits (tflag=1)
        tflag=0
        seg_flag=0
        do while(tflag.eq.0)
          albedo=0.d0
           do i =1,nlayer
              albedo=albedo+u_s_old(i)/(u_s_old(i)+u_a_old(i))*MASK(xcell,ycell,zcell,i)
           enddo
           call random_number(ran)
           if(ran.lt.albedo) then
              cost=nzp
              !************Scatter photon into new direction and update Stokes parameters

              call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
                   fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)
              scatter_count=scatter_count+1

           else

              EXIT
           endif

           !************Find next scattering location
            !if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
           call tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
                xcell,ycell,zcell,tflag,iseed,delta,u_s_old,u_a_old,b_wl,e,seg_flag,dnaval)
           !print*,xcell,ycell,zcell
           if (seg_flag.eq.1) EXIT
        end do



     end do ! end loop over nph photons

     print*, 'total number scatterings',scatter_count, 'per packet', real(scatter_count)/real(nphotons)

     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS


     print*,'photon count', pkt_count, nphotons, pkt_count/nphotons


     CALL PL_ESTIMATORS(pkt_count,nphotons,XMAX,YMAX,ZMAX,l,tot_irr,incident_spec_irr)


     print*,'Average number of scatterings = ',(scatter_count/nphotons)





endprogram mcpolar
!-------------------------------------------------------------------------------
