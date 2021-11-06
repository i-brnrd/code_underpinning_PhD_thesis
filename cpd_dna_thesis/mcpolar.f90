!************************************************************
program mcpolar
  include 'modules_full_list.txt'
  implicit none
  !*****Parameter declarations ****************************************

  integer ::  nphotons   !number of photon packets sent into simulation
  real*8 :: diffuse_fraction !diffuse fractions of incident light
  character*70 :: fname_incident_irradiation !filename of incident irradiation

  !initialise vars for spectra
  real*8 incident_spec_irr(nwl)
  real*8 tot_irr
  real*8 :: cdf(nwl)

  integer:: pkt_count,scatter_count   !Counts packets that actually travel through medium (not reflected)
  real*8 ran
  real*8 albedo

  !properties of packets
  real*8 hgg,g2  !Henyey Greenstein phase function variables
  integer b_wl

 ! (loop) indices, counters and seeds
 integer j  !
 INTEGER I
 real*8 kappa,pl,pc,sc
 integer r_flag


 logical :: reflect_flag, diffuse_flag, inside  !booleans for whether packet reflects/ is diffuse etc

 print*,'**********************'
 print*,'Simulation of MC-UVRT through upper layers of skin'

 nphotons = 10000000


 print*, 'model size(cm)', xmax,ymax,zmax
 print*, 'sim in',nlayer, ' layers'
 print*, 'using photons', nphotons
 print*,'**********************'

call iarray !Initialise the arrays shared by the modules to 0
! Initialise counter :
pkt_count=0


!Skin parameters/ Optical Properties parameters/ Tissue parameters
ns=1.38

kappa=1.
pl=0.
pc=0.
sc=1.

!INITIALISE GRID
call gridset(kappa)

!depth resolution

!INITIALISE OPTICAL PROPERTIES
if (nwl.gt.1) then !else jumpt to Jacques verification

  call optical_properties_init()

  !check loaded optical properties
  !print out n, g, all scatt & all absorb
  open(8,file='optical_properties_check.txt',status='unknown')
  do i=1,nwl
    write(8,*) l(i),u_s_all(1:5,i), u_a_all(1:5,i),g_all(i)
    print*,l(i)
  enddo
  close(8)

  fname_incident_irradiation='./spectra/sunbed_spec.txt'
  !diffuse_fraction=0.d0!
!  diffuse_fraction=0.13
  diffuse_fraction=1.d0

  call load_spec2(fname_incident_irradiation,incident_spec_irr,1.d0)
  !incident_spec_irr=1.0

  call get_cdf(cdf,l,incident_spec_irr,nwl)
  tot_irr=sum(incident_spec_irr)
!  incident_spec_irr=1.0

  print*,'Incident Spec',fname_incident_irradiation
  print*,'Total irradiance:',tot_irr
  print*,'Diffuse Fraction',diffuse_fraction


else

  call jacques_verification(diffuse_fraction, tot_irr,l,incident_spec_irr)

endif

!INITIALISE RGN
call random_seed()
!---------------------------------------------------------------------------------------------------------------
!MCRT Program
!---------------------------------------------------------------------------------------------------------------
scatter_count=0
print*, '**************'

do j=1,nphotons
  if(mod(j,100000).eq.0)then
    print *, j,' scattered photons completed'
  end if
!*****Release photon from source *******************************
  call sourceph(cdf,diffuse_fraction,b_wl,diffuse_flag)

  call n_interface(1.d0,ns,cost,diffuse_flag,reflect_flag) !fresnel reflection

  if (reflect_flag) CYCLE

  pkt_count=pkt_count+1 !counts actual packets reaching the medium
  n_pkt_wl(b_wl)= n_pkt_wl(b_wl) + 1  !add one to wavelength bin

  !obtain all the optical properties for the packet
  !assumes wavelength won't change

        u_a=u_a_all(:,b_wl)
        u_s=u_s_all(:,b_wl)
        ns=n_all(b_wl)
        hgg=g_all(b_wl)

        g2=hgg**2   ! Henyey-Greenstein parameter, hgg^2
        inside=.true.
        call tauint3(j,inside,b_wl)

        !********Photon scatters in grid until it exits
        !POSITIONALLY DEPENDANT OPTICAL PROPERTY: albedo
        do while(inside)
          albedo=0.d0
           do i =1,nlayer
              albedo=albedo+u_s(i)/(u_s(i)+u_a(i))*MASK(xcell,ycell,zcell,i)
           enddo
           call random_number(ran)
           if(ran.lt.albedo) then

              if(cost.ne.nzp) print*, 'cost & nzp are not the same'
              cost=nzp    !ensuring this

              !************Scatter photon into new direction and update Stokes parameters
              call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
                   fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi)
              scatter_count=scatter_count+1

           else

              !absorb... bin z position in the future?
              EXIT
           endif
           call tauint3(j,inside,b_wl)
        end do

     end do ! end loop over nph photons

     print*, 'total number scatterings',scatter_count, 'per packet', real(scatter_count)/real(nphotons)
     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS
     print*,'photon count', pkt_count, nphotons, pkt_count/nphotons
     CALL PL_ESTIMATORS(pkt_count,incident_spec_irr)
     print*,'Average number of scatterings = ',(real(scatter_count)/real(nphotons))
     print*,  Finished '

endprogram mcpolar
