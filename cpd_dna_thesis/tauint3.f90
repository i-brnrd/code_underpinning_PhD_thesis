 module tauint3_mod
   implicit none
   save
 contains
!integration routine
subroutine tauint3(j,inside,b_wl)
use optical_properties_mod, only: nlayer
use grid_mod, only: xmax,ymax,zmax,MASK,nxg,nyg,nzg,grid_max, PL_SUM
use packet_mod, only: xp,yp,zp,nxp,nyp,nzp,cost,xcell,ycell,zcell, u_a,u_s
implicit none
!INTENTS
integer, intent(in) :: j, b_wl
logical, intent(inout) :: inside

!mcrt, path lengths & optical
real*8:: ran !random random_number
real*8:: tau,taurun, taucell !optical depths
real*8:: d, d_cell  !actual distance travelled
real*8 :: rk_tot,u_a_tot !total optical properties over all species
real*8 :: pos(3),dir(3),faces(3,3)
integer :: cell(3)

logical :: face_flag(3),edge_flag(3),tau_done
integer :: i

!SWITCH COORDINATE SYSTEMS
pos(1)=xp + xmax
pos(2)=yp + ymax
pos(3)=zp + zmax

cell(1)=xcell
cell(2)=ycell
cell(3)=zcell

dir(1)=nxp
dir(2)=nyp
dir(3)=nzp

!Intitalise edge flags and cumulative optical depths/distances
edge_flag=(/.false.,.false.,.false./)
face_flag=(/.false.,.false.,.false./)

taurun=0.
d=0.
!get random optical depth
call random_number(ran)
tau=-log(ran)

do
  !find the current cell details, next cell wall and distance to it
  call find_next_cell(cell,pos,dir,d_cell,faces,face_flag)
  !set optical properties for current cell.
  !call get_optical_properties()
  rk_tot=0.
  u_a_tot=0.
  do i=1,nlayer
    rk_tot=rk_tot+ (MASK(cell(1),cell(2),cell(3),i)*(u_s(i) + u_a(i)))
    u_a_tot=u_a_tot+ (MASK(cell(1),cell(2),cell(3),i)*u_a(i))
  enddo
  taucell=d_cell*rk_tot
 !decide if packet keeps progrressing or reaches optical depth in current cell
  if (taurun + taucell .lt.tau) then
    !packet passes through current voxel
    !Update path length counters and continue
    taurun=taurun + taucell
    d=d+d_cell
    tau_done=.false.
    PL_SUM(b_wl,cell(1),cell(2),cell(3))=PL_SUM(b_wl,cell(1),cell(2),cell(3))+d_cell
    call update_position(tau_done,d_cell,pos,dir,cell,faces,face_flag,edge_flag)
    !if we are over any edge(outside domain)
    !then edge_flag should be triggered
    if (any(edge_flag)) then

      call boundary_behaviour(edge_flag,tau_done,cell,grid_max,pos,dir(3),inside)

      if (.not.inside) then
        EXIT
      endif
    endif

  else !within the voxel, we reach the interaction location
    !Update path length counters, and position, and return to main program
    d_cell=((tau-taurun)/rk_tot)
    d=d+d_cell
    tau_done=.true.
    PL_SUM(b_wl,cell(1),cell(2),cell(3))=PL_SUM(b_wl,cell(1),cell(2),cell(3))+d_cell
    call update_position(tau_done,d_cell,pos,dir,cell,faces,face_flag,edge_flag) !with an exit flag
    EXIT
  endif
enddo

if (inside.and.(.not.tau_done)) print*, 'Left integration INSIDE&&.not.TAU_DONE: hard stop TauInt3'
if (inside.and.(.not.tau_done)) stop

!Switch coordinate systems back to those used in main modules
xp=pos(1) -xmax
yp=pos(2) -ymax
zp=pos(3) -zmax

xcell=cell(1)
ycell=cell(2)
zcell=cell(3)

nzp =dir(3)
cost=nzp
end subroutine
end module


subroutine boundary_behaviour(edge_flag,tau_done,cell,grid_max,pos,z_dir,inside)
use optical_properties_mod, only: nlayer
use packet_mod, only: ns
use grid_mod, only: nxg,nyg,nzg,delta,MASK
use n_interface_mod
implicit none
real*8, intent(inout) :: pos(3), z_dir
real*8,intent(in) :: grid_max(3)
logical, intent(inout) :: edge_flag(3), tau_done,inside
real*8 ::deltas(3)
integer :: cell(3)
integer :: i,j
logical :: reflect_flag

deltas=delta

if (edge_flag(1)) then !reflect boundaries (tbh call...)
  call repeat_bounds(cell(1),deltas(1),grid_max(1),nxg, pos(1))
  call get_cells(pos, cell,edge_flag)

elseif (edge_flag(2)) then
  call repeat_bounds(cell(2),deltas(2),grid_max(2),nxg, pos(2))
  call get_cells(pos, cell,edge_flag)

elseif(edge_flag(3)) then
  if (cell(3).lt.1) then !hit bottom edge, leave
    inside=.false.
  elseif(cell(3).gt.nzg)then
    call n_interface(ns,1.d0,z_dir,.false.,reflect_flag)
    if (reflect_flag) then
      call reflect(z_dir)
      pos(3)=grid_max(3)-deltas(3)
      call get_cells(pos, cell,edge_flag)
    else
      inside=.false.
    endif
  else
    print*,'Error in boundary_behaviour',edge_flag,cell
  endif
 else
   print*,'Error in boundary_behaviour',edge_flag,cell
endif

if (any(edge_flag).and.inside) print*, 'leaving boundary_conds with inside&&edgeflag signalling: STOP'
if (any(edge_flag).and.inside) stop


end subroutine

subroutine repeat_bounds(cell,delta,grid_max,ncell_max,position)
implicit none
integer, intent(in) :: cell, ncell_max
real*8, intent(in) :: delta, grid_max
real*8, intent(out) :: position

if (cell.lt.1) then
  position=grid_max-delta
elseif (cell.gt.ncell_max) then
  position=delta
else
  print*,'Error in REPEAT BOUNDS'
endif

end subroutine

!(tau_done,d_cell,pos,dir,cell,faces,face_flag,edge_flag)
subroutine update_position(tau_done,d_cell,pos,dir,cell,faces,face_flag,edge_flag)
use grid_mod, only: delta
implicit none

logical :: face_flag(3),edge_flag(3)
real*8,intent(in) :: faces(3,3)
integer, intent(inout):: cell(3)
real*8, intent(inout) :: pos(3)
real*8, intent(in) :: dir(3),d_cell

real*8 :: deltas(3)
logical, intent(in) :: tau_done
integer::i

integer:: low, up

deltas=delta

call get_cells(pos,cell,edge_flag)

if(tau_done)then
          do i=1,3
          pos(i)=pos(i)+dir(i)*d_cell
          enddo
else
  do i=1,3
   if(face_flag(i))then!use faces
            if(dir(i).gt.0.) then
              pos(i)= faces(3,i) + deltas(i)
            elseif(dir(i).lt.0.) then
               pos(i) = faces(2,i) - deltas(i)
            else
               print*,'Error in update_pos', dir,face_flag
            end if
     else
     pos(i) = pos(i)+ dir(i)*d_cell
   endif
  enddo
endif

call get_cells(pos,cell,edge_flag)

end subroutine update_position

subroutine get_cells(current_position, cell,edge_flag)
!givem current position, updates the current voxels
!also returns edge_flag if out of bounds
use grid_mod, only : xface, yface, zface, nxg,nyg,nzg
use search_bisec_mod
implicit none

real*8,  intent(IN):: current_position(3)
integer, intent(OUT) :: cell(3)
logical, intent(OUT) :: edge_flag(3)
integer :: low,up

!assume in bounds
edge_flag=(/.false.,.false.,.false./)

call search_bisec(current_position(1),xface,low,up)
cell(1)=low + lbound(xface,1) -1

call search_bisec(current_position(2),xface,low,up)
cell(2)=low + lbound(yface,1) -1

call search_bisec(current_position(3),zface,low,up)
cell(3)=low + lbound(zface,1) -1

!if out of bounds
if ((cell(1).lt.1).or.(cell(1).gt.nxg)) edge_flag(1)=.true.
if ((cell(2).lt.1).or.(cell(2).gt.nyg)) edge_flag(2)=.true.
if ((cell(3).lt.1).or.(cell(3).gt.nzg)) edge_flag(3)=.true.

end subroutine get_cells

subroutine find_next_cell(cell,pos,dir,dist_min,faces,face_flag)
! given cell number, position and direction,
! finds the next cell wall that will be hit, and distance to it
! returns dist_min, face_flag & faces
use grid_mod, only: xface,yface,zface, grid_max
implicit none
real*8, intent(in) ::pos(3),dir(3)
integer, intent(in) :: cell(3)
real*8, intent(out) :: dist_min, faces(3,3)
logical, intent(out) :: face_flag(3)

real*8 :: distance(3), dist_to_cell_wall
integer :: i

!inititalise face)flag and distance
face_flag=(/.false.,.false.,.false./)

distance=0.d0
!extract faces surrounding given cell
faces(:,1)=(/xface(cell(1)-1),xface(cell(1)),xface(cell(1)+1)/)
faces(:,2)=(/yface(cell(2)-1),yface(cell(2)), yface(cell(2)+1)/)
faces(:,3)=(/zface(cell(3)-1),zface(cell(3)),zface(cell(3)+1)/)

do i=1,3
  distance(i)=dist_to_cell_wall(dir(i),pos(i),faces(:,i),grid_max(i))
enddo

dist_min=minval(distance(:))
if(dist_min < 0.)print*,'dcell < 0.0 warning!! ',dist_min,distance,dir
if(dist_min < 0.)stop

face_flag(minloc(distance(:)))=.true.
if (.not.any(face_flag))  print*,'Next Face ERROR'
if (.not.any(face_flag)) STOP

end subroutine find_next_cell


function dist_to_cell_wall(direction,current_position,face_position,grid_max)
  !returns the distance along a given direction to a cell wall
  !given the positions of both cell walls
 implicit none
  real*8 :: dist_to_cell_wall
   real*8,intent(in) :: direction, current_position, face_position(3),grid_max

   if(direction.gt.0.) then
      dist_to_cell_wall=(face_position(3)-current_position)/direction
   elseif(direction.lt.0.) then
      dist_to_cell_wall=(face_position(2)-current_position)/direction
  elseif(direction.eq.0.) then
      dist_to_cell_wall=1.e2*grid_max
   endif
 return
end function
