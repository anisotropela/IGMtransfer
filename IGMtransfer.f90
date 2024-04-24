!================================  Modules  ===================================

module PhysicalConstants
  implicit none
  real*8::           H_0                        !Hubble constant; cm/s/kpc
  real*8::           Omega_M                    !Matter density parameter
  real*8::           Omega_L                    !Dark energy density parameter
  real*8,parameter:: kpc      = 3.08567758131d21!Kiloparsec; cm
  real*8,parameter:: Mpc      = 3.08567758131d24!Kiloparsec; cm
  real*8,parameter:: c        = 2.99792458d10   !Speed of light; cm/s
  real*8,parameter:: m_H      = 1.673534d-24    !Hydrogen mass; g
  real*8,parameter:: m_e      = 9.1093897d-28   !Electron mass; g
  real*8,parameter:: e        = 4.803206d-10    !Electron charge; esu
  real*8,parameter:: h_Pl     = 6.6260755d-27   !Planck's constant; erg s
  real*8,parameter:: k_B      = 1.380658d-16    !Boltzmann's constant; erg/K
  real*8,parameter:: f_12     = 0.4162          !Oscillator strength
  real*8,parameter:: nu_0     = 2.46607d15      !Lya line-center frequency; Hz
  real*8,parameter:: lambda_0 = c / nu_0 * 1d8  !Line center wavelength; Angstrom
  real*8,parameter:: Dnu_L    = 9.936d7         !Lya natural line width; Hz
  real*8,parameter:: pi       = 3.14159265358979
end module PhysicalConstants

!------------------------------------------------------------------------------

module CellStructure
  implicit none
  integer, parameter::    posKind = kind(1.)    !Single precision should be fine
  type:: Cell
    real*4::              n_HI, Dnu_D, U_bulk(3)
    real(kind=posKind)::  C_pos(3)
    type(Cell), pointer:: Child(:,:,:)
# ifdef AMR
    integer(kind=1)::     Level                 !Max. 127 levels of refinement
    logical(kind=1)::     Refined
# endif
# ifdef dust
    real*4::              n_d
# endif
  end type Cell
end module CellStructure

!------------------------------------------------------------------------------

module DataArrays
  use CellStructure
  !Arrays for reading in data:
  real*4,    dimension(:), pointer:: n_HIString, TString,         &
                                     V_xString, V_yString, V_zString
# ifdef AMR
  integer*4, dimension(:), pointer:: LevelString
# endif
# ifdef dust
  real*4,    dimension(:), pointer:: ZString, n_HIIString
# endif
  !Mother grid with resolution ni,nj,nk:
  type(Cell), target::               BaseGrid
end module DataArrays

!------------------------------------------------------------------------------

module GlobalData
  integer::              ni,nj,nk               !Base grid resolution
  real*8::               dx0,dy0,dz0,dxTiny,dxWeeny!Base grid cell size / negligible distance
  real*8,allocatable,&   
         dimension(:)::  dx,dxH,dy,dyH,dz,dzH   !Cell size / Half cell size
  integer::              N_cells,i_cell         !Total # of cells; counter
  integer::              n_ph,n_lostot,n_los    !# of photons, # of sightlines
  integer::              l_lostot               !# of photons, # of sightlines
  integer::              L_max                  !Maximum level of refinement
  real*8::               D_box,R_box            !Phys. size/half size of Box
  real*8::               z                      !Redshift
  real*8::               d_tot,r_eff
  real*8::               offZet(3)
  real*4,allocatable::   Gollum(:,:)
  integer::              n_gal,n_write          !Total # of galaxies,LOSs between output
  real*4,allocatable::   maxx_gal(:),      &
                         maxy_gal(:),      &
                         maxz_gal(:),      &
                         maxr_vir(:),      &
                         maxv_x(:),        &
                         maxv_y(:),        &
                         maxv_z(:)
  real*4,allocatable::   x_gal(:),      &
                         y_gal(:),      &
                         z_gal(:),      &
                         r_vir(:),      &
                         v_x(:),        &
                         v_y(:),        &
                         v_z(:)
  real*8::               f_rvir                 !Fraction of r_vir to start los
  character(len=200)::   indir,subdir,CellData,GalData !Directories and names of input data
  character(len=200)::   outdir,Ioutfile        !Name of outfile
  character(len=3)::     DustType               !'LMC' or 'SMC'
  real*8::               f_ion                  !Dust destruction factor; determines to how large a degree ionized hydrogen contributes to dust
end module GlobalData

!------------------------------------------------------------------------------

module IJK
  integer:: i,j,k
end module IJK

!------------------------------------------------------------------------------

module ObservationalParameters
  implicit none
  real*8::             nu,lambda,BW(2)
  real*8,allocatable:: l_inj(:)
  integer::            SpecRes
end module ObservationalParameters

!------------------------------------------------------------------------------

module CellParameters
  real*8:: C_pos(3)                             !Position of cell center
  real*8:: n_HI                                 !Density of neutral hydrogen
  real*8:: n_d                                  !"Dust density", i.e. rescaled hydrogen density
  real*8:: Dnu_D                                !Doppler width of line
  real*8:: U_bulk(3)                            !Gas bulk velocity in Doppler widths
  real*8:: xp,xm,yp,ym,zp,zm,Faces(6)
end module CellParameters

!------------------------------------------------------------------------------

module INTERFACE_LocateHostCell
  interface
    subroutine LocateHostCell(X_pos,CurrentCell,i,j,k)
      use CellStructure
      implicit none
      real*8, intent(IN)::   X_pos(3)
      type(Cell), pointer::  CurrentCell        !intent(OUT)
      integer, intent(OUT):: i,j,k
    end subroutine LocateHostCell
  end interface
end module INTERFACE_LocateHostCell

!------------------------------------------------------------------------------

module INTERFACE_DigDeeper
  interface
    recursive subroutine DigDeeper(X_pos,CurrentCell)
      use CellStructure
      implicit none
      real*8, intent(IN)::   X_pos(3)
      type(Cell), pointer::  CurrentCell        !intent(INOUT)
    end subroutine DigDeeper
  end interface
end module INTERFACE_DigDeeper

!------------------------------------------------------------------------------

module INTERFACE_EveryoneShouldHaveACenter
  interface
    recursive subroutine EveryoneShouldHaveACenter(CurrentCell,i,j,k,ParentPos)
      use CellStructure
      use GlobalData
      implicit none
      type(Cell), target ::           CurrentCell
      integer,intent(in)::            i,j,k
      real(kind=posKind),intent(in):: ParentPos(3)
    end subroutine EveryoneShouldHaveACenter
  end interface
end module INTERFACE_EveryoneShouldHaveACenter

!------------------------------------------------------------------------------

module INTERFACE_BuildNestedUniverse
  interface
    recursive subroutine BuildNestedUniverse(CurrentCell,CurrentLevel)
      use CellStructure
      use GlobalData
      use DataArrays
      implicit none
      type(Cell), target::     CurrentCell
      integer, intent(IN)::    CurrentLevel
    end subroutine BuildNestedUniverse
  end interface
end module INTERFACE_BuildNestedUniverse

!------------------------------------------------------------------------------

module INTERFACE_GetCellData
  interface
    recursive subroutine GetCellData(HostCell,oldCell,DirVec)
      use CellStructure
      implicit none
      type(Cell), pointer::  HostCell,OldCell   !intent(in)
      integer, intent(IN)::  DirVec(3)
    end subroutine GetCellData
  end interface
end module INTERFACE_GetCellData

!------------------------------------------------------------------------------

module INTERFACE_EnterNeighbor
  interface
    recursive subroutine EnterNeighbor(CurrentCell,X_cut,DirVec,i,j,k)
      use CellStructure
      use GlobalData
      use INTERFACE_LocateHostCell
      implicit none
      type(Cell), pointer::   CurrentCell     !intent(INOUT)
      real*8,intent(IN)::     X_cut(3)
      integer, intent(IN)::   DirVec(3)
      integer,intent(INOUT):: i,j,k
    end subroutine EnterNeighbor
  end interface
end module INTERFACE_EnterNeighbor

!------------------------------------------------------------------------------

module INTERFACE_LorentzTransform
  interface
    subroutine LorentzTransform(SpecRes,x,nhat_i,oldCell,newCell)
      use CellStructure
      implicit none
      integer::               SpecRes
      real*8, intent(inout):: x(SpecRes)
      real*8, intent(in)::    nhat_i(3)
      type(Cell), pointer::   oldCell, newCell
    end subroutine LorentzTransform
  end interface
end module INTERFACE_LorentzTransform

!------------------------------------------------------------------------------

module INTERFACE_PrintCell
  interface
    subroutine PrintCell(CurrentCell)
      use CellStructure
      implicit none
      type(Cell), pointer::  CurrentCell              !intent(in)
    end subroutine PrintCell
  end interface
end module  

!------------------------------------------------------------------------------

!==============================================================================

program IGMtransfer

!-----------------  Declarations  ---------------------------------------------

use PhysicalConstants
use DataArrays
use GlobalData
use IJK
use ObservationalParameters
use CellParameters
use INTERFACE_LocateHostCell
use INTERFACE_DigDeeper
use INTERFACE_GetCellData
use INTERFACE_EnterNeighbor
use INTERFACE_LorentzTransform
use INTERFACE_PrintCell

implicit none
integer::             g,n,l
real*8::              r,dist
real*8,allocatable::  x(:),sigma_x(:),sigma_d(:),tau(:)
real*8::              norm,H
real*8,dimension(3):: X_seed,X_i,X_f,nhat_i,V_sys,U_sys
logical::             reflect
integer::             DirVec(3)
type(Cell), pointer:: HostCell => null(), oldCell => null()

!-----------------  Initialization  -------------------------------------------

call ReadInput
call ReadData
call ReadGalaxies
call ConvertParameters

r_eff    = r_eff * R_box                        !r_eff is converted from fraction to absolute distance
n_lostot = n_los * n_gal                        !TOTAL # of los for ALL galaxies

allocate(x(SpecRes))
allocate(sigma_x(SpecRes))
allocate(sigma_d(SpecRes)) ; sigma_d = 0d0
allocate(tau(SpecRes))
allocate(Gollum(n_lostot,SpecRes))
print*, "Final size of '", trim(Ioutfile)//"' is", kind(Gollum)*size(Gollum)/1.048576e6, 'Mb'

call CalcDX
call ConstructUniverse

d_tot = (lambda_0-BW(1))/BW(1) * c/H(z)         !Total distance to be integrated
print*, 'd_tot/Mpc =',real(d_tot/Mpc)
print*
print*, 'Go!'
print*

!-----------------  Fire! -----------------------------------------------------

do g = 1,n_gal                                  !Loop over galaxies
  if (mod(g,1) .eq. 0) print*, g, 'of', n_gal, 'galaxies'

  V_sys = (/v_x(g),v_y(g),v_z(g)/)              !Systemic velocity of galaxy

  do l = 1,n_los                                !Loop over sightlines for this galaxy
    if (mod(l,10) .eq. 0) print*,'  ', l, 'of', n_los, 'sightlines for this galaxy.'
    call InitPos(g,X_i)                         !Initial position
    call InitDir(g,X_i,nhat_i)                  !Initial direction
    call LocateHostCell(X_i,HostCell,i,j,k)     !Make pointer to host cell
    DirVec = (/7,9,13/)                         !DirVec will be a unit vector along the Cartesian axes; here it is set to anything else than that
    call GetCellData(HostCell,oldCell,DirVec)   !Physical parameters of cell ijk
    U_sys  = V_sys/c * nu_0/Dnu_D               !Systemic velocity of gal, expressed in terms of T in host cell
    call InitFreq(Dnu_D,U_bulk,U_sys,nhat_i,x)  !Inital frequency
    call sigma(SpecRes,Dnu_D,x,sigma_x)         !Lya cross section for each frequency bin
#   ifdef dust
      call XsecD(SpecRes,Dnu_D,x,sigma_d)       !Dust cross section for each frequency bin
#   endif
    tau     = 0.                                !Total tau encountered in each spectral bin for this sightline
    dist    = 0.                                !Total distance traversed for this sightline
    reflect = .false.
 
    do                                          !Bounce until exit
      X_f  = X_i + 2*D_box*nhat_i               !Somewhere outside box
      call CheckIntersect(X_i,X_f,nhat_i,DirVec,reflect) !Find exit point of cell
      r    = norm(X_f - X_i)                    !Dist. traveled in current step
      dist = dist + r                           !Cumulated r in this sightline
      tau  = tau  + r * (n_HI*sigma_x + n_d*sigma_d)!Cumulated tau in this bin
      if (dist .ge. d_tot) exit                 !Exit when spectrum is sufficiently redshifted

      if (.not. reflect) then                   !If sightline hasn't yet reached surface
        oldCell => HostCell
        call EnterNeighbor(HostCell,X_f,DirVec,i,j,k)           !Locate neighboring cell
        call GetCellData(HostCell,oldCell,DirVec)               !Get cell data of new cell
        call LorentzTransform(SpecRes,x,nhat_i,oldCell,HostCell)!Express frequency in reference frame of new cell
        call sigma(SpecRes,Dnu_D,x,sigma_x)                     !Re-calculate cross section for Lorentz transformed spectrum
#       ifdef dust
          call XsecD(SpecRes,Dnu_D,x,sigma_d)                   !   -"-
#       endif
      else
        call Bounce(X_f,nhat_i)                 !If sightline reaches surface
        reflect = .false.
      endif

      X_i = X_f                                 !Update position
    enddo

    l_lostot           = (g-1)*n_los + l        !"l", not "1"! Total # of sightlines so far
    Gollum(l_lostot,:) = real(tau(:))           !Store LAF for this sightline
  enddo

  if (mod(l_lostot,n_write).eq.0) call WriteData(l_lostot) !Write data obtained so far
enddo

call WriteData(l_lostot)

end program IGMtransfer

!==============================  Subroutines  =================================

subroutine Bounce(X_pos,nhat)

!When a line of sight reaches the edge of the computational volume,
!it "bounces", i.e. it continues in a random inward angle.
!This subroutine selects the new direction.
!If r_eff/R_box <= 1, the computational volume is an enclosed sphere.
!If r_eff/R_box >= sqrt(3), the computational volume is the full box
!Intermediate values are allowed, but is not really meaningful (and will use
!the "enclosed sphere" algorithm.
!To assure a more or less uniform sampling of the sphere, very wide angles
!(given by mu_min) are not allowed.

use PhysicalConstants
use GlobalData

implicit none
real*8,intent(in)::  X_pos(3)
real*8,intent(out):: nhat(3)
real*8,parameter::   mu_min = .3                !Max deflection from center in case of the computational volume being an enclosed sphere
real*8::             R,R1,R2,R3(3),n1,n2,mu,rsq,frevert,ToCent(3),norm,dmu
integer::            ii

dmu = 1d0 - mu_min                              !Range of allowed angles

if ((r_eff - R_box)/kpc .lt. 1e-6) then         !"Enclosed sphere"
  call random_number(R)
  mu = R*dmu + mu_min

  do
    call random_number(R1)
    call random_number(R2)
    n1 = 2*R1 - 1
    n2 = 2*R2 - 1

    rsq = n1**2 + n2**2
    if (rsq .le. 1d0) exit
  enddo

  frevert = sqrt((1 - mu**2) / rsq)
  n1      = frevert * n1
  n2      = frevert * n2

  offZet  = (/n1,n2,mu/)                        !nhat if X_pos were at (R_box,R_box,0)
  ToCent  = (R_box - X_pos) / norm(R_box-X_pos) !Direction toward center of box

  call RotateVector(ToCent,offZet,nhat)
else                                            !"Full box"
  do
    call random_number(R3)
    R3  = 2*R3 - 1
    rsq = sum(R3**2)
    if (rsq .le. 1d0) exit
  enddo
  nhat = R3 / sqrt(rsq)                         !Normalize to unit vector
  do ii = 1,3                                   !Make sure nhat points inward
    if (X_pos(ii)/kpc              .lt. 1d-13) nhat(ii) =  abs(nhat(ii))
    if (abs((X_pos(ii)-D_box)/kpc) .lt. 1d-13) nhat(ii) = -abs(nhat(ii))
  enddo
endif

end subroutine Bounce

!------------------------------------------------------------------------------

subroutine ConvertParameters

!Convert the following parameters:
!
! 1. Temperature is converted to "width of the associated Doppler broadened Lya
!    line in frequency", i.e. TString now contains Dnu_D.
!
! 2. Velocity is converted to "velocity in terms of 'width of the associated
!    Doppler broadened Lya line in velocity'", i.e. V_[xyz]String now contains
!    u_[xyz] = v_[xyz] / v_th.
!
! 3. Metallicity is converted to "dust density", i.e. a rescaled hydrogen
!    density such that multiplying this number by the dust cross section per
!    hydrogen nucleus (in the subroutine "XsecD") gives the optical depth of
!    dust per unit distance.

use GlobalData
use DataArrays

implicit none
integer:: i
real*8::  Z_0
real*4::  v2u,Doppler

do i = 1,N_cells
  TString(i)   = Doppler(TString(i))            !Convert T/K to Dnu_D/Hz
  V_xString(i) = v2u(V_xString(i),TString(i))   !Convert v/kms to u
  V_yString(i) = v2u(V_yString(i),TString(i))
  V_zString(i) = v2u(V_zString(i),TString(i))
enddo

# ifdef dust
  if (DustType .eq. 'SMC') then
    Z_0 = 10d0**(-.6)
  elseif (DustType .eq. 'LMC') then
    Z_0 = 10d0**(-.3)
  endif
  ZString = (n_HIString + f_ion*n_HIIString) * ZString/Z_0!Convert to "dust density"
  deallocate(n_HIIString)
# endif

end subroutine ConvertParameters

!------------------------------------------------------------------------------

subroutine ConstructUniverse

!Build the adaptively refined grid, on the basis of the levels of the cell
!(given in LevelString).

use GlobalData
use DataArrays
use CellStructure
use INTERFACE_BuildNestedUniverse
use INTERFACE_EveryoneShouldHaveACenter

implicit none
integer::            i,j,k
real(kind=posKind):: DummyPos(3)

print*, 'Constructing Universe:'

print*, ' - Initializing base grid...'
call InitializeBaseGrid

print*, ' - Building fully threaded structure...'
i_cell = 0
do i = 1,ni
  do j = 1,nj
    do k = 1,nk
      call BuildNestedUniverse(BaseGrid%Child(i,j,k),0)
    enddo
  enddo
enddo

print*, ' - Assigning cell centers...'
DummyPos = (/6,6,6/)
do i = 1,ni
  do j = 1,nj
    do k = 1,nk
      call EveryoneShouldHaveACenter(BaseGrid%Child(i,j,k),i,j,k,DummyPos)
    enddo
  enddo
enddo

# ifdef AMR
  deallocate(LevelString)
# endif
deallocate(n_HIString)
deallocate(TString)
deallocate(V_xString)
deallocate(V_yString)
deallocate(V_zString)
# ifdef dust
  deallocate(ZString)
# endif

end subroutine ConstructUniverse

!------------------------------------------------------------------------------

subroutine CalcDX

!Calculate the side length of cells of all levels.

use CellStructure
use GlobalData

implicit none
integer:: L

allocate(dx(0:L_max))
allocate(dy(0:L_max))
allocate(dz(0:L_max))
allocate(dxH(0:L_max))
allocate(dyH(0:L_max))
allocate(dzH(0:L_max))

dx0 = D_box / ni
dy0 = D_box / nj
dz0 = D_box / nk

do L = 0,L_max
  dx(L) = dx0 / 2d0**L
  dy(L) = dy0 / 2d0**L
  dz(L) = dz0 / 2d0**L
enddo

dxH = dx / 2d0
dyH = dy / 2d0
dzH = dz / 2d0

dxTiny  = 0.1d0 * dx(L_max)
dxWeeny = 1d-8  * dx(L_max)

end subroutine CalcDX

!------------------------------------------------------------------------------

subroutine CheckIntersect(X_i,X_f,nhat_i,DirVec,reflect)

!Determines the point of intersection between the path of the photon
!(sightline)and the face of the cell. If the photon escapes the cell,
!it is put at the face and a flag is set for updating parameters.
!If it stays in the the cell, a flag is set for scattering the photon.

use CellParameters
use IJK
use PhysicalConstants
use GlobalData

implicit none
real*8,intent(in)::    X_i(3),nhat_i(3)
real*8,intent(inout):: X_f(3)
integer, intent(out):: DirVec(3)
logical,intent(out)::  reflect
real*8::               X_cut(3),norm,d_f,d_rem

reflect = .true.

!Check if photon crosses yz-plane in positive direction
if (X_f(1) .ge. xp) then
  X_cut(1) = xp                                           !Intersection
  X_cut(2) = X_i(2) + (xp - X_i(1)) * nhat_i(2)/nhat_i(1) !with
  X_cut(3) = X_i(3) + (xp - X_i(1)) * nhat_i(3)/nhat_i(1) !plane

  if (X_cut(2) .lt. yp .and. &                  !If so, check if it happened
      X_cut(2) .ge. ym .and. &                  !through face of cell
      X_cut(3) .lt. zp .and. &
      X_cut(3) .ge. zm)        then
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/1,0,0/)                         !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif
endif

!Check if photon crosses yz-plane in negative direction
if (X_f(1) .lt. xm) then                        
  X_cut(1) = xm                                           !Intersection 
  X_cut(2) = X_i(2) + (xm - X_i(1)) * nhat_i(2)/nhat_i(1) !with
  X_cut(3) = X_i(3) + (xm - X_i(1)) * nhat_i(3)/nhat_i(1) !plane

  if (X_cut(2) .lt. yp .and. &                  !If so, check if it happened
      X_cut(2) .ge. ym .and. &                  !through face of cell
      X_cut(3) .lt. zp .and. &
      X_cut(3) .ge. zm)         then
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/-1,0,0/)                        !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif 
endif

!Check if photon crosses xz-plane in positive direction
if (X_f(2) .ge. yp) then                    
  X_cut(1) = X_i(1) + (yp - X_i(2)) * nhat_i(1)/nhat_i(2) !Intersection
  X_cut(2) = yp                                           !with
  X_cut(3) = X_i(3) + (yp - X_i(2)) * nhat_i(3)/nhat_i(2) !plane

  if (X_cut(1) .lt. xp .and. &                  !If so, check if it happened
      X_cut(1) .ge. xm .and. &                  !through face of cell
      X_cut(3) .lt. zp .and. & 
      X_cut(3) .ge. zm)        then 
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/0,1,0/)                         !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif 
endif

!Check if photon crosses xz-plane in negative direction   
if (X_f(2) .lt. ym) then                    
  X_cut(1) = X_i(1) + (ym - X_i(2)) * nhat_i(1)/nhat_i(2) !Intersection
  X_cut(2) = ym                                           !with
  X_cut(3) = X_i(3) + (ym - X_i(2)) * nhat_i(3)/nhat_i(2) !plane

  if (X_cut(1) .lt. xp .and. &                  !If so, check if it happened
      X_cut(1) .ge. xm .and. &                  !through face of cell
      X_cut(3) .lt. zp .and. &  
      X_cut(3) .ge. zm)        then
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/0,-1,0/)                        !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif
endif

!Check if photon crosses xy-plane in positive direction
if (X_f(3) .ge. zp) then                    
  X_cut(1) = X_i(1) + (zp - X_i(3)) * nhat_i(1)/nhat_i(3) !Intersection
  X_cut(2) = X_i(2) + (zp - X_i(3)) * nhat_i(2)/nhat_i(3) !with
  X_cut(3) = zp                                           !plane

  if (X_cut(1) .lt. xp .and. &                  !If so, check if it happened
      X_cut(1) .ge. xm .and. &                  !through face of cell
      X_cut(2) .lt. yp .and. &      
      X_cut(2) .ge. ym)        then
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/0,0,1/)                         !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif
endif

!Check if photon crosses xy-plane in negative direction
if (X_f(3) .lt. zm) then                   
  X_cut(1) = X_i(1) + (zm - X_i(3)) * nhat_i(1)/nhat_i(3) !Intersection
  X_cut(2) = X_i(2) + (zm - X_i(3)) * nhat_i(2)/nhat_i(3) !with
  X_cut(3) = zm                                           !plane

  if (X_cut(1) .lt. xp .and. &                  !If so, check if it happened
      X_cut(1) .ge. xm .and. &                  !through face of cell
      X_cut(2) .lt. yp .and. &       
      X_cut(2) .ge. ym)        then
    X_f     = X_cut                             !Place photon at face
    DirVec  = (/0,0,-1/)                        !Continue in neighboring cell
    reflect = .false.                           !...if edge is not reached
  endif 
endif

d_f = sqrt(sum((X_f - R_box)**2))               !Distance from center of box

if (d_f .ge. r_eff      .or. &
    any(X_f .le. 0.)    .or. &                  !Could probably be more efficient, e.g. by
    any(X_f .ge. D_box)) then                   !comparing mother cell indices which are integers, instead of positions which are reals
  reflect = .true.
  if (r_eff .le. R_box) call FindIntersect(X_i,nhat_i,d_rem,X_f) !If not, then X_f is simply that alreay found
endif

end subroutine CheckIntersect

!------------------------------------------------------------------------------

subroutine EnterNeighbor(CurrentCell,X_cut,DirVec,i,j,k)

!Locate the neighboring cell.

use CellStructure
use GlobalData
use INTERFACE_LocateHostCell
use DataArrays

implicit none
type(Cell), pointer::   CurrentCell             !intent(INOUT)
real*8,intent(IN)::     X_cut(3)
integer, intent(IN)::   DirVec(3)
integer,intent(INOUT):: i,j,k
integer::               NeighbInd(3)

! NeighbInd = CurrentCell%indexico + DirVec

! if (ALL(NeighbInd.ge.ones .and. NeighbInd.le.twos)) then
!   CurrentCell => CurrentCell%Parent%Child(NeighbInd(1), &
!                                           NeighbInd(2), &
!                                           NeighbInd(3))
! else
    call LocateHostCell(X_cut + DirVec*dxTiny,CurrentCell,i,j,k)
! endif

end subroutine EnterNeighbor

!------------------------------------------------------------------------------

subroutine FindIntersect(P1,nhat,u,P2)

!Given a sphere of radius r_eff centered at R_box, and a line specified by a
!point P1 and unit direction vector nhat (that points to interior of sphere),
!calculate the other point P2 of intersection as well as distance u to P2.
!For numerical reasons, positions are first divided and later multiplied by kpc.

use GlobalData

implicit none
real*8,intent(in)::  P1(3),nhat(3)
real*8,intent(out):: u,P2(3)
real*8,parameter::   kpc = 3.08567758131d21
real*8::             x1,y1,z1,x2,y2,z2,x3,y3,z3
real*8::             r,a,b,c,D

x1 = P1(1) / kpc                                !1st point on line
y1 = P1(2) / kpc
z1 = P1(3) / kpc

x2 = x1 + nhat(1)                               !2nd point on line
y2 = y1 + nhat(2)
z2 = z1 + nhat(3)

x3 = R_box / kpc                                !Center of sphere
y3 = R_box / kpc
z3 = R_box / kpc

r  = r_eff / kpc                                !Radius of sphere

a  = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2
b  = 2 * ((x2 - x1)*(x1 - x3) + (y2 - y1)*(y1 - y3) + (z2 - z1)*(z1 - z3))
c  = x3**2 + y3**2 + z3**2 &
   + x1**2 + y1**2 + z1**2 &
   - 2*(x3*x1 + y3*y1 + z3*z1) - r**2
D  = b**2 - 4*a*c

u  = (-b + sqrt(D)) / (2*a) * kpc
P2 = P1 + u*nhat

end subroutine FindIntersect

!------------------------------------------------------------------------------

subroutine GetCellData(HostCell,oldCell,DirVec)

!Store cell data and derived quantities in temporary variables.

use CellParameters
use DataArrays
use GlobalData
use PhysicalConstants

implicit none
type(cell), pointer:: HostCell,oldCell
integer, intent(IN):: DirVec(3)
integer::             L

n_HI   = dble(HostCell%n_HI)
Dnu_D  = dble(HostCell%Dnu_D)
U_bulk = dble(HostCell%U_bulk)

# ifdef AMR
  L    = HostCell%Level
# else
  L    = 0
# endif

#ifdef dust
  n_d  = dble(HostCell%n_d)
#else
  n_d  = 0.
#endif

if     (DirVec(1) .eq. 1)  then
  xm = xp
  xp = HostCell%C_pos(1) + dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(1) .eq. -1) then
  xp = xm
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(2) .eq. 1)  then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  ym = yp
  yp = HostCell%C_pos(2) + dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(2) .eq. -1) then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = ym
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(3) .eq. 1)  then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zm = zp
  zp = HostCell%C_pos(3) + dzH(L)
elseif (DirVec(3) .eq. -1) then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = zm
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(1) .eq. 7) then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
else
  stop "WTF?!"
endif

Faces  = (/xp,xm,yp,ym,zp,zm/)

end subroutine GetCellData

!------------------------------------------------------------------------------

subroutine InitializeBaseGrid

use CellStructure
use GlobalData
use DataArrays

implicit none
integer:: i,j,k

allocate(BaseGrid%Child(ni,nj,nk))

do k = 1, nk
  do j = 1, nj
    do i = 1, ni
#       ifdef AMR
        BaseGrid%Child(i,j,k)%Refined    = .false.
        BaseGrid%Child(i,j,k)%Level      = 0
#       endif
        BaseGrid%Child(i,j,k)%C_pos      = (/(i - .5d0), &
                                             (j - .5d0), &
                                             (k - .5d0)/) * dx(0)
        BaseGrid%Child(i,j,k)%n_HI       = 0e0
        BaseGrid%Child(i,j,k)%Dnu_D      = 0e0
        BaseGrid%Child(i,j,k)%U_bulk     = (/0e0,0e0,0e0/)
nullify(BaseGrid%Child(i,j,k)%Child)
    enddo
  enddo
enddo

end subroutine InitializeBaseGrid

!------------------------------------------------------------------------------

subroutine InitDir(g,X_pos,nhat)

!Given a position X_pos and a galaxy g, make a unit vector in the direction
!away from the center of the galaxy.

use PhysicalConstants
use GlobalData

implicit none
integer,intent(in):: g
real*8,intent(in)::  X_pos(3)
real*8,intent(out):: nhat(3)
real*8::             Cg2X(3),norm

Cg2X = X_pos - (/x_gal(g),y_gal(g),z_gal(g)/)
nhat = Cg2X / norm(Cg2X)

end subroutine InitDir

!------------------------------------------------------------------------------

subroutine InitFreq(Dnu_D,U_bulk,U_sys,nhat,x)

!Create an initial spectrum.

use PhysicalConstants
use ObservationalParameters

implicit none
real*8,intent(in)::  Dnu_D,U_bulk(3),U_sys(3),nhat(3)
real*8,intent(out):: x(SpecRes)
integer::            n
real*8::             dlambda,x_sys,U_rel(3),x_shift

dlambda = (BW(2) - BW(1)) / SpecRes

if (.not. allocated(l_inj)) allocate(l_inj(SpecRes))

do n = 1,SpecRes
  lambda   = BW(1) + (n-.5) * dlambda
  l_inj(n) = lambda
  nu       = c / (lambda * 1d-8)
  x_sys    = (nu - nu_0) / Dnu_D                !x in systemic frame
  U_rel    = U_bulk - U_sys                     !U of gas, relative to galaxy
  x_shift  = dot_product(nhat,U_rel)            !dx when leaving galaxy
  x(n)     = x_sys + x_shift
enddo

end subroutine InitFreq

!------------------------------------------------------------------------------

subroutine InitPos(g,X_seed)

!Find a random position X_seed on the surface of the sphere centered at a
!galaxy g, in the distance f_rvir virial radii from the center of the galaxy.

use PhysicalConstants
use GlobalData
use CellStructure
use DataArrays

implicit none
integer,intent(in):: g
real*8,intent(out):: X_seed(3)
real*8::             R(3),nhat(3),rsq,d,r_0

do
  do
    call random_number(R)
    nhat = 2*R - 1
    rsq  = sum(nhat**2)
    if (rsq .le. 1.) exit
  enddo

  nhat    = nhat / sqrt(rsq)
  r_0     = f_rvir * r_vir(g)
  X_seed  = (/x_gal(g),y_gal(g),z_gal(g)/) + r_0 * nhat

  d       = sqrt(sum((X_seed - R_box)**2))
  if (d .lt. r_eff .and. all(X_seed.gt.0) .and. all(X_seed.lt.D_box)) exit
enddo

end subroutine InitPos

!------------------------------------------------------------------------------

subroutine Int2Char(i,a,LZ)

!Convert integer to character.

implicit none
integer, intent(in)::           i
character(len=*), intent(out):: a
logical, intent(in)::           LZ              !Leading zeros?
integer, parameter::            zero = iachar('0')
integer::                       n, N_digits, length
character(len=len(a))::         b

a        = ''
b        = ''
length   = len(a)
N_digits = 0

do n = 1,length
  N_digits = N_digits + 1
  a = achar(zero + mod(i,10**n)/10**(n-1)) // trim(a)
  if (i .lt. 10**n) exit
enddo

if (LZ) then
  b = repeat('0',len(a) - N_digits)
  a = trim(b) // trim(a)
endif

end subroutine Int2Char

!------------------------------------------------------------------------------

subroutine LocateHostCell(X_pos,CurrentCell,i,j,k)

!Let 'CurrentCell' point to the cell hosting the point X_pos.

use CellStructure
use GlobalData
use DataArrays
use PhysicalConstants
use INTERFACE_DigDeeper

implicit none
real*8, intent(IN)::   X_pos(3)
type(Cell), pointer::  CurrentCell              !intent(out)
integer, intent(OUT):: i,j,k
real*8::               dist_c

i = ceiling(X_pos(1) / dx0)                     !
j = ceiling(X_pos(2) / dy0)                     !Indices of mother cell
k = ceiling(X_pos(3) / dz0)                     !

CurrentCell => BaseGrid%Child(i,j,k)
# ifdef AMR
  if (CurrentCell%Refined) call DigDeeper(X_pos,CurrentCell)
# endif

end subroutine LocateHostCell

!------------------------------------------------------------------------------

subroutine LorentzTransform(SpecRes,x,nhat_i,oldCell,newCell)

!Lorentz transform x according to temperature and velocity of new cell.

use CellStructure

implicit none
integer::               SpecRes
real*8, intent(inout):: x(SpecRes)
real*8, intent(in)::    nhat_i(3)
type(Cell), pointer::   oldCell, newCell
real*8::                oldD,oldU(3),newD,newU(3),u1n,u2n 

oldD = dble(oldCell%Dnu_D)
oldU = dble(oldCell%U_bulk)
newD = dble(newCell%Dnu_D)
newU = dble(newCell%U_bulk)

u1n  = dot_product(oldU,nhat_i)
u2n  = dot_product(newU,nhat_i)

x    = (x + u1n) * oldD/newD - u2n

end subroutine LorentzTransform

!------------------------------------------------------------------------------

subroutine PrintCell(CurrentCell)

!Print cell data (for debugging).

use CellStructure
use GlobalData
use PhysicalConstants

implicit none
type(Cell), pointer::  CurrentCell              !intent(in)

# ifdef AMR
  print*, 'Level:    ', CurrentCell%Level
# endif
  print*, 'n_HI:     ', CurrentCell%n_HI
# ifdef dust
  print*, 'n_d:      ', CurrentCell%n_d
# endif
  print*, 'Dnu_D:    ', CurrentCell%Dnu_D
  print*, 'T:        ', m_H/2/k_B * (CurrentCell%Dnu_D * c/nu_0)**2
  print*, 'U_bulk:   ', CurrentCell%U_bulk
  print*, 'V_bulk:   ', CurrentCell%U_bulk * CurrentCell%Dnu_D * c / nu_0 / 1e5
  print*, 'C_pos:    ', CurrentCell%C_pos / kpc

end subroutine PrintCell

!------------------------------------------------------------------------------

subroutine RandomDir(nhat_i)

!Randomly select a unit vector with isotropical distribution.

implicit none
real*8::             rsq
real*8::             R(3),n_i(3)
real*8,intent(out):: nhat_i(3)

do
  call random_number(R)
  n_i = 2*R - 1.
  rsq = sum(n_i**2)
  if (rsq .le. 1.) exit                         !Use r only if |r|<1 in order
enddo                                           !not to favor corners of box

nhat_i = n_i / sqrt(rsq)                        !Normalize to unit vector

end subroutine RandomDir

!------------------------------------------------------------------------------

subroutine ReadData

!Load cell data.
!
!Note the following difference from v1.0 in v1.1 and subsequent versions:
!
! 1. "TString" should read the real temperature; it will be converted to Dnu_D
!    later, i.e. the corresponding Doppler width.
!
! 2. "V_[xyz]String" should read the "real" velocities of the gas elements in
!    km/s, i.e. the sum of their peculiar motions and the Hubble flow, and with
!    respect to the center of the box (corresponding to the systemic velocities
!    of the galaxies in the GalData-file).
!    The velocities are later changed to "u", i.e. the velocities in term of the
!    thermal (Doppler) width.
!
! 3. If dust is included in the RT, then instead of supplying a "dust density",
!    a metallicity (in terms of Solar) and an ionized hydrogen density is
!    supplied.

use PhysicalConstants
use GlobalData
use DataArrays

implicit none
integer::           n,i,n_str
character(len=90):: infile
character(len=2)::  cn

infile = trim(indir) // '/' // trim(subdir) // '/' // trim(CellData)

open (14,file=infile,form='unformatted',status='old',action='read')

read(14) N_cells, D_box, ni,nj,nk

D_box = D_box * kpc                             !Convert kpc to cm
R_box = D_box / 2d0

allocate(n_HIString(N_cells))
allocate(TString(N_cells))
allocate(V_xString(N_cells))
allocate(V_yString(N_cells))
allocate(V_zString(N_cells))

n_str = 5                                       !Number of "strings"/records

# ifdef AMR
  allocate(LevelString(N_cells))
  n_str = n_str + 1
# endif
# ifdef dust
  allocate(ZString(N_cells))
  allocate(n_HIIString(N_cells))
  n_str = n_str + 2
# endif

call Int2Char(2*n_str+15,cn,.false.)

print*
write(*,'(a'//cn//')',advance='no') &
  'Loading data' // repeat(' ',n_str) // '|' // repeat(achar(8),n_str+1)
# ifdef AMR
  read(14) (LevelString(i),  i = 1,N_cells)
  write(*,'(a)',advance='no') '.'
# endif
read(14) (n_HIString(i),   i = 1,N_cells)
write(*,'(a)',advance='no') '.'
read(14) (TString(i),  i = 1,N_cells)
write(*,'(a)',advance='no') '.'
read(14) (V_xString(i),    i = 1,N_cells)
write(*,'(a)',advance='no') '.'
read(14) (V_yString(i),    i = 1,N_cells)
write(*,'(a)',advance='no') '.'
read(14) (V_zString(i),    i = 1,N_cells)
# ifdef dust
  write(*,'(a)',advance='no') '.'
  read(14) (ZString(i),    i = 1,N_cells)
  write(*,'(a)',advance='no') '.'
  read(14) (n_HIIString(i),    i = 1,N_cells)
# endif
write(*,'(a)',advance='yes') '.'

close(14)

where(n_HIString .lt. 1e-35) n_HIString = 1e-35 !Make sure there are no zeros

# ifdef AMR
  L_max = maxval(LevelString)                   !Maximum level of refinement
# else
  L_max = 0
# endif

end subroutine ReadData

!------------------------------------------------------------------------------

subroutine ReadGalaxies

!Load galactic data.

use PhysicalConstants
use GlobalData
use DataArrays

implicit none
integer::            i,ios
character(len=200):: header,infile

infile = trim(indir) // '/' // trim(subdir) // '/' // trim(GalData)
open (16,file=trim(infile),status='old',action='read')

do                                              !Skip lines beginning with '#'
  read(16,*) header
  if (index(adjustl(header),'#').eq. 0) exit    !Exit after first data line
enddo                                           !Pointer is now on 2nd data line
backspace(16)                                   !So go one back

n_gal = 0
do                                              !Count number of data line, i.e.
  read(16,*,iostat=ios) header                  !number of galaxies
  if (ios .lt. 0) exit
  n_gal = n_gal + 1
enddo
rewind(16)

allocate(x_gal(n_gal))
allocate(y_gal(n_gal))
allocate(z_gal(n_gal))
allocate(r_vir(n_gal))
allocate(v_x(n_gal))
allocate(v_y(n_gal))
allocate(v_z(n_gal))

do                                              !Read header one more time
  read(16,*) header
  if (index(adjustl(header),'#').eq. 0) exit
enddo        
backspace(16)

do i = 1,n_gal                                  !Read galaxy data
  read(16,*)  &
    x_gal(i), &
    y_gal(i), &
    z_gal(i), &
    r_vir(i), &
    v_x(i),   &
    v_y(i),   &
    v_z(i)
enddo

if (.not.any(abs(y_gal).gt.2.) .and. & !If GalData contains MoCaLaTA rather than
    .not.any(    x_gal .lt.0 ) .and. & !IGMtransfer format, x_gal will contain
    .not.any(    x_gal .gt.15)) then   !logM* and  y_gal will contain [O/H]
  print*, "Looks like '"// trim(GalData) //"' has the wrong format."
  stop
endif

close(16)

x_gal = x_gal*kpc + R_box                       !
y_gal = y_gal*kpc + R_box                       ! [-R/kpc, +R/kpc] -> [0, D/cm]
z_gal = z_gal*kpc + R_box                       !

r_vir = r_vir * kpc

v_x   = 1e5 * v_x                               !
v_y   = 1e5 * v_y                               ! km/s -> cm/s
v_z   = 1e5 * v_z                               !

end subroutine ReadGalaxies

!------------------------------------------------------------------------------

subroutine ReadInput

!Read input parameter file

use PhysicalConstants
use GlobalData
use ObservationalParameters

implicit none
character(len=200):: dummy

read*, indir    !Directory containing the data subdirectory
read*, subdir   !Subdirectory containing data to be read (and written)
read*, CellData !File containing cell parameters
read*, Galdata  !File containing galaxy parameters
read*, outdir   !Directory containing subdirectory for IGMtransfer's output
read*, Ioutfile !Output file from IGMtransfer / input for ProcessIGM
read*, dummy    !Output file from ProcessIGM
read*, n_write  !Number of sightlines traced between each write
read*, n_los    !Number of sightlines per galaxy
read*, SpecRes  !Spectral resolution in bins
read*, BW       !Lower and upper values of wavelength interval in Angstrom
read*, dummy    !Wavelength interval in which to calculate "T(z)"
read*, f_rvir   !Distance in virial radii from galaxy centers to start sightlines
read*, r_eff    !Fraction of total radius to be used for the IGM RT
read*, DustType !'SMC'- or 'LMC'-like dust
read*, f_ion    !Dust destruction factor (0 for max destruction, 1 for no destruction
read*, z        !Redshift of snapshot
read*, H_0      !Hubble constant in km/s/Mpc
read*, Omega_M  !Matter density fraction
read*, Omega_L  !Cosmological constant density fraction
                 
H_0 = H_0 * 1d2 !H_0 is converted from km/s/Mpc to cm/s/kpc

end subroutine ReadInput

!------------------------------------------------------------------------------

subroutine Real2Char(r,a)

!Convert real*8 number to character.
!For now, integral part must be 1 or 2 digits, while decimal fraction must be 2.
!To make delimiter, LZ, and TZ optional, make explicit interface

implicit none 
real*8, intent(in)::                     r
character(len=5), intent(out)::          a
character(len=1)::                       delimiter
integer::                                IntRepr,IntPart,DecFrac
character(len=1)::                       cIP
character(len=2)::                       cDF

if (kind(r) .ne. 8) stop 'Not real kind...' !DOES NOT WORK AS INTENDED!!!
!if (.not.present(delimiter)) delimiter = 'p'
!if (.not.present(LZ))        LZ        = .false.
!if (.not.present(TZ))        TZ        = .false.

delimiter = 'p'

IntRepr = nint(100*r)                           !3.25 -> 325
DecFrac = mod(IntRepr,100)                      !25
IntPart = (IntRepr - DecFrac) / 100             !3      

call Int2Char(IntPart,cIP,.true.)
call Int2Char(DecFrac,cDF,.true.)

a = cIP // delimiter // cDF

end subroutine Real2Char

!------------------------------------------------------------------------------

subroutine RotateVector(nhat_i,nhat,nhat_f)

!Rotate unit vector of outgoing photon so that it is given according to
!incident photon.
!nhat_i is direction of incoming photon, nhat is direction of deflection from
!from nhat_i, i.e. direction of outgoing photon if nhat_i lies along the global
!z-axis. nhat_f is true direction of outgoing photon in global coordinates.

implicit none
real*8,intent(in)::  nhat_i(3),nhat(3)
real*8,intent(out):: nhat_f(3)
real*8::             sq,r,ct,st,cp,sp
real*8::             xi,yi,zi,xf,yf,zf

xi = nhat_i(1)
yi = nhat_i(2)
zi = nhat_i(3)

ct = nhat(3)                                    !cos(theta)
st = sqrt(1d0 - nhat(3)**2)                     !sin(theta)
r  = sqrt(nhat(1)**2 + nhat(2)**2)
cp = nhat(1) / r                                !cos(phi)
sp = nhat(2) / r                                !sin(phi)

if (1d0 - abs(zi) .ge. 1d-4) then
  sq   = sqrt(1d0 - zi**2)
  xf = st/sq * (xi*zi*cp - yi*sp) + xi*ct
  yf = st/sq * (yi*zi*cp + xi*sp) + yi*ct
  zf = -st*cp*sq + zi*ct
else                                            !If nhat_i is close to z-axis
  xf = st*cp
  yf = st*sp
  zf = zi*ct
endif

nhat_f = (/xf,yf,zf/)

end subroutine RotateVector

!------------------------------------------------------------------------------

subroutine sigma(SpecRes,Dnu_D,x,sigma_x)

!Lyman alpha cross section.

use PhysicalConstants

implicit none
integer,intent(in)::        SpecRes
real*8,intent(in)::         x(SpecRes),Dnu_D
real*8,intent(out)::        sigma_x(SpecRes)
real*8,dimension(SpecRes):: xt,z,q,Voigt
real*8::                    a

a  = Dnu_L / (2.*Dnu_D)
xt = x**2
z  = (xt - .855) / (xt + 3.42)

where (z .le. 0.)
  q = 0.
elsewhere
  q = z*(1.+21./xt) * a/(pi*(xt+1.)) &
    * (.1117 + z * (4.421 + z*(-9.207 + 5.674*z)))
endwhere

Voigt   = sqrt(pi) * q + exp(-xt)
sigma_x = f_12 * sqrt(pi) * e**2 / (m_e * c * Dnu_D) * Voigt

end subroutine sigma

!------------------------------------------------------------------------------

subroutine WriteData(n_sofar)

!Write the results of the sightlines covered so far only.
!
!Since the data array ("Gollum") may be quite large, some systems may produce a
!segfault when writing. This may also depend on the compiler (e.g. a segfault
!occuring with ifort may not occur with gfortran), as well as the actual
!dimension sizes of Gollum (it may matter whether Gollum is, say, 1000*1000
!or 10000*100). For this reason Gollum is split up in several (n_rec) records,
!calculated from the maximum size of a record (recmaxMb) in bytes.
!If a segfault still occurs, decreasing recmaxMb should help.
!
!A possible, but not fully tested, alternative solution may be to type
!"ulimit -s unlimited" (csh) or "unlimit stack" (tcsh) in the terminal before
!execution.

use GlobalData
use ObservationalParameters

implicit none
integer::            n_sofar
character(len=200):: taudata
integer::            n,i,n_rec,ifrac
real*8::             recmaxMb

! real*8::             dlambda
! dlambda = (BW(2) - BW(1)) / SpecRes
! do i = 1, SpecRes
!   write (1,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(1,i))
!   write (2,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(2,i))
!   write (3,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(3,i))
!   write (4,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(4,i))
! ! XXXXX (5,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(5,i))
! ! XXXXX (6,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(6,i))
!   write (7,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(7,i))
!   write (8,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(8,i))
!   write (9,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(9,i))
!   write (10,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(10,i))
!   write (11,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(11,i))
!   write (12,*) BW(1)+(i-.5)*dlambda, exp(-Gollum(12,i))
! enddo

recmaxMb = 10.                                  !Maximum record size in Mb
n_rec    = size(Gollum(1:n_sofar,:)) / nint(1e6 * recmaxMb / kind(Gollum)) + 1
ifrac    = n_sofar / n_rec

print*, 'Writing data...'
taudata = trim(outdir)//'/'//trim(subdir)//'/'//trim(Ioutfile)
open(14,file=trim(taudata),form='unformatted',status='replace',action='write')
write(14) n_rec, n_sofar
do i = 1,n_rec
  if (i .lt. n_rec) then
    write(14) (Gollum(n,:), n = (i-1)*ifrac+1, i*ifrac)
  else
    write(14) (Gollum(n,:), n = (i-1)*ifrac+1, n_sofar)
  endif
enddo
close(14)

end subroutine WriteData

!------------------------------------------------------------------------------

subroutine XsecD(SpecRes,Dnu_D,x,sigma_d)

!"Gnedin/Pei/Weingartner & Draine" dust cross section.
!
!The cross section is expressed per hydrogen nucleus, such that the optical
!depth of dust, when traveling through a column of hydrogen atom N_H is
!N_H * XsecD, if the metallicity is equal to the metallicity of the reference
!galaxy given by the keyword 'DustType'.
!
!That is, for an SMC-like dust type, the "n_d" values should be n_H * Z/Z_SMC,
!while for LMC-like dust, it should be n_H * Z/Z_LMC. The reference metallicity
!can be done on an individual metal basis, or just using the average
!metallicity; the difference is not very large. The SMC (LMC) metallicity is
!on average 0.6 +/- 0.2 (0.3 +/- 0.2) dex below Solar (Welty et al. 1997/99).
!
!In Laursen et al. 2011, ApJ, 728, 52 we argue that rather than scaling with
!the total amount of hydrogen nuclei, dust destruction can be modeled by
!scaling with n_HI plus some fraction of n_HII, since regions with ionized
!hydrogen are often correlated with regions where dust may be destroyed.
!If this fraction is 0, dust is fully destroyed when hydrogen is ionized, while
!if it is 1, we're back at the original prescription. We further argue that
!observationally a value of 0.01 is fair, and show that using this value,
!the escape fraction of Lya from a galaxy lies approximately in the middle of
!what is obtained when using the two extrema.
!
!Bottom line: If you want to use this prescription, set the "dust density" in
!the i'th cell equal to
!   n_d(i) = (n_HI(i) + 0.01*n_HII(i)) * (Z(i)/Z_sun) / 10**(-0.6)
!for SMC dust. For LMC dust, simply swap the factor "0.6" with "0.3".
!
!Remember to set

use PhysicalConstants
use GlobalData

implicit none
integer,intent(in):: SpecRes
real*8,intent(in)::  Dnu_D, x(SpecRes)
real*8,intent(out):: sigma_d(SpecRes)

if (DustType .eq. 'SMC') then
  sigma_d = (.395083 + 1.7164865d-16 * Dnu_D * x) * 1d-21
elseif (DustType .eq. 'LMC') then
  sigma_d = (.723314 + 4.2156414d-16 * Dnu_D * x) * 1d-21
else
  print*, "Don't know dust type '" // DustType // "'"
  stop
endif

end subroutine XsecD

!------------------------------------------------------------------------------

recursive subroutine BuildNestedUniverse(CurrentCell,CurrentLevel)

use CellStructure
use GlobalData
use DataArrays

implicit none
type(Cell), target::     CurrentCell
integer, intent(IN)::    CurrentLevel
integer::                ii,jj,kk

i_cell = i_cell + 1      

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(i_cell)
    CurrentCell%n_HI   = n_HIString(i_cell)!6.583e-11
    CurrentCell%Dnu_D  = TString(i_cell)!1.0566362e+11!
    CurrentCell%U_bulk = (/V_xString(i_cell), &
                           V_yString(i_cell), &
                           V_zString(i_cell)/)
#   ifdef dust
    CurrentCell%n_d    = ZString(i_cell)
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel  
    allocate(CurrentCell%Child(2,2,2))
    
    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined  = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)   
          call BuildNestedUniverse(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop 
  endif
# else
    CurrentCell%n_HI   = n_HIString(i_cell)
    CurrentCell%Dnu_D  = TString(i_cell)
    CurrentCell%U_bulk = (/V_xString(i_cell), &
                           V_yString(i_cell), &
                           V_zString(i_cell)/)
#   ifdef dust
    CurrentCell%n_d    = ZString(i_cell)
#   endif
# endif

end subroutine BuildNestedUniverse

!------------------------------------------------------------------------------

recursive subroutine DigDeeper(X_pos,CurrentCell)

!Check if cell is refined. If not, this is the host cell. If so, go one level
!deeper.

use CellStructure
use GlobalData

implicit none
real*8, intent(IN)::  X_pos(3)
type(Cell), pointer:: CurrentCell               !intent(INOUT)
integer::             ii,jj,kk

# ifdef AMR
  ii = merge(1,2,X_pos(1) .lt. CurrentCell%C_pos(1))
  jj = merge(1,2,X_pos(2) .lt. CurrentCell%C_pos(2))
  kk = merge(1,2,X_pos(3) .lt. CurrentCell%C_pos(3))

  CurrentCell => CurrentCell%Child(ii,jj,kk)
    
  if (CurrentCell%Refined) call DigDeeper(X_pos,CurrentCell)
# endif

end subroutine DigDeeper

!------------------------------------------------------------------------------

recursive subroutine EveryoneShouldHaveACenter(CurrentCell,i,j,k,ParentPos)

!Assign center to all cells, even if refined.

use CellStructure
use GlobalData

implicit none
type(Cell), target ::           CurrentCell
integer,intent(in)::            i,j,k
real(kind=posKind),intent(in):: ParentPos(3)
integer::                       indexico(3),ii,jj,kk
real*8:: norm

indexico = (/i,j,k/)

# ifdef AMR
  if (CurrentCell%Level .eq. 0) then
    CurrentCell%C_pos = (indexico - .5d0) * dx(CurrentCell%Level)
  elseif (CurrentCell%Level .gt. 0) then
    CurrentCell%C_pos = ParentPos &
                      + 2*(indexico - 1.5d0) & ! -1 or +1
                      * dxH(CurrentCell%Level)
  else
    stop "Nej nej nej hvad sker der?!"
  endif

  if (CurrentCell%Refined) then
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          call EveryoneShouldHaveACenter(CurrentCell%Child(ii,jj,kk),ii,jj,kk,CurrentCell%C_pos)
        enddo
      enddo
    enddo
  endif
# else
    CurrentCell%C_pos = (indexico - .5d0) * dx(0)
# endif

! CurrentCell%U_bulk = (CurrentCell%C_pos-R_box)/kpc*0.028470507 * normalize?

end subroutine EveryoneShouldHaveACenter

!------------------------------------------------------------------------------

function Doppler(T)

!Doppler frequency width

use PhysicalConstants

implicit none
real*4:: Doppler,T

Doppler = sqrt(2.*k_B*T / m_H) * nu_0/c

end function 

!------------------------------------------------------------------------------

function H(z)

!Hubble parameter as a function of redshift in cm/s/cm (!)
!NOT THE SAME AS IN MoCaLaTA!!!

use PhysicalConstants

implicit none
real*8:: H,z

H = (1 + z) * sqrt(1 + Omega_M*z + Omega_L*(1/(1 + z)**2 - 1)) * H_0 / kpc

end function

!------------------------------------------------------------------------------

function norm(vec)

implicit none
real*8:: vec(3),norm

norm = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

end function 

!------------------------------------------------------------------------------

function v2u(v,Dnu_D)

!Converts a velocity in km/s to the dimensionless quantity "velocity in terms
!of velocity Doppler width".

use PhysicalConstants

implicit none
real*4:: v2u,v,Dnu_D

v2u = v*1e5 * nu_0 / (Dnu_D * c)

end function 

!------------------------------------------------------------------------------
