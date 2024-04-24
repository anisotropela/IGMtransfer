program dat2bin

implicit none
integer::               N_cells,ni,nj,nk,i
real*8::                D_box
integer*4,allocatable:: indexString(:), LevelString(:)
real*4,allocatable::    n_HIString(:), TString(:)
real*4,allocatable::    V_xString(:), V_yString(:), V_zString(:)
character(len=200)::    header

! Grid data -------------------------------------------------------------------
N_cells = 146 !Total number of cells
D_box   = 30. !Box side length in kpc
ni      = 3   !Base grid resolution in x-dir.
nj      = 3   !Base grid resolution in y-dir.
nk      = 3   !Base grid resolution in z-dir.

! Allocate arrays -------------------------------------------------------------
allocate(indexString(N_cells))
allocate(LevelString(N_cells))
allocate(n_HIString(N_cells))
allocate(TString(N_cells))
allocate(V_xString(N_cells))
allocate(V_yString(N_cells))
allocate(V_zString(N_cells))

! Read ASCII data -------------------------------------------------------------
open (1,file='CellData.dat',form='formatted',status='old',action='read')
do i = 1,4; read(1,*) header; enddo ! Read in first four lines
do i = 1,N_cells                    ! Read in consecutive lines
  read(1,*) indexString(i), LevelString(i), n_HIString(i), TString(i), &
    V_xString(i), V_yString(i), V_zString(i)
enddo
close(1)

! Write binary data -----------------------------------------------------------
open (2,file='CellData.bin',form='unformatted',status='replace',action='write')
write(2) N_cells,D_box,ni,nj,nk
write(2) (LevelString(i), i = 1,N_cells)
write(2) (n_HIString(i),  i = 1,N_cells)
write(2) (TString(i),     i = 1,N_cells)
write(2) (V_xString(i),   i = 1,N_cells)
write(2) (V_yString(i),   i = 1,N_cells)
write(2) (V_zString(i),   i = 1,N_cells)
close(2)

end program dat2bin
