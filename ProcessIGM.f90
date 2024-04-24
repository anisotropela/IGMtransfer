!================================  Modules  ===================================

module PhysicalConstants
  implicit none
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
end module PhysicalConstants

!------------------------------------------------------------------------------

module GlobalData
  integer::              N_cells,i_cell         !Total # of cells; counter
  integer::              n_los                  !# of sightlines
  real*8::               D_box,R_box            !Phys. size/half size of Box
  real*8::               z,BW(2),BW_stat(2),dlambda !Redshift; wavelength interval
  character(len=90)::    subdir,indir           !Home dir of data
  character(len=90)::    infile,Poutfile
  real*8::               d_tot,r_zoa,d_halfmax,mu_min,mu_max,dmu
  real*4,allocatable::   taudata(:,:)
  real*8,allocatable::   PDF_F(:),F_lam(:,:)
  real*4,allocatable::   Fi(:),Fn(:)
  integer::              n_gal                  !Total # of galaxies
  integer::              SpecRes,i1,i2,bintot
  real*8::               Tmean,T,Tlo,Thi        !Universal transmission
  real*8::               T1,T2,T3               !Universal transmission
end module GlobalData

!------------------------------------------------------------------------------

!==============================================================================

program ProcessIGM

use PhysicalConstants
use GlobalData

implicit none

call ReadInput
call ReadData
call CalcTransmission
call WriteData

end program ProcessIGM

!==============================  Subroutines  =================================

subroutine CalcTransmission

use GlobalData
use PhysicalConstants

implicit none
real*8::              lo,hi,loglo,loghi,DDF,n_dex,dF,aveF,sigF,skewF
integer::             i16,i50,i84
real*8::              p68,p34,p16,p84
real*8::              lambda,median,minus1sig,plus1sig
integer::             i,j,l,n,bpd,iave
integer,allocatable:: tauindex(:,:)
logical::             loggish

!Wavelength interval in AAngstroem in which to calculate the transmission
lo = BW_stat(1)                                 !Lower limit
hi = BW_stat(2)                                 !Upper limit

allocate(F_lam(SpecRes,5))                      !Flux & stddev as a func. of wavelength
allocate(Fi(n_los))                             !Unsorted fluxes at a given wavelength

p68 = 0.68268949
p34 = p68 / 2
p16 = .5 - p34
p84 = .5 + p34

i16 = nint(p16*n_los)
i50 = n_los / 2
i84 = nint(p84*n_los)

  write(*,'(a38)',advance='no') 'Processing...           |' // repeat(achar(8),12)
  do i = 1,SpecRes                         !
    Fi = exp(-taudata(:,i))                !   In each wavelangth bin, calc.
                                           !   the median of all los, along
    call Quicksort(n_los,Fi)               !   with 16 and 84 perc., and mean.
                                           !
    F_lam(i,1) = BW(1) + (i -.5) * dlambda !               ||
    F_lam(i,2) = Fi(i50)                   !               \/
    F_lam(i,3) = Fi(i16)                   !
    F_lam(i,4) = Fi(i84)                   !        THIS SHOULD BE THE
    F_lam(i,5) = sum(Fi) / n_los           !        CORRECT "F(lambda)"
    if (mod(i,SpecRes/10).eq.0) write(*,'(a)',advance='no') '.'
  enddo
  write(*,'(a)',advance='yes') '.'

  do i = 1,SpecRes                     !
    lambda = F_lam(i,1)                !
    if (lambda .ge. lo) then           !
      i1 = i                           !
      exit                             !    Find hi,lo-indices for interval
    endif                              !    in which to calculate T == <F_blue>
  enddo                                !
  do i = 1,SpecRes                     !
    lambda = F_lam(i,1)                !
    if (lambda .ge. hi) then           !
      i2 = i                           !
      exit                             !
    endif                              !
  enddo                                !
  bintot = i2 - i1 + 1                 !

  do n = 1,n_los                                !Calc. AVERAGE of all
    Fi(n) = sum(exp(-taudata(n,i1:i2))) / bintot!INTENSITIES in interval
  enddo                                         !

  call Quicksort(n_los,Fi)

  T1 = Fi(i50) ! MEDIAN of all AVERAGES
  T2 = Fi(i16)
  T3 = Fi(i84)

end subroutine CalcTransmission

!------------------------------------------------------------------------------

subroutine Quicksort(n,arr)

implicit none
integer,intent(in)::   n
real*4,intent(inout):: arr(n)
integer,parameter::    m=7, nstack=50
integer::              i,ir,j,jstack,k,l,istack(nstack)
real*4::               a,temp
logical::              flag

!This subroutine is a f90-ed version of the f77 subroutine given in Numerical
!Recipes in Fortran, [FIND VERSION]
!It sorts an array arr(1:n) into ascending numerical order using the Quicksort
!algorithm. n is input; arr is replaced on output by its sorted rearrangement.
!Parameters: m is the size of subarrays sorted by straight insertion and nstack
!is the required auxiliary storage. 

jstack = 0
l      = 1
ir     = n

do
  if (ir-l .lt. m) then     !Insertion sort when subarray small enough. 
    do j = l+1,ir
      a = arr(j)

      flag = .true.
      do i = j-1,l,-1
        if (arr(i) .le. a) then
          flag = .false.
          exit
        endif
        arr(i+1) = arr(i)
      enddo

      if (flag) i = l-1      !This flag shouldn't be necessary, if "i" always
      arr(i+1) = a           !has its final value when loop is exited. Does it?
    enddo

    if (jstack .eq. 0) exit

    ir     = istack(jstack)  !Pop stack and begin a new round of partitioning. 
    l      = istack(jstack-1)
    jstack = jstack - 2
  else
    k        = (l+ir)/2      !Choose median of left, center, and right elements
                             !as partitioning element a. Also rearrange so that
                             !a(l) ≤ a(l+1) ≤a(ir). 
    temp     = arr(k)
    arr(k)   = arr(l+1)
    arr(l+1) = temp

    if (arr(l) .gt. arr(ir)) then
      temp    = arr(l)
      arr(l)  = arr(ir)
      arr(ir) = temp
    endif

    if (arr(l+1) .gt. arr(ir)) then
      temp     = arr(l+1)
      arr(l+1) = arr(ir)
      arr(ir)  = temp
    endif

    if (arr(l) .gt. arr(l+1)) then
      temp     = arr(l)
      arr(l)   = arr(l+1)
      arr(l+1) = temp
    endif

    i = l + 1                !Initialize pointers for partitioning. 
    j = ir
    a = arr(l+1)             !Partitioning element. 

    do                       !Beginning of innermost loop. 
      i = i+1                !Scan up to find element > a. 
      if (arr(i) .lt. a) cycle

      do
        j = j - 1            !Scan down to find element < a. 
        if (arr(j) .le. a) exit
      enddo

      if (j .lt. i) exit     !Pointers crossed. Exit with partitioning complete. 

      temp   = arr(i)        !Exchange elements. 
      arr(i) = arr(j)
      arr(j) = temp
    enddo                    !End of innermost loop. 

    arr(l+1) = arr(j)        !Insert partitioning element. 
    arr(j)   = a
    jstack   = jstack + 2
                             !Push pointers to larger subarray on stack,
                             !process smaller subarray immediately. 
    if (jstack .gt. nstack) stop 'nstack too small in sort'

    if (ir-i+1 .ge. j-l) then
      istack(jstack)   = ir
      istack(jstack-1) = i
      ir               = j - 1
    else
      istack(jstack)   = j - 1
      istack(jstack-1) = l
      l                = i
    endif
  endif
enddo

end subroutine Quicksort

!------------------------------------------------------------------------------

subroutine ReadData

use GlobalData
use PhysicalConstants

implicit none
character(len=200):: infilefull
integer::            i,j,n,ifrac,n_rec

infilefull = trim(indir) //'/'// trim(subdir) //'/'// trim(infile)
print*, "Loading '", trim(infilefull), "'..."
open(14,file=trim(infilefull),form='unformatted',status='old',action='read')
read(14) n_rec, n_los

ifrac = n_los / n_rec

allocate(taudata(n_los,SpecRes))

do i = 1,n_rec
  if (i .lt. n_rec) then
    read(14) (taudata(n,:), n = (i-1)*ifrac+1, i*ifrac)
  else
    read(14) (taudata(n,:), n = (i-1)*ifrac+1, n_los)
  endif
enddo

close(14)

dlambda = (BW(2) - BW(1)) / SpecRes

!!Incomment this block to write ten spectra (#11-20) to files "fort.11-20"
! do i=1,SpecRes
!   do j=11,20
!     write(j,*) BW(1)+i*dlambda, exp(-taudata(j,i))
!   enddo
! enddo

end subroutine ReadData

!------------------------------------------------------------------------------

subroutine ReadInput

!Read input parameter file

use PhysicalConstants
use GlobalData

implicit none
character(len=200):: dummy

read*, dummy    !Directory containing the data subdirectory
read*, subdir   !Subdirectory containing data to be read (and written)
read*, dummy    !File containing cell parameters
read*, dummy    !File containing galaxy parameters
read*, indir    !Directory containing subdirectory for IGMtransfer's output
read*, infile   !Output file from IGMtransfer / input for ProcessIGM
read*, Poutfile !Output file from ProcessIGM
read*, dummy    !Number of sightlines traced between each write
read*, dummy    !Number of sightlines per galaxy
read*, SpecRes  !Spectral resolution in bins
read*, BW       !Lower and upper values of wavelength interval in Angstrom
read*, BW_stat  !Wavelength interval in which to calculate "T(z)"
read*, dummy    !Distance in virial radii from galaxy centers to start sightlines
read*, dummy    !Fraction of total radius to be used for the IGM RT
read*, dummy    !'SMC'- or 'LMC'-like dust
read*, dummy    !Dust destruction factor (0 for max destr., 1 for no destr.)
read*, z        !Redshift of snapshot
read*, dummy    !Hubble constant in km/s/Mpc
read*, dummy    !Matter density fraction
read*, dummy    !Cosmological constant density fraction

end subroutine ReadInput

!------------------------------------------------------------------------------

subroutine WriteData

use GlobalData

implicit none
integer::            i,n
character(len=200):: outfilefull

outfilefull = trim(indir) //'/'// trim(subdir) //'/'// trim(Poutfile)

open(14,file=trim(outfilefull),form='formatted',status='replace',action='write')

do i = 1,SpecRes
  write(14,*) real(F_lam(i,:))
enddo

print*, 'Transmission calculated in lambda = ', real((/F_lam(i1,1),F_lam(i2,1)/))
print*
print*, ' Transmitted fraction at z =', real(z), 'is:'
print*, ' median          16%             84%'
print*, real((/T1, T2, T3/))

close(14)

end subroutine WriteData

!------------------------------------------------------------------------------
