; Example tool for visualizing the output from ProcessIGM.
; Input file is the same as used for IGMtransfer and ProcessIGM.
; Output file will be put in the same directory as input file.

pro F_lam, infile, outfile=outfile, xrange=xrange, yrange=yrange

if n_elements(outfile) eq 0 then outfile = 'F_lam.eps'
if n_elements(xrange)  eq 0 then xrange  = [1213., 1218.]
if n_elements(yrange)  eq 0 then yrange  = [0., 1.1]

; Read infile

n_in  = 20                                      ;# of lines in .in-file
l_dir = 4                                       ;Line # of mother-dir (1st is 0)
l_sub = 1                                       ;Line # of subdir
l_dat = 6                                       ;Line # of ProcessIGM-output-file
close,1
openr,1,infile
indata = sindgen(n_in)
readf,1,indata
close,1

; Determine name of indata-file

dir = indata[l_dir]
ih  = strpos(dir,'#')
if ih ne -1 then dir = strmid(dir,0,ih)
is  = strpos(dir,' ')
if is ne -1 then dir = strmid(dir,0,is)
if strmid(dir,0,1) eq "'" or strmid(dir,0,1) eq '"' then $
  dir = strmid(dir,1,strlen(dir)-2)

sub = indata[l_sub]
ih  = strpos(sub,'#')
if ih ne -1 then sub = strmid(sub,0,ih)
is  = strpos(sub,' ')
if is ne -1 then sub = strmid(sub,0,is)
if strmid(sub,0,1) eq "'" or strmid(sub,0,1) eq '"' then $
  sub = strmid(sub,1,strlen(sub)-2)

dat = indata[l_dat]
ih  = strpos(dat,'#')
if ih ne -1 then dat = strmid(dat,0,ih)
is  = strpos(dat,' ')
if is ne -1 then dat = strmid(dat,0,is)
if strmid(dat,0,1) eq "'" or strmid(dat,0,1) eq '"' then $
  dat = strmid(dat,1,strlen(dat)-2)

datafile = dir +'/'+ sub +'/'+ dat

; Read indata

close,1
openr,1,datafile
spawn, 'wc -l ' + datafile, cn
n    = long(total(long(cn)))                    ; = SpecRes
data = findgen(5,n)
readf,1,data
close,1

; Define arrays

lam = findgen(n)
med = findgen(n)
p16 = findgen(n)
p84 = findgen(n)

for i = 0,n-1 do begin    ; necessary so that arrays are [n] and not [1,n]
  lam[i] = data(0,i)
  med[i] = data(1,i)
  p16[i] = data(2,i)
  p84[i] = data(3,i)
endfor

; Find lower and upper indices

dlam = lam[1] - lam[0]

for i = 0,n-1 do begin
  l  = lam[0] + i*dlam
  i1 = i-1
  if l gt xrange[0] then break
endfor
for i = 0,n-1 do begin
  l  = lam[0] + i*dlam
  i2 = i
  if l gt xrange[1] then break
endfor

; Plot 
; print, dir +'/'+ sub +'/'+ outfile

!p.multi=[0,1,1]
set_plot,'ps'
device,/cmyk,/encap,$
       filename=dir +'/'+ sub +'/'+ outfile, $
       /color,xsize=14, ysize=12, xoffset=1., yoffset=1.
loadct,0,/silent

lam0 = 1215.67                                  ;Line center
c    = 3e5                                      ;Speed of light
gray = 220                                      ;Color index of conf. interval
th   =  4                                       ;Line thickness
sz   =  1.2                                     ;Character size
AA   = string("305B);"                          ;Angstrom symbol

plot, lam, data[1,*], xr=xrange,yr=yrange,xst=5,/yst, $
      xtit='!4k!3!Drest!N/'+AA, ytit='F',$
      position=[.15,.15,.85,.85], $
      xth=th,yth=th,th=th,charth=th,chars=sz,/nodata

; Draw confidence interval
  polyx = [lam[i1:i2]   >xrange[0]<xrange[1],lam[i2]<xrange[1], $
   reverse(lam[i1:i2-1])>xrange[0]<xrange[1],lam[i1]>xrange[0]]
  polyy = [p16[i1:i2]   >yrange[0]<yrange[1],$
   reverse(p84[i1:i2])  >yrange[0]<yrange[1]]
  polyfill,polyx,polyy, col=gray

; Redraw axes + make vel-axis
  axis,xaxis=0,xr=xrange,xst=1, xth=th, chars=sz, charth=th, xtit='!4k!3!Drest!N/'+AA
  axis,xaxis=1,xr=(xrange-lam0)/lam0*c,xst=1, xth=th, chars=sz, charth=th, $
       xtit='!4D!3v/km s!U-1!N'
  axis,yaxis=0,yr=yrange,yst=1, yth=th, chars=0.001
  axis,yaxis=1,yr=yrange,yst=1, yth=th, chars=0.001

; Median + line center
  oplot, lam, med, th=th
  oplot, [lam0,lam0], yrange, lines=2, th=th

device,/close
set_plot,'x'

end
