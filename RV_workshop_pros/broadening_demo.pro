; Create some plots to illustrate the rotational broadening vs. LSF


function make_lsf_kernel, LSF_0=lsf_0, OVERSAMP=oversamp, FNUM=fnum

if n_elements(lsf_0) eq 0 then lsf_0 = 3000
if n_elements(oversamp) eq 0 then oversamp = 7
if n_elements(fnum) eq 0 then fnum = 0 ; can only choose up to 9 right now

nplsf = 20L*oversamp+1 ; npix lsf kernel

; Read in LSF coefficients
napofiles = 10
lsf_str = readin_apo(nfiles = napofiles, hdu = 8, /flux)

lsf_coeff_arr = dblarr(28, napofiles)
for i = 0, napofiles-1 do begin
    lsf_coeff_arr[0,i] = *(lsf_str[i].output)
endfor

; make lsf
lx = (dindgen(nplsf)-(nplsf/2))/oversamp + lsf_0
lsf_xnorm = lx - lsf_0
lsf_coeff = lsf_coeff_arr[*,fnum]
lsf = lsf_gh(lx, lsf_0, lsf_coeff)/oversamp

; make xvector in velocity units
bigc = 3d5 ;km/s
dlogwl = 6d-6
; change in velocity equal to 1 oversampled logwl pixel
deltav = ((10d0^(dlogwl))-1d0)*bigc
lsf_xvel = lsf_xnorm * deltav

return, {lsf:lsf, lsf_x:lx, lsf_xnorm:lsf_xnorm, lsf_xvel:lsf_xvel}

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 


function make_broadening_kernel, vsini, OVERSAMP=oversamp
;;; ROTATION
if n_elements(oversamp) eq 0 then oversamp=77
bigc = 3d5 ;km/s
dlogwl = 6d-6
; change in velocity equal to 1 oversampled logwl pixel
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc

help, deltav

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid)
;help, velgrid
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))
;help, velgrid_kernel_x
; Make a similar vector which matches the logwl-grid
nlsfrot = n_elements(velgrid)+6 ;pad to make sure we get the whole kernel
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
lsf_x = (dindgen(nlsfrot)-(nlsfrot/2))/oversamp
;help, lsf_rot_kernel_x
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x, /spline)

; Make sure interpolation doesn't introduce negatives
rot_neg_idx = where(lsf_rot lt 0d0, ncnt)
if ncnt gt 0 then lsf_rot[rot_neg_idx]=0d0
; Normalize

lsf_rot = lsf_rot/total(lsf_rot, /double)

; Recover velgrid
nvelgrid = n_elements(velgrid) 
lsf_velgrid_x = (dindgen(nvelgrid)-(nvelgrid/2))/oversamp
velgrid_out = interpol(velgrid, lsf_velgrid_x, lsf_x)

return, {lsf_rot:lsf_rot, kernel_logwl:lsf_rot_kernel_x, kernel_x:lsf_x, velgrid:velgrid_out}

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro broadening_demo

oversamp = 99

; Get the LSF vectors
lsfstr = make_lsf_kernel(oversamp=oversamp, fnum=1)
lsf = lsfstr.lsf
lsf_xnorm = lsfstr.lsf_xnorm
lsf_xvel = lsfstr.lsf_xvel

;Set up the plot

;load rainbow color table with discrete color tags	      
loadct,39,/silent
setcolors,/system_variables,/silent,decomposed=0  ;requires setcolors.pro

;plotting parameters to set axes and fonts
!p.background=0    ;default background color
!p.color=255
!p.charsize=1.7		;text default size
!p.charthick=6		;text default thickness
!x.thick=5		;thicken x-axis 
!y.thick=5		;thicken y-axis
!p.font=1		;set default font
!p.thick=5		;set default plotting line thickness
!x.margin=[7,4]

;set psym=8 to circle by default
circind=findgen(40) * (!pi*2/39.)
defsysv,'!circ',transpose( [[cos(circind)],[sin(circind)]])  ;user symbol vertices
usersym,!circ(0,*),!circ(1,*),/fill  ;circle plot

psopen, 'broadening_demo2.eps', /encapsulated, xs=10, ys=10, /inches, /color

!p.multi=[0,2,2]

; Get the rotational broadening
vsinivec = [1d0, 5d0, 10d0, 20d0]

foreach vsini, vsinivec, vidx do begin

    lsfrot = make_broadening_kernel(vsini, oversamp=oversamp)
    lsf_rot = lsfrot.lsf_rot
    lsf_rotx = lsfrot.kernel_x
    vel_x = lsfrot.velgrid

    ;Plot
;    if vidx le 1 then begin
        plot, lsf_xnorm, lsf, $
          xr=[min(lsf_xnorm)<min(lsf_rotx), max(lsf_xnorm)>max(lsf_rotx)], $
          yr=[0,max(lsf)>max(lsf_rot)], $
          charsize=2, linest=2, xstyle=8, xtit="Pixel"
          ;title="Vsini = "+strtrim(fix(vsini),2)+" km/s", charsize=2, linest=2, $
          
    
        oplot, lsf_rotx, lsf_rot
        xtickvalue = [-40,-30,-20,-10,10,20,30,40]
        help, lsf_xnorm, lsf_xvel
        print, minmax(lsf_xnorm)
        print, minmax(lsf_xvel)
;        stop
        tickval_map = interpol(lsf_xnorm, lsf_xvel, xtickvalue)
        axis, xaxis=1, xticks = n_elements(xtickvalue), xtickv=tickval_map, xtickn=strtrim(xtickvalue,2), charthick=6, charsize=1.7, xtit="Delta V (km/s)"
        
    ; endif else begin

;         plot, lsf_xnorm, lsf, $
;           xr=[min(lsf_xnorm)<min(lsf_rotx), max(lsf_xnorm)>max(lsf_rotx)], $
;           yr=[0,max(lsf)>max(lsf_rot)], $
;           title="Vsini = "+strtrim(fix(vsini),2)+" km/s", charsize=2, linest=2, $
;           xtit="Delta V (km/s)"
    
;         oplot, lsf_rotx, lsf_rot
;         axis, xaxis=1, xtickv=[-20,-15,-10,-5,0,5,10,15,20], xtickn=strtrim(xtickv,2), charthick=6, charsize=2
;     endelse

endforeach

psclose


stop


end
