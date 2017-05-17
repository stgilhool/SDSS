; Routine to plot vsini results vs ASPCAP Teff
pro show_vsini_results

; Load ASPCAP data
a = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)

; Load latest vsini data

r = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_logg5_compare/rfile_master.fits',1)
;r = mrdfits('/home/stgilhool/APOGEE/vsini_results/initial_run/rfile_master.fits',1)
r2 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_logg5_compare/rfile_master.fits',2)

t = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)

gd_idx = where(r2.valid_flag eq 1, ngood)

teff = a.teff
vsini = r.vsini_best_45
;vsini = r.vsini_best
v2 = t.vsini_vfit

teff = teff[gd_idx]
vsini = vsini[gd_idx]
v2 = v2[gd_idx]


;;; New plot
;load rainbow color table with discrete color tags	      
loadct,39,/silent
setcolors,/system_variables,decomposed=0,/silent  ;requires setcolors.pro

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


PSOPEN,'vsini_results.eps',encapsulated=1,/INCHES,XSIZE=10,YSIZE=6,/COLOR

errstring_teff = textoidl('\sigma_{Teff}')
errstring_vsini = textoidl('\sigma_{Vsini}')
errstring = textoidl('\sigma')
xtitstr = textoidl('T_{eff} (K)')

det_idx = where(vsini gt 5, ndet, complement=ndet_idx, ncomplement=nndet)

plot, teff, vsini, xr=[4000,2600], ps=8, /nodata, $
  tit = "Vsini results from PHOENIX fitting", $
  xtit = xtitstr, $
  ytit = "Vsini (km/s)"
oplot, teff[det_idx], vsini[det_idx], ps=8, symsize=1.5
oplot, teff[det_idx], vsini[det_idx], ps=8, symsize=0.7, color=!blue
plotsym, 5, /fill; upper limit plot
oplot, teff[ndet_idx], replicate(5d0,nndet), ps=8, symsize=1.5
oplot, teff[ndet_idx], replicate(5d0,nndet), ps=8, symsize=0.7, color=!red
usersym,!circ(0,*),!circ(1,*),/fill  ;circle plot

oploterror, [2700], [70], replicate(68d0,2), replicate(3,2), /nohat, ps=8, symsize=1.5
oploterror, [2700], [70], replicate(68d0,2), replicate(3,2), /nohat, ps=8, errcolor=!blue

;xyouts, [2690], [61], errstring_teff
;xyouts, [2690], [57], errstring_vsini
xyouts, [2675], [67], errstring, charsize = 2
oplot, [2600,2800],[65,65]
oplot, [2600,2800],[75,75]
oplot, [2600,2600],[65,75]
oplot, [2800,2800],[65,75]

psclose

stop

end
