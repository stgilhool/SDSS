; Quick thing to make a plot of the vsini's from the epsilon test

pro epsilon_test_display

a6 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_with_fit_full/rfile_master.fits',1)
b6 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_with_fit_full/rfile_master.fits',2)

a4 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_eptest_4/rfile_master.fits',1)
b4 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_eptest_4/rfile_master.fits',2)

a8 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_eptest_8/rfile_master.fits',1)
b8 = mrdfits('/home/stgilhool/APOGEE/vsini_results/logg45_eptest_8/rfile_master.fits',2)

vi = where(b4.valid_flag eq 1, ngood)

if ngood ne 688 then message, "oops"

; sort by difference between 4 and 8 guys
a6 = a6[vi]
a4 = a4[vi]
a8 = a8[vi]
b6 = b6[vi]
b4 = b4[vi]
b8 = b8[vi]

res48 = a4.vsini_best-a8.vsini_best

dsort = sort(abs(res48))

v6 = a6[dsort].vsini_best
v4 = a4[dsort].vsini_best
v8 = a8[dsort].vsini_best

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

psopen, 'limb_darkening_test.eps', /encapsulated, xs=10, ys=10, /inches, /color

plot, v8, /nodata, ytit = "Vsini", xr=[375,700], title='Results of Vsini fit using different LD parameters'
oplot, v8, ps=8, color = !red
oplot, v6, ps=8, color=!yellow
oplot, v4, ps=8, color=!blue

al_legend, ['0.8','0.6','0.4'], colors=[!red,!yellow,!blue], psym=8

psclose

stop

end
