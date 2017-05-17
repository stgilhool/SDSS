;example plotting routine
;requires user defined symbol (psym=8) to be already defined, setcolors 'active'
PRO plot_example


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

;dummy arrays for plotting
npts = 100
seed = 4363254325d
x = DINDGEN(npts) + 50d
y = RANDOMU(seed,npts)*10d + 100d
err = SQRT(y)/2d

;axis ranges, used for all subsequent plot commands
xr = [MIN(x) - 5d,MAX(x) + 5d]
yr = [MIN(y) * 0.95d, MAX(y) * 1.05d]

;save eps file
PSOPEN,'plotting_example.eps',encapsulated=1,/INCHES,XSIZE=10,YSIZE=6,/COLOR

;plot grayed out axes with gridded lines first (no data involved)
PLOT,x,y,xrange=xr,yrange=yr,/XSTY,/YSTY,title='!4Plot title',xtitle='!4X value',ytitle='!4Y value', $
	xgridstyle=1,ygridstyle=1,xticklen=1,yticklen=1,co=!gray,/NODATA

;replot the same axis in black, without the gridded lines
PLOT,x,y,xrange=xr,yrange=yr,/XSTY,/YSTY,title='!4Plot title',xtitle='!4X value',ytitle='!4Y value', $
	/NOERASE,/NODATA

;plot error bars in gray first, if desired
FOR i = 0, N_ELEMENTS(y) - 1 DO OPLOT,[x[i],x[i]],[y[i] - err[i]/2d, y[i] + err[i]/2d],co=!gray,thick=7

;plot x and y values in black with large symbols
OPLOT,x,y,psym = 8,symsize=1.3

;replot x and y values in red with slightly smaller symbols (yields a small black outline on each symbol)
OPLOT,x,y,psym = 8,symsize=0.8,co=!red

PSCLOSE

END
