; Create and use a widget to measure EW of spectral lines
; in PHOENIX spectra and SDSS spectra
pro pho_initialize

common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp

; FIX Change this later
;restore, 'phoenix_spectra.sav'
restore, 'phoenix_spectra_debug.sav'

sdss_info = mrdfits('/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits',1)

; sdss_temp = []
; for tspec = 0, n_elements(sdss_info)-1 do begin
;     tempstr = sdss_info[tspec]
;     newstr = create_struct(tempstr, 'STATUS', 0)
;     sdss_temp = [sdss_temp, newstr]
; endfor

; sdss_info = sdss_temp
nsdss = n_elements(sdss_info) 
; Initialize 
npix_sdss = 8575L
pix_vec = dindgen(npix_sdss)
; Read headers to get wl info
head_str = readin_apo(hdu=0, nfiles=1)
head = *(head_str.header)
wl0_sdss = fxpar(head, 'CRVAL1')
wl1_sdss = fxpar(head, 'CDELT1')
log_wl_sdss = poly(pix_vec, [wl0_sdss,wl1_sdss])
wl_grid_sdss = 10d0^log_wl_sdss

; also get snr and hmag
sdss_head = mrdfits('/home/stgilhool/APOGEE/apg_head_info.fits',1)
snr = sdss_head.snr
hmag = sdss_head.h
nvisits = sdss_head.nvisits

sdss_supp = {status:replicate(1, nsdss), $
             snr:snr, $
             hmag:hmag, $
             nvisits:nvisits, $
             wl_grid_sdss:wl_grid_sdss}



nphospec = n_elements(pho_str)
pspec_status = replicate(1, nphospec)

; Read in linelist
readcol, linelist, wl_ctr, species, format = 'D,A', delimiter = "'", nlines=nlines
; truncate
line_idx = where(wl_ctr gt wl0 and wl_ctr lt wl1, nlines)
wl_ctr = wl_ctr[line_idx]
species = species[line_idx]

;fe_idx = where(species eq 'Fe 1', nfe)
;if nfe gt 0 then begin
;    wl_ctr = wl_ctr[fe_idx]
;    species = species[fe_idx]
;endif else message, "no iron lines found"

; PLOT FIRST LINE
;MAKE WL SECTION
;make feature wl section
;make window wl section
wl_ctr_0 = wl_ctr[0]
species_0 = species[0]
feature_hwhm = 0.4d0
feature_range = wl_ctr_0 + [-1d0*feature_hwhm, feature_hwhm]
feature_idx = value_locate(pho_wl_grid, feature_range)
feature_idx[0] = (feature_idx[0]-1) > 0
feature_wl = pho_wl_grid[feature_idx[0]:feature_idx[1]]


window_hw = 2d0
window_range = wl_ctr_0 + [-1d0*window_hw, window_hw]
window_idx = value_locate(pho_wl_grid, window_range)
window_idx[0] = (window_idx[0]-1) > 0
window_wl = pho_wl_grid[window_idx[0]:window_idx[1]]

continuum_left = pho_wl_grid[feature_idx[0]-3:feature_idx[0]-1]
continuum_right = pho_wl_grid[feature_idx[-1]+1:feature_idx[-1]+3]

pho_yr = [0.94, 1.02]
pho_xr = [window_wl[0], window_wl[1]]

teff_yr = [1900, 4300]
feh_yr = [-4.5, 1.0]
logg_yr = [4.0, 5.5]

;;; Calculate EW in the PHOENIX spectra by summing over the feature range
pho_depth = d_wl*(1d0-pho_arr)
pho_ew = total(pho_depth[feature_idx[0]:feature_idx[1],*], 1, /double)
;;; Calculate Pearson correlations
pho_ew = reform(pho_ew)
teff_corr = find_correlation(pho_ew, pho_teff)
logg_corr = find_correlation(pho_ew, pho_logg)
feh_corr  = find_correlation(pho_ew, pho_feh)

;Save
pho_spec_info = {phoenix_spectra:pho_str, $
                 d_wl:d_wl, $
                 pspec_status:pspec_status, $
                 pteff_status:pspec_status, $
                 pfeh_status:pspec_status, $
                 plogg_status:pspec_status, $
                 pho_wl_grid:pho_wl_grid, $
                 feature_idx:feature_idx, $
                 window_idx:window_idx, $
                 feature_hwhm:feature_hwhm, $
                 continuum_left:continuum_left, $
                 continuum_right:continuum_right, $
                 pho_yr:pho_yr, $
                 pho_xr:pho_xr, $
                 pho_ew:pho_ew, $
                 pho_teff:pho_teff, $
                 pho_feh:pho_feh, $
                 pho_logg:pho_logg, $
                 teff_corr:teff_corr, $
                 feh_corr:feh_corr, $
                 logg_corr:logg_corr, $
                 teff_yr:teff_yr, $
                 feh_yr:feh_yr, $
                 logg_yr:logg_yr}

line_info = {wl_ctr:wl_ctr, $
             species:species, $
             wl_ctr_updated:wl_ctr, $
             feature_widths:replicate(pho_spec_info.feature_hwhm, nlines), $
             line_status:replicate(0, nlines), $
             wl_ctr_current:wl_ctr_0, $
             species_current:species_0, $
             line_current:0, $
             nlines:nlines}


end

pro update_pspec_idx
common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp



end

pro line_switch
common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp

; change line info fields 
line_current = line_info.line_current
vald_wl = line_info.wl_ctr[line_current]
wl_ctr_current = line_info.wl_ctr_updated[line_info.line_current]
line_info.wl_ctr_current = wl_ctr_current
species_current = line_info.species[line_info.line_current]
line_info.species_current = species_current
pho_wl_grid = pho_spec_info.pho_wl_grid

; change feature idx and window idx

;feature_hwhm = pho_spec_info.feature_hwhm
feature_hwhm = line_info.feature_widths[line_info.line_current]
feature_range = wl_ctr_current[0] + [-1d0*feature_hwhm, feature_hwhm]
feature_idx = value_locate(pho_wl_grid, feature_range)
feature_idx[0] = (feature_idx[0]-1) > 0
feature_wl = pho_wl_grid[feature_idx[0]:feature_idx[1]]

window_hw = 2d0
window_range = wl_ctr_current[0] + [-1d0*window_hw, window_hw]
window_idx = value_locate(pho_wl_grid, window_range)
window_idx[0] = (window_idx[0]-1) > 0
window_wl = pho_wl_grid[window_idx[0]:window_idx[1]]

continuum_left = pho_wl_grid[feature_idx[0]-3:feature_idx[0]-1]
continuum_right = pho_wl_grid[feature_idx[-1]+1:feature_idx[-1]+3]

pho_yr = [0.94, 1.02]
pho_xr = [window_wl[0], window_wl[1]]

; save
pho_spec_info.window_idx = window_idx
pho_spec_info.feature_idx = feature_idx
pho_spec_info.continuum_left = continuum_left
pho_spec_info.continuum_right = continuum_right

;;; Calculate EW in the PHOENIX spectra by summing over the feature
;;; range
pspec_status = pho_spec_info.pspec_status
pspec_idx = where(pspec_status eq 1, npspec)

;pho_arr =
;pho_spec_info.phoenix_spectra[pspec_idx].spec;[window_idx[0]:window_idx[1]]
pho_arr = pho_spec_info.phoenix_spectra.spec;[window_idx[0]:window_idx[1]]

pho_depth = pho_spec_info.d_wl*(1d0-pho_arr)

pho_ew = total(pho_depth[feature_idx[0]:feature_idx[1],*], 1, /double)

;;; Calculate Pearson correlations
pho_ew = reform(pho_ew)

teff_corr = find_correlation(pho_ew, pho_spec_info.pho_teff)
logg_corr = find_correlation(pho_ew, pho_spec_info.pho_logg)
feh_corr  = find_correlation(pho_ew, pho_spec_info.pho_feh)

; save
pho_spec_info.pho_ew = pho_ew
pho_spec_info.teff_corr = teff_corr
pho_spec_info.feh_corr = feh_corr
pho_spec_info.logg_corr = logg_corr


widget_control, status.feature_width, set_value=feature_hwhm
          
widget_control, status.central_wl, set_value=wl_ctr_current
widget_control, status.vald_wl, set_value='VALD WL: '+strtrim(vald_wl,2)
widget_control, status.species_txt, set_value='Species: '+line_info.species_current
widget_control, status.line_sel, /update
          

end

pro ew_plot
common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp



; what is needed?
; spectral vectors
; wl vector (window)
; wl vector (feature)
; mask vector (sdss)
; mask on/off
; nmasked (sdss)
; err vector (sdss)
; err on/off
; snr (sdss)
; continuum region selections
; ew width calculations
; correlation coefficients with relevant parameters
; plot yr
; plot xr
  ; pho
  ; sdss
  ; corr teff
  ; corr logg
  ; corr feh
; feature info
  ; species
  ; wl ctr
  ; median gaussian parameters?
  ; median ew (numerical integration)?
  ; IDEA: plot distribution of gaussian parameters?




; PLOT PHOENIX SPECTRUM
; get winid
wset, status.pho_draw_id
; get vectors
pho_wl_grid = pho_spec_info.pho_wl_grid
window_idx = pho_spec_info.window_idx
feature_idx = pho_spec_info.feature_idx
window_wl = pho_wl_grid[window_idx[0]:window_idx[1]]
feature_wl = pho_wl_grid[feature_idx[0]:feature_idx[1]]
pho_yr = pho_spec_info.pho_yr
pspec_status = pho_spec_info.pspec_status
pspec_idx = where(pspec_status eq 1, npspec)
if npspec lt 1 then message, "Must select at least one spectrum"
pho_spec = pho_spec_info.phoenix_spectra[pspec_idx].spec[window_idx[0]:window_idx[1]]

;plot, rebin(window_wl, n_elements(window_wl), npspec), pho_spec, yr =
;pho_yr, ps=3, /xs
plot, window_wl, pho_spec[*,0], yr = pho_yr, /xs
for i = 1, npspec-1 do oplot, window_wl, pho_spec[*,i]
vline, feature_wl[0]
vline, feature_wl[-1]

;SDSS
;Choose only plx and irtf_feh ones
sdss_sel_idx = where(finite(sdss_info.feh_irtf) and finite(sdss_info.plx_dtm),nsel)
if nsel gt 0 then begin
    sdss_supp.status = 0
    sdss_supp.status[sdss_sel_idx]=1
endif
;print, nsel


sdss_spec = sdss_info[sdss_sel_idx].spec
wl_sdss = sdss_supp.wl_grid_sdss
sdss_phospec = []
for i = 0, nsel-1 do begin
    sdss_phogrid = interpol(sdss_spec[*,i], wl_sdss, pho_wl_grid)
    sdss_phospec = [[sdss_phospec], [sdss_phogrid]]
endfor

sdss_drawspec = sdss_phospec[window_idx[0]:window_idx[1],*]
; get winid
wset, status.sdss_draw_id
plot, rebin(window_wl, n_elements(window_wl), nsel), sdss_drawspec, yr=[0.7,1.1], ps=3, /xs

snr_hist = histogram(sdss_supp.snr[sdss_sel_idx], binsize = 100, omin =om)
hist_x = lindgen(n_elements(snr_hist))*100 + om
hmag_hist = histogram(sdss_supp.hmag[sdss_sel_idx], binsize = 0.5, omin = om)
hhist_x = lindgen(n_elements(hmag_hist))*0.5 + om 
plot, hhist_x, hmag_hist, ps=10

;print, [reform(sdss_supp.snr[sdss_sel_idx],1,nsel),$
crit_arr= [reform(sdss_supp.snr[sdss_sel_idx],1,nsel),$
        reform(sdss_supp.hmag[sdss_sel_idx],1,nsel),$
        reform(sdss_supp.nvisits[sdss_sel_idx],1,nsel)]

crit_idx = where(crit_arr[0,*] ge 400 and crit_arr[2,*] gt 3, ncrit)
;print, ncrit
;print, sdss_sel_idx[crit_idx]
;print, 'feh range = ', minmax(sdss_info[sdss_sel_idx[crit_idx]].feh_irtf)
;help, sdss_info, /str
;stop

;plot, rebin(window_wl, n_elements(window_wl), ncrit), sdss_drawspec[*,crit_idx], yr=[0.7,1.1], /xs;, ps=3
plot, window_wl, sdss_drawspec[*,crit_idx[0]], yr=[0.6,1.1], /xs
for i = 1, ncrit-1 do begin
    oplot, window_wl, sdss_drawspec[*,crit_idx[i]]
endfor
vline, feature_wl[0]
vline, feature_wl[-1]
;for i = 0, nsel-1 do begin
;for w = 0, ncrit-1 do begin
;i = crit_idx[w]
;plot, wl_sdss, sdss_spec[*,i], title = 'SNR: '+strtrim(sdss_supp.snr[sdss_sel_idx[i]],2)+' | Nvisits: '+strtrim(sdss_supp.nvisits[sdss_sel_idx[i]],2)+' | HMag: '+strtrim(sdss_supp.hmag[sdss_sel_idx[i]],2)
;dxstop
;endfor


;ALSO plot neighboring lines




; PLOT CORRELATIONS
pho_ew = pho_spec_info.pho_ew[pspec_idx]
teff_corr = pho_spec_info.teff_corr
pho_teff = pho_spec_info.pho_teff[pspec_idx]
feh_corr = pho_spec_info.feh_corr
pho_feh = pho_spec_info.pho_feh[pspec_idx]
logg_corr = pho_spec_info.logg_corr
pho_logg = pho_spec_info.pho_logg[pspec_idx]
; get winid
wset, status.teff_draw_id
; plot teff_corr
plot, pho_ew, pho_teff, ps=3, title='Teff Correlation: '+strtrim(teff_corr,2)

; get winid
wset, status.feh_draw_id
; plot feh_corr
plot, pho_ew, pho_feh, ps=3, title='[Fe/H] Correlation: '+strtrim(feh_corr,2)

; get winid
wset, status.logg_draw_id
; plot logg_corr
plot, pho_ew, pho_logg, ps=3, title='logg Correlation: '+strtrim(logg_corr,2)


;stop
end


PRO pho_widget_event, ev

common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp

  widget_control, ev.id, get_uvalue = uvalue
  
  ;nteff = n_elements(stellar_params.teff_set) 
  ;nfeh = n_elements(stellar_params.feh_set) 
  ;nlogg = n_elements(stellar_params.logg_set) 

  pspec_status = pho_spec_info.pspec_status
  pteff_status = pho_spec_info.pteff_status
  pfeh_status = pho_spec_info.pfeh_status
  plogg_status = pho_spec_info.plogg_status
  pho_teff_long = pho_spec_info.phoenix_spectra.teff
  pho_feh_long = pho_spec_info.phoenix_spectra.feh
  pho_logg_long = pho_spec_info.phoenix_spectra.logg

  case uvalue of
      
      'teff_sel': begin
          stellar_params.teff_set[ev.value] = ev.select
          teff_val = stellar_params.teff_grid[ev.value]
          
          change_idx = where(pho_teff_long eq teff_val, nchange)
          
          if nchange gt 0 then begin
              pteff_status[change_idx] = ev.select
              pho_spec_info.pteff_status = pteff_status
              pho_spec_info.pspec_status = $
                pteff_status and pfeh_status and plogg_status
          endif else message, "Problem with teff selection"
          
          ew_plot
      end
      
      'feh_sel': begin
          stellar_params.feh_set[ev.value] = ev.select
          
          feh_val = stellar_params.feh_grid[ev.value]
          
          change_idx = where(float(pho_feh_long) eq feh_val, nchange)
          if nchange gt 0 then begin
              pfeh_status[change_idx] = ev.select
              pho_spec_info.pfeh_status = pfeh_status
              pho_spec_info.pspec_status = $
                pteff_status and pfeh_status and plogg_status
          endif else message, "Problem with feh selection"
          
          ew_plot
      end
      
      'logg_sel': begin
          stellar_params.logg_set[ev.value] = ev.select
          logg_val = stellar_params.logg_grid[ev.value]
          
          change_idx = where(pho_logg_long eq logg_val, nchange)
          if nchange gt 0 then begin
              plogg_status[change_idx] = ev.select
              pho_spec_info.plogg_status = plogg_status
              pho_spec_info.pspec_status = $
                pteff_status and pfeh_status and plogg_status
          endif else message, "Problem with logg selection"
          
          ew_plot
      end
      
      'line_sel': begin
          line_num = ev.value
          line_info.line_current = line_num
          line_switch
          
          ew_plot
      end
      
      'first_btn': begin
          current_line = line_info.line_current
          current_line = 0
          line_info.line_current = current_line
          widget_control, status.line_sel, set_value=current_line
          widget_control, status.line_sel, send_event={id:0l,top:0l,handler:0l,value:current_line}
      end
      
      'prev_btn': begin
          current_line = line_info.line_current
          current_line = (current_line-1) > 0
          line_info.line_current = current_line
          widget_control, status.line_sel, set_value=current_line
          widget_control, status.line_sel, send_event={id:0l,top:0l,handler:0l,value:current_line}
      end
      'next_btn': begin
          current_line = line_info.line_current
          current_line = (current_line+1) < (line_info.nlines-1)
          line_info.line_current = current_line
          widget_control, status.line_sel, set_value=current_line
          widget_control, status.line_sel, send_event={id:0l,top:0l,handler:0l,value:current_line}
          
      end
      'last_btn': begin
          current_line = line_info.line_current
          current_line = line_info.nlines-1
          line_info.line_current = current_line
          widget_control, status.line_sel, set_value=current_line
          widget_control, status.line_sel, send_event={id:0l,top:0l,handler:0l,value:current_line}
      end

      'feature_width': begin
          feature_hwhm = ev.value
          ;pho_spec_info.feature_hwhm = feature_hwhm
          line_info.feature_widths[line_info.line_current] = feature_hwhm
          wl_ctr_current = line_info.wl_ctr_current
          pho_wl_grid = pho_spec_info.pho_wl_grid
          
          feature_range = wl_ctr_current[0] + [-1d0*feature_hwhm, feature_hwhm]
          feature_idx = value_locate(pho_wl_grid, feature_range)
          feature_idx[0] = (feature_idx[0]-1) > 0
          pho_spec_info.feature_idx = feature_idx
          line_switch
          ew_plot
      end
      
      'central_wl': begin
          wl_ctr_updated = ev.value
          line_info.wl_ctr_updated[line_info.line_current] = wl_ctr_updated
          
          line_switch
          ew_plot
      end 

  endcase
  
END


 
PRO pho_ew_widget

common plot_info, status, stellar_params, pho_spec_info, line_info, sdss_info, sdss_supp

; Initialize
pho_initialize

;;; Create widget and populate its fields
base = WIDGET_BASE(/COLUMN)
spec_base = widget_base(base, /row)
corr_base = widget_base(base, /row)
corr_btn_base = widget_base(base, /row)


;Teff
nteff = 23
teff_grid = lindgen(nteff)*100 + 2000L
teff_names = strtrim(teff_grid,2)
teff_set = replicate(1, nteff)
;FeH
feh_names = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', $
             '-1.0', '-0.5', '-0.0', '+0.3', '+0.5'] 
nfeh = n_elements(feh_names) 

feh_vals = float(feh_names)
feh_vals[where(feh_vals eq 0)] = 0.0

feh_grid = feh_vals
feh_set = replicate(1,nfeh)
;logg
logg_names = ['4.5', '5.0']
nlogg = n_elements(logg_names) 
logg_set = [1,1]


stellar_params = {teff_set:teff_set, $
                  teff_grid:teff_grid, $
                  feh_set:feh_set, $
                  feh_grid:feh_names, $
                  logg_set:logg_set, $
                  logg_grid:logg_names}

teff_sel = cw_bgroup(spec_base, teff_names, label_top='Teff', /frame,$
                     uvalue='teff_sel', column=2, $
                     /nonexclusive, set_value=teff_set)
feh_sel = cw_bgroup(spec_base, feh_names, label_top='Fe/H', /frame, $
                    uvalue='feh_sel', /column, $
                    /nonexclusive, set_value=feh_set)
logg_sel = cw_bgroup(spec_base, logg_names, label_top='logg', /frame, $
                     uvalue='logg_sel', /column, $
                     /nonexclusive, set_value=logg_set)

winxs = 500
winys = 400

pho_draw = WIDGET_draw(spec_base, uvalue='pspec_vis', xs = winxs, ys = winys)
sdss_draw = WIDGET_draw(spec_base, uvalue='sdss_vis', xs = winxs, ys = winys)

field_base = widget_base(spec_base, /column)

teff_draw = WIDGET_draw(corr_base, uvalue='teff_vis', xs = winxs, ys = winys)
feh_draw = WIDGET_draw(corr_base, uvalue='feh_vis', xs = winxs, ys = winys)
logg_draw = WIDGET_draw(corr_base, uvalue='logg_vis', xs = winxs, ys = winys)



WIDGET_CONTROL, base, /REALIZE
widget_control, pho_draw, get_value=pho_draw_id 
widget_control, sdss_draw, get_value=sdss_draw_id 
widget_control, teff_draw, get_value=teff_draw_id 
widget_control, feh_draw, get_value=feh_draw_id 
widget_control, logg_draw, get_value=logg_draw_id 


first_btn = widget_button(corr_btn_base, value="First", uvalue='first_btn')
prev_btn = widget_button(corr_btn_base, value="Previous", uvalue='prev_btn')

line_sel = widget_slider(corr_btn_base, minimum = 0, maximum = line_info.nlines-1, $
                         title='Line Number', uvalue='line_sel', /frame, $
                         value = 0)

next_btn = widget_button(corr_btn_base, value="Next", uvalue='next_btn')
last_btn = widget_button(corr_btn_base, value="Last", uvalue='last_btn')


line_current = line_info.line_current

feature_width = cw_field(field_base, uvalue='feature_width', title='Feature width',$
                            value = pho_spec_info.feature_hwhm, $
                            /return_events, /floating)



central_wl = cw_field(field_base, uvalue='central_wl', title='Central WL',$
                            value = line_info.wl_ctr_current[line_current], $
                            /return_events, /floating)

vald_wl = widget_text(field_base, uvalue='vald_wl', $
                      value='VALD WL: '+strtrim(line_info.wl_ctr[line_current],2))

species_txt = widget_text(field_base, uvalue='species_txt', $
                      value='Species: '+line_info.species[line_current])

status = {pho_draw_id:pho_draw_id, $
          sdss_draw_id:sdss_draw_id, $
          teff_draw_id:teff_draw_id, $
          feh_draw_id:feh_draw_id, $
          logg_draw_id:logg_draw_id, $
          line_sel:line_sel, $
          feature_width:feature_width, $
          central_wl:central_wl, $
          vald_wl:vald_wl, $
          species_txt:species_txt}


ew_plot

XMANAGER, 'pho_widget', base


stop

END
