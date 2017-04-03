;+
; NAME:
;
;
;
; PURPOSE:
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
pro vline, val,_extra=extra, min=min, max=max
if !y.type eq 1 then yrange = 10^!y.crange else yrange = !y.crange
nv = n_elements(val)
for i = 0, nv-1 do oplot,fltarr(2)+val[i],yrange,_extra=extra
end

function feh_fit, p, ew=ew, feh=feh, err=err, vis=vis

common outputs, model, res

feh_col = reform(feh, 1, n_elements(feh))

model = ew ## p

model = reform(model)

data = feh

res = data - model

dev = res/err

chi2 = sqrt(total(dev^2, /double))

sorti = sort(data)
sortd = data[sorti]
sortm = model[sorti]
sortr = res[sorti]

if vis eq 1 then begin
    ;plot, data, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chi2,2), /xs
    ;oplot, model, ps = 6, color = 200
    ;plot, res, ps = 6, title = "Residuals", /xs
    plot, sortd, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chi2,2), /xs
    oplot, sortm, ps = 6, color = 200
    plot, sortr, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(res),2), /xs
endif

return, dev

end




pro apo_feh_redo

; Constants
bigc = 299792.458 ;km/s

; Options
vis = 1

; Paths
root_dir = '/home/stgilhool/APOGEE/'
data_path = root_dir + 'APOGEE_data/'

listfile = data_path + 'mdwarflist.txt'

; Readin irtf data
irtf = read_irtf_data()

irtf_id = strtrim(irtf.APG_APOGEE_ID,2)
feh_afk = irtf.AFK
feh_afke= irtf.AFKE

; Constants
nx = 8575L
n_apg = 1350
n_irtf_all = n_elements(irtf)

; Readin APG info
;;;;
goto, loadskip

apg = readin_apo(hdu = 1, nfiles = n_apg)

save, /variables, filename='readinapotemp.sav'

loadskip:
restore, 'readinapotemp.sav'
;;;;

spectra = apg.spectra

;window, 0, xsize = 1400, ysize=1000

; Loop through and save IRTF matches
matchidx = []
apg_id = []
for i = 0, n_apg-1 do begin

    ; Find the matching IRTF file
    head = *(apg[i].header)
    apg_id_i = fxpar(head, 'OBJID')
    apg_id_i = strtrim(apg_id_i,2)
    apg_id = [apg_id, apg_id_i]

    irtfmatch_idx = where(irtf_id eq apg_id_i, matchcnt)
    if matchcnt eq 1 then begin
        
        matchidx = [[matchidx], [i,irtfmatch_idx[0]]]
        
    endif else matchidx = [[matchidx], [i,!values.d_nan]]
    
endfor


; for ir = 0, n_irtf-1 do begin
;     irtfmatch_idx = where(apg_id eq irtf_id[ir], matchcnt)

;     if matchcnt eq 1 then begin
        
;         matchidx = [matchidx, irtfmatch_idx]

;     endif else matchidx = [matchidx, !values.l_nan]

; endfor

mwrfits, matchidx, '/home/stgilhool/APOGEE/apg_irtf_match_idx.fits', /create 

apgidx = matchidx[0,*]
irtfidx = matchidx[1,*]

bothidx = where(finite(irtfidx),fincnt)
n_irtf = fincnt


apgfinidx = apgidx[bothidx]
irtffinidx = irtfidx[bothidx]

irtfsort = sort(irtffinidx)
irtf_idx = irtffinidx[irtfsort]
apg_idx = apgfinidx[irtfsort]
objid = apg_id[bothidx]
apgid = objid[irtfsort]

irtf_id_test = irtf_id[irtf_idx]
feh_afk = feh_afk[irtf_idx]
feh_afke = feh_afke[irtf_idx]


; Sum up chunks

npix_win = 25L

depth = 1d0-spectra


; This gives average value in window centered at each pixel
tophat = replicate(1d0, npix_win)
ew_arr_all = convol(depth, tophat)
ew_arr = ew_arr_all[*,apg_idx]
    
; Fit the metalicity of each chunk
; Initialize

; Make median bspec for display
spec_med = median(spectra, dimension=2)


save, /variables, filename='correlation_redo.sav'
test:
restore, 'correlation_redo.sav'


!p.multi=0
if vis eq 1 then begin
    window, 0, xsize=1300, ysize=1000
    ;!p.multi=[0,1,2]
endif

; Line centers
corr_x0b = [720L, 2000L] + 403L
    
corr_x0g = [300, 400, 510, 1090, 1170, 1290, 1350, 1390, 1580, $
            1610, 1790, 1810, 2080, 2110, 2130, 2200]+3740L

corr_x0r = [12, 100, 120, 150, 890, 910, 950, 1050, 1080]+6496L

ew_corr_x0 = [corr_x0b, corr_x0g, corr_x0r]
ew_corr_y0 = ew_arr[ew_corr_x0,*]
  
n_ewc = n_elements(ew_corr_x0)

ew_corr_fp = ew_corr_x0 - (npix_win/2)
ew_corr_lp = ew_corr_x0 + (npix_win/2)

; Add column for constant term
ew_corr_y = [reform(replicate(1d0,n_irtf),1,n_irtf),ew_corr_y0]
help, ew_corr_y

!p.multi=0
if vis eq 1 then begin
    window, 0, xsize=1300, ysize=1000
    !p.multi=[0,1,2]
endif
    
;;; Fit the metalicity of each chunk

common outputs, model, res

functargs = {ew:ew_corr_y, $
             feh:feh_afk, $
             err:feh_afke, $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_ewc+1)


r = mpfit('feh_fit', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2)


print, "result:"
print, r

outlier = where(abs(res) ge 0.2, ocount)
print, ocount
print, outlier

ew_out = ew_arr_all[ew_corr_x0, *]
ew_out = [reform(replicate(1d0,n_apg),1,n_apg),ew_out]

feh_out = ew_out ## r
feh_out = reform(feh_out)

mwrfits, feh_out, '/home/stgilhool/APOGEE/feh_results/feh_results_redo.fits', /create


outstr = {result:r, $
          model:model, $
          data:feh_afk, $
          res:res, $
          err:feh_afke, $
          nlines:n_ewc, $
          x0:ew_corr_x0, $
          fp:ew_corr_fp, $
          lp:ew_corr_lp, $
          width:npix_win, $
          status:status, $
          irtf_idx:matchidx, $
          chi2:chi2, $
          chi2dof:chi2/dof $
          }
         
          
outpath = '/home/stgilhool/APOGEE/feh_results/'
outfile = outpath + 'feh_fit_redo.fits'
mwrfits, outstr, outfile, /create


;reg = regress(ew_corr_y[1:*,*], feh_afk, measure_errors=feh_afke, sigma=parsigma, const=const, chisq=chisq, /double, correlation=correlation, mcorrelation=mcorrelation, ftest=ftest, status=regstatus, yfit=regfit)

;window, 2, xsize=1300, ysize=1000
;!p.multi = [0,1,2]

;sortfeh = sort(feh_afk)
;sfeh = feh_afk[sortfeh]
;sregfit = regfit[sortfeh]
;sres = sfeh-sregfit

;plot, sfeh, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chisq,2), /xs
;oplot, sregfit, ps = 6, color = 200
;plot, sres, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(sres),2), /xs


;print, "IDL REGRESS"
;print, "Result: ", [const, reform(reg)]
;print, "Correlation: "
;print, correlation

;print, v_src

stop


end
