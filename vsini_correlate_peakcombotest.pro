
; Testing to see if I can combine the fft's to get at the real peak
function fft_log, yvec

fft_vec = fft(yvec, /double)

fft_amp = sqrt(real_part(fft_vec)^2 + imaginary(fft_vec)^2)

log_fft = alog10(fft_amp)

return, log_fft

end

function zero_wings, yvec, min_idx

; get lower and upper idx
lmin = min(min_idx, max=umin)

rvec = yvec

; zero entries below lmin
rvec[0:lmin-5] = 0d0
; zero entries about umin
rvec[umin+5:*] = 0d0

; raise kernel to have min at zero
minyvec = min(rvec)

rvec[lmin-4:umin+4] = rvec[lmin-4:umin+4] - minyvec

; Smooth
rvec_cp = rvec

rvec = gauss_smooth(rvec_cp, 11)

return, rvec

end

function get_min_idx, xvec, yvec

ysmooth = gauss_smooth(yvec, 9)
; !p.multi = [0,1,2]
; plot, xvec, yvec
; plot, xvec, ysmooth
; !p.multi=0
; dxstop

min_idx = extremes(ysmooth, -1)
min_idx_diff = min_idx - (n_elements(xvec)/2)
min_idx_abs = abs(min_idx_diff)
min_idx_sort = sort(min_idx_abs)
min_idx_idx = min_idx[min_idx_sort]
min_idx_diff_sort = min_idx_diff[min_idx_sort]
;min_idx_fin = min_idx_idx[0:1]
min_idx_neg = where(min_idx_diff_sort lt 0, nneg)
min_idx_pos = where(min_idx_diff_sort gt 0, npos)

min_idx_fin = [min_idx_idx[min_idx_neg[0]], min_idx_idx[min_idx_pos[0]]]


return, min_idx_fin

end


function vct_find_zeroes_2, xvec, yvec


show=0

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

vsini_sigvec = k1/xvec
; define number of pixels to search
null = where(vsini_sigvec gt 3, npix_search)
fft_domain = where(vsini_sigvec le 100 and vsini_sigvec gt 3, npix_domain)
first_pix = min(fft_domain)
last_pix = max(fft_domain)

; Slide a window along and find local minima
ysmooth = gauss_smooth(yvec, 1.5)

window_width = 21L
window_idx_x = lindgen(window_width)-(window_width/2)

minidx_vec = dblarr(npix_search)
max_vec = dblarr(npix_search)

for i = first_pix, last_pix do begin

    window_idx = window_idx_x + i

    ; Get the minima
    minval = min(ysmooth[window_idx], minidx)
    min_idx = window_idx[minidx]
    ; Take the max on either side of the window and average
    maxval_minus = max(ysmooth[window_idx[0:(window_width/2)]])
    maxval_plus = max(ysmooth[window_idx[(window_width/2)+1:*]])
    maxval_avg = mean([maxval_minus,maxval_plus])
    
    minidx_vec[i] = min_idx
    max_vec[i] = maxval_avg

endfor

hist = histogram(minidx_vec, min=0, reverse_indices=ri)
minima_idx_all = where(hist gt (window_width/2)+1, nminima_all)
minima_idx_idx = where(minima_idx_all gt 10, nminima)
minima_idx = minima_idx_all[minima_idx_idx]

;dxstop, 'minima_idx_all', 'minima_idx'


if nminima gt 0 then begin
   
    minima = ysmooth[minima_idx]
    maxima = max_vec[minima_idx]
    
    ; Calculate the "depth"
    peak_depth = (maxima-minima)/abs(maxima)

    ;Take only peaks that are 10%
    if show then begin
        vs_tickvec = [60,30,25,20,15,10,8,5]	
        vs_tickval = (k1/vs_tickvec)*1e3

        plot, xvec*1d3, ysmooth+10, xr=[0,150], ps=-8, xs=8
        
        axis, xaxis=1, xticks=n_elements(vs_tickvec), $
          xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
        
        oplot, xvec[0:n_elements(hist)-1]*1d3, hist, ps=10, $
          thick=2, co=!red
        for w = 0, (5<(n_elements(minima_idx)-1)) do begin
            oplot, replicate(xvec[minima_idx[w]],2)*1d3, [-100,100], thick=1, co=!blue
            print, peak_depth[w]
            dxstop
        endfor
        
    endif
    real_minima_idx = where(peak_depth gt 0.06, nreal)
    if nreal gt 0 then begin
        sigpeak_vec = xvec[minima_idx[real_minima_idx]]
        first_peak_idx = minima_idx[real_minima_idx[0]]
        sigpeak = xvec[first_peak_idx]
        vsini_vec = k1/sigpeak_vec
        vsini = k1/sigpeak
        ;print, vsini
        ;dxstop

        if show then begin
            oplot, replicate(sigpeak,2)*1d3, [-100,100], thick=3, co=!orange
            xyouts, sigpeak*1d3 + 1, 8, strtrim(vsini,2)
            dxstop
        endif
    endif else begin
        sigpeak_vec = [-99999]
        sigpeak = -99999
        vsini = 0d0
    endelse
    
endif

sigpeak_str = {minima_idx:minima_idx,$
               minima_depth:peak_depth ,$
               minima_sigval:xvec[minima_idx], $
               minima_vsini:k1/xvec[minima_idx], $
               sigvec:xvec[fft_domain],$
               fftvec:yvec[fft_domain]}

return, sigpeak_str

end


pro vsini_correlate_peakcombotest

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

;files = [840,431,793,326,953,391,286,1076,1306,289,740,1142,60,43,622,479,523,867,520]

;files = lindgen(50)

;files = [431, 793,326,397,289,1142,60,622]
;files = [431, 793, 289]
;files = [289]
;files=[793]
mi = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)

vsini_vec_all = mi.vsini_vfit

fil_idx = where(vsini_vec_all gt 10, nfilidx)

files = fil_idx
;files = fil_idx[[0,1,3,4,5,6,7,8,9,10]]

teff_vec = mi[files].teff_vfit
feh_vec = mi[files].feh_vfit
vsini_known_vec = mi[files].vsini_vfit


; Run the correlation procedure
vsini_correlate_wrapper, 350, 1400, b1, files=files
vsini_correlate_wrapper, 1400,3000, b2, files=files
;vsini_correlate_wrapper, 3700,5800, g1, files=files
vsini_correlate_wrapper, 4186,5532, g1, files=files
;vsini_correlate_wrapper, 6500,8100, r1, files=files
vsini_correlate_wrapper, 7020,7886, r1, files=files

nfiles = n_elements(files) 
result_vec = []

for fil = 0, nfiles-1 do begin

    bstr1 = *b1[fil]
    bstr2 = *b2[fil]
    gstr1 = *g1[fil]
    rstr1 = *r1[fil]

    ; make common velgrid
    vx = gstr1.vel_x_over

    ; make common sigvec
    svec = gstr1.sigvec

    ; interpolate all onto common grids
    ; FFTs
    b1fft = interpol(bstr1.fftvec, bstr1.sigvec, svec)
    b2fft = interpol(bstr2.fftvec, bstr2.sigvec, svec)
    g1fft = gstr1.fftvec
    r1fft = interpol(rstr1.fftvec, rstr1.sigvec, svec)
    ; Kernels
    b1gout = interpol(bstr1.kernel, bstr1.vel_x_over, vx)
    b2gout = interpol(bstr2.kernel, bstr2.vel_x_over, vx)
    g1gout = gstr1.kernel
    r1gout = interpol(rstr1.kernel, rstr1.vel_x_over, vx)

    ; ID the local minima around the kernel
    b1mins = get_min_idx(vx, b1gout)
    b2mins = get_min_idx(vx, b2gout)
    g1mins = get_min_idx(vx, g1gout)
    r1mins = get_min_idx(vx, r1gout)
    
    b1fix = zero_wings(b1gout, b1mins)
    b2fix = zero_wings(b2gout, b2mins)
    g1fix = zero_wings(g1gout, g1mins)
    r1fix = zero_wings(r1gout, r1mins)

    ; Scale kernels
    max_kern = max(b1fix)>max(b2fix)>max(g1fix)>max(r1fix)
    b1gscale = b1fix/max(b1fix) * max_kern
    b2gscale = b2fix/max(b2fix) * max_kern
    g1gscale = g1fix/max(g1fix) * max_kern
    r1gscale = r1fix/max(r1fix) * max_kern

    
    ; Take means and medians
    fftarr = [[b1fft],[b2fft],[g1fft],[r1fft]]
    goutarr = [[b1gout],[b2gout],[g1gout],[r1gout]]
    gscalearr = [[b1gscale],[b2gscale],[g1gscale],[r1gscale]]

    meanfft = mean(fftarr, dimension=2)
    medianfft = median(fftarr, dimension=2)

    meangout = mean(goutarr, dimension=2)
    mediangout = median(goutarr, dimension=2)

    meangscale = mean(gscalearr, dim=2)
    mediangscale = median(gscalearr, dim=2)


    ; NEW FFTS
    b1f = fft_log(b1gscale)
    b2f = fft_log(b2gscale)
    g1f = fft_log(g1gscale)
    r1f = fft_log(r1gscale)

    medianf = fft_log(mediangscale)
    meanf = fft_log(meangscale)
    
    ; Plot
    ; Kernels
    plot, vx, b1gscale, xr = [-100,100], thick=1
    oplot, vx, b2gscale, color=!blue, thick=1
    oplot, vx, g1gscale, color=!green, thick=1
    oplot, vx, r1gscale, color=!red, thick=1

    oplot, vx, meangscale, color=!forest, thick=3
    oplot, vx, mediangscale, color=!orange, thick=3

    ;dxstop




    ; FFTs
    plot, svec*1d3, b1f, thick=1, xs=8, xr=[0,150]
    oplot, svec*1d3, b2f, color=!blue, thick=1
    oplot, svec*1d3, g1f, color=!green, thick=1
    oplot, svec*1d3, r1f, color=!red, thick=1

    oplot, svec*1d3, meanf, color=!forest, thick=3
    oplot, svec*1d3, medianf, color=!orange, thick=3

    ; Put corresponding Vsini values on top
    epsilon = 0.6d0
    k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

    vs_tickvec = [60,30,25,20,15,10,8,5]	
    vs_tickval = (k1/vs_tickvec)*1d3
    axis, xaxis=1, xticks=n_elements(vs_tickvec), $
      xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
        

    ;dxstop

    b1v_str = vct_find_zeroes_2(svec, b1f)
    b2v_str = vct_find_zeroes_2(svec, b2f)
    g1v_str = vct_find_zeroes_2(svec, g1f)
    r1v_str = vct_find_zeroes_2(svec, r1f)
    mean_str = vct_find_zeroes_2(svec, meanf)
    median_str = vct_find_zeroes_2(svec, medianf)

    b1v = b1v_str.minima_vsini[0]
    b2v = b2v_str.minima_vsini[0]
    g1v = g1v_str.minima_vsini[0]
    r1v = r1v_str.minima_vsini[0]
    meanv = mean_str.minima_vsini[0]
    medianv = median_str.minima_vsini[0]

    vvec = [b1v, b2v, g1v, r1v, meanv, medianv]

    result_vec = [[result_vec], [vvec]]    

endfor

; Plot
one2one = lindgen(fix(max(vsini_known_vec)*1.2))
plot, one2one, one2one, linest=2, xtit='Vsini from fitting', ytit='Vsini from Cross Correlation', charsize=1.5

for fil = 0, nfiles-1 do begin & $

    vels = result_vec[*,fil] & $
    vknow = vsini_known_vec[fil] & $

    oplot, [vknow], [vels[0]], co=!blue, ps=8, symsize=0.8 & $
    oplot, [vknow], [vels[1]], co=!purple, ps=8, symsize=0.8 & $
    oplot, [vknow], [vels[2]], co=!green, ps=8, symsize=0.8 & $
    oplot, [vknow], [vels[3]], co=!red, ps=8, symsize=0.8 & $
    oplot, [vknow], [vels[4]], co=!orange, ps=8, symsize=1 & $
    oplot, [vknow], [vels[5]], co=!yellow, ps=8, symsize=1 & $

endfor
    
stop


end

