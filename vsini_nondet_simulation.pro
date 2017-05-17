; Simulate vsini's from a V_eq distribution, and get a handle of how
; many non-detections one would expect from N observations

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sini_probability, sample
; Give the probability for each element in sample

sini_prob = sample/sqrt(1d0-sample^2)

return, sini_prob

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function veq_probability, sample
; Give the probability for each element in sample
common pdf_info, veq_pdf, veq_vec

veq_prob = interpol(veq_pdf, veq_vec, sample)

return, veq_prob

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro vsini_sim_input, output
; Get the input vector/pdf
; this file has young stars only
input = mrdfits('vsini_simulation_tempfile.fits',1)

; define teff range to use
idx = where(input.teff le 2800 and input.teff ge 2600, nidx)

vector = input.veq[idx]

output = vector

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function vsini_nondet_simulate, nsamples=nsamples, niter=niter

common pdf_info, veq_pdf, veq_vec

if n_elements(nsamples) eq 0 then nsamples=100
if n_elements(niter) eq 0 then niter=1d3

non_detections = lonarr(niter)
non_det_hist = lonarr(nsamples)

for iter = 0, niter - 1 do begin
    ; Draw V_eq and sini for N stars from pdf
    ; draw sini
    sample_pop_sini = create_sample_distribution('sini_probability',[0d0,1d0],nsamples)
    ; draw V_eq
    sample_pop_veq = create_sample_distribution('veq_probability',$
                                                minmax(veq_vec),nsamples)

    ; count the number of non-detections
    sample_vsini = sample_pop_veq * sample_pop_sini
    ndet_idx = where(sample_vsini le 5d0, n_nondet)

    ; store the number of non-detections
    non_detections[iter] = n_nondet

    non_det_hist[n_nondet]++

    ;plot, non_det_hist, ps=10, title = "Number of non-detections distribution"
    ;wait, 0.05

endfor

return, non_det_hist

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro vsini_nondet_simulation

common pdf_info, veq_pdf, veq_vec

; input vector or pdf
vsini_sim_input, veq_vector
; Create veq_pdf
binsize = 0.5d0
binmin = 0d0
veq_hist = histogram(veq_vector, binsize = binsize, min=binmin)
nbins = n_elements(veq_hist) 

veq_vec = dindgen(nbins)*binsize + binmin + (binsize/2d0)
veq_pdf = veq_hist/total(veq_hist, /double)

; Settings
nsamples = 10
niter = 1d2

non_det_hist10 = vsini_nondet_simulate(niter=niter, nsamples=nsamples)

nsamples = 20

non_det_hist20 = vsini_nondet_simulate(niter=niter, nsamples=nsamples)

nsamples = 100

non_det_hist100 = vsini_nondet_simulate(niter=niter, nsamples=nsamples)

; ask: what is the probability of X non-detections, given N stars?
nd10 = non_det_hist10/total(non_det_hist10, /double)
nd20 = non_det_hist20/total(non_det_hist20, /double)
nd100 = non_det_hist100/total(non_det_hist100, /double)

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

 set_plot, 'ps'
 device, filename = "vsini_nondet_simulation1.eps"
 device, /color
device, xs=10,ys=6, /inches
; loadct,13

!p.multi = [0,1,3]

plot, nd10, ps=10, title="!4Prob. of X non-detections given 10 stars",xtit="!4# Non-detections", xr = [0,5], /xs
plot, nd20, ps=10, title="!4Prob. of X non-detections given 20 stars",xtit="!4# Non-detections", xr = [0,5], /xs
plot, nd100, ps=10, title="!4Prob. of X non-detections given 100 stars",xtit="!4# Non-detections", xr = [0,15], /xs

device, /close_file

print, "Simulation finished"

stop


end
