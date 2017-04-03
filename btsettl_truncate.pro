pro btsettl_truncate

ppath = '/home/stgilhool/PHOENIX/'

pfile = ppath + 'lte0'+strtrim(teff,2)+'-4.5'+feh+'.BT-Settl.7'

if file_test(pfile) eq 0 then begin
    message, pfile+' not found'
endif
df = -8d0
; Read in PHOENIX model
readcol, pfile, pwl, c1, c2, c3, c4, f1, f2, f3, f4, f5, f6, c5, c6, c7, c8, format = "D,D,D,D,F,F,F,F,F,F,D,D,D,D"
pws = sort(pwl)
pwl = pwl[pws]
c1 = c1[pws]
c2 = c2[pws]
c3 = c3[pws]
c4 = c4[pws]
f1 = f1[pws]
f2 = f2[pws]
f3 = f3[pws]
f4 = f4[pws]
f5 = f5[pws]
f6 = f6[pws]
c5 = c5[pws]
c6 = c6[pws]
c7 = c7[pws]
c8 = c8[pws]
window, 1
;plot, pwl, 10^(c1+df)
plot, pwl, c1, yr=[0,30]
window, 2
;plot, pwl, 10^(c2+df)
plot, pwl, c2, yr=[0,30]
window, 3
plot, pwl, c3,yr=[0,30]
window, 4
plot, pwl, c4, yr=[0,30]

stop


end
