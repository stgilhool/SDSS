pro nlines_check

; Look at target directory
path = '/home/stgilhool/APOGEE/results/'

filenames = file_search(path+'*', count=nfiles)

r = ptrarr(nfiles, /allocate_heap)


rmsres = dblarr(nfiles)
thresh = dblarr(nfiles)
nlines = dblarr(nfiles)
chi2dof= dblarr(nfiles)

; Read in data
foreach file, filenames, fidx do begin
    
    ; Read in file
    instr = mrdfits(file, 1)

    ; Make array of structures from input files
    *(r[fidx]) = instr
    
    rmsres[fidx] = stddev(instr.res)
    thresh[fidx] = instr.pthresh
    nlines[fidx] = instr.nlines
    chi2dof[fidx]= instr.chi2dof

endforeach

; Look at residuals as a function of correlation threshhold

plot, thresh, chi2dof, ps=6

stop

end

    
    
