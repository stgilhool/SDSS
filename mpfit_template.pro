pro mpfit_template

parinfo_str = {VALUE:0,$
               FIXED:0, $
               LIMITED:[0,0], $
               LIMITS:[0,0], $
               PARNAME:'', $
               STEP:0, $
               RELSTEP:0, $
               MPSIDE:0, $
               MPDERIV_DEBUG:0, $
               MPDERIV_ABSTOL:0, $
               MPDERIV_RELTOL:0, $
               MPMAXSTEP:0, $
               TIED:0, $
               MPPRINT:0, $
               MPFORMAT:''}               



result = MPFIT(MYFUNCT, $
               start_parms, $
               FUNCTARGS=fcnargs, $
               PARINFO=parinfo, $
               MAXITER=maxiter, $
               STATUS=status, $
               BESTNORM=bestnorm, $
               NFEV=nfev, $
               NITER=niter, $
               NPRINT=nprint, $
               QUIET=quiet, $
               PERROR=perror, $
               ERRMSG=errmsg, $
               FTOL=ftol, $
               XTOL=xtol, $
               GTOL=gtol, $
               ITERPROC=iterpronc, $
               ITERARGS=iterargs, $
               COVAR=covar)
               
            
   
;    Constraining Parameter Values With The Parinfo Keyword

;   The behavior of MPFIT can be modified with respect to each 
;   parameter to be fitted. A parameter value can be fixed; simple 
;   boundary constraints can be imposed; limitations on the parameter 
;   changes can be imposed; properties of the automatic derivative can 
;   be modified; and parameters can be tied to one another. 
;   These properties are governed by the PARINFO structure, which is 
;   passed as a keyword parameter to MPFIT. 
;   PARINFO should be an array of structures, one for each parameter. 
;   Each parameter is associated with one element of the array, in 
;   numerical order. The structure can have the following entries 
;   (none are required): 
  
;     .VALUE - the starting parameter value (but see the START_PARAMS 
;               parameter for more information). 
  
;     .FIXED - a boolean value, whether the parameter is to be held 
;               fixed or not. Fixed parameters are not varied by 
;               MPFIT, but are passed on to MYFUNCT for evaluation. 
  
;     .LIMITED - a two-element boolean array. If the first/second 
;                 element is set, then the parameter is bounded on the 
;                 lower/upper side. A parameter can be bounded on both 
;                 sides. Both LIMITED and LIMITS must be given 
;                 together. 
  
;     .LIMITS - a two-element float or double array. Gives the 
;               parameter limits on the lower and upper sides, 
;               respectively. Zero, one or two of these values can be 
;               set, depending on the values of LIMITED. Both LIMITED 
;               and LIMITS must be given together. 
  
;     .PARNAME - a string, giving the name of the parameter. The 
;                 fitting code of MPFIT does not use this tag in any 
;                 way. However, the default ITERPROC will print the 
;                 parameter name if available. 
  
;     .STEP - the step size to be used in calculating the numerical 
;             derivatives. If set to zero, then the step size is 
;             computed automatically. Ignored when AUTODERIVATIVE=0. 
;             This value is superceded by the RELSTEP value. 
;     .RELSTEP - the *relative* step size to be used in calculating 
;                 the numerical derivatives. This number is the 
;                 fractional size of the step, compared to the 
;                 parameter value. This value supercedes the STEP 
;                 setting. If the parameter is zero, then a default 
;                 step size is chosen. 
;     .MPSIDE - selector for type of derivative calculation. This 
;               field can take one of five possible values: 
;                   0 - one-sided derivative computed automatically 
;                   1 - one-sided derivative (f(x+h) - f(x) )/h 
;                 -1 - one-sided derivative (f(x) - f(x-h))/h 
;                   2 - two-sided derivative (f(x+h) - f(x-h))/(2*h) 
;                   3 - explicit derivative used for this parameter 
;               In the first four cases, the derivative is approximated 
;               numerically by finite difference, with step size 
;               H=STEP, where the STEP parameter is defined above. The 
;               last case, MPSIDE=3, indicates to allow the user 
;               function to compute the derivative explicitly (see 
;               section on "EXPLICIT DERIVATIVES"). AUTODERIVATIVE=0 
;               overrides this setting for all parameters, and is 
;               equivalent to MPSIDE=3 for all parameters. For 
;               MPSIDE=0, the "automatic" one-sided derivative method 
;               will chose a direction for the finite difference which 
;               does not violate any constraints. The other methods 
;               (MPSIDE=-1 or MPSIDE=1) do not perform this check. The 
;               two-sided method is in principle more precise, but 
;               requires twice as many function evaluations. Default: 
;               0. 
;     .MPDERIV_DEBUG - set this value to 1 to enable debugging of 
;               user-supplied explicit derivatives (see "TESTING and 
;               DEBUGGING" section above). In addition, the 
;               user must enable calculation of explicit derivatives by 
;               either setting AUTODERIVATIVE=0, or MPSIDE=3 for the 
;               desired parameters. When this option is enabled, a 
;               report may be printed to the console, depending on the 
;               MPDERIV_ABSTOL and MPDERIV_RELTOL settings. 
;               Default: 0 (no debugging) 
    
;     .MPDERIV_ABSTOL, .MPDERIV_RELTOL - tolerance settings for 
;               print-out of debugging information, for each parameter 
;               where debugging is enabled. See "TESTING and 
;               DEBUGGING" section above for the meanings of these two 
;               fields. 
;     .MPMAXSTEP - the maximum change to be made in the parameter 
;                   value. During the fitting process, the parameter 
;                   will never be changed by more than this value in 
;                   one iteration. 
;                   A value of 0 indicates no maximum. Default: 0. 
  
;     .TIED - a string expression which "ties" the parameter to other 
;             free or fixed parameters as an equality constraint. Any 
;             expression involving constants and the parameter array P 
;             are permitted. 
;             Example: if parameter 2 is always to be twice parameter 
;             1 then use the following: parinfo[2].tied = '2 * P[1]'. 
;             Since they are totally constrained, tied parameters are 
;             considered to be fixed; no errors are computed for them, 
;             and any LIMITS are not obeyed. 
;             [ NOTE: the PARNAME can't be used in a TIED expression. ] 
;     .MPPRINT - if set to 1, then the default ITERPROC will print the 
;                 parameter value. If set to 0, the parameter value 
;                 will not be printed. This tag can be used to 
;                 selectively print only a few parameter values out of 
;                 many. Default: 1 (all parameters printed) 
;     .MPFORMAT - IDL format string to print the parameter within 
;                 ITERPROC. Default: '(G20.6)' (An empty string will 
;                 also use the default.) 
;   Future modifications to the PARINFO structure, if any, will involve 
;   adding structure tags beginning with the two letters "MP". 
;   Therefore programmers are urged to avoid using tags starting with 
;   "MP", but otherwise they are free to include their own fields 
;   within the PARINFO structure, which will be ignored by MPFIT. 
  
;   PARINFO Example: 
;   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $ 
;                       limits:[0.D,0]}, 5) 
;   parinfo[0].fixed = 1 
;   parinfo[4].limited[0] = 1 
;   parinfo[4].limits[0] = 50.D 
;   parinfo[*].value = [5.7D, 2.2, 500., 1.5, 2000.] 
  
;   A total of 5 parameters, with starting values of 5.7, 
;   2.2, 500, 1.5, and 2000 are given. The first parameter 
;   is fixed at a value of 5.7, and the last parameter is 
;   constrained to be above 50.
            
end
