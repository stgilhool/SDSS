pro moment_of_inertia


mass = dindgen(12)*0.1d0 + 0.1

i_cz = [0.0003d0, $
        0.002d0, $
        0.0055d0, $
        0.009d0, $
        0.011, $
        0.012, $
        0.013, $
        0.014, $
        0.013, $
        0.009, $
        0.0055, $
        0.0021]

i_rc = [!values.d_nan, $
        !values.d_nan, $
        !values.d_nan, $
        0.0025d0, $
        0.009, $
        0.018, $
        0.027, $
        0.037, $
        0.05, $
        0.063, $
        0.08, $
        0.09]

window, 0, xs = 1300, ys = 900
!p.multi = [0,1,2]

plot, mass, i_cz, ytit = 'I_cz', ps=6, charsize = 2, title = "Moments of Inertia vs. Mass for 1 Gyr star", /ylog, xr = [0,1.3], /xs

plot, mass, i_rc, ytit = 'I_rad', xtit = 'Mass (M_sun)', charsize = 2, ps=6, /ylog, xr = [0,1.3], /xs

stop

end
