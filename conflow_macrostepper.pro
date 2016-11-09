; this routine will advance the state by a designated time step by
; performing courant substeps until the necessary time is acheived.

; st_struc contains the fluid state
; cs_struc contains the coordinate system
; t_step is the desired time step
; RK_safe is the Rhunga-kutta factor
; CT_safe is the courant time factor

pro ConFlow_macrostepper, SV_struc=SV_struc, cs_struc=cs_struc, dt_macro=dt_macro, RK_frac=RK_frac, CT_safe=CT_safe, vis=vis, loud=loud, zero_beta=zero_beta, debug=debug

if keyword_set(CT_safe) eq 0 then CT_safe = 2.5
if keyword_set(RK_frac) eq 0 then RK_frac = 0.51

; we'll set up a counter and keep track of time elapsed

T_tot = 0
N_tot = 0
dt_mean = 0

field_grad, CS_struc.f, df, temp, lbc='open', rbc='open', tbc='open',bbc='open'
field_grad, CS_struc.g, temp, dg, lbc='open', rbc='open', tbc='open',bbc='open'
temp = 0

while T_tot lt dt_macro do begin

   field_grad, SV_struc.phi, dphi_df, dphi_dg, xx = CS_struc.f, yy = CS_struc.g
   mag_b_sqr = SV_struc.phi_not^2. * CS_struc.met^2. * ((1 + dphi_df)^2. + dphi_dg^2.)

; we'll express the wave characteristic speeds in cartesian
; coordinats first. 
; For viscocity, recall that df ~ s dx
   
; The viscocity length scale will always go as at least the grid scale for
; fastest growing modes but may be higher for strong gradients in the velocity.

   field_grad, alog(abs(SV_struc.uf*CS_struc.met) + 0.001), dx=df, dy=dg, kff, kfg
   field_grad, alog(abs(SV_struc.ug*CS_struc.met) + 0.001), dx=df, dy=dg, kgf, kgg
   kf = sqrt(kff^2. + kfg^2.) > 1. / df
   kg = sqrt(kgf^2. + kgg^2.) > 1. / dg

   Vf_visc = SV_struc.mu * CS_struc.met * kf / sv_struc.rho
   Vg_visc = SV_struc.mu * CS_struc.met * kg / sv_struc.rho

   Va = sqrt( mag_b_sqr  / ( SV_struc.rho * 4. * !pi ) )
   Cs = sqrt( SV_struc.p / SV_struc.rho )


; pick off the biggest wave speed

   Vw_f = Va > Cs > Vf_visc
   Vw_g = Va > Cs > Vg_visc

; The wave crossing time is the bulk speed + wave speed / grid size.
; The bulk speed in cartesian is s * bulk speed in conformal basis

   ct_f = ( df / CS_struc.met ) / ( abs( SV_struc.uf * CS_struc.met ) + Vw_f)
   ct_g = ( dg / CS_struc.met ) / ( abs( SV_struc.ug * CS_struc.met ) + Vw_g)

   cour_time = min([ct_f[where(finite(ct_f) eq 1 and ct_f ne 0)],ct_g[where(finite(ct_g) eq 1 and ct_g ne 0)]])


   dt_micro = (cour_time / CT_safe) < (dt_macro - T_tot)

   if dt_micro le 0 then begin
      print, 'courant time too short, aborting'
      stop
      GOTO, jump_abort
   endif

   ConFlow_microstepper, SV_struc=SV_struc, CS_struc=CS_struc, dt_micro=dt_micro, RK_frac=RK_frac, zero_beta=zero_beta, debug=debug

   T_tot = T_tot + dt_micro
   dt_mean = (dt_mean * float(N_tot) + dt_micro)/(N_tot + 1.)
   N_tot += 1

endwhile

if keyword_set(vis) then begin
   if min(CS_struc.f) lt 0 and max(CS_struc.f) gt 0 then begin
      sf = where(CS_struc.f[*,0] eq min(abs(CS_struc.f)) > 0)
      !p.multi = [0,2,1]
      display, SV_struc.rho, CS_struc.f[*,0], CS_struc.g[0,*], aspect = 1
      contour, /over, levels = indgen(21) - 10 , SV_struc.phi + CS_struc.f, CS_struc.f, CS_struc.g
      !p.multi = [1,2,1]
      plot, sv_struc.ug[sf,*]*cs_struc.met[sf,*], cs_struc.g[sf,*], linestyle = 0, xr = [0,2.*(SV_struc.ug[0,0]>1)]
      oplot, sv_struc.ug[sf,*], cs_struc.g[sf,*], linestyle = 2
      oplot, sv_struc.rho[sf,*], cs_struc.g[sf,*], linestyle = 3
      oplot, cs_struc.met[sf,*], cs_struc.g[sf,*], linestyle = 4
      legend, ['u_mag', 'u_g', 'rho', 'm'], linestyle = [0,2,3,4]
      !p.multi = 0
   endif else begin
      sf = 0
      !p.multi = [0,2,1]
      display, SV_struc.rho, CS_struc.f[*,0], CS_struc.g[0,*], aspect = 1
      contour, /over, levels = indgen(21) - 10 , SV_struc.phi + CS_struc.f, CS_struc.f, CS_struc.g
      !p.multi = [1,2,1]
      plot, sv_struc.ug[sf,*]*cs_struc.met[sf,*], cs_struc.g[sf,*], linestyle = 0, xr = [0,2.*(SV_struc.ug[0,0]>1)]
      oplot, sv_struc.ug[sf,*], cs_struc.g[sf,*], linestyle = 2
      oplot, sv_struc.rho[sf,*], cs_struc.g[sf,*], linestyle = 3
      oplot, cs_struc.met[sf,*], cs_struc.g[sf,*], linestyle = 4
      legend, ['u_mag', 'u_g', 'rho', 'm'], linestyle = [0,2,3,4]
      !p.multi = 0
   endelse
endif

if keyword_set(loud) then begin
   print, 'Macro-step acheived in '+trim(string(N_tot))+' micro-steps.'
   print, 'Mean microstep time = '+trim(string(dt_mean))+' seconds.'
endif
   
jump_abort:

if keyword_set(debug) then stop

end
