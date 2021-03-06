; This is a collection of routines for finding the change in each
; state variable. Each variable gets its own routine. The state
; variables are (for now):
; phi : the magnetic flux function
; uf : the projection of velocity onto grad f
; ug : the projection of velocity onto grad g
; rho : the density
; p : the pressure.


function ConFlow_SV_derivs, SV_struc = SV_struc, CS_struc = CS_struc, DSV_struc = DSV_struc, debug=debug
  
; First we need to set up a velocity direction array that will be used
; for upwind differencing.

f_vel_sgn = sign(SV_struc.uf*0+1., SV_struc.uf)
g_vel_sgn = sign(SV_struc.ug*0+1., SV_struc.ug)

; In order to deal with the intrusion, we'll set the velocity
; direction to be right_to_left on the right side and visa-versa.
; This avoids taking derivatives across the intrusion.

ss_f_p = where(CS_struc.f gt 0)
if n_elements(ss_f_p) ne 1 then begin
   f_p_min = min(CS_struc.f[ss_f_p])
   ss_p_r = where(CS_struc.f eq f_p_min and abs(CS_struc.g) le 2)
   f_vel_sgn[ss_p_r] = -1.
endif

ss_f_n = where(CS_struc.f lt 0)   
if n_elements(ss_f_n) ne 1 then begin
   f_n_max = max(CS_struc.f[ss_f_n])
   ss_n_r = where(CS_struc.f eq f_n_max and abs(CS_struc.g) le 2)
   f_vel_sgn[ss_n_r] = 1.
endif







;;;;;;;;;;;;;;
; Continuity ;
;;;;;;;;;;;;;;




; SV_struc.rho is the density
; SV_struc.uf and SV_struc.ug and the components of velocity
; CS_struc.met is the metric

; In the conformal space, the continuity equation is:

; drho_dt = - s^2 ( d(SV_struc.rho SV_struc.uf)/df + d(SV_struc.rho SV_struc.ug)/dg)

field_grad, SV_struc.rho * SV_struc.uf, df_ruf, unimportant, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
field_grad, SV_struc.rho * SV_struc.ug, unimportant, dg_rug, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
udef = size(temporary(unimportant))

drho_dt = - CS_struc.met^2. * (df_ruf + dg_rug)









;;;;;;;;;;;;;
; Induction ;
;;;;;;;;;;;;;



; SV_struc.phi is the magnetic SV_struc variable
   ; The actual field is phi + f
   ; df_field = 1 + df_phi
   ; dg_field = dg_phi
   ; lap_field = lap_phi
; CS_struc.f and CS_struc.g are the CS_strucinate grid
; SV_struc.uf and SV_struc.ug are the velocity salars
; CS_struc.met is the metric of the space.

; In the conformal space, the induction equation is:

; dphi_dt = - s^2 (SV_struc.ug dphi_dg + SV_struc.uf dphi_df)

field_grad, SV_struc.phi, dphi_df, dphi_dg, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff

dphi_dt = - CS_struc.met^2. * (SV_struc.uf * (1 + dphi_df) + SV_struc.ug * dphi_dg)
















;;;;;;;;;;
; Energy ;
;;;;;;;;;;




; In the conformal space the adiabatic energy equation is
; dp_dt = - s^2 ( uf*dp/df + ug*dp/dg + gamma p duf/df + gamma p dug/dg)
; dp_dt = - s^2 (d(p uf)/df + d(p ug)/dg + (gamma - 1) * p * (duf/df + dug/dg) )


field_grad, SV_struc.uf, df_uf, dg_uf, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
field_grad, SV_struc.ug, df_ug, dg_ug, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
field_grad, SV_struc.p * SV_struc.uf, df_puf, dg_puf, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
field_grad, SV_struc.p * SV_struc.ug, df_pug, dg_pug, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
udef = size(temporary(df_pug))
udef = size(temporary(dg_puf))

dp_dt = - CS_struc.met^2. * (df_puf + dg_pug + (CS_struc.gamma - 1) * SV_struc.p * (df_uf + dg_ug) )






;;;;;;;;;;;;
; Momentum ;
;;;;;;;;;;;;

; several derivatives are common to both forms of the momentum
; equation. 

;We'll be using centered differencing in every case but we'll also
;have to evaluate the uw diff for each value in order to replace the
;interior boundary.

; some of the uw stuff we've already done
df_phi_uw = temporary(dphi_df)
dg_phi_uw = temporary(dphi_dg)
df_uf_uw  = temporary(df_uf)
df_ug_uw  = temporary(df_ug)
dg_uf_uw  = temporary(dg_uf)
dg_ug_uw  = temporary(dg_ug)

; several terms, however, must be done fresh
field_grad, CS_struc.met, df_met_uw, dg_met_uw, xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff
field_grad, SV_struc.p,   df_p_uw,   dg_p_uw,   xx = CS_struc.f, yy = CS_struc.g, uxb = f_vel_sgn, uyb = g_vel_sgn, /uw_diff

; Now we get all these for centered differencing
field_grad, SV_struc.phi, df_phi, dg_phi, xx=CS_struc.f, yy=CS_struc.g
field_grad, SV_struc.p,   df_p,   dg_p,   xx=CS_struc.f, yy=CS_struc.g
field_grad, SV_struc.uf,  df_uf,  dg_uf,  xx=CS_struc.f, yy=CS_struc.g
field_grad, SV_struc.ug,  df_ug,  dg_ug,  xx=CS_struc.f, yy=CS_struc.g
field_grad, CS_struc.met, df_met, dg_met, xx=CS_struc.f, yy=CS_struc.g

; And now we put them together
if n_elements(ss_p_r) eq 0 then ss_r = ss_n_r
if n_elements(ss_n_r) eq 0 then ss_r = ss_p_r
if n_elements(ss_n_r) ne 0 and n_elements(ss_p_r) ne 0 then ss_r = [ss_n_r,ss_p_r]

df_phi[ss_r] = df_phi_uw[ss_r]
dg_phi[ss_r] = dg_phi_uw[ss_r]

df_p[ss_r]   = df_p_uw[ss_r]
dg_p[ss_r]   = dg_p_uw[ss_r]

df_uf[ss_r]  = df_uf_uw[ss_r]
dg_uf[ss_r]  = dg_uf_uw[ss_r]
df_ug[ss_r]  = df_ug_uw[ss_r]
dg_ug[ss_r]  = dg_ug_uw[ss_r]

df_met[ss_r] = df_met_uw[ss_r]
dg_met[ss_r] = dg_met_uw[ss_r]

; The momentum equation has some complicated terms.


;;;;;;;;;;;;;;;;
; Momentum - f ;
;;;;;;;;;;;;;;;;

; For now, we'll just use the easiest possible viscocity.

; first term: - ( u dot grad ) u
; fmom_t1: = - s^2 ( SV_struc.uf*df_uf + SV_struc.ug*dg_uf) - s (SV_struc.uf^2 + SV_struc.ug^2) df_s

fmom_t1 = - CS_struc.met^2. * ( SV_struc.uf * df_uf + SV_struc.ug * dg_uf) - CS_struc.met * (SV_struc.uf^2. + SV_struc.ug^2.)*df_met

; second term = - ( grad p ) / SV_struc.rho
; fmom_t2 = - df_p / SV_struc.rho

fmom_t2 = - df_p / SV_struc.rho

; thrid term = (1/4 pi) ( ( curl b) cross b) / SV_struc.rho
; fmom_t3 = - (1/4pi) * (s^2 / SV_struc.rho) * df_phi * lap_phi

field_lap, SV_struc.phi, lap_phi, xx = CS_struc.f, yy = CS_struc.g, lbc = 'zero', rbc = 'zero'

fmom_t3 = - ( CS_struc.met^2. / ( 4.*!pi*SV_struc.rho) ) * (1 + df_phi) * lap_phi * SV_struc.phi_not^2.

; mom_t4 = mu s^2 * lap U / SV_struc.rho ::::::: this will be updated later.

field_lap, CS_struc.met * SV_struc.uf, lap_muf, xx = CS_struc.f, yy = CS_struc.g

fmom_t4 = SV_struc.mu * CS_struc.met * lap_muf / sv_struc.rho

; So we put it all together

duf_dt = fmom_t1 + fmom_t2 + fmom_t3 + fmom_t4




;;;;;;;;;;;;;;;;
; Momentum - g ;
;;;;;;;;;;;;;;;;

; For now, we'll just use the easiest possible viscocity.

; first term: - ( u dot grad ) u
; gmom_t1: = - s^2 ( SV_struc.uf*df_ug + SV_struc.ug*dg_ug) - s (SV_struc.uf^2 + SV_struc.ug^2) df_s

gmom_t1 = - CS_struc.met^2. * ( SV_struc.uf * df_ug + SV_struc.ug * dg_ug) - CS_struc.met * (SV_struc.uf^2. + SV_struc.ug^2.)*dg_met

; second term = - ( grad p ) / SV_struc.rho
; fmom_t2 = - df_p / SV_struc.rho

gmom_t2 = - dg_p / SV_struc.rho

; thrid term = (1/4 pi) ( ( curl b) cross b) / SV_struc.rho
; fmom_t3 = (1/4pi) * (s^2 / SV_struc.rho) * df_phi * lap_phi

gmom_t3 = ( CS_struc.met^2. / ( 4.*!pi*SV_struc.rho) ) * dg_phi * lap_phi * SV_struc.phi_not^2.

; mom_t4 = mu s^2 * lap U / SV_struc.rho ::::::: this will be updated later.

field_lap, CS_struc.met * SV_struc.ug, lap_mug, xx = CS_struc.f, yy = CS_struc.g

gmom_t4 = SV_struc.mu * CS_struc.met * lap_mug / sv_struc.rho

; So we put it all together

dug_dt = gmom_t1 + gmom_t2 + gmom_t3 + gmom_t4
















DSV_struc = {dphi_dt:dphi_dt, drho_dt:drho_dt, dp_dt:dp_dt, duf_dt:duf_dt, dug_dt:dug_dt}

if keyword_set(debug) then stop

return, DSV_struc

end
