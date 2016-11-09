; This routine will advance the system state by a single courant time
; step.
; The step is acheived in two substeps using a RK midpoint of RK_frac*dt

; SV_struc contains State Variables
; CS_struc contains Coordinate System variables
; dt_micro is the time jump
; RK_frac is the RK midpoint fraction

pro ConFlow_microstepper, SV_struc=SV_struc, CS_struc=CS_struc, dt_micro=dt_micro, RK_frac=RK_frac, zero_beta=zero_beta, debug=debug

; we have to select the grid that most closely hugs the intrusion.
ss_pos = where(CS_struc.f ge 0)
ss_neg = where(CS_struc.f le 0)
f_min_p = 0
f_max_m = 0
if n_elements(ss_pos) ne 1 then f_min_p = min(CS_struc.f[ss_pos])
if n_elements(ss_neg) ne 1 then f_max_m = max(CS_struc.f[ss_neg])
ss_r = where(CS_struc.f eq f_min_p and abs(CS_struc.g) le 2.)
ss_l = where(CS_struc.f eq f_max_m and abs(CS_struc.g) le 2.)

DSV_struc = ConFlow_SV_derivs(SV_struc=SV_struc, CS_struc=CS_struc, debug=debug)

; Now we'll calculate the state at a time RK_frac * dt_micro
; into the future

dt_RK = RK_frac * dt_micro
SV_struc_RK = SV_struc

if keyword_set(zero_beta) then begin
   DSV_struc.dphi_dt = 0.*DSV_struc.dphi_dt
   DSV_struc.duf_dt = 0.*DSV_struc.duf_dt
endif

SV_struc_RK.time = SV_struc.time + dt_RK
SV_struc_RK.phi  = SV_struc.phi  + dt_RK * DSV_struc.dphi_dt
SV_struc_RK.rho  = SV_struc.rho  + dt_RK * DSV_struc.drho_dt
SV_struc_RK.p    = SV_struc.p    + dt_RK * DSV_struc.dp_dt
SV_struc_RK.uf   = SV_struc.uf   + dt_RK * DSV_struc.duf_dt
SV_struc_RK.ug   = SV_struc.ug   + dt_RK * DSV_struc.dug_dt

;density will be held fixed at the bottom boundary but flow on the
;sides and top
SV_struc_RK.rho[*, 0] = SV_struc.rho[*, 0]

SV_struc_RK.rho[*,-1] = SV_struc_RK.rho[*,-2]
SV_struc_RK.rho[ 0,*] = SV_struc_RK.rho[ 1,*]
SV_struc_RK.rho[-1,*] = SV_struc_RK.rho[-2,*]

;flux function will be fixed on the sides but flow on the bottom
SV_struc_RK.phi[*, 0] = SV_struc_RK.phi[*, 1]
SV_struc_RK.phi[*,-1] = SV_struc_RK.phi[*,-2]

SV_struc_RK.phi[ 0,*] = SV_struc.phi[ 0,*]
SV_struc_RK.phi[-1,*] = SV_struc.phi[-1,*]

;flux function must also respect the interior boundary
if n_elements(ss_l) ne 1 then SV_struc_RK.phi[ss_l] = SV_struc_RK.phi[ss_l] < 0
if n_elements(ss_r) ne 1 then SV_struc_RK.phi[ss_r] = SV_struc_RK.phi[ss_r] > 0

;pressure will mimic density 
SV_struc_RK.p[*, 0] = SV_struc.p[*, 0]

SV_struc_RK.p[*,-1] = SV_struc_RK.p[*,-2]
SV_struc_RK.p[ 0,*] = SV_struc_RK.p[ 1,*]
SV_struc_RK.p[-1,*] = SV_struc_RK.p[-2,*]

;uf will be fixed to zero on the sides
SV_struc_RK.uf[*, 0] = SV_struc_RK.uf[*, 1]
SV_struc_RK.uf[*,-1] = SV_struc_RK.uf[*,-2]

SV_struc_RK.uf[ 0,*] = 0
SV_struc_RK.uf[-1,*] = 0

;uf will also be zero along the intrusion
if n_elements(ss_l) ne 1 then SV_struc_RK.uf[ss_l] = 0
if n_elements(ss_r) ne 1 then SV_struc_RK.uf[ss_r] = 0

;ug will be fixed on bottom
SV_struc_RK.ug[*, 0] = SV_struc.ug[*, 0]

SV_struc_RK.ug[*,-1] = SV_struc_RK.ug[*,-2]
SV_struc_RK.ug[ 0,*] = SV_struc_RK.ug[ 1,*]
SV_struc_RK.ug[-1,*] = SV_struc_RK.ug[-2,*]



; Now we get the full time step using the RK midpoint to estimate the derivatives

DSV_struc = ConFlow_SV_derivs(SV_struc = SV_struc_RK, CS_struc=CS_struc, debug=debug)

SV_struc_next = temporary(SV_struc_RK)

; Here we'll build in a stop-check against negative densities
if min(SV_struc.rho + dt_micro * DSV_struc.drho_dt) le 0 then begin
   print, 'step size leads to unphysical density'
   print, 'step size '+trim(string(dt_micro))+' reduced to '+trim(string(0.75*dt_micro))
   dt_micro = 0.75*dt_micro
endif

if keyword_set(zero_beta) then begin
   DSV_struc.dphi_dt = 0.*DSV_struc.dphi_dt
   DSV_struc.duf_dt = 0.*DSV_struc.duf_dt
endif

SV_struc_next.time = SV_struc.time + dt_micro
SV_struc_next.phi  = SV_struc.phi  + dt_micro * DSV_struc.dphi_dt
SV_struc_next.rho  = SV_struc.rho  + dt_micro * DSV_struc.drho_dt
SV_struc_next.p    = SV_struc.p    + dt_micro * DSV_struc.dp_dt
SV_struc_next.uf   = SV_struc.uf   + dt_micro * DSV_struc.duf_dt
SV_struc_next.ug   = SV_struc.ug   + dt_micro * DSV_struc.dug_dt


;density will be held fixed at the bottom boundary but flow on the
;sides and top
SV_struc_next.rho[*, 0] = SV_struc.rho[*, 0]

SV_struc_next.rho[*,-1] = SV_struc_next.rho[*,-2]
SV_struc_next.rho[ 0,*] = SV_struc_next.rho[ 1,*]
SV_struc_next.rho[-1,*] = SV_struc_next.rho[-2,*]

;flux function will be fixed on the sides but flow on the bottom
SV_struc_next.phi[*, 0] = SV_struc_next.phi[*, 1]
SV_struc_next.phi[*,-1] = SV_struc_next.phi[*,-2]

SV_struc_next.phi[ 0,*] = SV_struc.phi[ 0,*]
SV_struc_next.phi[-1,*] = SV_struc.phi[-1,*]

;flux function must also respect the interior boundary.
SV_struc_next.phi[ss_l] = SV_struc_next.phi[ss_l] < 0
SV_struc_next.phi[ss_r] = SV_struc_next.phi[ss_r] > 0

;pressure will mimic density
SV_struc_next.p[*, 0] = SV_struc.p[*, 0]

SV_struc_next.p[*,-1] = SV_struc_next.p[*,-2]
SV_struc_next.p[ 0,*] = SV_struc_next.p[ 1,*]
SV_struc_next.p[-1,*] = SV_struc_next.p[-2,*]

;uf will be fixed to zero on the sides
SV_struc_next.uf[*, 0] = SV_struc_next.uf[*, 1]
SV_struc_next.uf[*,-1] = SV_struc_next.uf[*,-2]

SV_struc_next.uf[ 0,*] = 0
SV_struc_next.uf[-1,*] = 0

;uf will also be zero along the intrusion
if n_elements(ss_l) ne 1 then SV_struc_next.uf[ss_l] = 0
if n_elements(ss_r) ne 1 then SV_struc_next.uf[ss_r] = 0

;ug will be fixed on bottom
SV_struc_next.ug[*, 0] = SV_struc.ug[*, 0]

SV_struc_next.ug[*,-1] = SV_struc_next.ug[*,-2]
SV_struc_next.ug[ 0,*] = SV_struc_next.ug[ 1,*]
SV_struc_next.ug[-1,*] = SV_struc_next.ug[-2,*]

; And finally we advance the system be replacing the state with the
; advanced state

if keyword_set(debug) then stop

SV_struc = temporary(SV_struc_next)

end
