; This is the main level function for the ConFlow code set
; It defines the coordinate system and the general simulation grid

; CS_struc is the coordinate system structure
; SV_struc_arr is the state variables (it is an array of structures)

function ConFlow_main, f_range=f_range, g_range=g_range, f_num=f_num, g_num=g_num, phi_not=phi_not, t_num=t_num, t_range=t_range, mu=mu, CT_safe=CT_safe, uf_not=uf_not, ug_not=ug_not, conflow_struc_init = conflow_struc_init, cloan_prev = cloan_prev, zero_beta = zero_beta, debug=debug, rho_init = rho_init, prs_init = prs_init, gamma = gamma, vis=vis, loud=loud, mach_a = mach_a, beta_not=beta_not

if n_elements(conflow_struc_init) ne 0 then goto, jump_resume

if keyword_set(gamma) eq 0 then gamma = 1.

if keyword_set(zero_beta) then begin
   uf_not = 0.
   phi_not = -1.
endif

if n_elements(ug_not)  eq 0 then ug_not  = 1.
if n_elements(uf_not)  eq 0 then uf_not  = 0.
if n_elements(phi_not) eq 0 then phi_not = sqrt(8.*!pi)
if n_elements(mach_a)  eq 1 then phi_not = sqrt(4.*!pi)*(ug_not / mach_a)
if n_elements(beta_not)eq 1 then phi_not = sqrt(8.*!pi / beta_not)
if n_elements(mu)      eq 0 then mu = 0.1

if n_elements(f_range) eq 0 then f_range = [0,10]
if n_elements(g_range) eq 0 then g_range = [-20,20]
if n_elements(t_range) eq 0 then t_range = [0,20]
if n_elements(f_num)   eq 0 then f_num   = 200
if n_elements(g_num)   eq 0 then g_num   = 800
if n_elements(t_num)   eq 0 then t_num   = 400

f_range = float(f_range)
g_range = float(g_range)
t_range = float(t_range)


f = min(f_range) + (max(f_range) - min(f_range)) * (findgen(f_num) + 0.5) / f_num
f = f # replicate(1, g_num)

g = min(g_range) + (max(g_range) - min(g_range)) * (findgen(g_num) + 0.5) / g_num
g = g # replicate(1, f_num)
g = transpose(g)

xp = 0.5 * ( f + Real_Part( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
xm = 0.5 * ( f - Real_Part( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
yp = 0.5 * ( g + Imaginary( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
ym = 0.5 * ( g - Imaginary( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )

ssp = where(f ge 0)
ssm = where(f le 0)
x = xp
y = yp
if n_elements(ssm) ne 1 then begin
   x[ssm] = xm[ssm]
   y[ssm] = ym[ssm]
endif

met = sqrt(1 + ((1 + 2.*(x^2. - y^2.)) / ((x^2. + y^2.)^2.)))
; let's filter this to get rid of the extreme behavior. 
ss_sm = where(met le 1)

met[ss_sm] = 0.5*(met[ss_sm]^2. + 1.)
for i = 0, 4 do met = smooth(met, [1,7])


phi = f*0.
; note that the magnetic field will be phi_not * curl z (phi + f)
rho = met
uf = f*0. + uf_not
ug = ug_not / rho
p = rho
if keyword_set(rho_init) then rho = rho_init
if keyword_set(prs_init) then p   = prs_init

SV_struc = {phi_not:phi_not, phi:phi, uf:uf, ug:ug, rho:rho, p:p, mu:mu, time:min(t_range)}

CS_struc = {f:f, g:g, met:met, gamma:gamma}

SV_struc_arr = replicate(SV_struc, t_num + 1)

dt_macro = (max(t_range) - min(t_range)) / float(t_num)
nt_supp = 0

jump_resume: if n_elements(conflwo_struc_init) ne 0 then print, 'resuming from previous structure'

if n_elements(conflow_struc_init) ne 0 then begin

   if keyword_set(t_num) then print, 'cannot specify number of frames with provided previous structure'
   if keyword_set(t_range) eq 0 then print, 'cannot resume without specified time range'

   nt_supp = n_elements(conflow_struc_init.SV_struc_arr) - 1
   tr_supp = [min(conflow_struc_init.SV_struc_arr.time), max(conflow_struc_init.SV_struc_arr.time)]

   if t_range[0] ne tr_supp[0] then print, 'resetting start time to match supplied structure'
   t_range[0] = tr_supp[0]
   t_num = round( float(nt_supp) * (max(t_range) - min(t_range)) / (max(tr_supp) - min(tr_supp)))
   dt_macro = ( max(t_range) - min(t_range) ) / float(t_num)

   CS_struc = conflow_struc_init.CS_struc
   SV_struc_arr = replicate(conflow_struc_init.SV_struc_arr[nt_supp] , t_num - nt_supp + 1)

endif


for i = 1, (t_num - nt_supp) do begin
   SV_struc_temp = SV_struc_arr[i-1]
   ConFlow_macrostepper, SV_struc=SV_struc_temp, CS_struc=CS_struc, dt_macro=dt_macro, CT_safe=CT_safe, vis=vis, loud=loud, zero_beta=zero_beta, debug=debug
   SV_struc_arr[i] = temporary(SV_struc_temp)
   print, 'finished time step '+trim(string(i + nt_supp))+' of '+trim(string(t_num))
endfor

if keyword_set(cloan_prev) then begin
   SV_struc_arr_temp = replicate(SV_struc_arr[0], t_num + 1)
   SV_struc_arr_temp[0:nt_supp] = conflow_struc_init.SV_struc_arr
   SV_struc_arr_temp[nt_supp + 1: t_num] = SV_struc_arr[1: t_num - nt_supp]
   SV_struc_arr = temporary(SV_struc_arr_temp)
endif

sys_struc = {SV_struc_arr:SV_struc_arr, CS_struc:CS_struc}

if keyword_set(debug) then stop

return, sys_struc

end
