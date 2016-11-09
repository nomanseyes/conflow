; This routine remaps from the conformal coordinates back to cartesian
; for visualization.

pro ConFlow_remap, ConFlow_struc=ConFlow_struc, CarFlow_struc=CarFlow_struc, double = double

f = ConFlow_struc.CS_struc.f
g = ConFlow_struc.CS_struc.g
m = ConFlow_struc.CS_struc.met

if keyword_set(double) then begin
   f = [-reverse(f,1), f]
   g = [reverse(g,1), g]
   m = [reverse(m,1), m]
endif

xp = 0.5 * ( f + Real_Part( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
xm = 0.5 * ( f - Real_Part( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
yp = 0.5 * ( g + Imaginary( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )
ym = 0.5 * ( g - Imaginary( sqrt( f^2. + 2. * f * g * complex(0,1) - g^2. + 4 ) ) )

x_conf = xp
x_conf[where(f le 0)] = xm[where(f le 0)]
y_conf = yp
y_conf[where(f le 0)] = ym[where(f le 0)]

nx = n_elements(x_conf[*,0])
ny = n_elements(x_conf[0,*])
nt = n_elements(ConFlow_struc.SV_struc_arr)

x_cart = ( min(x_conf) + (max(x_conf) - min(x_conf)) * indgen(nx)/(nx-1.) ) # replicate(1, ny)
y_cart = replicate(1, nx) # ( min(y_conf) + (max(y_conf) - min(y_conf)) * indgen(ny)/(ny-1.) )
dx = x_cart[1,0] - x_cart[0,0]
dy = y_cart[0,1] - y_cart[0,0]
limits = [min(x_cart), min(y_cart), max(x_cart), max(y_cart)]
ss_r = where(x_cart^2 + y_cart^2 le 1.01)

triangulate, x_conf, y_conf, triangles
x_tg = trigrid( x_conf, y_conf, x_conf, triangles, nx=nx, ny=ny)
y_tg = trigrid( x_conf, y_conf, y_conf, triangles, nx=nx, ny=ny)
f_tg = trigrid( x_conf, y_conf, f     , triangles, nx=nx, ny=ny)
g_tg = trigrid( x_conf, y_conf, g     , triangles, nx=nx, ny=ny)
m_tg = trigrid( x_conf, y_conf, m     , triangles, nx=nx, ny=ny)

phi_tg = ConFlow_struc.SV_struc_arr.phi
uf_tg = ConFlow_struc.SV_struc_arr.uf
ug_tg = ConFlow_struc.SV_struc_arr.ug
rho_tg = ConFlow_struc.SV_struc_arr.rho
p_tg = ConFlow_struc.SV_struc_arr.p

if keyword_set(double) then begin
   temp = [-reverse(phi_tg, 1), phi_tg]
   phi_tg = temporary(temp)
   temp = [-reverse(uf_tg, 1), uf_tg]
   uf_tg = temporary(temp)
   temp = [reverse(ug_tg, 1), ug_tg]
   ug_tg = temporary(temp)
   temp = [reverse(rho_tg, 1), rho_tg]
   rho_tg = temporary(temp)
   temp = [reverse(p_tg,1), p_tg]
   p_tg = temporary(temp)
endif


; these values will be replaced by their triangulated equivalents

for i = 0, nt - 1 do begin
   temp_arr = trigrid( x_conf, y_conf, phi_tg[*,*,i], triangles, nx=nx, ny=ny)
   temp_arr[ss_r] = 0
   phi_tg[0,0,i] = temp_arr

   temp_arr = trigrid( x_conf, y_conf,  uf_tg[*,*,i], triangles, nx=nx, ny=ny)
   temp_arr[ss_r] = 0
   uf_tg[0,0,i]  = temp_arr

   temp_arr = trigrid( x_conf, y_conf,  ug_tg[*,*,i], triangles, nx=nx, ny=ny)
   temp_arr[ss_r] = ug_tg[0,0,0]
   ug_tg[0,0,i]  = temp_arr

   temp_arr = trigrid( x_conf, y_conf, rho_tg[*,*,i], triangles, nx=nx, ny=ny)
   temp_arr[ss_r] = 0
   rho_tg[0,0,i] = temp_arr

   temp_arr = trigrid( x_conf, y_conf,   p_tg[*,*,i], triangles, nx=nx, ny=ny)
   temp_arr[ss_r] = 0
   p_tg[0,0,i]   = temp_arr

   print, 'Finished remapping frame '+trim(string(i+1))+' of '+trim(string(nt))
   
endfor

carflow_struc = {phi_tg:phi_tg, uf_tg:uf_tg, ug_tg:ug_tg, rho_tg:rho_tg, p_tg:p_tg, x_cart:x_cart, y_cart:y_cart, f_tg:f_tg, g_tg:g_tg, m_tg:m_tg}

end
