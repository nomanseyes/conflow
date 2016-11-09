; this procedure takes an array of input filenames and then generates
; a visualization of the fluid state for every time step in every
; file.

pro vis_gen, files, double = double

mov_dirs = file_dirname(files)+'/'+file_basename(files, '.sav')+'_mov/'

for j = 0, n_elements(files) - 1 do begin
   if file_test(mov_dirs[j]) eq 0 then spawn, 'mkdir '+mov_dirs[j]
   restore, files[j], /v

   if keyword_set(double) then begin
      f_double = [-reverse(con_struc.cs_struc.f,1),   con_struc.cs_struc.f]
      g_double = [+reverse(con_struc.cs_struc.g,1),   con_struc.cs_struc.g]
      m_double = [+reverse(con_struc.cs_struc.met,1), con_struc.cs_struc.met]
   endif

   conflow_remap, conflow_struc = con_struc, carflow_struc = car_struc, double=double

   for i = 0, n_elements(con_struc.sv_struc_arr) - 1 do begin

      if keyword_set(double) then begin
         rho_double = [+reverse(con_struc.sv_struc_arr[i].rho,1), con_struc.sv_struc_arr[i].rho]
         phi_double = [-reverse(con_struc.sv_struc_arr[i].phi,1), con_struc.sv_struc_arr[i].phi]
      endif

      beta = 8.*!pi / con_struc.sv_struc_arr[0].phi_not^2.
      if con_struc.sv_struc_arr[0].phi_not eq -1 then beta = 0
      ms = con_struc.sv_struc_arr[0].ug[0,0]
      ma = sqrt(ms^2 * beta / 2)

      set_plot, 'ps'
      device, filename = mov_dirs[j]+trim(string(i, format = '(I3.3)'))+'.eps', xs = 10, ys = 3, /inches, /encaps, /color, bits = 12
      !p.font = 0
      !p.multi = [0,3,1]
      plot, con_struc.sv_struc_arr[i].ug[5,*]*con_struc.cs_struc.met[5,*], con_struc.cs_struc.g[0,*], xr = [0,5], xstyle = 1, yr = [min(con_struc.cs_struc.g),max(con_struc.cs_struc.g)], ystyle = 1, ytitle = 'height', thick = 3
      oplot, con_struc.sv_struc_arr[i].rho[5,*], con_struc.cs_struc.g[0,*], linestyle = 3, thick = 3
      !p.multi = [2,3,1]
      !y.tickname = ' '

      loadct, 3

      if keyword_set(double) then begin
         display, alog10(rho_double^2.), f_double[*,0], g_double[0,*], min = -0.5, max = 0.5, aspect = 1, /noerase, $
               title = 'T = '+trim(string(con_struc.sv_struc_arr[i].time, format = '(f4.1)'))+'R/Cs'
         contour, /over, /noerase, phi_double + f_double, f_double, g_double, levels = indgen(21) - 10
      endif else begin
         display, alog10(con_struc.sv_struc_arr[i].rho^2.), con_struc.cs_struc.f[*,0], con_struc.cs_struc.g[0,*], min = -0.5, max = 0.5, aspect = 1, /noerase, $
               title = 'T = '+trim(string(con_struc.sv_struc_arr[i].time, format = '(f4.1)'))+'R/Cs'
         contour, /over, /noerase, con_struc.sv_struc_arr[i].phi + con_struc.cs_struc.f, con_struc.cs_struc.f, con_struc.cs_struc.g, levels = ((indgen(21) - 10.)/10.)*max(con_struc.cs_struc.f)
      endelse


      !p.multi = [1,3,1]
      
      if keyword_set(double) then begin
         display, alog10(car_struc.rho_tg[*,*,i]^2), f_double[*,0], g_double[0,*], min = -0.5, max = 0.5, aspect = 1, /noerase, $
                  title = 'beta = '+trim(string(beta, format = '(f4.2)'))+', Ms = '+trim(string(ms, format = '(f4.2)'))+', Ma = '+trim(string(Ma, format = '(f4.2)'))
         contour, /over, /noerase, car_struc.phi_tg[*,*,i] + car_struc.f_tg, f_double, g_double, levels = indgen(21) - 10
      endif else begin
         display, alog10(car_struc.rho_tg[*,*,i]^2), con_struc.cs_struc.f[*,0], con_struc.cs_struc.g[0,*], min = -0.5, max = 0.5, aspect = 1, /noerase, $
                  title = 'beta = '+trim(string(beta, format = '(f4.2)'))+', Ms = '+trim(string(ms, format = '(f4.2)'))+', Ma = '+trim(string(Ma, format = '(f4.2)'))
         contour, /over, /noerase, car_struc.phi_tg[*,*,i] + car_struc.f_tg, con_struc.cs_struc.f, con_struc.cs_struc.g, levels = ((indgen(21) - 10.)/10.)*max(con_struc.cs_struc.f)
      endelse
      !p.multi = 0
      !p.font = -1
      colorbar, range = [-0.5, 0.5], pos = [0.32, 0.1, 0.35, 0.9], /vertical, /normal, format = '(f4.1)', divisions = 10
      loadct, 0
      device, /close
      set_plot, 'x'
      !y.tickname = ''

   endfor
endfor
end
