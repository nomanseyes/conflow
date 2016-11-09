;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

function mf_shock_angle, vel, v_s, v_a

; all possible angles of shock front wrt flow
theta = findgen(1000)*(!pi / 2) / 999

; assuming the waves propagate in the direciton of the shock normal,
; this is the fms wave speed
vf_not = sqrt(v_s^2 + v_a^2)
vf_theta = sqrt(0.5*vf_not^2. + 0.5*sqrt(vf_not^4 - 4. * v_s^2 * v_a^2 * cos(theta)^2.))

; this is the velocity of the flow projected to the shock normal
vel_projected = vel * cos(theta)

; the mach number of that projected velocity wrt the angular fms wave speed
mach_projected = vel_projected / vf_theta

; extracting the angle that gives mach = 1
theta_mach1 = interpol(theta, mach_projected,1)

; the complement of which is the opening angle
theta_cone_angle = !pi / 2 - theta_mach1

return, theta_cone_angle

end
