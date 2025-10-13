import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

aspect_ratio_vp = 2.0

d_vp = 1 # ft
theta_vp = np.deg2rad(60.0)

w_vp = 2*d_vp*np.tan(theta_vp)
h_vp = w_vp/aspect_ratio_vp

x_c_vp = np.array([d_vp, d_vp, d_vp, d_vp])
y_c_vp = 0.5*np.array([-w_vp, -w_vp, w_vp, w_vp])
z_c_vp = 0.5*np.array([-h_vp, h_vp, h_vp, -h_vp])

print(x_c_vp)
print(y_c_vp)
print(z_c_vp)

P_fc = np.array([-20.0, 0.0, -5.0])

euler = np.deg2rad(np.array([5.0, -10.0, 0.0]))
c_euler = np.cos(euler)
s_euler = np.sin(euler)
c_phi = c_euler[0]
c_theta = c_euler[1]
c_psi = c_euler[2]
s_phi = s_euler[0]
s_theta = s_euler[1]
s_psi = s_euler[2]

# Rotation matrix
R = np.zeros((3,3))
R[0,0] = c_theta*c_psi
R[0,1] = c_theta*s_psi
R[0,2] = -s_theta
R[1,0] = s_phi*s_theta*c_psi - c_phi*s_psi
R[1,1] = s_phi*s_theta*s_psi + c_phi*c_psi
R[1,2] = s_phi*c_theta
R[2,0] = c_phi*s_theta*c_psi + s_phi*s_psi
R[2,1] = c_phi*s_theta*s_psi - s_phi*c_psi
R[2,2] = c_phi*c_theta

R_inv = sp.linalg.inv(R)

x_f_vp = np.zeros_like(x_c_vp)
y_f_vp = np.zeros_like(y_c_vp)
z_f_vp = np.zeros_like(z_c_vp)

for i in range(len(x_c_vp)):
    temp = np.matmul(R_inv, [x_c_vp[i], y_c_vp[i], z_c_vp[i]]) + P_fc
    x_f_vp[i] = temp[0]
    y_f_vp[i] = temp[1]
    z_f_vp[i] = temp[2]

print(x_f_vp)
print(y_f_vp)
print(z_f_vp)

fig = plt.figure(figsize=(aspect_ratio_vp*5.0, 5.0))
ax = fig.add_subplot(111)
plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
#plt.axis('off')
ax.axes.set_xlim( y_c_vp[0],y_c_vp[2])
ax.axes.set_ylim( -z_c_vp[1], -z_c_vp[0])
#ax.axes.xaxis.set_ticklabels([])
#ax.axes.yaxis.set_ticklabels([])
#ax.set_xticks([])
#ax.set_yticks([])
ax.axes.set_aspect('equal')
fig.canvas.draw()
plt.show()

