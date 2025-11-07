import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import scipy as sp

class Camera:
    def __init__(self, **kwargs):
        self.d_vp = kwargs['depth']
        self.theta_vp = kwargs['angle']
        self.AR = kwargs['aspect_ratio']

        self.theta_vp = np.radians(self.theta_vp)
        self.w_vp = 2*self.d_vp*np.tan(self.theta_vp)
        self.h_vp = self.w_vp/self.AR

        self.r_c_vp = np.zeros((3,4))

        self.r_c_vp[0,:] = self.d_vp
        self.r_c_vp[1,:] = 0.5*self.w_vp
        self.r_c_vp[1,0] = -self.r_c_vp[1,0]
        self.r_c_vp[1,1] = -self.r_c_vp[1,1]
        self.r_c_vp[2,:] = 0.5*self.h_vp
        self.r_c_vp[2,0] = -self.r_c_vp[2,0]
        self.r_c_vp[2,3] = -self.r_c_vp[2,3]
                
    
    def update(self, **kwargs):
        self.pos = np.array(kwargs['position'])
        self.euler = kwargs['euler']
        
        R = np.zeros((3,3))

        c_euler = np.cos(np.radians(self.euler))
        s_euler = np.sin(np.radians(self.euler))
        c_phi = c_euler[0]
        c_theta = c_euler[1]
        c_psi = c_euler[2]
        s_phi = s_euler[0]
        s_theta = s_euler[1]
        s_psi = s_euler[2]

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

        self.r_f_vp = np.matmul(R_inv, self.r_c_vp) + self.pos[:,np.newaxis]


if __name__ == '__main__':

    # Check for interactive mode
    if "-i" in sys.argv:
        pass
    else:
        # Get input filename from command line arguments
        input_filename = sys.argv[-1]

        # Check for valid input
        if ".json" not in input_filename:
            raise IOError("Please specify a .json input file (got {0}).".format(input_filename))

    # Load json file
    with open(input_filename) as json_file_handle:
        input_dict = json.load(json_file_handle)

    cam_input = input_dict["camera"]
    

    aspect_ratio_vp = 2.0

    d_vp = 1.0 # ft
    theta_vp = np.deg2rad(60.0)

    print("Starting graphics...")
    print("Initializing Camera...")
    cam = Camera(depth=1.0, angle=60.0, aspect_ratio=2.0)
    print("Done.")
    print("Updating position...")
    cam.update(position=[-20.0, 0.0, -5.0], euler=[5.0, -10.0, 0.0])
    print("Done.")

    print(cam.r_c_vp)
    print(cam.r_f_vp)

    points = np.zeros((3,20))

    points[0,:] = [-10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,-10.0,-5.0,-5.0,0.0,0.0,5.0,5.0,10.0,10.0]
    points[1,:] = [-10.0,-10.0,-5.0,-5.0,0.0,0.0,5.0,5.0,10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,10.0,-10.0,10.0]
    points[2,:] = 0.0

    print(points)

    #P_fc = np.array([-20.0, 0.0, -5.0])

    #euler = np.radians(np.array([5.0, -10.0, 0.0]))

    #cam.update(position=P_fc, euler=euler)

    #c_euler = np.cos(euler)
    #s_euler = np.sin(euler)
    #c_phi = c_euler[0]
    #c_theta = c_euler[1]
    #c_psi = c_euler[2]
    #s_phi = s_euler[0]
    #s_theta = s_euler[1]
    #s_psi = s_euler[2]

    #phi   = np.deg2rad(5.0)
    #theta = np.deg2rad(-10.0)
    #psi   = np.deg2rad(0.0)

    ## Rotation matrix
    #R = np.zeros((3,3))
    ##R[0,0] = np.cos(theta)*np.cos(psi)
    ##R[0,1] = np.cos(theta)*np.sin(psi)
    ##R[0,2] = -np.sin(theta)
    ##R[1,0] = np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi)
    ##R[1,1] = np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi)
    ##R[1,2] = np.sin(phi)*np.cos(theta)
    ##R[2,0] = np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)
    ##R[2,1] = np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi)
    ##R[2,2] = np.cos(phi)*np.cos(theta)

    #R[0,0] = c_theta*c_psi
    #R[0,1] = c_theta*s_psi
    #R[0,2] = -s_theta
    #R[1,0] = s_phi*s_theta*c_psi - c_phi*s_psi
    #R[1,1] = s_phi*s_theta*s_psi + c_phi*c_psi
    #R[1,2] = s_phi*c_theta
    #R[2,0] = c_phi*s_theta*c_psi + s_phi*s_psi
    #R[2,1] = c_phi*s_theta*s_psi - s_phi*c_psi
    #R[2,2] = c_phi*c_theta

    #R_inv = sp.linalg.inv(R)

    #x_f_vp = np.zeros_like(x_c_vp)
    #y_f_vp = np.zeros_like(y_c_vp)
    #z_f_vp = np.zeros_like(z_c_vp)

    #for i in range(len(x_c_vp)):
        #temp = np.matmul(R_inv, np.array([x_c_vp[i], y_c_vp[i], z_c_vp[i]])) + P_fc
        #print(temp)
        #x_f_vp[i] = temp[0]
        #y_f_vp[i] = temp[1]
        #z_f_vp[i] = temp[2]

    fig = plt.figure(figsize=(aspect_ratio_vp*5.0, 5.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
    #plt.axis('off')
    ax.axes.set_xlim( cam.r_c_vp[1,1], cam.r_c_vp[1,2])
    ax.axes.set_ylim( -cam.r_c_vp[2,1], -cam.r_c_vp[2,0])
    #ax.axes.xaxis.set_ticklabels([])
    #ax.axes.yaxis.set_ticklabels([])
    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.axes.set_aspect('equal')
    fig.canvas.draw()
    #plt.show()

