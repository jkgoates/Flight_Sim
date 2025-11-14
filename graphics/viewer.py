import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import scipy as sp
import time
from connection_m import connection

class Camera:
    def __init__(self, **kwargs):
        self.d_vp = kwargs['depth']
        self.theta_vp = kwargs['angle']
        self.AR = kwargs['aspect_ratio']

    #self.theta_vp = np.radians(self.theta_vp)
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

        print(self.r_c_vp)

    
    def update(self, **kwargs):
        self.pos = np.array(kwargs['position'])
        self.euler = kwargs['euler']
        
        self.R = np.zeros((3,3))

        c_euler = np.cos(np.radians(self.euler))
        s_euler = np.sin(np.radians(self.euler))
        c_phi = c_euler[0]
        c_theta = c_euler[1]
        c_psi = c_euler[2]
        s_phi = s_euler[0]
        s_theta = s_euler[1]
        s_psi = s_euler[2]

        self.R[0,0] = c_theta*c_psi
        self.R[0,1] = c_theta*s_psi
        self.R[0,2] = -s_theta
        self.R[1,0] = s_phi*s_theta*c_psi - c_phi*s_psi
        self.R[1,1] = s_phi*s_theta*s_psi + c_phi*c_psi
        self.R[1,2] = s_phi*c_theta
        self.R[2,0] = c_phi*s_theta*c_psi + s_phi*s_psi
        self.R[2,1] = c_phi*s_theta*s_psi - s_phi*c_psi
        self.R[2,2] = c_phi*c_theta

        R_inv = sp.linalg.inv(self.R)

        self.r_f_vp = np.matmul(R_inv, self.r_c_vp) + self.pos[:,np.newaxis]

        self.P_0 = np.zeros(3)
        self.P_0[0] = np.average(self.r_f_vp[0,:])
        self.P_0[1] = np.average(self.r_f_vp[1,:])
        self.P_0[2] = np.average(self.r_f_vp[2,:])
        self.P_1 = self.r_f_vp[:,3]
        self.P_2 = self.r_f_vp[:,0]
        self.n_vp = np.cross(self.P_1-self.P_0, self.P_2-self.P_0)



class Grid:
    def __init__(self, dict, ax):

        self.ax = ax

        # Parse json
        color = dict["color"]
        ground_altitude = dict['altitude[ft]']
        grid_number = dict['grid_number']
        grid_scale = dict['grid_scale[ft]']

        # Initialize grid
        nxlines = grid_number*2 + 1
        self.nlines = 2*nxlines
        self.points_f = np.zeros((4*nxlines,3))
        self.lines_f = np.zeros((2*nxlines,2), dtype= int)


        for i in range(nxlines):
            self.points_f[2*i, :] = [-grid_number*grid_scale, (i-grid_number)*grid_scale, -ground_altitude]
            self.points_f[2*i+1, :] = [grid_number*grid_scale, (i-grid_number)*grid_scale, -ground_altitude]
            self.lines_f[i,0] = 2*i
            self.lines_f[i,1] = 2*i+1

        for i in range(nxlines):
            self.points_f[2*i+2*nxlines, :] =   [(i-grid_number)*grid_scale, -grid_number*grid_scale, -ground_altitude]
            self.points_f[2*i+2*nxlines+1, :] = [(i-grid_number)*grid_scale, grid_number*grid_scale, -ground_altitude]
            self.lines_f[i+nxlines,0] = 2*i+2*nxlines
            self.lines_f[i+nxlines,1] = 2*i+2*nxlines+1

        self.npoints = len(self.points_f)
        self.nlines = len(self.lines_f)


        self.ax, = ax.plot([],[],ls='-',color=color)

        self.I_ca = np.zeros((self.npoints,3))
        self.lamb = np.zeros(self.npoints)
        self.P_b_c = np.zeros((self.npoints,3))

        self.points_vp = np.zeros((self.npoints, 2))
        self.lines_vp = np.full((self.nlines*3,2), None, dtype=object)

    def draw(self, cam):

        self.I_ca = self.points_f - cam.pos[np.newaxis,:]

        for i in range(self.npoints):
            self.lamb[i] = np.dot(cam.P_0-cam.pos,cam.n_vp)/np.dot(self.I_ca[i,:],cam.n_vp)
            self.P_b_c[i,:] = self.lamb[i]*self.I_ca[i,:]

        temp = np.transpose(np.matmul((cam.R), np.transpose(self.P_b_c)))
        self.points_vp[:,0] = temp[:,1]
        self.points_vp[:,1] = -temp[:,2]


        for i in range(self.nlines):
            i0 = self.lines_f[i,0]
            i1 = self.lines_f[i,1]
            
            # both points are behind camera
            if (self.lamb[i0] < 0.0 and self.lamb[i1] < 0.0):
                self.lines_vp[3*i,   :] = None
                self.lines_vp[3*i+1, :] = None
            # only one behind camera
            elif (self.lamb[i1] < 0.0):
                l10_f = self.points_f[i1,:] - self.points_f[i0,:]
                lamb_temp = np.dot(cam.P_0-self.points_f[i0,:], cam.n_vp)/np.dot(l10_f, cam.n_vp)
                point_f = self.points_f[i0,:] + lamb_temp*l10_f - cam.pos
                temp = np.matmul((cam.R),point_f)
                point = np.zeros(2)
                point[0] = temp[1]
                point[1] = -temp[2]
                self.lines_vp[3*i, :] =   self.points_vp[i0,:]
                self.lines_vp[3*i+1, :] = point
            elif (self.lamb[i0] < 0.0):
                l10_f = self.points_f[i0,:] - self.points_f[i1,:]
                lamb_temp = np.dot(cam.P_0-self.points_f[i1,:], cam.n_vp)/np.dot(l10_f, cam.n_vp)
                point_f = self.points_f[i1,:] + lamb_temp*l10_f -cam.pos
                temp = np.matmul((cam.R),point_f)
                point = np.zeros(2)
                point[0] = temp[1]
                point[1] = -temp[2]
                self.lines_vp[3*i, :] =   point
                self.lines_vp[3*i+1, :] = self.points_vp[i1,:] 
            else:
                self.lines_vp[3*i, :] =   self.points_vp[i0,:]
                self.lines_vp[3*i+1, :] = self.points_vp[i1,:]

        self.ax.set_data(self.lines_vp[:,0],self.lines_vp[:,1])


def on_close(event):
    global run
    run = False

def on_move(event):
    global aileron
    global elevator
    if event.inaxes:
        print(f'data coords {event.xdata}, {event.ydata}')


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


    graphics_conn = connection(input_dict["connections"]["receive_states"])


    d_vp = input_dict['camera']['view_plane']['distance[ft]']
    aspect_ratio_vp = input_dict['camera']['view_plane']['aspect_ratio']
    theta_vp = np.deg2rad(input_dict['camera']['view_plane']['angle[deg]'])

    position = input_dict['camera']['location[ft]']
    euler = input_dict['camera']['orientation[deg]']

    print("Starting graphics...")
    print("Initializing Camera...")
    cam = Camera(depth=d_vp, angle=theta_vp, aspect_ratio=aspect_ratio_vp)
    print("Done.")
    print("Updating position...")
    cam.update(position=position, euler=euler)
    print("Done.")


    fig = plt.figure(figsize=(aspect_ratio_vp*5.0, 5.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
    plt.axis('off')
    ax.axes.set_xlim( cam.r_c_vp[1,1], cam.r_c_vp[1,2])
    ax.axes.set_ylim( -cam.r_c_vp[2,1], -cam.r_c_vp[2,0])
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axes.set_aspect('equal')
    fig.canvas.draw() 
    plt.show(block=False)

    # Set up close_event
    fig.canvas.mpl_connect('close_event', on_close)
    fig.canvas.mpl_connect('motion_notify_event', on_move)
    
    ground = Grid(input_dict["scene"]["ground"], ax)


    run = True

    states = np.zeros(7)
    cnt = 0

    while run:
        if cnt == 0:
            t_i = time.time()
        # Receive states
        states = graphics_conn.recv()
        cam.update(position=states[1:4], euler=np.degrees(states[4:7]))

        ground.draw(cam)

        fig.canvas.draw()
        fig.canvas.flush_events()

        cnt += 1
        if cnt == 100:
            t_o = time.time()
            fps = 100/(t_o-t_i)
            print("graphics rate [hz] = ", fps)
            cnt = 0

        #run = False



