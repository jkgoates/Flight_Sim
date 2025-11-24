import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import scipy as sp
import time
from connection_m import connection

class Camera:
    def __init__(self, dict):
        self.d_vp = dict['view_plane']['depth[ft]']
        self.theta_vp = np.radians(dict['view_plane']['angle[deg]'])
        self.AR = dict['view_plane']['aspect_ratio']

        if 'type' in dict:
            self.type = dict['type']
        else:
            self.type = 'follow'

        if self.type == 'follow':
            self.distance = dict['follow_distance[ft]']

        

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
        self.pos = np.array(kwargs['pos'])
        self.euler = kwargs['euler']
        self.vel = np.array(kwargs['vel'])
        
        self.R = np.zeros((3,3))

        c_euler = np.cos((self.euler))
        s_euler = np.sin((self.euler))
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

        if self.type == 'follow':
            self.pos = self.pos + self.distance*np.matmul(R_inv, self.vel/np.linalg.norm(self.vel))

        self.r_f_vp = np.matmul(R_inv, self.r_c_vp) + self.pos[:,np.newaxis]

        self.P_0 = np.zeros(3)
        self.P_0[0] = np.average(self.r_f_vp[0,:])
        self.P_0[1] = np.average(self.r_f_vp[1,:])
        self.P_0[2] = np.average(self.r_f_vp[2,:])
        self.P_1 = self.r_f_vp[:,3]
        self.P_2 = self.r_f_vp[:,0]
        self.n_vp = np.cross(self.P_1-self.P_0, self.P_2-self.P_0)



class LinesObject:
    def __init__(self, dict, ax):

        self.ax = ax

        # Parse json
        color = dict["color"]
        self.clipping = False

        if "type" in dict:
            self.type = dict["type"]
        else:
            self.type = 'vtk'

        if self.type == 'vtk':
            vtkfile = dict['filename']

            with open(vtkfile, 'r') as reader:
                content = reader.readlines()

            points = []
            lines = []

            started_points = False
            started_lines = False

            for line in content:
                if 'POINTS' in line:
                    started_points = True
                    npoints = int(line.split()[1])
                    continue

                if 'LINES' in line:
                    started_lines = True
                    nlines = int(line.split()[1])
                    continue

                if started_points and not started_lines:
                    coords = line.strip().split()
                    if len(coords) == 3:
                        points.append([float(coords[0]), float(coords[1]), float(coords[2])])

                if started_lines:
                    indices = line.strip().split()
                    if len(indices) == 3 and indices[0] == '2':
                        lines.append([int(indices[1]), int(indices[2])])

            self.points = np.array(points)
            self.lines = np.array(lines)


        if self.type == 'grid':
            self.clipping = True
            ground_altitude = dict['altitude[ft]']
            grid_number = dict['grid_number']
            grid_scale = dict['grid_scale[ft]']

            # Initialize grid
            nxlines = grid_number*2 + 1
            self.nlines = 2*nxlines
            self.points = np.zeros((4*nxlines,3))
            self.lines = np.zeros((2*nxlines,2), dtype= int)


            for i in range(nxlines):
                self.points[2*i, :] = [-grid_number*grid_scale, (i-grid_number)*grid_scale, -ground_altitude]
                self.points[2*i+1, :] = [grid_number*grid_scale, (i-grid_number)*grid_scale, -ground_altitude]
                self.lines[i,0] = 2*i
                self.lines[i,1] = 2*i+1

            for i in range(nxlines):
                self.points[2*i+2*nxlines, :] =   [(i-grid_number)*grid_scale, -grid_number*grid_scale, -ground_altitude]
                self.points[2*i+2*nxlines+1, :] = [(i-grid_number)*grid_scale, grid_number*grid_scale, -ground_altitude]
                self.lines[i+nxlines,0] = 2*i+2*nxlines
                self.lines[i+nxlines,1] = 2*i+2*nxlines+1


        self.npoints = len(self.points)
        self.nlines = len(self.lines)

        self.ax, = ax.plot([],[],ls='-',color=color)

        self.points3D = np.zeros((self.npoints,3))
        self.update([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

        self.I_ca = np.zeros((self.npoints,3))
        self.lamb = np.zeros(self.npoints)
        self.P_b_c = np.zeros((self.npoints,3))

        self.points_vp = np.zeros((self.npoints, 2))
        self.lines_vp = np.full((self.nlines*3,2), None, dtype=object)

    def update(self, pos, euler):

        pos = np.array(pos)
        euler = np.array(euler)

        self.R = np.zeros((3,3))

        c_euler = np.cos((euler))
        s_euler = np.sin((euler))
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

        self.points3D = np.transpose(np.matmul(R_inv, np.transpose(self.points))) + pos[np.newaxis,:]


    def draw(self, cam):

        self.I_ca = self.points3D - cam.pos[np.newaxis,:]

        for i in range(self.npoints):
            self.lamb[i] = np.dot(cam.P_0-cam.pos,cam.n_vp)/np.dot(self.I_ca[i,:],cam.n_vp)
            self.P_b_c[i,:] = self.lamb[i]*self.I_ca[i,:]

        temp = np.transpose(np.matmul((cam.R), np.transpose(self.P_b_c)))
        self.points_vp[:,0] = temp[:,1]
        self.points_vp[:,1] = -temp[:,2]


        for i in range(self.nlines):
            i0 = self.lines[i,0]
            i1 = self.lines[i,1]
            
            # both points are behind camera
            if (self.lamb[i0] < 0.0 and self.lamb[i1] < 0.0):
                self.lines_vp[3*i,   :] = None
                self.lines_vp[3*i+1, :] = None
            # only one behind camera
            elif (self.lamb[i1] < 0.0):
                l10_f = self.points3D[i1,:] - self.points3D[i0,:]
                lamb_temp = np.dot(cam.P_0-self.points3D[i0,:], cam.n_vp)/np.dot(l10_f, cam.n_vp)
                point_f = self.points3D[i0,:] + lamb_temp*l10_f - cam.pos
                temp = np.matmul((cam.R),point_f)
                point = np.zeros(2)
                point[0] = temp[1]
                point[1] = -temp[2]
                self.lines_vp[3*i, :] =   self.points_vp[i0,:]
                self.lines_vp[3*i+1, :] = point
            elif (self.lamb[i0] < 0.0):
                l10_f = self.points3D[i0,:] - self.points3D[i1,:]
                lamb_temp = np.dot(cam.P_0-self.points3D[i1,:], cam.n_vp)/np.dot(l10_f, cam.n_vp)
                point_f = self.points3D[i1,:] + lamb_temp*l10_f -cam.pos
                temp = np.matmul((cam.R),point_f)
                point = np.zeros(2)
                point[0] = temp[1]
                point[1] = -temp[2]
                self.lines_vp[3*i, :] =   point
                self.lines_vp[3*i+1, :] = self.points_vp[i1,:] 
            # Both points in front
            else:
                self.lines_vp[3*i, :] =   self.points_vp[i0,:]
                self.lines_vp[3*i+1, :] = self.points_vp[i1,:]


        self.ax.set_data(self.lines_vp[:,0],self.lines_vp[:,1])

class TickerTape:
    def _init_(self, ax, color, orientation, lim, num, step, pos, width):
        
        pass

    def update(self, val):
        pass


class HUD:
    def __init__(self,json_data,ax,camera):
        color = json_data["color"]
        box_background_color = 'lightgrey'
        dx = camera.dx
        dy = camera.dy

        # Altitude
        self.altitude_minor = TickerTape(ax,color,'vertical',0.4*dy, 10, 100, 0.4*dx, -0.02*dx)
        self.altitude_major = TickerTape(ax,color,'vertical',0.4*dy,1,1000,0.4*dx, -0.05*dx, True, 0.01*dx - 0.02*dy)
        ax.fill([0.4*dx, 0.42*dx, 0.5*dx, 0.5*dx, 0.42*dx, 0.4*dx],[0.0, 0.05*dy, 0.05*dy, -0.05*dy, -0.05*dy, 0.0], facecolor=box_background_color, edgecolor = color, linewidth=1, zorder=100)
        self.altitude_box = ax.text(0.415*dx, -0.02*dy, str("{:0.0f}".format(0.0)), color=color, zorder=101)
        
        # Heading
        self.heading_minor = TickerTape(ax, color, 'horizontal', 0.2*dx, 4,5,-0.48*dy, 0.02*dy)
        self.heading_major = TickerTape(ax, color, 'horizontal', 0.2*dx, 2, 10, -0.48*dy, 0.05*dy, True, -0.03*dx, -0.03*dy)
        ax.fill([0.0, -0.04*dx, -0.04*dx, 0.04*dx, 0.04*dx, 0.0],[-0.48*dy, -0.5*dy, -0.57*dy, -0.57*dy, -0.5*dy, -0.48*dy], facecolor=box_background_color, edgecolor=color, linewidth=1, zorder=100)
        self.heading_box = ax.text(-0.04*dx, -0.55*dy, str("{:0.0f}".format(0.0)), color=color, zorder=101)

    def draw(self, camera, state):
        dx = camera.dx
        dy = camera.dy

        # altitude ticker
        self.altitude_minor.update(-state.location[2])
        self.altitude_major.update(-state.location[2])
        self.altitude_box.set_text(str("{:0.0f}".format(-state.location[2])))

        # heading ticker
        self.heading_minor.update(state.euler[2]*180.0/np.pi)
        self.heading_major.update(state.euler[2]*180.0/np.pi)
        self.heading_box.set_text(str("{:0.0f}".format(state.euler[2]*180.0/np.pi)))



# MATPLOTLIB Events
def on_close(event):
    global run
    run = False

def on_move(event):
    global aileron
    global elevator
    #if event.inaxes:
        #print(f'data coords {event.xdata}, {event.ydata}')

    aileron = -np.radians(30*event.xdata/1.7)
    elevator = np.radians(30*event.ydata/0.85)

def on_keypress(event):
    global throttle
    if (event.key == 'pageup'):
        throttle += 0.02
    elif (event.key == 'pagedown'):
        throttle -= 0.02

    if throttle > 1.0 : 
        throttle = 1.0
    elif throttle < 0.0:
        throttle = 0.0

    print('throttle adjusted to: ', throttle)


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

    print("Starting graphics...")
    print("Initializing Camera...")
    cam = Camera(cam_input)
    print("Done.")
    print("Updating position...")
    cam.update(pos=[0.0, 0.0, 0.0], euler=[0.0, 0.0, 0.0], vel=[1.0, 0.0, 0.0])
    print("Done.")

    graphics_conn = connection(input_dict["connections"]["receive_states"])
    controls_conn = connection(input_dict["connections"]["send_controls"])

    fig = plt.figure(figsize=(cam.AR*5.0, 5.0))
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

    # Set up matplotlib events 
    fig.canvas.mpl_connect('close_event', on_close)
    fig.canvas.mpl_connect('motion_notify_event', on_move)
    fig.canvas.mpl_connect('key_press_event', on_keypress)
    
    ground = LinesObject(input_dict["scene"]["ground"], ax)
    vehicle = LinesObject(input_dict["scene"]["vehicle"], ax)

    vehicle.update([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    run = True

    states = np.zeros(13)
    cnt = 0

    # Initialize controls
    controls = np.zeros(4)
    aileron = 0.0
    rudder = 0.0
    elevator = 0.0
    throttle = 0.5
    

    while run:
        if cnt == 0:
            t_i = time.time()
        # Receive states
        states = graphics_conn.recv()
        vehicle.update(states[7:10], states[10:13])
        cam.update(pos=states[7:10], euler=states[10:13], vel=states[1:4])

        ground.draw(cam)
        vehicle.draw(cam)

        fig.canvas.draw()
        fig.canvas.flush_events()

        controls = [aileron, elevator, rudder, throttle]
        controls_conn.send(controls)

        cnt += 1
        if cnt == 50:
            t_o = time.time()
            fps = 50/(t_o-t_i)
            print("graphics rate [hz] = ", fps)
            #print('controls:', controls)
            cnt = 0




