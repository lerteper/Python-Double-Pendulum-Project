import numpy as np
import os
import matplotlib.pyplot as plt

import streamlit as st

fig, ax = plt.subplots()

line, = ax.plot([], [], lw = 2)
ax.set_aspect(1, adjustable='box')




class simple_pendulum:
    def __init__(self, radius, starting_angle, pivot_coords = [0, 0], time = 0, mass = 1):
        self.mass = mass
        self.radius = radius
        self.starting_angle = starting_angle
        self.pivot_coords = pivot_coords
        self.time = time
        self.period = 2 * np.pi * np.sqrt(radius / 9.81)
        self.angle = starting_angle
        self.ang_accel = 0
        self.ang_velocity = 0

    @property
    def bob_coords(self):
        bob_coords = [0,0]
        bob_coords[0] = self.pivot_coords[0] + self.radius * np.sin(self.angle)
        bob_coords[1] = self.pivot_coords[1] - self.radius * np.cos(self.angle)

        return bob_coords
    def set_physics(self, angle, ang_accel, ang_velocity):
        self.angle = angle
        self.ang_accel = ang_accel
        self.ang_velocity = ang_velocity
        return
    def reset_physics(self):
        self.set_physics(self.starting_angle, 0, 0)
        return

#class pendulum_pair:
#    def __init__(self, pendulum1, pendulum2, solutiontype):
#        self.pendulum1 = pendulum1
#        self.pendulum2 = pendulum2
#        self.pair = [pendulum1, pendulum2]
#        self.solutiontype = solutiontype
        
    
def pairing(list):
    if len(list) % 2 == 1:
        return TypeError
    res = []
    for i in range(int(len(list)/2)):
        res.append((list[2*i], list[2*i + 1]))
    return res


st.title("Double Pendulum Simulation")
""
""
""
pendulum_pair_creator = st.columns(2)
with pendulum_pair_creator[0]:
    st.title("Pendulum 1")
    pendulum1_details = {}
    pendulum1_details["radius"] = st.number_input("Length", key = "L1",  min_value = 0.01, value = 1.00)
    pendulum1_details["mass"] = st.number_input("Mass", key = "M1",  min_value = 0.01, value = 1.00)
    pendulum1_details["starting angle"] = st.number_input("Starting Angle", key = "A1")
    if st.checkbox("Radian",True, key = "CB1") == False:
        pendulum1_details["starting angle"] = pendulum1_details["starting angle"]/180 * np.pi
with pendulum_pair_creator[1]:
    st.title("Pendulum 2")
    pendulum2_details = {}
    pendulum2_details["radius"] = st.number_input("Length", key = "L2", min_value = 0.01, value = 1.00)
    pendulum2_details["mass"] = st.number_input("Mass", key = "M2", min_value = 0.01, value = 1.00)
    pendulum2_details["starting angle"] = st.number_input("Starting Angle", key = "A2")
    if st.checkbox("Radian",True, key = "CB2") == False:
        pendulum2_details["starting angle"] = pendulum2_details["starting angle"]/180 * np.pi

if "pendulums" not in st.session_state:
    st.session_state.pendulums = []

def create_pendulum_pair(details1, details2, solutiontype):
    R1 = details1["radius"]
    M1 = details1["mass"]
    A1 = details1["starting angle"]

    pendulum1 = simple_pendulum(radius = R1, mass = M1, starting_angle = A1)
    st.session_state.pendulums.append(pendulum1)

    R2 = details2["radius"]
    M2 = details2["mass"]
    A2 = details2["starting angle"]

    pendulum2 = simple_pendulum(radius = R2, mass = M2, starting_angle = A2)
    st.session_state.pendulums.append(pendulum2)


solutiontype = st.selectbox("", ["Euler's Method", "Rungne-Kratta 4"])
if st.button("Create Pendulum Pair"):
    create_pendulum_pair(pendulum1_details, pendulum2_details, solutiontype)

st.write(st.session_state.pendulums)
lines = []
for _ in range(len(st.session_state.pendulums)): 
    new_line, = ax.plot([], [], lw=2)
    lines.append(new_line)

pendulum_pairs = pairing(st.session_state.pendulums)
line_pairs = pairing(lines)

plotsize = 0
for pendulum_pair in pendulum_pairs:
    max_length = pendulum_pair[0].radius + pendulum_pair[1].radius
    if max_length > plotsize:
        plotsize = max_length
ax.set_xlim(-plotsize, plotsize)
ax.set_ylim(-plotsize, plotsize)


def init():
    for line in lines:
        line.set_data([], [])
    return lines

def ang_accel_solver(pendulum1, pendulum2):
    A = (pendulum1.mass + pendulum2.mass) * pendulum1.radius
    B = pendulum2.mass * pendulum2.radius * np.cos(pendulum1.angle - pendulum2.angle)
    C = - pendulum2.mass * pendulum2.radius * (pendulum2.ang_velocity ** 2) * np.sin(pendulum1.angle - pendulum2.angle) - ((pendulum1.mass + pendulum2.mass) * 9.81 * np.sin(pendulum1.angle))
    
    D = pendulum1.radius * np.cos(pendulum1.angle - pendulum2.angle)
    E = pendulum2.radius
    F = pendulum1.radius * (pendulum1.ang_velocity ** 2) * np.sin(pendulum1.angle - pendulum2.angle) - 9.81 * np.sin(pendulum2.angle)

    LHS = np.array([[A, B],
                    [D, E]])
    RHS = np.array([C, F])
    soln = np.linalg.solve(LHS, RHS)
    return soln

def euler_method(pendulum1, pendulum2):
    soln = ang_accel_solver(pendulum1, pendulum2)
    ang_accel1 = soln[0]
    ang_accel2 = soln[1]

    ang_velocity1 = pendulum1.ang_velocity + pendulum1.ang_accel * 0.02
    ang_velocity2 = pendulum2.ang_velocity + pendulum2.ang_accel * 0.02
    angle1 = pendulum1.angle + pendulum1.ang_velocity * 0.02
    angle2 = pendulum2.angle + pendulum2.ang_velocity * 0.02

    pendulum1.set_physics(angle1, ang_accel1, ang_velocity1)
    pendulum2.set_physics(angle2, ang_accel2, ang_velocity2)
    
    return [[angle1, ang_accel1, ang_velocity1], [angle2, ang_accel2, ang_velocity2]]

def update(frame):
    
    
    for i in range(int(len(pendulum_pairs))):
        
        line1 = line_pairs[i][0]
        line2 = line_pairs[i][1]
        pendulum1 = pendulum_pairs[i][0]
        pendulum2 = pendulum_pairs[i][1]

        euler_method(pendulum1, pendulum2)
        
        pendulum2.pivot_coords = pendulum1.bob_coords
        line1.set_data([pendulum1.pivot_coords[0], pendulum1.bob_coords[0]],[pendulum1.pivot_coords[1], pendulum1.bob_coords[1]])
        line2.set_data([pendulum2.pivot_coords[0], pendulum2.bob_coords[0]],[pendulum2.pivot_coords[1], pendulum2.bob_coords[1]])
    return lines




if len(st.session_state.pendulums) == 0:
    st.info("Create a pendulum pair, then run the simulation.")
    st.stop()          


pendulum_pairs = pairing(st.session_state.pendulums)
init()
update(0)
st.write("update(0) succeeded")
st.write(pendulum_pairs)

 
if st.button("Run simulation", "start_sim"):
    init()
    for pendulum in st.session_state.pendulums:
        pendulum.reset_physics()
    placeholder = st.empty()          
    frame = 0
    stop = st.button("Stop simulation", "stop_sim")

    while stop == False:
        update(frame)
        placeholder.pyplot(fig)
        frame += 1  
    

if st.button('Clear Pendulums', "reset lists"):
    st.session_state.pendulums = []
    pendulum_pairs = []
    lines = []