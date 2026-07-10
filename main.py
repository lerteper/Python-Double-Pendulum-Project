import numpy as np
import os
import matplotlib.pyplot as plt
import time
import streamlit as st
import pandas as pd

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
    def create_pendulum(details):
        R = details["radius"]
        M = details["mass"]
        A = details["starting angle"]

        pendulum = simple_pendulum(radius = R, mass = M, starting_angle = A)
        return pendulum
    def set_physics(self, angle, ang_accel, ang_velocity):
        self.angle = angle
        self.ang_accel = ang_accel
        self.ang_velocity = ang_velocity
        return
    def get_physics(self):
        physics_data = {
            "angle" : self.angle,
            "ang_accel" : self.ang_accel,
            "ang_velocity" : self.ang_velocity,
            "bob_coords" : self.bob_coords
        }

        return physics_data
    def reset_physics(self):
        self.set_physics(self.starting_angle, 0, 0)
        return
    

class pendulum_pair:
    def __init__(self, pendulum1, pendulum2, solutiontype):
        self.pendulum1 = pendulum1
        self.pendulum2 = pendulum2
        self.pair = [pendulum1, pendulum2]
        self.solutiontype = solutiontype

    def get_physics(self):
        return {"pendulum1" : self.pendulum1.get_physics(), "pendulum2" : self.pendulum2.get_physics()}
    
    def reset_physics(self):
        self.pendulum1.reset_physics()
        self.pendulum2.reset_physics()
    
    def ang_accel_solver(self):
        pendulum1 = self.pendulum1
        pendulum2 = self.pendulum2

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
    
    def euler_method(self, time, ang_accel1, ang_accel2):
        pendulum1 = self.pendulum1
        pendulum2 = self.pendulum2

        
        ang_velocity1 = pendulum1.ang_velocity + ang_accel1 * time
        ang_velocity2 = pendulum2.ang_velocity + ang_accel2 * time
        angle1 = pendulum1.angle + ang_velocity1 * time
        angle2 = pendulum2.angle + ang_velocity2 * time

    
        return [[angle1, ang_velocity1], [angle2, ang_velocity2]]
    
    def step(self, method, time = 0.02):
        
        self.pendulum2.pivot_coords = self.pendulum1.bob_coords
        coords_set_1 = [self.pendulum1.pivot_coords, self.pendulum1.bob_coords]
        coords_set_2 = [self.pendulum2.pivot_coords, self.pendulum2.bob_coords]

        new_ang_accel = self.ang_accel_solver()
        ang_accel1 = new_ang_accel[0]
        ang_accel2 = new_ang_accel[1]


        soln = self.euler_method(time, ang_accel1, ang_accel2)
        angle1 = soln[0][0]
        angle2 = soln[1][0]
        ang_velocity1 = soln[0][1]
        ang_velocity2 = soln[1][1]

        self.pendulum1.set_physics(angle1, ang_accel1, ang_velocity1)
        self.pendulum2.set_physics(angle2, ang_accel2, ang_velocity2)

        return [coords_set_1, coords_set_2]

class pendulum_pair_with_drag(pendulum_pair):
    def __init__(self, pendulum1, pendulum2, solutiontype, k):
        super().__init__(pendulum1, pendulum2, solutiontype)
        self.k = k
        # since D = rho * v^2 * A * Cd / 2,
        # D is proportional to v^2
        # D = kv^2
        # D = k(r * omega)^2
        # ma = k(r * omega)^2
        # a = k((r * omega)^2) / m
        # alpha * r = k((r * omega)^2) / m
        # alpha = k * r * omega^2 / m 
    def apply_drag(self, matrix):
        
        for i, pendulum in enumerate(self.pair):
            matrix[i] -= self.k * pendulum.radius * pendulum.ang_velocity * abs(pendulum.ang_velocity) / pendulum.mass
        return matrix
    def ang_accel_solver(self):
        return self.apply_drag(super().ang_accel_solver())

def pairing(list):
    if len(list) % 2 == 1:
        return TypeError
    res = []
    for i in range(int(len(list)/2)):
        res.append((list[2*i], list[2*i + 1]))
    return res


st.title("Double Pendulum Simulation")
""
"The Double Pendulum System is an unpredictable but deterministic system."
"This simulation assumes that each pendulum is made of a rigid, massless rod and a point mass at the bottom edge of the rod."
"Every frame, the physics is recalculated (using Euler's or RK4) based on the changes in angle, angular velocity and angular acceleration."
"A pendulum's angular acceleration at any point in time is calculated through the simultaneous equations: "
""
""
pendulum_pair_creator = st.columns(2)
with pendulum_pair_creator[0]:
    st.title("Pendulum 1")
    pendulum1_details = {}
    pendulum1_details["radius"] = st.number_input("Length (metres)", key = "L1",  min_value = 0.01, value = 1.00)
    pendulum1_details["mass"] = st.number_input("Mass (kg)", key = "M1",  min_value = 0.01, value = 1.00)
    pendulum1_details["starting angle"] = st.number_input("Starting Angle (degrees/radian)", key = "A1", value = 90.00)
    if st.checkbox("Radian", key = "CB1") == False:
        pendulum1_details["starting angle"] = pendulum1_details["starting angle"]/180 * np.pi
with pendulum_pair_creator[1]:
    st.title("Pendulum 2")
    pendulum2_details = {}
    pendulum2_details["radius"] = st.number_input("Length (metres)", key = "L2", min_value = 0.01, value = 1.00)
    pendulum2_details["mass"] = st.number_input("Mass (kg)", key = "M2", min_value = 0.01, value = 1.00)
    pendulum2_details["starting angle"] = st.number_input("Starting Angle (degrees/radian)", key = "A2", value = 90.00)
    if st.checkbox("Radian", key = "CB2") == False:
        pendulum2_details["starting angle"] = pendulum2_details["starting angle"]/180 * np.pi

if "pendulum_pairs" not in st.session_state:
    st.session_state.pendulum_pairs = []
if "overall_physics_data" not in st.session_state:
    st.session_state.overall_physics_data = []





solutiontype = st.selectbox("", ["Euler's Method", "Rungne-Kutta 4"])
drag = st.checkbox("Air Resistance")
if drag:
    drag_explanation = r"since D = $\rho * v^2 * A * C_d / 2,$" + "\n\n" + r"$D \propto v^2$" + "\n\n" + r"$D = kv^2$" + "\n\n" + r"$D = k(r * \omega)^2$" + "\n\n" + r"$ma = k(r * \omega)^2$" + "\n\n" + r"$a = k((r * \omega)^2) / m$" + "\n\n" + r"$r * \alpha = k((r * \omega)^2) / m$" + "\n\n" + r"$\alpha = k(r * \omega^2) / m$"
    k = st.number_input("Coefficient k", min_value = 0.00, value = 0.10, help = drag_explanation)
if st.button("Create Pendulum Pair"):

    new_pendulum1 = simple_pendulum.create_pendulum(pendulum1_details)
    new_pendulum2 = simple_pendulum.create_pendulum(pendulum2_details)
    if not drag:
        st.session_state.pendulum_pairs.append(pendulum_pair(new_pendulum1, new_pendulum2, "euler"))
    else:
        st.session_state.pendulum_pairs.append(pendulum_pair_with_drag(new_pendulum1, new_pendulum2, "euler", k))

st.write(st.session_state.pendulum_pairs)
lines = []
for _ in range(2 * len(st.session_state.pendulum_pairs)): 
    new_line, = ax.plot([], [], lw=2)
    lines.append(new_line)

line_pairs = pairing(lines)

plotsize = 0
for pair in st.session_state.pendulum_pairs:
    max_length = pair.pendulum1.radius + pair.pendulum2.radius
    if max_length > plotsize:
        plotsize = max_length
ax.set_xlim(-plotsize, plotsize)
ax.set_ylim(-plotsize, plotsize)


def init():
    for line in lines:
        line.set_data([], [])
    return lines


def update(frame, seconds_per_frame):
    
    physics_data = []

    for i, pair in enumerate(st.session_state.pendulum_pairs):
        
        line1 = line_pairs[i][0]
        line2 = line_pairs[i][1]

        
        line_data = pair.step("euler", time = seconds_per_frame)

        ## line_data[0] : [[bob_coord_x_1, bob_coord_y_1], [pivot_coord_x_1, pivot_coord_y_1]]
        ## line_data[1] : [[bob_coord_x_2, bob_coord_y_2], [pivot_coord_x_2, pivot_coord_y_2]]
        
        physics_data.append(pair.get_physics())
        line1.set_data([line_data[0][0][0], line_data[0][1][0]],[line_data[0][0][1], line_data[0][1][1]])
        line2.set_data([line_data[1][0][0], line_data[1][1][0]],[line_data[1][0][1], line_data[1][1][1]])
    
    return physics_data




if len(st.session_state.pendulum_pairs) == 0:
    st.info("Create a pendulum pair, then run the simulation.")
    st.stop()          




seconds_per_frame = st.number_input("Seconds per frame", 0.01, 1.00, 0.02, help = "Physics calculations run every frame, Seconds per frame determines the time interval between each calculation.")

table_placeholder = st.empty() 

if st.button("Run simulation", "start_sim"):
    init()
    for pair in st.session_state.pendulum_pairs:
        pair.reset_physics()
    placeholder = st.empty() 
    
    frame = 0
    stop = st.button("Stop simulation", "stop_sim")

    while stop == False:
        physics_data = update(frame, seconds_per_frame)
        placeholder.pyplot(fig)
        frame += 1  
        st.session_state.overall_physics_data.append(physics_data)
        time.sleep(seconds_per_frame)
    
table_placeholder.dataframe(st.session_state.overall_physics_data[0][0])
if st.button('Clear Pendulums', "reset lists"):
    st.session_state.pendulum_pairs = []
    st.session_state.overall_physics_data = []
    lines = []
    st.rerun()