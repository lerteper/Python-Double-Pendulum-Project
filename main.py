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
            "bob_x_coords" : self.bob_coords[0],
            "bob_y_coords" : self.bob_coords[1]
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
        self.name = f"Pendulum Pair: {pendulum1.mass} kg, {pendulum1.radius} m, {pendulum1.starting_angle} rad, {pendulum2.mass} kg, {pendulum2.radius} m, {pendulum2.starting_angle} rad "

    def get_physics(self):
        res = {"pendulum1" : self.pendulum1.get_physics(), "pendulum2" : self.pendulum2.get_physics()}
        energies = self.energies()
        res["pendulum1"]["KE"] = energies["pendulum1"]["KE"]
        res["pendulum2"]["KE"] = energies["pendulum2"]["KE"]
        res["pendulum1"]["GPE"] = energies["pendulum1"]["GPE"]
        res["pendulum2"]["GPE"] = energies["pendulum2"]["GPE"]
        return res
    def reset_physics(self):
        self.pendulum1.reset_physics()
        self.pendulum2.reset_physics()
    
    def ang_accel_solver(self, m1, m2, r1, r2, w1, w2, th1, th2):
        #pendulum1 = self.pendulum1
        #pendulum2 = self.pendulum2
        #m1, m2 = pendulum1.mass, pendulum2.mass
        #r1, r2 = pendulum1.radius, pendulum2.radius
        #w1, w2 = pendulum1.ang_velocity, pendulum2.ang_velocity
        #th1, TH2 = pendulum1.angle, pendulum2.angle



        A = (m1 + m2) * r1
        B = m2 * r2 * np.cos(th1 - th2)
        C = - m2 * r2 * (w2 ** 2) * np.sin(th1 - th2) - ((m1 + m2) * 9.81 * np.sin(th1))
    
        D = r1 * np.cos(th1 - th2)
        E = r2
        F = r1 * (w1 ** 2) * np.sin(th1 - th2) - 9.81 * np.sin(th2)

        LHS = np.array([[A, B],
                    [D, E]])
        RHS = np.array([C, F])
        soln = np.linalg.solve(LHS, RHS)
        return soln    
    
    def symplectic_euler(self, time):
        pendulum1 = self.pendulum1
        pendulum2 = self.pendulum2

        new_ang_accel = self.ang_accel_solver(self.pendulum1.mass, self.pendulum2.mass, self.pendulum1.radius, self.pendulum2.radius, self.pendulum1.ang_velocity, self.pendulum2.ang_velocity, self.pendulum1.angle, self.pendulum2.angle)
        ang_accel1 = new_ang_accel[0]
        ang_accel2 = new_ang_accel[1]
        
        ang_velocity1 = pendulum1.ang_velocity + ang_accel1 * time
        ang_velocity2 = pendulum2.ang_velocity + ang_accel2 * time
        angle1 = pendulum1.angle + ang_velocity1 * time
        angle2 = pendulum2.angle + ang_velocity2 * time

    
        return [[angle1, ang_velocity1, ang_accel1], [angle2, ang_velocity2, ang_accel2]]
    
    def derivatives(self, state):
        th1, w1, th2, w2 = state
        a1, a2 = self.ang_accel_solver(
            self.pendulum1.mass, self.pendulum2.mass,
            self.pendulum1.radius, self.pendulum2.radius,
            w1, w2, th1, th2
        )
        return np.array([w1, a1, w2, a2])

    def runge_kutta_4(self, time):
        h = time
        y = np.array([self.pendulum1.angle, self.pendulum1.ang_velocity,
                      self.pendulum2.angle, self.pendulum2.ang_velocity])

        k1 = self.derivatives(y)
        k2 = self.derivatives(y + 0.5 * h * k1)
        k3 = self.derivatives(y + 0.5 * h * k2)
        k4 = self.derivatives(y + h * k3)

        th1, w1, th2, w2 = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        
        return [[th1, w1, k1[1]], [th2, w2, k1[3]]]
    def step(self, method, time = 0.02):
        
        self.pendulum2.pivot_coords = self.pendulum1.bob_coords
        coords_set_1 = [self.pendulum1.pivot_coords, self.pendulum1.bob_coords]
        coords_set_2 = [self.pendulum2.pivot_coords, self.pendulum2.bob_coords]


        #m1, m2, r1, r2, w1, w2, th1, th2 =
        

        if method == "Symplectic Euler":
            soln = self.symplectic_euler(time)
            angle1 = soln[0][0]
            angle2 = soln[1][0]
            ang_velocity1 = soln[0][1]
            ang_velocity2 = soln[1][1]
            ang_accel1 = soln[0][2]
            ang_accel2 = soln[1][2]
        if method == "Runge-Kutta 4":
            soln = self.runge_kutta_4(time)
            angle1 = soln[0][0]
            angle2 = soln[1][0]
            ang_velocity1 = soln[0][1]
            ang_velocity2 = soln[1][1]
            ang_accel1 = soln[0][2]
            ang_accel2 = soln[1][2]

        self.pendulum1.set_physics(angle1, ang_accel1, ang_velocity1)
        self.pendulum2.set_physics(angle2, ang_accel2, ang_velocity2)

        return [coords_set_1, coords_set_2]
    
    def energies(self):
        pendulum1 = self.pendulum1
        pendulum2 = self.pendulum2

        KE1 = 0.5 * pendulum1.mass * (pendulum1.radius * pendulum1.ang_velocity) ** 2
        GPE1 = - pendulum1.mass * 9.81 * pendulum1.radius * np.cos(pendulum1.angle)
        KE2 = 0.5 * pendulum2.mass * (((pendulum1.radius * pendulum1.ang_velocity) ** 2) + ((pendulum2.radius * pendulum2.ang_velocity) ** 2) + 2 * pendulum1.radius * pendulum2.radius * pendulum1.ang_velocity * pendulum2.ang_velocity * np.cos(pendulum1.angle - pendulum2.angle))
        GPE2 = - pendulum2.mass * 9.81 * (pendulum1.radius * np.cos(pendulum1.angle) + pendulum2.radius * np.cos(pendulum2.angle))
        return {"pendulum1" : {"KE" : KE1, "GPE" : GPE1}, "pendulum2" : {"KE" : KE2, "GPE" : GPE2}}

class pendulum_pair_with_drag(pendulum_pair):
    def __init__(self, pendulum1, pendulum2, solutiontype, k):
        super().__init__(pendulum1, pendulum2, solutiontype)
        self.k = k
        self.name = f"Pendulum Pair (k = {k}): {pendulum1.mass} kg, {pendulum1.radius} m, {pendulum1.starting_angle} rad, {pendulum2.mass} kg, {pendulum2.radius} m, {pendulum2.starting_angle} rad "
        # since D = rho * v^2 * A * Cd / 2,
        # D is proportional to v^2
        # D = kv^2
        # D = k(r * omega)^2
        # ma = k(r * omega)^2
        # a = k((r * omega)^2) / m
        # alpha * r = k((r * omega)^2) / m
        # alpha = k * r * omega^2 / m 
    def apply_drag(self, matrix, w1, w2):
        velocities = [w1, w2]
        for i, pendulum in enumerate(self.pair):
            matrix[i] -= self.k * pendulum.radius * velocities[i] * abs(velocities[i]) / pendulum.mass
        return matrix
    def ang_accel_solver(self, m1, m2, r1, r2, w1, w2, th1, th2):
        return self.apply_drag(super().ang_accel_solver(m1, m2, r1, r2, w1, w2, th1, th2), w1, w2)

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
"Every frame, the physics is recalculated based on the changes in angle, angular velocity and angular acceleration."
"A double pendulum's angular accelerations at any point in time is calculated through the simultaneous equations: "
st.latex(r"(m_1+m_2)L_1\ddot\theta_1 + m_2 L_2\ddot\theta_2\cos\Delta + m_2 L_2\dot\theta_2^{\,2}\sin\Delta + (m_1+m_2)g\sin\theta_1 = 0")
st.latex(r"L_2\ddot\theta_2 + L_1\ddot\theta_1\cos\Delta - L_1\dot\theta_1^{\,2}\sin\Delta + g\sin\theta_2 = 0")
st.latex(r"\Delta = \theta_1 - \theta_2")
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
    pendulum2_details["starting angle"] = st.number_input("Starting Angle (degrees/radian)", key = "A2", value = 150.00)
    if st.checkbox("Radian", key = "CB2") == False:
        pendulum2_details["starting angle"] = pendulum2_details["starting angle"]/180 * np.pi

if "pendulum_pairs" not in st.session_state:
    st.session_state.pendulum_pairs = []
if "overall_physics_data" not in st.session_state:
    st.session_state.overall_physics_data = []





solutiontype = st.selectbox("", ["Symplectic Euler", "Runge-Kutta 4"], index = 1)
drag = st.checkbox("Air Resistance")
if drag:
    drag_explanation = r"since D = $\rho * v^2 * A * C_d / 2,$" + "\n\n" + r"$D \propto v^2$" + "\n\n" + r"$D = kv^2$" + "\n\n" + r"$D = k(r * \omega)^2$" + "\n\n" + r"$ma = k(r * \omega)^2$" + "\n\n" + r"$a = k((r * \omega)^2) / m$" + "\n\n" + r"$r * \alpha = k((r * \omega)^2) / m$" + "\n\n" + r"$\alpha = k(r * \omega^2) / m$"
    k = st.number_input("Coefficient k", min_value = 0.00, value = 0.10, help = drag_explanation)
if st.button("Create Pendulum Pair"):

    new_pendulum1 = simple_pendulum.create_pendulum(pendulum1_details)
    new_pendulum2 = simple_pendulum.create_pendulum(pendulum2_details)
    if not drag:
        st.session_state.pendulum_pairs.append(pendulum_pair(new_pendulum1, new_pendulum2, solutiontype))
    else:
        st.session_state.pendulum_pairs.append(pendulum_pair_with_drag(new_pendulum1, new_pendulum2, solutiontype, k))

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

        
        line_data = pair.step(pair.solutiontype, time = seconds_per_frame)

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



if st.button("Run simulation", "start_sim"):
    init()
    st.session_state.overall_physics_data = []
    for pair in st.session_state.pendulum_pairs:
        pair.reset_physics()
    placeholder = st.empty() 
    
    frame = 0
    stop = st.button("Stop simulation", "stop_sim")
    while stop == False:
        start_time = time.perf_counter()
        physics_data = update(frame, seconds_per_frame)
        placeholder.pyplot(fig)
        frame += 1  
        st.session_state.overall_physics_data.append(physics_data)
        #print(str(time.perf_counter() - start_time))
        time.sleep(max(seconds_per_frame - time.perf_counter() - start_time, 0.0))


table_placeholder, graph_placeholder = [st.empty().container() for _ in range(2)]

pair_selection_options = {}
for number, pair in enumerate(st.session_state.pendulum_pairs):
    pair_selection_options[pair.name] = number

if len(st.session_state.overall_physics_data) > 0 and len(st.session_state.pendulum_pairs) > 0:
    table_placeholder.title("Analysis")
    pair_selection = table_placeholder.selectbox("Select pair to view physics data", pair_selection_options)
    pair_index = pair_selection_options[pair_selection]
    
    table_data = []
    for frame, frame_data in enumerate(st.session_state.overall_physics_data):
        table_data.append([])
        pair_data = frame_data[pair_index]
        table_data[frame].append((frame + 1) * seconds_per_frame)
        for value in pair_data["pendulum1"].values():
            table_data[frame].append(value)
        for value in pair_data["pendulum2"].values():
            table_data[frame].append(value)
    dataframe_column_names = ["time elapsed"]
    for key in pair_data["pendulum1"].keys():
        dataframe_column_names.append(key + "1")
    for key in pair_data["pendulum2"].keys():
        dataframe_column_names.append(key + "2")
    refined_data = pd.DataFrame(data = table_data, columns = dataframe_column_names)
    
    refined_data["totalKE"] = refined_data["KE1"] + refined_data["KE2"]
    refined_data["totalGPE"] = refined_data["GPE1"] + refined_data["GPE2"]
    refined_data["total energy"] = refined_data["totalKE"] + refined_data["totalGPE"]
    table_placeholder.dataframe(refined_data)
    
    with graph_placeholder:
        graph_axis = st.selectbox("x-axis",refined_data.columns), st.selectbox("y-axis", refined_data.columns)
        st.scatter_chart(refined_data, x = graph_axis[0], y = graph_axis[1])

if st.button('Clear Pendulums', "reset lists"):
    st.session_state.pendulum_pairs = []
    st.session_state.overall_physics_data = []
    lines = []
    st.rerun()