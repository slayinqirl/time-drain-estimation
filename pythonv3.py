
# dt = 1
# h = 
# final_h =
# tolerance = 1e-6


# def compute_v():
    # pass

# def computer_f():
"""
check if it is laminar, trnasitions, turbulent
    in the if statement calculate f for each case using different equations 
"""
    # pass

# def simulate_drainage():
"""
iterate f values for each incremental height

while h > hf
    f = 0.02 -> initial guess

    for i in range(1000)
        v = computer_v(f, h, L)
        f_new, Re = compute_Re+f(v)

        if abs(new_f - f) < tol:
            break
        f = new_f

    dhdt = (Ap / At) * v
    h = dhdt * dt
    h = max(h, hf)

    time += dt -> in seconds
    increments += 1
    veolcity_sum += v
    f_sum += f
    Re_sum += Re

avrage_v = veolcity_sum/steps
average_f_sum = f_sum/steps
average_Re_sum = Re_sum/steps

timestring = time but converted into minutes and second....

return time_string, avrage_v, average_f_sum, average_Re_sum
        
"""
    # pass

import math

# -----------------------------
# Constants
# -----------------------------
g = 9.81  # gravity (m/s^2)
rho = 1000  # density of water (kg/m^3)
mu = 0.001  # dynamic viscosity of water (Pa·s)
D = 0.00794  # diameter of the discharge tube (m)
Ap = math.pi * (D / 2) ** 2  # pipe cross-sectional area (m^2)
At = 0.1  # assumed tank cross-sectional area (m^2) — adjust if known
tol = 1e-6  # tolerance for friction factor convergence
dt = 1.0  # time step in seconds
# hf = 0  # final water height (m)
# h0 = 0.08  # initial water height (m) — 8 cm
k = 0.3 # coefficent for contraction
rel_pipe_roughness = 0.00015/0.00794 #relative roughness (epislon/diameter)

# -----------------------------
# Velocity Calculation
# -----------------------------
def compute_v(f, h, L):
    """Calculate outlet velocity using modified Torricelli with friction loss."""
    veloctiy = math.sqrt((2 * g * h) / (1 + k + f * L / D))
    print(f"v: {veloctiy}")
    return veloctiy

def colebrook_newton(Re, D, eps, tol=1e-6, max_iter=1000):
    """
    Solves the Colebrook equation for turbulent flow using Newton-Raphson method.
    Parameters:
        Re      : Reynolds number (dimensionless)
        D       : pipe diameter (m)
        eps     : roughness height (m)
        tol     : tolerance for convergence
        max_iter: maximum number of iterations
    Returns:
        f       : Darcy friction factor
    """

    if Re <= 0:
        raise ValueError("Reynolds number must be positive.")

    # Initial guess
    f = 0.02

    for i in range(max_iter):
        sqrt_f = math.sqrt(f)
        lhs = 1 / sqrt_f
        rhs = -2.0 * math.log10((eps / D) / 3.7 + 2.51 / (Re * sqrt_f))
        func = lhs - rhs

        # Derivative of the function wrt f (from chain rule)
        dfdx = (-0.5 / f ** 1.5) + \
               (2.0 * 2.51) / (math.log(10) * Re * ((eps / D) / 3.7 + 2.51 / (Re * sqrt_f)) * f ** 1.5)

        # Newton-Raphson update
        f_new = f - func / dfdx

        if abs(f_new - f) < tol:
            return f_new

        f = f_new

    # If it did not converge
    raise RuntimeError("Colebrook-Newton solver did not converge.")

# -----------------------------
# Friction Factor Calculation
# -----------------------------
def compute_f_and_re(v):
    """Calculate Reynolds number and corresponding friction factor."""
    Re = (rho * v * D) / mu

    if Re < 2300:
        print("\n This is laminar flow")
        f = 64 / Re  # Laminar
    elif 2300 <= Re < 4000:

        print("\n This is transitional flow")
        # Linear interpolation between laminar and turbulent
        f_laminar = 64 / Re

        # Haaland approximation (approximate Colebrook)
        # f_turblent = (-1.8 * math.log10((6.9 / Re))) ** -2

        # s-jain
        # f_turbulent = 0.25 / (math.log10(rel_pipe_roughness / (3.7) + 5.74 / Re**0.9))**2
        f_turblent = colebrook_newton(Re, D, rel_pipe_roughness)
        # scale factor between 0 (at Re=2300) and 1 (at Re=4000)
        scale = (Re - 2300) / (4000 - 2300)
        f = f_laminar + scale * (f_turblent - f_laminar)
    else:
        print("\n This is turbulent flow")
        # Turbulent flow — use approximate Colebrook-White via Haaland's equation

        # Haaland approximation (approximate Colebrook)
        # f = (-1.8 * math.log10((6.9 / Re))) ** -2

        # s-jain
        # f = 0.25 / (math.log10(rel_pipe_roughness / (3.7) + 5.74 / Re**0.9))**2
        f = colebrook_newton(Re, D, rel_pipe_roughness)
    return f, Re

# -----------------------------
# Simulation Function
# -----------------------------
def simulate_drainage(L_cm):
    L = L_cm / 100  # convert to meters
    hf = 0.02 + L/150  # final water height (m)
    h0 = 0.1 + L/150  # initial water height (m) — 8 cm
    h = h0
    t = 0

    # To calculate averages later
    velocity_sum = 0
    f_sum = 0
    Re_sum = 0
    steps = 0

    while h > hf:
        f = 0.02  # initial guess

        for _ in range(1000):
            v = compute_v(f, h, L)
            f_new, Re = compute_f_and_re(v)
            # print(f"f_new: {f_new}; Re: {Re}")

            if abs(f_new - f) < tol:
                break
            f = f_new

        dhdt = -(Ap / At) * v
        dh = dhdt * dt
        h += dh
        h = max(h, hf)

        t += dt
        velocity_sum += v
        f_sum += f
        Re_sum += Re
        steps += 1

    avg_v = velocity_sum / steps
    avg_f = f_sum / steps
    avg_Re = Re_sum / steps

    minutes = int(t // 60)
    seconds = int(t % 60)
    time_string = f"{minutes}:{seconds:02d}"

    return {
        "tube_length_cm": L_cm,
        "time_string": time_string,
        "avg_velocity": avg_v,
        "avg_friction_factor": avg_f,
        "avg_reynolds": avg_Re,
        "steps": steps,
        "total_seconds": t
    }

# -----------------------------
# Run Simulations
# -----------------------------
results = []

for length in [20, 30, 40, 60]:
    results.append(simulate_drainage(length))

# Display results
for r in results:
    print(f"\nTube Length: {r['tube_length_cm']} cm")
    print(f"Time to Drain: {r['time_string']} (or {r['total_seconds']} s)")
    print(f"Average Velocity: {r['avg_velocity']:.4f} m/s")
    print(f"Average Friction Factor: {r['avg_friction_factor']:.5f}")
    print(f"Average Reynolds Number: {r['avg_reynolds']:.0f}")
    print(f"Steps Taken: {r['steps']}")




# jane sent me this code
# f = 0.02
# for  in range(max_iter):
#     # Compute pipe velocity from Bernoulli with current f
#     denominator = 1 + f * (L / D) + K_sum
#     vp = math.sqrt((2 * g * H) / denominator)

#     # Compute Reynolds number
#     Re = (rho * vp * D) / mu

#     # Update friction factor (laminar vs turbulent)
#     if Re < 2000:
#         f_new = 64 / Re
#     else:
#         # Colebrook iteration
#         try:
#             lhs = -2.0 * math.log10((epsilon / D) / 3.7 + 2.51 / (Re * math.sqrt(f)))
#             f_new = 1 / (lhs ** 2)
#         except:
#             f_new = f  # fallback if log breaks

#     if abs(f_new - f) < 1e-6:
#         break
#     f = f_new

# # 2. Use converged f and vp to compute dH/dt
# dHdt = - (A_pipe / A_tank) * vp

# # 3. Update height
# H += dHdt * dt

# # 4. Increment time
# t += dt
