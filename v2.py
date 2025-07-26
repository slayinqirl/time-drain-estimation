import math
import matplotlib.pyplot as plt

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
k = 0.5 # coefficent for entrance
eps = 1.5e-6  # roughness height (m) — assume plastic pipe
from scipy.optimize import fsolve
rel_pipe_roughness = eps/D #relative roughness (epislon/diameter)

# -----------------------------
# Velocity Calculation
# -----------------------------
def compute_v(f, h, L):
    """Calculate outlet velocity using modified Torricelli with friction loss."""
    veloctiy = math.sqrt((2 * g * h) / (1 + k + ((f * L) / D)))
    #print(f"v: {veloctiy}")
    return veloctiy

def colebrook(f, Re, eps, D):
    return 1.0 / math.sqrt(f) + 2.0 * math.log10((eps / D) / 3.7 + 2.51 / (Re * math.sqrt(f)))

# -----------------------------
# Friction Factor Calculation
# -----------------------------
def compute_f_and_re(v,L):
    """Calculate Reynolds number and corresponding friction factor."""
    Re = (rho * v * D) / mu

    if Re < 2300:
        f = 64 / Re  # Laminar


    elif 2300 <= Re < 4000:

        #print("\n This is transitional flow", Re,  "for length", L)
        # Linear interpolation between laminar and turbulent
        f_laminar = 64 / Re



        #Newton approximation using fsolve
        f_guess = 0.02   
        f_turbulent = fsolve(colebrook, f_guess, args=(Re, eps, D))
        f_turbulent_float = float(f_turbulent[0])


        scale = (Re - 2300) / (4000 - 2300)
        f = f_laminar + scale * (f_turbulent_float- f_laminar)
    else:
        #print("\n This is turbulent flow", Re,  "for length", L)
        # Turbulent flow — use approximate Colebrook-White via Haaland's equation

        #Newton approximation using fsolve
        f_guess = 0.02   
        f_solution = fsolve(colebrook, f_guess, args=(Re, eps, D))
        f = float(f_solution[0])
      


        # s-jain
        # f = 0.25 / (math.log10(rel_pipe_roughness / (3.7) + 5.74 / Re**0.9))**2
        # f= colebrook_newton(Re, D, rel_pipe_roughness)

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

    # To store time and velocity points for plottin
    time_points = []
    velocity_points = []

    while h > hf:
        f = 0.02  # initial guess

        # Lists to store time & velocity for plotting
      

        for _ in range(1000):
            v = compute_v(f, h, L)
            f_new, Re = compute_f_and_re(v,L)
            # print(f"f_new: {f_new}; Re: {Re}")

            if abs(f_new - f) < tol:
                break
            f = f_new

        print(f"For length = {L:2f} t = {t:.1f} s | v = {v:.4f} m/s | Re = {Re:.0f} | f = {f:.5f}")
        
        dhdt = -(Ap / At) * v
        dh = dhdt * dt
        h += dh
        h = max(h, hf)


        # Record time and velocity
        time_points.append(t)
        velocity_points.append(v)

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
        "total_seconds": t,
         "time_points": time_points,
        "velocity_points": velocity_points
    }

# -----------------------------
# Run Simulations
# -----------------------------
results = []

for length in [20, 30, 40, 60]:
    results.append(simulate_drainage(length))
    print("\n")

# Display results
for r in results:
    print(f"\nTube Length: {r['tube_length_cm']} cm")
    print(f"Time to Drain: {r['time_string']} (or {r['total_seconds']} s)")
    print(f"Average Velocity: {r['avg_velocity']:.4f} m/s")
    print(f"Average Friction Factor: {r['avg_friction_factor']:.5f}")
    print(f"Average Reynolds Number: {r['avg_reynolds']:.0f}")
    print(f"Steps Taken: {r['steps']}")


# -----------------------------
# Experimental Data (convert min:sec to total seconds)
# -----------------------------
exp_lengths = [20, 30, 40, 60]
exp_times = [
    3 * 60 + 19,  # 3:19 = 199 seconds
    3 * 60 + 34,  # 3:34 = 214 seconds
    4 * 60 + 26,  # 4:26 = 266 seconds
    4 * 60 + 48   # 4:48 = 288 seconds
]

# Simulation results
sim_lengths = [r['tube_length_cm'] for r in results]
sim_times = [r['total_seconds'] for r in results]

# -----------------------------
# Plotting
# -----------------------------
plt.figure(figsize=(8, 5))
plt.plot(sim_lengths, sim_times, 'o-', label='Simulation', color='blue')
plt.plot(exp_lengths, exp_times, 's--', label="Professor's Data", color='red')

plt.xlabel('Tube Length (cm)')
plt.ylabel('Drain Time (s)')
plt.title('Drain Time vs Tube Length')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



plt.figure(figsize=(8, 5))
for r in results:
    plt.plot(r["time_points"], r["velocity_points"], label=f"L = {r['tube_length_cm']} cm")

plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity vs Time for Different Pipe Lengths")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()