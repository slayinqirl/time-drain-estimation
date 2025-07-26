import math

# -----------------------------
# Constants
# -----------------------------
g = 9.81  # gravity (m/s^2)
rho = 1000  # density of water (kg/m^3)
mu = 0.001  # dynamic viscosity of water (Pa·s)
D = 0.00794  # diameter of the discharge tube (m)
D_tjoint = 0.0111125 # diameter of the T joint (m)

Ap = math.pi * (D / 2) ** 2  # pipe cross-sectional area 1 (m^2)
A_tjoint = math.pi * (D_tjoint/2)**2 # t joint cross-sectional area 2 (m^2)


At = 0.1  # assumed tank cross-sectional area (m^2) — adjust if known
tol = 1e-6  # tolerance for friction factor convergence
dt = 1.0  # time step in seconds
k = 0.5 # coefficent for contraction
k_tjoint = 2.4 # coefficent for minor loss in T joint
eps = 1.5e-6  # roughness height (m) — assume plastic pipe
from scipy.optimize import fsolve
rel_pipe_roughness = eps/D #relative roughness (epislon/diameter)

# -----------------------------
# Velocity Calculation for pipe 1
# -----------------------------
def compute_v(f, h, L):
    """Calculate outlet velocity using modified Torricelli with friction loss."""
    denominator = ((Ap/2*A_tjoint)**2)*(1+k_tjoint)+(f*L/D)+k
    veloctiy = math.sqrt((2 * g * h) / denominator)
    #print(f"v: {veloctiy}")
    return veloctiy

def colebrook(f, Re, eps, D):
    return 1.0 / math.sqrt(f) + 2.0 * math.log10((eps / D) / 3.7 + 2.51 / (Re * math.sqrt(f)))

# -----------------------------
# Friction Factor Calculation
# -----------------------------
def compute_f_and_re(v):
    """Calculate Reynolds number and corresponding friction factor."""
    Re = (rho * v * D) / mu

    if Re < 2300:
        print("\n This is laminar flow with ",Re)
        f = 64 / Re  # Laminar
    elif 2300 <= Re < 4000:

        print("\n This is transitional flow", Re)
        # Linear interpolation between laminar and turbulent
        f_laminar = 64 / Re

        # Haaland approximation (approximate Colebrook)


        #Newton approximation using fsolve
        f_guess = 0.02   
        f_turbulent = fsolve(colebrook, f_guess, args=(Re, eps, D))
        f_turbulent_float = float(f_turbulent[0])

        # s-jain
        # f_turbulent = 0.25 / (math.log10(rel_pipe_roughness / (3.7) + 5.74 / Re**0.9))**2
        #f_turblent = colebrook_newton(Re, D, rel_pipe_roughness)
        # scale factor between 0 (at Re=2300) and 1 (at Re=4000)
        scale = (Re - 2300) / (4000 - 2300)
        f = f_laminar + scale * (f_turbulent_float- f_laminar)
    else:
        #print("\n This is turbulent flow", Re)
        # Turbulent flow — use approximate Colebrook-White via Haaland's equation

        # Haaland approximation (approximate Colebrook)
        #f = (-1.8 * math.log10((6.9 / Re))) ** -2

        #Newton approximation using fsolve
        f_guess = 0.02   
        f_solution = fsolve(colebrook, f_guess, args=(Re, eps, D))
        f = float(f_solution[0])
        #print(f"Friction factor: {f:.6f}")


        # s-jain
        # f = 0.25 / (math.log10(rel_pipe_roughness / (3.7) + 5.74 / Re**0.9))**2
        # f= colebrook_newton(Re, D, rel_pipe_roughness)

    print(f"This is friction: {f:.6f}")  
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

        for _ in range(10):
            v = compute_v(f, h, L)
            f_new, Re = compute_f_and_re(v)
            # print(f"f_new: {f_new}; Re: {Re}")

            if abs(f_new - f) < tol:
                break
            f = f_new

        #compute tjoint velocity again using pipe and tjoint flow rate
        v_tjoint = (Ap / 2*A_tjoint) * v

        dhdt = -(A_tjoint * v_tjoint * 2)/At
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

for length in [40]:
    results.append(simulate_drainage(length))

# Display results
for r in results:
    print(f"\nTube Length: {r['tube_length_cm']} cm")
    print(f"Time to Drain: {r['time_string']} (or {r['total_seconds']} s)")
    print(f"Average Velocity: {r['avg_velocity']:.4f} m/s")
    print(f"Average Friction Factor: {r['avg_friction_factor']:.5f}")
    print(f"Average Reynolds Number: {r['avg_reynolds']:.0f}")
    print(f"Steps Taken: {r['steps']}")