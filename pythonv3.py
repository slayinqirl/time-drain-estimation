
dt = 1
# h = 
# final_h =
tolerance = 1e-6


def compute_v():
    pass

def computer_f():
    """
    check if it is laminar, trnasitions, turbulent
        in the if statement calculate f for each case using different equations 
    """
    pass

def simulate_drainage():
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
    pass



results = []

results.append(simulate_drainage())

# jane sent me this code
f = 0.02
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
