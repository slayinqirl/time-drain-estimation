import math 

#tube lengths
L1 = 0.20
L2 = 0.40
L3 = 0.60
Current_tube_length = 0.20


#tube and tank constants
G = 9.81 #m/s^2
D = 0.00794 #m
Area_pipe = math.pi*(D**2)/4
Area_tank=0.0832 #m^2
rel_pipe_roughness=0.00015/0.00794 #relative roughness (epislon/diameter)
kinematic_viscosity= 1.004e-6  #how did she find this value?
#dt= 0.01 #time step
tolerance = 0.000001 # convergence tolerance for friction
max_iter=100 #mAX number of iterations for friction
water_density = 1000 #kg/m^3
K_entrance = 0.5 #entrnace loss coefficient
#K = D/(10+(Current_tube_length*math.sin(1/150))) #entrance loss

#tank heights
Initial_height = 0.08 + 0.02 + (Current_tube_length*math.sin(1/150))
Final_height =  0.02 +  (Current_tube_length*math.sin(1/150))
Average_height = (Initial_height+Final_height)/2

def find_time_to_drain (length):

    l = length

    #Call the function to compute the friction factor
    friction = compute_friction_factor(l)

    print("The area of pipe is:}", Area_pipe, "Area of tank is:", Area_tank)

    print("\nThis is final guessed friction factor", friction)

    #Compute the velocity in the pipe
    numerator = 1 + friction * (length / D) + K_entrance
    denominator = (2 * G)
    time_to_drain = (Area_tank/Area_pipe)*2*(math.sqrt(10+(length*math.sin(1/150)))-math.sqrt(2+(length*math.sin(1/150))))*(math.sqrt(numerator/denominator))
    return time_to_drain



def compute_friction_factor(length):
    l = length
    friction = 0.02
    new_friction = friction

    # Step 1: calculate the velocity in the pipe
    denominator = 1 + (friction * (l / D)) + K_entrance
    vp = math.sqrt((2 * G * Average_height) / denominator)

    print("This is the velocity in the pipe:", vp)

    # Step 2: calculate the Reynolds number
    Re = (vp * D) / kinematic_viscosity
    print("This is the Reynolds number based on the frictoin guess:", Re)

    # Step 3: determine flow type and friction factor
    if Re < 2300:
        new_friction = 64 / Re
        return new_friction
    # translational case
    # elif 
    # turblent case
    else:
        # If you use Colebrook, perform iteration for turbulent flow
        print("This is turbulent flow, we are using Colebrook equation")
        lhs = -2.0 * math.log10((rel_pipe_roughness / (3.7)) + 2.51 / (Re * math.sqrt(friction)))
        new_friction = 1 / (lhs ** 2)

        while abs(new_friction - friction) > tolerance:
            friction = new_friction
            print("Iterating with friction factor:", friction)
            lhs = -2.0 * math.log10((rel_pipe_roughness / (3.7)) + 2.51 / (Re * math.sqrt(friction)))
            new_friction = 1 / (lhs ** 2)

        #Return the final friction factor after convergence
        return new_friction
    

    
time = find_time_to_drain(Current_tube_length)
print(f"The time is{time}")