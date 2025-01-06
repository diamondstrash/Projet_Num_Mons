# Library to import
from biblio import *

# Data loading
r_data, q_data = np.loadtxt('q_structure.txt', comments='%', unpack=True)
q_data *= 0.5  # Scale the source term

# Interpolate q(r) for use on the grid
q_interp = e.interp1d(r_data, q_data, fill_value="extrapolate")

# Unit test

x_min = 0
x_max = 0.02
num_points = 500

##################
# constants values
##################

k_acier = 44.5    # W/mK
k_air = 0.05      # W/mK
k_cuivre = 400    # W/mK

T_amb = 20        # °C
h_air = 20        # W/m²·K

nodes_acier = 4097 # number of nodes for core magnetic side (2^m + 1), with m an interger
nodes_air = 1024 # number of nodes for the air side (2^p), with p an interger
nodes_cuivre = 256 # number of nodes for core magnetic side (2^x), with x an interger

n_tot = nodes_acier + nodes_air + nodes_cuivre # notal number of nodes

length_acier = 0.015 # thinkness of core magnetic part
length_air = 0.004 # thinkness of air part
length_coil = 0.001 # thinkness of coil part

total_length = length_acier + length_air + length_coil # total length (in radial direction)

# Generate non-uniform mesh
r_acier = np.linspace(0, length_acier, nodes_acier, endpoint=False)
r_air = np.linspace(length_acier, length_acier + length_air, nodes_air, endpoint=False)
r_cuivre = np.linspace(length_acier + length_air, total_length, nodes_cuivre)
r = np.concatenate([r_acier, r_air, r_cuivre])

k = np.zeros(n_tot) # matrix used to store every thermal conductivity

# Assign thermal conductivity per material
k[:nodes_acier] = k_acier
k[nodes_acier:nodes_acier+nodes_air] = k_air
k[nodes_acier+nodes_air:n_tot] = k_cuivre

# Interpolated source term
q = q_interp(r)

# creation of matrix
T=np.zeros(n_tot)  
b=np.zeros(n_tot)
A=np.zeros((n_tot,n_tot))

# Boundary condition 

#####################################################
# Neumman boundary condition dT/dr = 0 (no flux term)
#####################################################

##########################################################################
# Compute average radius i.e. r_(i+1/2) or r_(i-1/2)
# And also average heat conductivity i.e. r_(i+1/2) or r_(i-1/2)
##########################################################################

r0_plus_half = return_moy(r[0], r[1])

V0 = return_volume(r0_plus_half, 0) # half volume considered

A[0,0] = 1
A[0,1] = -1

b[0] = q[0] * V0 # No external source term at the first node (due to zero flux)

#####################################################
# Fourier boundary condition -k(dT/dr) = h(T_f - T_n) 
#   with T_f, the ambiant air temperature
#   with h the convection air coefficient
#####################################################

##########################################################################
# Compute average radius i.e. r_(i-1/2)
# And also average heat conductivity i.e. k_(i-1/2)
##########################################################################

r_not_minus_half = return_moy(r[n_tot-2], r[n_tot-1])
ki_not_minus_half = return_moy(k[n_tot-2], k[n_tot-1])

##########################################################################
# Compute the 2 surface area to compute :
# S_i and S_{i-1/2}
##########################################################################

S_ntot_1 = return_surface(r[n_tot-1])  # Surface area for convection at the last node (outermost radius)
S_ntot_minus_half = return_surface(r_not_minus_half)  # Surface area for conduction between second-to-last and last node

##########################################################################
# Compute the volume : V_i 
##########################################################################

V_ntot = return_volume(r[n_tot-1], r_not_minus_half) # half of a volume is considered

##########################################################################
# Compute the three coefficient i.e. a_{i,i} and a_{i,i-1} 
# Compute also the source term coefficient i.e. b_{i}
##########################################################################

# Applying the Fourier boundary condition at the last node
dr_last = r[n_tot-1] - r[n_tot-2]  # Radial step for the last node
A[n_tot-1, n_tot-2] = - ( ki_not_minus_half / dr_last ) * S_ntot_minus_half 
A[n_tot-1, n_tot-1] = ( ki_not_minus_half / dr_last ) * S_ntot_minus_half + h_air * S_ntot_1 

# Source term due to convection at the last node
b[n_tot-1] = q[n_tot-1] * V_ntot + h_air * T_amb * S_ntot_1 # At boundary S{i+1/2} = S{i}

for i in range (1,n_tot-1):
    
    ######################################################################################
    # Compute average radius i.e. r_(i+1/2) or r_(i-1/2)
    # And also average heat conductivity i.e. r_(i+1/2) or r_(i-1/2)
    ######################################################################################
    
    ri_plus_half = return_moy(r[i], r[i+1])    
    ri_minus_half = return_moy(r[i-1], r[i])
    
    ki_plus_half = k[i+1]
    ki_minus_half = k[i]
    
    ######################################################################################
    # Compute the radial step (difference in radius) for neighboring nodes
    # dr_minus_one = r[i] - r[i-1] corresponds to the radial step between node i and i-1
    # dr_plus_one = r[i+1] - r[i] corresponds to the radial step between node i and i+1
    ######################################################################################

    dr_minus_one = r[i] - r[i-1]
    dr_plus_one = r[i+1] - r[i]
    
    ############################################################################################################
    # Compute the surface areas at the average radii:
    # Si_plus_half corresponds to the surface area at the average radius (r_(i+1/2)) between nodes i and i+1.
    # Si_minus_half corresponds to the surface area at the average radius (r_(i-1/2)) between nodes i-1 and i.
    ############################################################################################################

    Si_plus_half = return_surface(ri_plus_half)
    Si_minus_half = return_surface(ri_minus_half)

    ############################################################################################################
    # Compute the volumes at the average radii:
    # Vi corresponds to the volume between the average radii r_(i-1/2) and r_(i+1/2).
    ############################################################################################################

    Vi = return_volume(ri_plus_half, ri_minus_half)
        
    ##########################################################################
    # Compute the three coefficient i.e. a_{i,i-1} ; a_{i,i} and a_{i,i+1}
    # Compute also the source term coefficient i.e. b_{i}
    ##########################################################################
    
    A[i,i-1] = - ( ki_minus_half / dr_minus_one ) * Si_minus_half
    A[i,i] = ( ki_minus_half / dr_minus_one ) * Si_minus_half + ( ki_plus_half / dr_plus_one ) * Si_plus_half
    A[i,i+1] = - ( ki_plus_half / dr_plus_one ) * Si_plus_half
    
    b[i] = q[i] * Vi

T = np.linalg.solve(A,b)
plot_temperature_vs_radius(r, T)

# Create a plot to visualize the mesh
plt.figure(figsize=(8, 6))

# Plot the radial positions of the mesh (r)
plt.plot(r, np.zeros_like(r), 'bo', label='Mesh Nodes')

# Highlight different regions with different colors
plt.axvspan(0, length_acier, color='gray', alpha=0.3, label='Acier (Core Magnetic)')
plt.axvspan(length_acier, length_acier + length_air, color='yellow', alpha=0.3, label='Air')
plt.axvspan(length_acier + length_air, total_length, color='red', alpha=0.3, label='Cuivre (Coil)')

# Add labels and title
plt.xlabel('Radial Position (r)')
plt.title('Radial Mesh Distribution')
plt.legend()

# Show the plot
plt.show()