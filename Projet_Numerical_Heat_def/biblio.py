import numpy as np
import scipy.interpolate as e
import matplotlib.pyplot as plt
import math
import os
import sys

def return_moy(r_1, r_2):
    """
    Return the average radius, i.e., r_(i+1/2) or r_(i-1/2).
    Return also the average heat conductivity i.e. k_(i+1/2) or k_(i-1/2).

    Parameters:
    r_1 (float): The first radius (e.g., r_i or r_{i-1} or k).
    r_2 (float): The second radius (e.g., r_{i+1} or r_i or).
    

    Returns:
    float: The average radius.
    """
    return (r_1 + r_2) / 2

def return_surface(r):
    """
    Return the surface for a specific radius.
    
    The surface is calculated as:
        Surface = 2 * pi * r
    
    Parameters:
    r (float): The radius at a specific point (e.g., r_i, r_{i+1/2}, r_{i-1/2}).

    Returns:
    float: The computed surface.
    """
    return 2 * math.pi * r

def return_volume(r_after, r_before):
    """
    Calculate the volume using the average of squared radii.

    Parameters:
    r_after (float): Radius after (e.g., r_i+1 or r_i).
    r_before (float): Radius before (e.g., r_i or r_i-1).

    Returns:
    float: The computed volume.
    """
    return math.pi * (r_after**2 - r_before**2)

def plot_temperature_map(r, T, cmap='coolwarm', colorbar_label="Temperature (°C)", title="Structure Temperature Color Map", save_path="Figure/heat_map_temperature.png"):
    """
    Plots a 2D temperature map with specified separating radii and annotations.
    
    Parameters:
    - r (1D array): Array of radial positions.
    - T (1D array): Array of temperature values corresponding to r.
    - separating_radii (list): List of radii to separate regions with dashed lines.
    - cmap (str): Colormap for the temperature plot (default: 'coolwarm').
    - colorbar_label (str): Label for the color bar.
    - title (str): Title for the plot.
    """
    
    # z (1D array): Array of axial positions (length along z-axis).
    z = np.linspace(0, 0.01, 100)  # Length in z, from 0 to 10 mm (100 points)
    
    # Define the radii separating the 3 regions
    separating_radii = [15e-3, 19e-3]  # 15 mm and 19 mm
    
    # Create a 2D matrix for the temperature (assume constant along z)
    T_2D = np.tile(T, (len(z), 1))
    
    # Find the min and max temperature values
    min_temp = np.min(T)
    max_temp = np.max(T)
    
    # Ensure the directory exists
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    # Plot the color map
    plt.figure(figsize=(8, 6))
    mesh = plt.pcolormesh(r, z, T_2D, shading='auto', cmap=cmap)
    
    # Add the color bar
    cbar = plt.colorbar(mesh, label=colorbar_label)
    
    # Add min and max values on the color bar
    cbar.ax.text(0.5, -0.05, f"{min_temp:.2f} °C", ha='center', va='center', transform=cbar.ax.transAxes, color='blue', fontsize=10)
    cbar.ax.text(0.5, 1.05, f"{max_temp:.2f} °C", ha='center', va='center', transform=cbar.ax.transAxes, color='red', fontsize=10)
    
    # Add dashed lines for the regions
    for radius in separating_radii:
        plt.axvline(x=radius, color='black', linestyle='--', linewidth=1.5, label=f'r = {radius:.3f} m')
    
    # Add labels and title
    plt.xlabel("Radius (m)")
    plt.ylabel("Length (m)")
    plt.title(title)
    plt.legend(loc='upper right')  # Legend for the separating lines
    plt.tight_layout()  # Adjust layout for text visibility
    plt.grid(True)
    
    # Save the figure
    plt.savefig(save_path)
    print(f"Figure saved as {save_path}")
    plt.show()

def plot_temperature_vs_radius(r, T, save_path="Figure/temperature_vs_radius.png"):
    """
    Plots the temperature as a function of the radius and saves the plot to a file.

    Parameters:
    r (array): Array of radius values.
    T (array): Array of temperature values corresponding to the radius values.
    save_path (str): The path where the plot should be saved.
    """
    # Ensure the directory exists
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(r, T, label="Temperature (T)", color="red")
    plt.title("Temperature as a Function of Radius")
    plt.xlabel("Radius (r) [m]")
    plt.ylabel("Temperature (T) [°C]")
    plt.legend()
    plt.grid(True)
    
    # Save the figure
    plt.savefig(save_path)
    print(f"Figure saved as {save_path}")
    plt.show()

def return_fourier_number(k, rho, cp, dr, dt):
    """
    Returns the Fourier number (array or float) for multiple/single materials.
    
    Parameters:
    
        k (float or np.ndarray): Thermal conductivity of the material(s) [W/(m·K)].
        rho (float or np.ndarray): Density of the material(s) [kg/m³].
        cp (float or np.ndarray): Specific heat capacity of the material(s) [J/(kg·K)].
        dr (float): Spatial step size (discretization length) [m].
        dt (float or np.ndarray): Time step size (scalar for all materials or array for different materials) [s].
        
    Returns:
        Fo (np.ndarray): Fourier number for each material.
    """
    
    # Ensure dt is an array for element-wise operations
    dt = np.asarray(dt)
    
    # Compute the Fourier number (vectorized)
    Fo = (k / (rho * cp)) * (dt / dr**2)
    
    return Fo

def stability_condition(Fo):
    """
    Checks if the explicit method stability condition is satisfied.
    
    Parameters:
        Fo (np.ndarray): Fourier number array.
        
    Returns:
        bool: True if max(Fo) <= 0.5, otherwise False.
    """
    return np.max(Fo) <= 0.5