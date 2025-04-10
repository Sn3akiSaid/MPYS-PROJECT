#!/usr/bin/env python3
"""
Electronic Band Structure 3D Visualization
------------------------------------------
This script creates a 3D visualization of electronic band structures from a data file.
It plots four energy bands as separate surfaces, with options for customization.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

def load_data(filename):
    """
    Load data from file with columns: kx, ky, E1, E2, E3, E4
    where E1-E4 are the energy values for the four bands.
    """
    try:
        data = np.loadtxt(filename)
        kx = data[:, 0]
        ky = data[:, 1]
        bands = [data[:, i] for i in range(2, 6)]  # E1, E2, E3, E4
        return kx, ky, bands
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None, None

def create_grid(kx, ky, resolution=100):
    """
    Create a regular grid for interpolation.
    """
    kx_min, kx_max = kx.min(), kx.max()
    ky_min, ky_max = ky.min(), ky.max()
    
    # Create grid for interpolation
    xi = np.linspace(kx_min, kx_max, resolution)
    yi = np.linspace(ky_min, ky_max, resolution)
    xi_grid, yi_grid = np.meshgrid(xi, yi)
    
    return xi_grid, yi_grid

def interpolate_band(kx, ky, band_data, xi_grid, yi_grid, method='cubic'):
    """
    Interpolate the scattered band data onto a regular grid.
    """
    # Perform interpolation
    zi_grid = griddata((kx, ky), band_data, (xi_grid, yi_grid), method=method)
    return zi_grid

def plot_bands(kx, ky, bands, output_file=None, interpolation_method='cubic', 
               resolution=100, alpha=0.7, colormap='viridis', figsize=(12, 10),
               elev=30, azim=45, scatter_points=True):
    """
    Create a 3D visualization of electronic bands.
    
    Parameters:
    -----------
    kx, ky : arrays
        k-point coordinates
    bands : list of arrays
        Energy values for each band
    output_file : str, optional
        If provided, save the figure to this file
    interpolation_method : str, default 'cubic'
        Interpolation method ('linear', 'cubic', 'nearest')
    resolution : int, default 100
        Resolution of the interpolated grid
    alpha : float, default 0.7
        Transparency of surfaces
    colormap : str, default 'viridis'
        Colormap to use
    figsize : tuple, default (12, 10)
        Figure size
    elev, azim : float, default 30, 45
        Viewing angle elevation and azimuth
    scatter_points : bool, default True
        Whether to also plot the original data points
    """
    # Create grid for interpolation
    xi_grid, yi_grid = create_grid(kx, ky, resolution)
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Color maps for different bands
    colormaps = ['viridis', 'plasma', 'inferno', 'magma']
    
    # Plot each band
    for i, band_data in enumerate(bands):
        # Interpolate data to grid
        zi_grid = interpolate_band(kx, ky, band_data, xi_grid, yi_grid, method=interpolation_method)
        
        # Plot surface
        cmap = plt.get_cmap(colormaps[i % len(colormaps)])
        surf = ax.plot_surface(xi_grid, yi_grid, zi_grid, cmap=cmap, 
                              alpha=alpha, antialiased=True, linewidth=0.5,
                              label=f'Band {i+1}')
        
        # Add original points if requested
        if scatter_points:
            norm = plt.Normalize(band_data.min(), band_data.max())
            ax.scatter(kx, ky, band_data, c=band_data, cmap=cmap, s=20, alpha=0.8, 
                      norm=norm, marker='o', edgecolor='k', linewidth=0.3)
    
    # Set labels and title
    ax.set_xlabel('$k_x$', fontsize=14, labelpad=10)
    ax.set_ylabel('$k_y$', fontsize=14, labelpad=10)
    ax.set_zlabel('Energy (eV)', fontsize=14, labelpad=10)
    ax.set_title('Electronic Band Structure', fontsize=16, pad=20)
    
    # Set viewing angle
    ax.view_init(elev=elev, azim=azim)
    
    # Add a color bar for each band
    for i, band_data in enumerate(bands):
        norm = plt.Normalize(band_data.min(), band_data.max())
        sm = plt.cm.ScalarMappable(cmap=colormaps[i % len(colormaps)], norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.05, aspect=30, shrink=0.7)
        cbar.set_label(f'Band {i+1} Energy (eV)', fontsize=12)
    
    # Enhance the appearance
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save the figure if an output file is specified
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_file}")
    
    return fig, ax

def main():
    # Configuration
    data_file = "NewDataWSMUnperturbed1.dat"  # Update with your file name
    output_file = "band_structure_3d.png"  # Output file name
    
    # Load data
    kx, ky, bands = load_data(data_file)
    
    if kx is None:
        return
    
    # Plot bands
    fig, ax = plot_bands(
        kx, ky, bands,
        output_file=output_file,
        interpolation_method='cubic',
        resolution=100,
        alpha=0.7,
        colormap='viridis',
        figsize=(12, 10),
        elev=30,
        azim=45,
        scatter_points=True
    )
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    main()