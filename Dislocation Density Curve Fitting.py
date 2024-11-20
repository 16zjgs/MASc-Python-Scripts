#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[13]:


# Dose values in dpa
D = np.array([0, 1, 2])  # Dose in dpa
D_2 = np.array([0, 0.05, 0.5, 5])

# Dislocation density measurements (in m^-2)
# Replace these with your actual data
rho_32000_HPT = np.array([301.3, 50.23, 47.35])
rho_600 = np.array([0.64, 0.8, 1.12, 1.2])


rho_32000_HPT = np.array([301.3, 50.23, 47.35])
rho_SA_HPT = np.array([286.9, 81.1, 77.6])
rho_957_HPT = np.array([180, 19.75, 17.35])
rho_957 = np.array([0.5, 18.35, 22.44])

rho_32000 = np.array([0.01, 2.54, 3.91])


# In[14]:


class Sample:
    def __init__(self, doses, rho_values):
        """
        Initialize the Sample object with doses and rho values.

        Parameters:
        - doses: Array of dose values (dpa).
        - rho_values: Array of dislocation density values corresponding to the doses.
        """
        self.doses = doses
        self.rho_values = rho_values
        
        # Calculate rho_ss and k using curve fitting
        self.rho_ss, self.k = self.fit_dislocation_density()

    def dislocation_density_model(self, D, rho_ss, k, increasing=True):
        """Model for dislocation density as a function of dose."""
        rho_0 = self.rho_values[0]  # Initial dislocation density at D = 0
        
        if increasing:
            # Model for increasing dislocation density
            return rho_0 + (rho_ss - rho_0) * (1 - np.exp(-k * D))
        else:
            # Model for decreasing dislocation density (original model)
            return rho_ss + (rho_0 - rho_ss) * np.exp(-k * D)

    def fit_dislocation_density(self):
        """Fit the dislocation density model to the given data to find rho_ss and k."""
        
        # Determine if the dislocation density is increasing or decreasing
        if self.rho_values[-1] > self.rho_values[0]:
            increasing = True  # Dislocation densities are increasing
        else:
            increasing = False  # Dislocation densities are decreasing
        
        initial_guess = [self.rho_values[-1], 0.5]  # Initial guess for rho_ss and k
        
        # Curve fitting with respect to increasing or decreasing model
        popt, _ = curve_fit(lambda D, rho_ss, k: self.dislocation_density_model(D, rho_ss, k, increasing),
                            self.doses, self.rho_values, p0=initial_guess, maxfev=2000)
        
        rho_ss_fitted, k_fitted = popt
        return rho_ss_fitted, k_fitted

    def plot_dislocation_density(self):
        """Plot the experimental dislocation density data and the fitted model."""
        # Generate doses for plotting the fitted curve
        D_plot = np.linspace(min(self.doses), max(self.doses), 100)
        
        # Check if the dislocation density is increasing or decreasing
        if self.rho_values[-1] > self.rho_values[0]:
            increasing = True
        else:
            increasing = False
        
        # Calculate the fitted dislocation densities
        rho_fitted = self.dislocation_density_model(D_plot, self.rho_ss, self.k, increasing=increasing)

        # Create the plot
        plt.figure(figsize=(8, 6))
        
        # Plot the experimental data (smaller markers)
        plt.plot(self.doses, self.rho_values, 'o', label='Data', markersize=5)

        # Plot the fitted model (less thick line)
        plt.plot(D_plot, rho_fitted, '-', label=f'Fit\n$k={self.k:.2f}$', linewidth=1)

        # Plot the steady-state dislocation density line and add the label for rho_ss
        plt.axhline(y=self.rho_ss, color='r', linestyle='--', linewidth=1)
        plt.text(max(self.doses)*0.8, self.rho_ss * 1.05, f'$\\rho_{{ss}}={self.rho_ss:.2e}$', color='r')
        
        # Remove gridlines
        plt.grid(False)

        # Adjust x and y axis labels
        plt.xlabel('Dose (dpa)')
        plt.ylabel('Dislocation Density (m$^{-2}$)')

        # Set the y-axis to logarithmic scale
        plt.yscale('log')

        # Update the plot title
        plt.title('Dislocation Density vs. Dose')

        # Add a legend
        plt.legend()

        # Show the plot
        plt.show()


# In[15]:


sample_957_HPT = Sample(D, rho_957_HPT)
sample_957 = Sample(D, rho_957)


# In[16]:


sample_32000_HPT = Sample(D, rho_32000_HPT)
sample_32000 = Sample(D, rho_32000)


# In[17]:


samples = [sample_32000_HPT, sample_32000]


# In[18]:


def plot_multiple_samples(samples, title, asymptote_threshold=1e-3):
    """
    Plots multiple samples on the same figure, with data points, fitted models, and steady-state lines.
    
    Parameters:
    - samples: List of Sample objects to plot.
    - title: Title for the plot.
    - asymptote_threshold: A threshold value used to extend the x-axis to show the asymptote.
    """
    
    plt.figure(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(samples)))  # Generate different colors for each sample
    
    for idx, sample in enumerate(samples):
        color = colors[idx]  # Assign a color to each sample

        # Initial dislocation density
        rho_0 = sample.rho_values[0]
        
        # Generate a fine range of dose values for plotting the fitted curve
        max_dose = max(sample.doses)
        extended_doses = np.linspace(0, max_dose * 2, 500)
        
        # Calculate the fitted dislocation densities
        rho_fitted = sample.dislocation_density_model(extended_doses, sample.rho_ss, sample.k)
        
        # Plot the experimental data points
        plt.plot(sample.doses, sample.rho_values, 'o', label=f'Data {idx+1}', color=color, markersize=5)

        # Plot the fitted curve
        plt.plot(extended_doses, rho_fitted, '-', label=f'Fit {idx+1}', color=color, linewidth=1)
        
        # Plot the steady-state dislocation density line
        plt.axhline(y=sample.rho_ss, color=color, linestyle='--', linewidth=1)
        
        # Add rho_ss label to the plot
        plt.text(max_dose * 1.5, sample.rho_ss * 1.05, f'$\\rho_{{ss}}={sample.rho_ss:.2e}$', color=color)

    # Set plot title and labels
    plt.title(title)
    plt.xlabel('Dose (dpa)')
    plt.ylabel('Dislocation Density (m$^{-2}$)')
    plt.yscale('log')  # Log scale for the y-axis
    plt.grid(False)  # No gridlines
    plt.legend()
    
    # Show the plot
    plt.show()


# In[19]:


def plot_multiple_samples(samples, sample_names, title, asymptote_threshold=1e-3):
    """
    Plots multiple samples on the same figure, with data points, fitted models, and steady-state lines.
    
    Parameters:
    - samples: List of Sample objects to plot.
    - sample_names: List of sample variable names (as strings), where each name corresponds to a sample.
    - title: Title for the plot.
    - asymptote_threshold: A threshold value used to extend the x-axis to show the asymptote.
    """
    
    plt.figure(figsize=(16, 9))
    
    # Define specific colors
    colors = ['red', 'blue', 'green']
    
    for idx, sample in enumerate(samples):
        color = colors[idx % len(colors)]  # Cycle through colors if more than 3 samples
        
        # Extract sample name without "_sample"
        sample_name = sample_names[idx].replace('_sample', '')

        # Initial dislocation density
        rho_0 = sample.rho_values[0]
        
        # Generate a fine range of dose values for plotting the fitted curve
        max_dose = max(sample.doses)
        extended_doses = np.linspace(0, max_dose * 2, 500)
        
        # Calculate the fitted dislocation densities
        rho_fitted = sample.dislocation_density_model(extended_doses, sample.rho_ss, sample.k)
        
        # Plot the experimental data points with hollow circles of the same color
        plt.plot(sample.doses, sample.rho_values, 'o', label=f'{sample_name} meas', color=color, markersize=7, markerfacecolor='none')
        
        # Plot the fitted curve
        plt.plot(extended_doses, rho_fitted, '-', label=f'{sample_name} fit', color=color, linewidth=1.5)
        
        # Plot the steady-state dislocation density line
        plt.axhline(y=sample.rho_ss, color=color, linestyle='--', linewidth=1.5)
        
        # Add rho_ss label to the plot
        plt.text(max_dose * 1.5, sample.rho_ss * 1.05, f'$\\rho_{{ss}}={sample.rho_ss:.4f}$', color=color)


    # Set plot title and labels
    plt.title(title)
    plt.xlabel('Dose [dpa]')
    plt.ylabel('Dislocation Density (1e-14) [m$^{-2}$]')
    plt.yscale('log')  # Log scale for the y-axis
    plt.grid(False)  # No gridlines
    plt.legend()
    
    # Show the plot
    plt.show()


# In[20]:


plot_multiple_samples(samples, ['617_32000_HPT', '617_32000'],'Dislocation Density vs dpa', asymptote_threshold=1e-3)


# In[21]:


def simultaneous_fit(sample1, sample2):
    """
    Perform a simultaneous fit on two Sample instances to find a common steady-state dislocation density (rho_ss)
    while allowing for different k values for each dataset (one increasing, one decreasing).

    Parameters:
    - sample1: First instance of the Sample class (usually descending).
    - sample2: Second instance of the Sample class (usually ascending).

    Returns:
    - rho_ss_fitted: Fitted steady-state dislocation density.
    - k1_fitted: Fitted k value for sample1 (descending).
    - k2_fitted: Fitted k value for sample2 (ascending).
    """
    
    def combined_model(D_combined, rho_ss, k1, k2):
        """Combined model for fitting both datasets with a common rho_ss."""
        D1 = D_combined[:len(sample1.doses)]  # Doses for sample1
        D2 = D_combined[len(sample1.doses):]  # Doses for sample2
        
        # Model for sample1 (typically decreasing)
        rho_0_1 = sample1.rho_values[0]
        model1 = rho_ss + (rho_0_1 - rho_ss) * np.exp(-k1 * D1)
        
        # Model for sample2 (typically increasing)
        rho_0_2 = sample2.rho_values[0]
        model2 = rho_0_2 + (rho_ss - rho_0_2) * (1 - np.exp(-k2 * D2))
        
        # Combine the two models
        return np.concatenate((model1, model2))

    # Combine dose and rho values from both samples
    combined_doses = np.concatenate((sample1.doses, sample2.doses))
    combined_rho_values = np.concatenate((sample1.rho_values, sample2.rho_values))
    
    # Initial guesses: use rho_ss as the average of final rho values, k1 and k2 are guesses for the decay constants
    initial_guess = [np.mean([sample1.rho_values[-1], sample2.rho_values[-1]]), 0.5, 0.5]
    
    # Perform the simultaneous curve fit
    popt, _ = curve_fit(combined_model, combined_doses, combined_rho_values, p0=initial_guess, maxfev=5000)
    
    # Extract fitted parameters
    rho_ss_fitted, k1_fitted, k2_fitted = popt
    
    # Print or return the results
    print(f"Fitted steady-state dislocation density (rho_ss): {rho_ss_fitted:.4e}")
    print(f"Fitted k1 (Sample 1 - decreasing): {k1_fitted:.4e}")
    print(f"Fitted k2 (Sample 2 - increasing): {k2_fitted:.4e}")
    
    return rho_ss_fitted, k1_fitted, k2_fitted


# In[22]:


rho_ss_fitted, k1_fitted_, k2_fitted = simultaneous_fit(sample_32000_HPT, sample_32000)


# In[23]:


(20.0174 + 12.714)/2


# In[24]:


def plot_multiple_samples(samples, sample_names, title, rho_ss_fitted, asymptote_threshold=1e-3):
    """
    Plots multiple samples on the same figure, with data points, fitted models, and a common steady-state line.
    
    Parameters:
    - samples: List of Sample objects to plot.
    - sample_names: List of sample variable names (as strings), where each name corresponds to a sample.
    - title: Title for the plot.
    - rho_ss_fitted: The fitted steady-state dislocation density (common rho_ss across samples).
    - asymptote_threshold: A threshold value used to extend the x-axis to show the asymptote.
    """
    
    plt.figure(figsize=(16, 9))
    
    # Define specific colors
    colors = ['red', 'blue', 'green']
    
    for idx, sample in enumerate(samples):
        color = colors[idx % len(colors)]  # Cycle through colors if more than 3 samples
        
        # Extract sample name without "_sample"
        sample_name = sample_names[idx].replace('_sample', '')

        # Initial dislocation density
        rho_0 = sample.rho_values[0]
        
        # Generate a fine range of dose values for plotting the fitted curve
        max_dose = max(sample.doses)
        extended_doses = np.linspace(0, max_dose * 2, 500)
        
        # Calculate the fitted dislocation densities using the provided rho_ss_fitted
        rho_fitted = sample.dislocation_density_model(extended_doses, rho_ss_fitted, sample.k)
        
        # Plot the experimental data points with hollow circles of the same color
        plt.plot(sample.doses, sample.rho_values, 'o', label=f'{sample_name} meas', color=color, markersize=7, markerfacecolor='none')
        
        # Plot the fitted curve
        plt.plot(extended_doses, rho_fitted, '-', label=f'{sample_name} fit', color=color, linewidth=1.5)

    # Plot the common steady-state dislocation density (rho_ss_fitted) as a green dashed line
    plt.axhline(y=rho_ss_fitted, color='green', linestyle='--', linewidth=2, label=f'Common $\\rho_{{ss}}={rho_ss_fitted:.4f}$')
    
    # Set plot title and labels
    plt.title(title)
    plt.xlabel('Dose [dpa]')
    plt.ylabel('Dislocation Density (1e-14) [m$^{-2}$]')
    plt.yscale('log')  # Log scale for the y-axis
    plt.grid(False)  # No gridlines
    plt.legend()
    
    # Show the plot
    plt.show()


# In[25]:


plot_multiple_samples(samples, ['617_32000_HPT', '617_32000'],'Dislocation Density vs dpa', rho_ss_fitted, asymptote_threshold=1e-3)


# # 
