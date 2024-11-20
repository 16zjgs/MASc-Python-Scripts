#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import math
from scipy.optimize import minimize
from sklearn.metrics import r2_score
import mystic as my
from mystic.solvers import diffev2
from mystic.monitors import VerboseMonitor
from scipy.optimize import minimize, differential_evolution
from lmfit import Minimizer, Parameters
import sys
import io
import xlsxwriter
from matplotlib.ticker import MultipleLocator
from scipy.spatial import distance


# In[2]:


class Sample:
    def __init__(self, file_path, crystal_structure, sample_name, a, C_h00, q_avg, a_i, b_i, c_i):
        
        
        
        self.file_path = file_path
        self.crystal_structure = crystal_structure
        self.sample_name = sample_name
        self.a = a
        self.C_h00 = C_h00
        self.q_avg = q_avg
        self.a_i = a_i
        self.b_i = b_i
        self.c_i = c_i

        # Initialize attributes for optimal values, means, stds, and SEs
        self.eps_test = []
        self.rho_test = []
        self.q_test = []
        self.eps_mean = None
        self.eps_std = None
        self.rho_mean = None
        self.rho_std = None
        self.q_mean = None
        self.q_std = None
        self.eps_se = None
        self.rho_se = None
        self.q_se = None

        # Call the methods to process the data
        self.k_delta_k = self.obtain_k_delta_k()
        if self.k_delta_k[0] is not None:  # Check if k_delta_k was properly obtained
            self.peaks = self.determine_peaks()
            self.k_delta_k_hkl = self.append_hkl()
            self.k_delta_k_hkl_C = self.append_c_hkl()
            self.k_delta_k_hkl_C_isub = self.subtract_instrumental()
        
        # Store the line fit coefficients
        self.line_fit_coefficients = self.calculate_mWH_coefficients()
        
        self.perform_optimizations()
        self.perform_bootstrap_se()
    

    def obtain_k_delta_k(self):
        centers = []
        fwhms = []

        if os.path.isfile(self.file_path):
            try:
                with open(self.file_path, 'r') as file:
                    for line in file:
                        # Process lines that start with '%_' and contain 'Pearson7'
                        if line.startswith('%_') and "Pearson7" in line:
                            # Split the line using '\t' as the delimiter
                            columns = line[5:].split('\t')

                            # Extract specific values: center (2nd value) and fwhm (5th value)
                            try:
                                center = float(columns[1])  # Center value (second column)
                                fwhm = float(columns[4])    # FWHM value (fifth column)
                                centers.append(center)
                                fwhms.append(fwhm)
                            except ValueError:
                                # Skip lines that don't contain valid float data
                                continue

                if centers and fwhms:
                    return np.array(centers), np.array(fwhms)
                else:
                    print("No valid data found.")
                    return None, None

            except Exception as e:
                print(f"Error processing file: {e}")
        else:
            print(f"The file {self.file_path} does not exist.")

        return None, None

    def extract_information(self):
        centers = []
        fwhms = []

        if os.path.isfile(self.file_path):
            with open(self.file_path, 'r') as file:
                lines = file.readlines()

                for line in lines:
                    # Process lines starting with '%_' and containing 'Pearson7'
                    if line.startswith('%_') and "Pearson7" in line:
                        # Remove '%_' and split by '\t'
                        columns = line[5:].split('\t')

                        try:
                            # Extract center and fwhm values from the split data
                            center = float(columns[1])  # 2nd column
                            fwhm = float(columns[4])    # 5th column
                            centers.append(center)
                            fwhms.append(fwhm)
                        except ValueError:
                            continue

            if centers and fwhms:
                return centers, fwhms
            else:
                return None, None

        else:
            print(f"The file {self.file_path} does not exist.")
            return None, None
        
    def determine_peaks(self):
        if self.crystal_structure == 'fcc':
            # removing 222 peak from both because they give weird values
            hkl_list = ['111', '200', '220', '311', '222' ,'400', '331', '420', '422', '440']
        elif self.crystal_structure == 'bcc':
            hkl_list = ['110', '200', '211', '220', '310', '222', '321', '400', '420', '332']
        else:
            raise ValueError("Invalid crystal structure specified. Choose 'fcc' or 'bcc'.")
        
        K_list = [((int(i[0])**2 + int(i[1])**2 + int(i[2])**2)**0.5) / self.a for i in hkl_list]
        return [K_list, hkl_list]

    def append_hkl(self, err=0.08):
        theoretical_peaks, hkl_list = self.peaks
        k_delta_k = [[self.k_delta_k[0][i], self.k_delta_k[1][i]] for i in range(len(self.k_delta_k[0]))]

        hkl_indices = []
        for i in range(len(k_delta_k)):
            center = k_delta_k[i][0]
            match_found = False
            for theoretical_peak, hkl in zip(theoretical_peaks, hkl_list):
                if abs(center - theoretical_peak) <= err:
                    hkl_indices.append(hkl)
                    match_found = True
                    break
            if not match_found:
                hkl_indices.append('nan')

        filtered_centers = []
        filtered_fwhms = []
        filtered_hkl_indices = []
        for i in range(len(k_delta_k)):
            if hkl_indices[i] != 'nan':
                filtered_centers.append(k_delta_k[i][0])
                filtered_fwhms.append(k_delta_k[i][1])
                filtered_hkl_indices.append(hkl_indices[i])
        return [filtered_centers, filtered_fwhms, filtered_hkl_indices]

    def big_h(self, h, k, l):
        denominator = h**2 + k**2 + l**2
        if denominator == 0:
            return 0
        return ((h**2 * k**2) + (h**2 * k**2) + (k**2 * l**2)) / (denominator ** 2)

    def append_c_hkl(self):
        hkl = self.k_delta_k_hkl[2]
        contrast_factors = []
        for i in hkl:
            h, k, l = int(i[0]), int(i[1]), int(i[2])
            H = self.big_h(h, k, l)
            C = self.C_h00 * (1 - self.q_avg * H**2)
            contrast_factors.append(C)
        self.k_delta_k_hkl.append(contrast_factors)
        return self.k_delta_k_hkl

    def subtract_instrumental(self):
        k = self.k_delta_k_hkl_C[0]
        delta_k = self.k_delta_k_hkl_C[1]
        hkl = self.k_delta_k_hkl_C[2]
        C = self.k_delta_k_hkl_C[3]

        ib = [self.a_i * j**2 + self.b_i * j + self.c_i for j in k]

        delta_k_new = [(delta_k[j]**2 - ib[j]**2)**0.5 for j in range(len(delta_k))]
        return [k, delta_k_new, hkl, C]
    
    def plot_WH(self):
        # k values and delta_k values from the data after instrumental subtraction
        k = self.k_delta_k_hkl_C_isub[0]
        delta_k = self.k_delta_k_hkl_C_isub[1]
        
        # Create the scatter plot
        plt.scatter(k, delta_k, color='g')  # 'g' for green, can be changed
        plt.grid(True)
        plt.xlabel('k')
        plt.ylabel('delta_k')
        plt.title(f'WH Plot {self.sample_name}')
        plt.show()
        
    def calculate_mWH_coefficients(self):
        k = self.k_delta_k_hkl_C_isub[0]
        delta_k = self.k_delta_k_hkl_C_isub[1]
        c_hkl = self.k_delta_k_hkl_C_isub[3]

        # Calculate k * sqrt(C_hkl)
        rootC = [i ** 0.5 for i in c_hkl]
        k_rootC = [k[i] * rootC[i] for i in range(len(k))]

        # Perform linear fit (slope and intercept)
        slope, intercept = np.polyfit(k_rootC, delta_k, 1)

        # Return the tuple (slope, intercept)
        return slope, intercept
    
    def plot_mWH(self):
        k = self.k_delta_k_hkl_C_isub[0]
        delta_k = self.k_delta_k_hkl_C_isub[1]
        hkl = self.k_delta_k_hkl_C_isub[2]
        c_hkl = self.k_delta_k_hkl_C_isub[3]

        # Calculate k * sqrt(C_hkl)
        rootC = [i ** 0.5 for i in c_hkl]
        k_rootC = [k[i] * rootC[i] for i in range(len(k))]

        # Plot the measured data
        plt.scatter(k_rootC, delta_k, color='b', label='Measured data')
        plt.grid(True)

        # Plot the line of best fit using the slope and intercept
        slope, intercept = self.line_fit_coefficients
        best_fit_line = [slope * x + intercept for x in k_rootC]
        plt.plot(k_rootC, best_fit_line, 'r--', label=f'Best fit line: y = {slope:.2f}x + {intercept:.2f}')

        # Add labels and title
        plt.xlabel('kC^0.5')
        plt.ylabel('delta_k')
        plt.title(f'mWH Plot {self.sample_name}')
        plt.legend()

        # Display the plot
        plt.show()
    
        
    def delta_k(self, eps, big_m, b_, rho, k_, c_hkl):
        return (0.9 / eps) + ((np.pi * (big_m**2) * (b_**2) / 2) ** 0.5) * ((rho / 10**18) ** 0.5) * k_ * (c_hkl ** 0.5)

    # Optimization function
    def optimization_function_v2(self, params, data):
        eps = params['eps'].value
        rho = params['rho'].value
        q = params['q'].value

        hkl = data[2]
        c_hkl = data[3]

        obj_result = 0
        for i in range(len(data[0])):
            if self.crystal_structure == 'bcc':
                obj_result += (data[1][i] - self.delta_k(eps, 1, self.a * 3**0.5 / 2, rho, data[0][i], c_hkl[i])) ** 2
            else:
                obj_result += (data[1][i] - self.delta_k(eps, 1, self.a * 2**0.5 / 2, rho, data[0][i], c_hkl[i])) ** 2

        return obj_result

    def find_optimal_values(self, data):
        params = Parameters()
        params.add('eps', value=75, min=0, max=1000)
        params.add('rho', value=1e14, min=1e12, max=1e17)
        params.add('q', value=self.q_avg, min=1.416, max=2.3)

        minimizer = Minimizer(self.optimization_function_v2, params, fcn_args=(data,))
        
        original_stdout = sys.stdout
        sys.stdout = io.StringIO()
        
        try:
            result_nm = minimizer.minimize(method='Nelder-Mead')
            if result_nm.success:
                result_de = minimizer.minimize(method='differential_evolution', params=result_nm.params)
                if result_de.success:
                    optimal_eps = result_de.params['eps'].value
                    optimal_rho = result_de.params['rho'].value
                    optimal_q = result_de.params['q'].value
                    return optimal_eps, optimal_rho / 10**14, optimal_q
                else:
                    return None, None, None
            else:
                return None, None, None
        finally:
            sys.stdout = original_stdout

    # Bootstrapping method for standard error
    def bootstrap_se(self, data, n_bootstrap=100):
        eps_values = []
        rho_values = []
        q_values = []

        for i in range(n_bootstrap):
            random_indices = np.random.randint(0, len(data[0]), len(data[0]))
            resampled_data = [[data[j][i] for i in random_indices] for j in range(4)]

            result = self.find_optimal_values(resampled_data)

            if result[0] is not None:
                eps_values.append(result[0])
                rho_values.append(result[1])
                q_values.append(result[2])

        # Filter out None values
        rho_values = [x for x in rho_values if x is not None]
        q_values = [x for x in q_values if x is not None]
        eps_values = [x for x in eps_values if x is not None]

        # Calculate standard errors
        rho_se = np.std(rho_values) if rho_values else "Not Calculated"
        q_se = np.std(q_values) if q_values else "Not Calculated"
        eps_se = np.std(eps_values) if eps_values else "Not Calculated"

        return (eps_se / np.mean(eps_values)) * 100 if eps_values else "Not Calculated", \
               (rho_se / np.mean(rho_values)) * 100 if rho_values else "Not Calculated", \
               (q_se / np.mean(q_values)) * 100 if q_values else "Not Calculated"
    
    # Method to calculate optimal values and their statistics (means and stds)
    def perform_optimizations(self, iterations=10):
        # Reset the test lists
        self.eps_test = []
        self.rho_test = []
        self.q_test = []

        # Perform multiple iterations of the optimization
        for _ in range(iterations):
            eps, rho, q = self.find_optimal_values(self.k_delta_k_hkl_C_isub)
            if eps is not None:
                self.eps_test.append(eps)
                self.rho_test.append(rho)
                self.q_test.append(q)

        # Calculate the means and standard deviations
        self.eps_mean = np.mean(self.eps_test)
        self.eps_std = np.std(self.eps_test)
        self.rho_mean = np.mean(self.rho_test)
        self.rho_std = np.std(self.rho_test)
        self.q_mean = np.mean(self.q_test)
        self.q_std = np.std(self.q_test)

    # Bootstrapping method for standard error
    def perform_bootstrap_se(self, n_bootstrap=75):
        eps_values = []
        rho_values = []
        q_values = []

        for i in range(n_bootstrap):
            random_indices = np.random.randint(0, len(self.k_delta_k_hkl_C_isub[0]), len(self.k_delta_k_hkl_C_isub[0]))
            resampled_data = [[self.k_delta_k_hkl_C_isub[j][i] for i in random_indices] for j in range(4)]

            result = self.find_optimal_values(resampled_data)

            if result[0] is not None:
                eps_values.append(result[0])
                rho_values.append(result[1])
                q_values.append(result[2])

        # Filter out None values
        rho_values = [x for x in rho_values if x is not None]
        q_values = [x for x in q_values if x is not None]
        eps_values = [x for x in eps_values if x is not None]

        # Calculate standard errors
        self.eps_se = np.std(eps_values) if eps_values else "Not Calculated"
        self.rho_se = np.std(rho_values) if rho_values else "Not Calculated"
        self.q_se = np.std(q_values) if q_values else "Not Calculated"

        # Convert standard errors to percentages of the mean optimized values
        self.eps_se = (self.eps_se / np.mean(eps_values)) * 100 if eps_values else "Not Calculated"
        self.rho_se = (self.rho_se / np.mean(rho_values)) * 100 if rho_values else "Not Calculated"
        self.q_se = (self.q_se / np.mean(q_values)) * 100 if q_values else "Not Calculated"


# In[38]:


def find_peak_files(root_dir):
    peak_files = []

    # Walk through the directory tree
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # If the current directory has no subdirectories (base level)
        if not dirnames:
            # Loop through the filenames and check for '_peaks'
            for file in filenames:
                if '_peaks' in file or '_params' in file:
                    peak_files.append(os.path.join(dirpath, file))
    
    return peak_files


# In[39]:


root_dir = "C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\Winter 2024\\New LXRD Data\\Analysis v3"


# In[40]:


peak_files = find_peak_files(root_dir)


# In[41]:


peak_files


# In[7]:


a_617 = 0.359
a_957 = 0.289
a_754 = 0.355

C_617 = 0.2646
C_957 = 0.2895
C_754 = 0.2781

q_617 = (1.416+2.3)/2
q_957 = (1.278+2.6)/2
q_754 = (1.04+1.95)/2


# In[42]:


input_params = [['fcc', '617_32000_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_32000_HPT_1dpa', a_617, C_617, q_617, 0.0001, -0.002, 0.0154],
                ['fcc', '617_32000_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_32000_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_SA_HPT_1dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', 'MA754_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_c_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_e_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_m_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['bcc', 'MA957_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_2dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_1dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_c_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_e_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_m_0dpa', a_957, C_957, q_957, 0,0,0],
               ]
                
## SWITCHED 957 hpt 1 AND 2 DPA


# In[43]:


input_params = [['fcc', '617_32000_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_32000_HPT_1dpa', a_617, C_617, q_617, 0.0001, -0.002, 0.0154],
                ['fcc', '617_32000_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_32000_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_SA_HPT_1dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', 'MA754_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_c_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_e_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_m_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['bcc', 'MA957_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_2dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_1dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_c_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_e_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_m_0dpa', a_957, C_957, q_957, 0,0,0],
               ]
                
## SWITCHED 957 hpt 1 AND 2 DPA


# In[87]:


input_params = [['fcc', '617_32000_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_32000_HPT_1dpa', a_617, C_617, q_617, 0.0001, -0.002, 0.0154],
                ['fcc', '617_32000_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_32000_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_32000_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_0dpa', a_617, C_617, q_617, 0, 0, 0],
                ['fcc', '617_SA_HPT_1dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_2dpa', a_617, C_617, q_617, -0.00007, 0.0008, 0.0085],
                ['fcc', '617_SA_HPT_c_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_m_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', '617_SA_HPT_e_0dpa', a_617, C_617, q_617, 0,0,0],
                ['fcc', 'MA754_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_c_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_e_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['fcc', 'MA754_HPT_m_0dpa', a_754, C_754, q_754, 0, 0, 0],
                ['bcc', 'MA957_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_1dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_2dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_1dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_2dpa', a_957, C_957, q_957, -0.00007, 0.0008, 0.0085],
                ['bcc', 'MA957_HPT_c_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_e_0dpa', a_957, C_957, q_957, 0,0,0],
                ['bcc', 'MA957_HPT_m_0dpa', a_957, C_957, q_957, 0,0,0],
               ]
                
## SWITCHED 957 hpt 1 AND 2 DPA


# In[88]:


peak_files[19]


# In[89]:


sample_objects = []


# In[90]:


for i in range(len(peak_files)):
    sample_objects.append(Sample(peak_files[i], input_params[i][0], input_params[i][1], input_params[i][2], input_params[i][3], input_params[i][4], input_params[i][5], input_params[i][6], input_params[i][7]))


# In[91]:


def export_mWH_data(samples, root_directory):
    # Ensure the root directory exists
    if not os.path.exists(root_directory):
        os.makedirs(root_directory)

    # Create lists to store data for all samples
    coefficients_list = []
    sample_statistics_list = []

    for sample in samples:
        # Prepare the data for this sample
        k = sample.k_delta_k_hkl_C_isub[0]
        delta_k = sample.k_delta_k_hkl_C_isub[1]
        hkl = sample.k_delta_k_hkl_C_isub[2]
        c_hkl = sample.k_delta_k_hkl_C_isub[3]
        slope, intercept = sample.line_fit_coefficients

        # Calculate k * sqrt(C_hkl) for the x-axis
        rootC = [i ** 0.5 for i in c_hkl]
        k_rootC = [k[i] * rootC[i] for i in range(len(k))]

        # Prepare data in the format for saving to a CSV
        data_dict = {
            'k (nm^-1)': k,
            'delta_k (nm^-1)': delta_k,
            'hkl': hkl,
            'C_hkl': c_hkl,
            'k * sqrt(C_hkl)': k_rootC
        }

        # Convert the data into a pandas DataFrame
        df = pd.DataFrame(data_dict)

        # Create the filename for the CSV file
        csv_filename = f"{sample.sample_name}_mWH_data.csv"
        csv_path = os.path.join(root_directory, csv_filename)

        # Write the data to a CSV file
        df.to_csv(csv_path, index=False)

        # Store the coefficients in the coefficients list
        coefficients_list.append({
            'sample_name': sample.sample_name,
            'slope': slope,
            'intercept': intercept
        })

        # Collect sample statistics for the new CSV file
        sample_statistics_list.append({
            'sample_name': sample.sample_name,
            'eps_mean': sample.eps_mean,
            'eps_se': sample.eps_se,
            'rho_mean': sample.rho_mean,
            'rho_se': sample.rho_se,
            'q_mean': sample.q_mean,
            'q_se': sample.q_se
        })

    # Write the coefficients to a CSV file in the root directory
    coefficients_df = pd.DataFrame(coefficients_list)
    coefficients_filename = os.path.join(root_directory, "mWH_coefficients.csv")
    coefficients_df.to_csv(coefficients_filename, index=False)

    # Optionally, you can also write the coefficients to a text file
    coefficients_txt_path = os.path.join(root_directory, "mWH_coefficients.txt")
    with open(coefficients_txt_path, 'w') as f:
        for coeff in coefficients_list:
            f.write(f"Sample: {coeff['sample_name']}, Slope: {coeff['slope']}, Intercept: {coeff['intercept']}\n")

    # Write the sample statistics to a CSV file
    statistics_df = pd.DataFrame(sample_statistics_list)
    statistics_filename = os.path.join(root_directory, "mWH_statistics.csv")
    statistics_df.to_csv(statistics_filename, index=False)

    print(f"Data and statistics exported to {root_directory} successfully.")


# In[92]:


mWH_data_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\MASc Thesis\\Results\\Results for paper\\mWH plots for paper\\New 2'


# In[93]:


export_mWH_data(sample_objects, mWH_data_dir)


# In[94]:


def plot_mWH_samples(samples, output_directory, plot_title):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Colors for different data sets
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    
    # Collect all sample names to use in the filename
    sample_names = "_".join([sample.sample_name for sample in samples])

    # Create a filename based on concatenated sample names
    plot_filename = os.path.join(output_directory, f'{sample_names}_mWH.png')

    # Create the figure for the plot
    plt.figure(figsize=(16, 9))  # Set figure size to 16:9

    # Initialize variables to track global max values for x and y
    global_max_x = 0
    global_max_y = 0

    # This will store labeled points to prevent repeated hkl labeling
    labeled_points = []

    # Define a small error threshold to check if points are close to each other
    error_threshold = 0.05

    # First loop through to find the global max values for k_rootC and delta_k
    for sample in samples:
        k = sample.k_delta_k_hkl_C_isub[0]
        delta_k = sample.k_delta_k_hkl_C_isub[1]
        c_hkl = sample.k_delta_k_hkl_C_isub[3]

        # Calculate k * sqrt(C_hkl)
        rootC = [i ** 0.5 for i in c_hkl]
        k_rootC = [k[i] * rootC[i] for i in range(len(k))]

        # Update global maximums
        global_max_x = max(global_max_x, max(k_rootC))
        global_max_y = max(global_max_y, max(delta_k))

    # Add a little buffer (e.g., 10%) to the max values
    global_max_x += global_max_x * 0.1
    global_max_y += global_max_y * 0.1

    # Loop through each sample object to plot
    for idx, sample in enumerate(samples):
        # Prepare the data for this sample
        k = sample.k_delta_k_hkl_C_isub[0]
        delta_k = sample.k_delta_k_hkl_C_isub[1]
        hkl = sample.k_delta_k_hkl_C_isub[2]  # hkl indices for labeling
        c_hkl = sample.k_delta_k_hkl_C_isub[3]

        # Calculate k * sqrt(C_hkl)
        rootC = [i ** 0.5 for i in c_hkl]
        k_rootC = [k[i] * rootC[i] for i in range(len(k))]

        # Plot the data (Remove the word "data" from the label)
        plt.scatter(k_rootC, delta_k, color=colors[idx % len(colors)], s=100, label=f'{sample.sample_name}', alpha=0.7)

        # Add the line of best fit (dashed), but make it cover the entire global x-range (0 to global_max_x)
        slope, intercept = sample.line_fit_coefficients
        best_fit_x = np.linspace(0, global_max_x, 100)  # Generate 100 points from 0 to global_max_x
        best_fit_y = slope * best_fit_x + intercept
        # Change "best fit" to "calc" in the label
        plt.plot(best_fit_x, best_fit_y, linestyle='--', color=colors[idx % len(colors)], alpha=0.7, label=f'{sample.sample_name} calc')

        # Now label the hkl indices, ensuring no duplicates and adjusting label position
        for i in range(len(k_rootC)):
            k_val = k_rootC[i]
            delta_k_val = delta_k[i]
            hkl_val = hkl[i]

            # Check if a point within a small error range is already labeled
            already_labeled = any(np.sqrt((k_val - x) ** 2 + (delta_k_val - y) ** 2) < error_threshold for x, y in labeled_points)

            if not already_labeled:
                # Identify clusters of points (i.e., points that are close to each other in the x-direction)
                cluster_indices = [j for j in range(len(k_rootC)) if abs(k_rootC[j] - k_val) < error_threshold]
                
                if len(cluster_indices) > 1:
                    # This is a cluster of points
                    cluster_delta_k_vals = [delta_k[j] for j in cluster_indices]
                    topmost_index = cluster_indices[np.argmax(cluster_delta_k_vals)]
                    bottommost_index = cluster_indices[np.argmin(cluster_delta_k_vals)]
                    k_top = k_rootC[topmost_index]
                    k_bottom = k_rootC[bottommost_index]
                    delta_k_top = delta_k[topmost_index]
                    delta_k_bottom = delta_k[bottommost_index]

                    # Calculate the trendline values at the top and bottom of the cluster
                    y_trendline_top = slope * k_top + intercept
                    y_trendline_bottom = slope * k_bottom + intercept

                    if y_trendline_top < delta_k_bottom:
                        # If the trendline is below the cluster, place the label above the topmost point
                        y_pos = delta_k_top + 0.05 * global_max_y
                    else:
                        # If the trendline is above the cluster, place the label below the bottommost point
                        y_pos = delta_k_bottom - 0.05 * global_max_y

                else:
                    # Single point, avoid trendline and label in available space
                    y_trendline = slope * k_val + intercept
                    if abs(delta_k_val - y_trendline) < 0.02 * global_max_y or delta_k_val > global_max_y / 2:
                        # If it's near the trendline or high in the plot, place the label below
                        y_pos = delta_k_val - 0.05 * global_max_y
                    else:
                        # Otherwise, place the label above the point
                        y_pos = delta_k_val + 0.05 * global_max_y

                # Add the annotation with black text
                plt.text(k_val, y_pos, hkl_val, fontsize=10, ha='center', color='black')

                # Add this point to the list of labeled points
                labeled_points.append((k_val, delta_k_val))

    # Set axis labels and titles
    plt.xlabel(r'$k \cdot \sqrt{C_{hkl}}$', fontsize=14, fontname='Arial')
    plt.ylabel(r'$\Delta k$', fontsize=14, fontname='Arial')

    # Set axis limits based on the global max values
    plt.xlim(0, global_max_x)
    plt.ylim(0, global_max_y)

    # Set x-axis major tick intervals to 0.5
    plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))

    # Remove gridlines
    plt.grid(False)

    # Set tick label font size
    plt.xticks(fontsize=12, fontname='Arial')
    plt.yticks(fontsize=12, fontname='Arial')

    # Add a legend outside of the plot to the right
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

    # Set plot title using the provided input
    plt.title(plot_title, fontsize=16, fontname='Arial')

    # Adjust layout to fit the legend
    plt.tight_layout()

    # Save the plot as a PNG file with the concatenated sample names in the filename
    plt.savefig(plot_filename, format='png')

    # Close the plot to free memory
    plt.close()

    print(f"Plot saved as {plot_filename}")


# In[95]:


mWH_plots_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\MASc Thesis\\Results\\Results for paper\\mWH plots for paper\\New 2'


# In[ ]:


index_title_combos = [([3,4,5],'mWH Plot 617_32000_HPT Various Radii'),
                      ([9,10,11],'mWH Plot 617_SA_HPT Various Radii'),
                      ([13,14,15], 'mWH Plot MA754_HPT Various Radii'),
                      ([19,20,21], 'mWH Plot MA957_HPT Various Radii'),
                      ([1,2,5], 'mWH Plot 617_32000_HPT Various Doses'),
                      ([7,8,10], 'mWH Plot 617_SA_HPT Various Doses'),
                      ([17,18,19], 'mWH Plot MA957_HPT Various Doses')
                     ]


# In[ ]:


index_title_combos = [([16, 17, 18], 'mWH Plot MA957 Various Doses')]


# In[ ]:


for i in index_title_combos:
    plot_mWH_samples([sample_objects[x] for x in i[0]], mWH_plots_dir, i[1])


# In[102]:


def plot_WH_samples(samples, output_directory, plot_title):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Colors for different data sets
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    
    # Collect all sample names to use in the filename
    sample_names = "_".join([sample.sample_name for sample in samples])

    # Create a filename based on concatenated sample names
    plot_filename = os.path.join(output_directory, f'{sample_names}_WH.png')

    # Create the figure for the plot
    plt.figure(figsize=(16, 9))  # Set figure size to 16:9

    # Loop through each sample object to plot
    for idx, sample in enumerate(samples):
        # Prepare the data for this sample
        k = sample.k_delta_k_hkl_C_isub[0]
        delta_k = sample.k_delta_k_hkl_C_isub[1]
        hkl = sample.k_delta_k_hkl_C_isub[2]  # hkl indices for labeling

        # Sort the data by k values to connect points in the correct order
        sorted_indices = np.argsort(k)
        k_sorted = np.array(k)[sorted_indices]
        delta_k_sorted = np.array(delta_k)[sorted_indices]
        hkl_sorted = np.array(hkl)[sorted_indices]

        # Plot the data points and connect them with a solid line (sorted by k)
        plt.plot(k_sorted, delta_k_sorted, color=colors[idx % len(colors)], marker='o', linestyle='-', label=f'{sample.sample_name}', alpha=0.7)

        # Annotate each point with its hkl value
        for i in range(len(k_sorted)):
            plt.text(k_sorted[i], delta_k_sorted[i], hkl_sorted[i], fontsize=10, ha='center', color='black')

    # Set x and y axis limits based on data with a buffer
    all_k_values = np.concatenate([sample.k_delta_k_hkl_C_isub[0] for sample in samples])
    all_delta_k_values = np.concatenate([sample.k_delta_k_hkl_C_isub[1] for sample in samples])

    # Calculate buffer size (10% buffer)
    x_buffer = 0.1 * (max(all_k_values) - min(all_k_values))
    y_buffer = 0.1 * (max(all_delta_k_values) - min(all_delta_k_values))

    # Set axis limits with buffer, and ensure y-axis has extra space below the x-axis
    plt.xlim(min(all_k_values) - x_buffer, max(all_k_values) + x_buffer)
    plt.ylim(min(all_delta_k_values) - y_buffer - 0.05 * abs(min(all_delta_k_values)), max(all_delta_k_values) + y_buffer)

    # Set x-axis major tick intervals to 0.5
    plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))

    # Remove gridlines
    plt.grid(False)

    # Set tick label font size
    plt.xticks(fontsize=12, fontname='Arial')
    plt.yticks(fontsize=12, fontname='Arial')

    # Set axis labels with the desired units
    plt.xlabel(r'$k \, [\mathrm{nm}^{-1}]$', fontsize=14, fontname='Arial')
    plt.ylabel(r'$\Delta k \, [\mathrm{nm}^{-1}]$', fontsize=14, fontname='Arial')

    # Add a legend outside of the plot to the right
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

    # Set plot title using the provided input
    plt.title(plot_title, fontsize=16, fontname='Arial')

    # Adjust layout to fit the legend
    plt.tight_layout()

    # Save the plot as a PNG file with the concatenated sample names in the filename
    plt.savefig(plot_filename, format='png')

    # Close the plot to free memory
    plt.close()

    print(f"Plot saved as {plot_filename}")


# In[103]:


def plot_WH_samples(samples, output_directory, plot_title):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # Colors for different data sets
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    
    # Collect all sample names to use in the filename
    sample_names = "_".join([sample.sample_name for sample in samples])

    # Create a filename based on concatenated sample names
    plot_filename = os.path.join(output_directory, f'{sample_names}_WH.png')

    # Create the figure for the plot
    plt.figure(figsize=(16, 9))  # Set figure size to 16:9

    # Loop through each sample object to plot
    for idx, sample in enumerate(samples):
        # Prepare the data for this sample
        k = sample.k_delta_k_hkl_C_isub[0]
        delta_k = sample.k_delta_k_hkl_C_isub[1]
        hkl = sample.k_delta_k_hkl_C_isub[2]  # hkl indices for labeling

        # Sort the data by k values to connect points in the correct order
        sorted_indices = np.argsort(k)
        k_sorted = np.array(k)[sorted_indices]
        delta_k_sorted = np.array(delta_k)[sorted_indices]
        hkl_sorted = np.array(hkl)[sorted_indices]

        # Plot the data points and connect them with a solid line (sorted by k)
        plt.plot(k_sorted, delta_k_sorted, color=colors[idx % len(colors)], marker='o', linestyle='-', label=f'{sample.sample_name}', alpha=0.7)

        # Annotate only one point per cluster
        cluster_distance_threshold = 0.2  # Define a threshold for clustering points
        annotated_points = []

        for i in range(len(k_sorted)):
            current_point = (k_sorted[i], delta_k_sorted[i])

            # Calculate the distance from current point to already annotated points
            distances = [distance.euclidean(current_point, ann) for ann in annotated_points]

            # Annotate only if the point is far enough from already annotated points
            if not distances or min(distances) > cluster_distance_threshold:
                plt.text(k_sorted[i], delta_k_sorted[i], hkl_sorted[i], fontsize=12, ha='center', color='black')
                annotated_points.append(current_point)

    # Set x and y axis limits based on data with a buffer
    all_k_values = np.concatenate([sample.k_delta_k_hkl_C_isub[0] for sample in samples])
    all_delta_k_values = np.concatenate([sample.k_delta_k_hkl_C_isub[1] for sample in samples])

    # Calculate buffer size (10% buffer)
    x_buffer = 0.1 * (max(all_k_values) - min(all_k_values))
    y_buffer = 0.1 * (max(all_delta_k_values) - min(all_delta_k_values))

    # Set axis limits with buffer, and ensure y-axis has extra space below the x-axis
    plt.xlim(min(all_k_values) - x_buffer, max(all_k_values) + x_buffer)
    plt.ylim(min(all_delta_k_values) - y_buffer - 0.05 * abs(min(all_delta_k_values)), max(all_delta_k_values) + y_buffer)

    # Set x-axis major tick intervals to 1
    plt.gca().xaxis.set_major_locator(MultipleLocator(1))  # Changed from 0.5 to 1

    # Remove gridlines
    plt.grid(False)

    # Set tick label font size
    plt.xticks(fontsize=18, fontname='Arial')  # Increased font size
    plt.yticks(fontsize=18, fontname='Arial')  # Increased font size

    # Set axis labels with the desired units and increased font size
    plt.xlabel(r'$k \, [\mathrm{nm}^{-1}]$', fontsize=22, fontname='Arial')  # Larger axis label
    plt.ylabel(r'$\Delta k \, [\mathrm{nm}^{-1}]$', fontsize=22, fontname='Arial')  # Larger axis label

    # Add a legend at the bottom right corner
    plt.legend(loc='lower right', fontsize=16)  # Moved legend to bottom right

    # Set plot title using the provided input and increased font size
    plt.title(plot_title, fontsize=24, fontname='Arial')  # Larger plot title

    # Adjust layout to fit the plot and legend
    plt.tight_layout()

    # Save the plot as a PNG file with the concatenated sample names in the filename
    plt.savefig(plot_filename, format='png')

    # Close the plot to free memory
    plt.close()

    print(f"Plot saved as {plot_filename}")


# In[104]:


WH_plots_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Conferences\\MS&T 2024'


# In[105]:


index_title_combos_WH = [([3,4,5],' Plot 617_32000_HPT Various Radii'),
                      ([9,10,11],'WH Plot 617_SA_HPT Various Radii'),
                      ([13,14,15], 'WH Plot MA754_HPT Various Radii'),
                      ([19,20], 'WH Plot MA957_HPT Various Radii'),
                      ([1,2,5], 'WH Plot 617_32000_HPT Various Doses'),
                      ([7,8,10], 'WH Plot 617_SA_HPT Various Doses'),
                      ([17,18,19], 'WH Plot MA957_HPT Various Doses'),
                      ([0, 6, 12, 16], 'WH Plots Undeformed Samples')
                     ]


# In[106]:


index_title_combos_WH = []


# In[111]:


index_title_combos_WH = [([16,18,17], 'WH Plot MA957 Various Doses')]


# In[ ]:





# In[112]:


for i in index_title_combos_WH:
    plot_WH_samples([sample_objects[x] for x in i[0]], WH_plots_dir, i[1])


# In[ ]:





# In[ ]:




