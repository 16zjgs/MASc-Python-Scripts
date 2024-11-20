#!/usr/bin/env python
# coding: utf-8

# In[8]:


import os
import numpy as np
from scipy.interpolate import CubicSpline


# In[9]:


root_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\MASc Thesis\\Results\\CMWP_BG_plots_Sept_2024\\CMWP Outputs'


# In[10]:


def find_spline_dat_files(directory):
    spline_files = {}
    
    # Loop through all subdirectories and files
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Check if the file is a .dat file and contains 'spline' in the name
            if file.endswith('.dat') and 'spline' in file.lower():
                file_path = os.path.join(root, file)
                
                # Read the file and store its content as tuples (x, y)
                with open(file_path, 'r') as f:
                    file_values = []
                    for line in f:
                        # Split each line into x and y values and convert to float
                        if line.strip():  # Ensure line is not empty
                            x, y = map(float, line.strip().split())
                            file_values.append((x, y))
                
                # Store the file values in the dictionary with the filename as the key
                spline_files[file] = file_values
    
    return spline_files


# In[11]:


spline_points = find_spline_dat_files(root_dir)


# In[12]:


def compute_cubic_coefficients(directory):
    spline_files = find_spline_dat_files(directory)  # Use the previous function to get spline points
    cubic_coefficients_with_domain = {}
    
    for file, points in spline_files.items():
        # Unpack points into separate x and y lists
        x_values, y_values = zip(*points)
        
        # Fit a cubic spline to the x and y data
        cs = CubicSpline(x_values, y_values)
        
        # Extract the cubic spline coefficients for each segment and include the domain
        file_coefficients = []
        for i in range(len(cs.x) - 1):
            # Coefficients for segment i in the form of (a, b, c, d)
            a, b, c, d = cs.c[:, i]
            # Domain for this segment (x_start, x_end)
            x_start, x_end = cs.x[i], cs.x[i+1]
            # Append tuple containing the coefficients and the x-range for the segment
            file_coefficients.append(((a, b, c, d), (x_start, x_end)))
        
        cubic_coefficients_with_domain[file] = file_coefficients
    
    return cubic_coefficients_with_domain


# In[13]:


spline_coeffs = compute_cubic_coefficients(root_dir)


# In[14]:


spline_coeffs.keys()


# In[17]:


def rename_keys(original_dict, new_key_names):
    """
    Rename the keys of the original dictionary based on a provided list of new key names.
    
    Parameters:
    original_dict (dict): The original dictionary.
    new_key_names (list): A list of new key names, in the same order as the original dictionary.
    
    Returns:
    dict: A new dictionary with renamed keys.
    """
    if len(new_key_names) != len(original_dict):
        raise ValueError("The length of the new key names list must match the number of keys in the original dictionary.")
    
    renamed_dict = {}
    
    # Loop through original dict and the new key names simultaneously
    for (old_key, new_key) in zip(original_dict.keys(), new_key_names):
        renamed_dict[new_key] = original_dict[old_key]
    
    return renamed_dict


# In[18]:


new_keys = ['617_32000_HPT',
            '617_SA_HPT',
            'MA754_HPT',
            'MA957_HPT'
           ]


# In[19]:


spline_coeffs_named = rename_keys(spline_coeffs, new_keys)


# In[35]:


def apply_spline_subtraction(spline_coeffs, directory):
    """
    Subtract the spline y-values from the corresponding file y-values in the directory.
    
    Parameters:
    spline_coeffs (dict): A dictionary where the keys are folder names and the values are spline coefficients with x-ranges.
    directory (str): The directory containing the folders and files to process.
    """
    
    for key, segments in spline_coeffs.items():
        # Find the folder in the directory that contains the key
        folder_path = None
        for root, dirs, files in os.walk(directory):
            for dir_name in dirs:
                if key in dir_name:
                    folder_path = os.path.join(root, dir_name)
                    break
            if folder_path:
                break
        
        if not folder_path:
            print(f"Folder containing '{key}' not found.")
            continue
        
        # Find the file ending with "_k_converted"
        target_file = None
        for file in os.listdir(folder_path):
            if file.endswith('_k_converted.txt') or file.endswith('_k_converted.xye'):
                target_file = file
                break
        
        if not target_file:
            print(f"No '_k_converted' file found in folder '{folder_path}'.")
            continue
        
        target_file_path = os.path.join(folder_path, target_file)
        
        # Read the x and y values from the target file
        data = np.loadtxt(target_file_path)
        x_values, y_values = data[:, 0], data[:, 1]
        
        # Prepare to store the new y-values (with background subtraction)
        new_y_values = np.copy(y_values)
        
        # For each segment, calculate the spline y-values and subtract them from the actual y-values
        for (coefficients, (x_start, x_end)) in segments:
            # Find indices within the x-range
            indices_in_range = np.where((x_values >= x_start) & (x_values <= x_end))[0]
            
            if len(indices_in_range) == 0:
                continue  # Skip if no x-values are within this range
            
            a, b, c, d = coefficients  # Unpack coefficients
            
            for idx in indices_in_range:
                x = x_values[idx]  # Actual x-value
                # Calculate the spline y-value using the cubic polynomial
                spline_y = a + b * (x - x_start) + c * (x - x_start)**2 + d * (x - x_start)**3
                # Subtract the spline y-value from the actual y-value
                new_y_values[idx] -= spline_y
        
        # Write the new x-y data to a new file with '_bgsub' added to the filename
        new_file_name = target_file.replace('_k_converted', '_bgsub')
        new_file_path = os.path.join(folder_path, new_file_name)
        
        # Save the new x-y data
        np.savetxt(new_file_path, np.column_stack((x_values, new_y_values)), fmt='%.6f', delimiter=' ')
        
        print(f"Background-subtracted file saved as: {new_file_path}")


# In[30]:


def apply_spline_subtraction(spline_coeffs, directory):
    """
    Subtract the spline y-values from the corresponding file y-values in all valid folders in the directory.
    
    Parameters:
    spline_coeffs (dict): A dictionary where the keys are folder names and the values are spline coefficients with x-ranges.
    directory (str): The directory containing the folders and files to process.
    """
    
    for key, segments in spline_coeffs.items():
        # Search all folders and files in the directory for matching criteria
        for root, dirs, files in os.walk(directory):
            for dir_name in dirs:
                # Check if the directory contains the key and does not contain "1dpa" or "2dpa"
                if key in dir_name and "1dpa" not in dir_name and "2dpa" not in dir_name:
                    folder_path = os.path.join(root, dir_name)
                    
                    # Iterate through files in the valid folder
                    for file in os.listdir(folder_path):
                        # Check if the file does NOT end with "_k_converted", "_bgsub", and does not contain "peaks"
                        if (
                            (file.endswith('.txt') or file.endswith('.xye')) 
                            and not (file.endswith('_k_converted.txt') or file.endswith('_bgsub.txt') 
                                     or file.endswith('_k_converted.xye') or file.endswith('_bgsub.xye')) 
                            and 'peaks' not in file
                        ):
                            target_file_path = os.path.join(folder_path, file)
                            
                            # Read the x and y values from the target file
                            data = np.loadtxt(target_file_path)
                            x_values, y_values = data[:, 0], data[:, 1]
                            
                            # Prepare to store the new y-values (with background subtraction)
                            new_y_values = np.copy(y_values)
                            
                            # For each segment, calculate the spline y-values and subtract them from the actual y-values
                            for (coefficients, (x_start, x_end)) in segments:
                                # Find indices within the x-range
                                indices_in_range = np.where((x_values >= x_start) & (x_values <= x_end))[0]
                                
                                if len(indices_in_range) == 0:
                                    continue  # Skip if no x-values are within this range
                                
                                a, b, c, d = coefficients  # Unpack coefficients
                                
                                for idx in indices_in_range:
                                    x = x_values[idx]  # Actual x-value
                                    # Calculate the spline y-value using the cubic polynomial
                                    spline_y = a + b * (x - x_start) + c * (x - x_start)**2 + d * (x - x_start)**3
                                    # Subtract the spline y-value from the actual y-value
                                    new_y_values[idx] -= spline_y
                            
                            # Write the new x-y data to a new file with '_bgsub' added to the filename
                            new_file_name = file.replace('.txt', '_bgsub.txt').replace('.xye', '_bgsub.xye')
                            new_file_path = os.path.join(folder_path, new_file_name)
                            
                            # Save the new x-y data
                            np.savetxt(new_file_path, np.column_stack((x_values, new_y_values)), fmt='%.6f', delimiter=' ')
                            
                            print(f"Background-subtracted file saved as: {new_file_path}")


# In[33]:


process_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\MASc Thesis\\Results\\Results for paper\\Peak Narrowing\\MA957_HPT\\2 dpa.xye'


# In[32]:


apply_spline_subtraction(spline_coeffs_named, process_dir)


# In[ ]:




