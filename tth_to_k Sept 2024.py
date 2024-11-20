#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import math
import pandas as pd
import numpy as np


# In[2]:


def tth_to_k(file_path, wavelength_nm=0.154056): 
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"The file {file_path} does not exist.")
        return

    # Prepare the new file name
    directory, filename = os.path.split(file_path)
    new_filename = filename.rsplit('.', 1)[0] + '_k_converted.' + filename.rsplit('.', 1)[1]
    new_file_path = os.path.join(directory, new_filename)

    try:
        with open(file_path, 'r') as file, open(new_file_path, 'w') as new_file:
            for line in file:
                # Split the line into columns based on spaces
                columns = line.strip().split()
                
                # Try to process only lines that contain at least two numeric columns
                try:
                    if len(columns) >= 2:
                        # Check if the first two columns are numeric
                        tth_degrees = float(columns[0])
                        intensity = float(columns[1])  # Second column (intensity or similar)

                        # Proceed with conversion for 2theta to k
                        th_degrees = tth_degrees / 2
                        theta_radians = math.radians(th_degrees)
                        k = (2 * math.sin(theta_radians)) / wavelength_nm

                        # Write only the first two columns (k and intensity) to the new file
                        new_file.write(f"{k} {intensity}\n")
                    else:
                        print(f"Line skipped due to insufficient columns: {line.strip()}")
                except ValueError:
                    # If conversion to float fails, skip the line (likely non-numeric content)
                    continue

        # Remove the original file after successful conversion
        os.remove(file_path)
        print(f"Original file {file_path} deleted.")
        print(f"File saved successfully as {new_file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")
def process_directory(directory):
    for root, dirs, files in os.walk(directory):
        # Only process if there are no more subdirectories (bottom-level directory)
        if not dirs:
            for file in files:
                # Check if the file has the .txt or .xye extension and does not contain '_k_converted'
                if (file.endswith('.txt') or file.endswith('.xye')) and '_k_converted' not in file:
                    file_path = os.path.join(root, file)
                    # Apply the tth_to_k function to the file
                    tth_to_k(file_path)


# In[3]:


def tth_to_k(file_path, wavelength_nm=0.067): 
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"The file {file_path} does not exist.")
        return

    # Prepare the new file name
    directory, filename = os.path.split(file_path)
    new_filename = filename.rsplit('.', 1)[0] + '_k_converted.' + filename.rsplit('.', 1)[1]
    new_file_path = os.path.join(directory, new_filename)

    try:
        with open(file_path, 'r') as file, open(new_file_path, 'w') as new_file:
            for line in file:
                # Split the line into columns based on spaces
                columns = line.strip().split()
                
                # Try to process only lines that contain at least two numeric columns
                try:
                    if len(columns) >= 2:
                        # Check if the first two columns are numeric
                        tth_degrees = float(columns[0])
                        intensity = float(columns[1])  # Second column (intensity or similar)

                        # Proceed with conversion for 2theta to k
                        th_degrees = tth_degrees / 2
                        theta_radians = math.radians(th_degrees)
                        k = (2 * math.sin(theta_radians)) / wavelength_nm

                        # Write only the first two columns (k and intensity) to the new file
                        new_file.write(f"{k} {intensity}\n")
                    else:
                        print(f"Line skipped due to insufficient columns: {line.strip()}")
                except ValueError:
                    # If conversion to float fails, skip the line (likely non-numeric content)
                    continue

        # Remove the original file after successful conversion
        os.remove(file_path)
        print(f"Original file {file_path} deleted.")
        print(f"File saved successfully as {new_file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

def process_directory(directory, exclude_folders=None):
    """
    Process the files in a directory and exclude certain folders from being processed.

    Parameters:
    directory (str): The root directory to search for files.
    exclude_folders (list): A list of folder names to exclude from processing.
    """
    if exclude_folders is None:
        exclude_folders = []

    for root, dirs, files in os.walk(directory):
        # Skip directories that are in the exclude_folders list
        if any(excluded in root for excluded in exclude_folders):
            print(f"Skipping folder: {root}")
            continue

        # Only process if there are no more subdirectories (bottom-level directory)
        if not dirs:
            for file in files:
                # Check if the file ends with '_bgsub.txt' or '_bgsub.xye'
                if (file.endswith('_bgsub.txt') or file.endswith('_bgsub.xye')):
                    file_path = os.path.join(root, file)
                    # Apply the tth_to_k function to the file
                    tth_to_k(file_path)


# In[4]:


def tth_to_k(file_path, wavelength_nm=0.154056): 
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"The file {file_path} does not exist.")
        return

    # Prepare the new file name
    directory, filename = os.path.split(file_path)
    new_filename = filename.rsplit('.', 1)[0] + '_k_converted.' + filename.rsplit('.', 1)[1]
    new_file_path = os.path.join(directory, new_filename)

    try:
        with open(file_path, 'r') as file, open(new_file_path, 'w') as new_file:
            for line in file:
                # Split the line into columns based on spaces
                columns = line.strip().split()
                
                # Try to process only lines that contain at least two numeric columns
                try:
                    if len(columns) >= 2:
                        # Check if the first two columns are numeric
                        tth_degrees = float(columns[0])
                        intensity = float(columns[1])  # Second column (intensity or similar)

                        # Proceed with conversion for 2theta to k
                        th_degrees = tth_degrees / 2
                        theta_radians = math.radians(th_degrees)
                        k = (2 * math.sin(theta_radians)) / wavelength_nm

                        # Write only the first two columns (k and intensity) to the new file
                        new_file.write(f"{k} {intensity}\n")
                    else:
                        print(f"Line skipped due to insufficient columns: {line.strip()}")
                except ValueError:
                    # If conversion to float fails, skip the line (likely non-numeric content)
                    continue

        # Remove the original file after successful conversion
        os.remove(file_path)
        print(f"Original file {file_path} deleted.")
        print(f"File saved successfully as {new_file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


def process_directory(directory, exclude_folders=None):
    """
    Process the files in a directory and exclude certain folders from being processed.

    Parameters:
    directory (str): The root directory to search for files.
    exclude_folders (list): A list of folder names to exclude from processing.
    """
    if exclude_folders is None:
        exclude_folders = []

    for root, dirs, files in os.walk(directory):
        # Skip directories that are in the exclude_folders list
        if any(excluded in root for excluded in exclude_folders):
            print(f"Skipping folder: {root}")
            continue

        # Only process if there are no more subdirectories (bottom-level directory)
        if not dirs:
            for file in files:
                # Check if the file ends with '_bgsub.txt' or '_bgsub.xye'
                if file:
                    file_path = os.path.join(root, file)
                    # Apply the tth_to_k function to the file
                    tth_to_k(file_path)


# In[5]:


def tth_to_k(file_path, wavelength_nm=0.154056):
    """
    Convert 2theta (tth) to k and write the result to a new file with _k_converted appended to the filename.
    
    Args:
    file_path (str): The path to the file (.txt or .xye) that contains 2theta and intensity data.
    wavelength_nm (float): The wavelength in nm (default is 0.067 nm).
    
    Returns:
    None
    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"The file {file_path} does not exist.")
        return

    # Prepare the new file name
    directory, filename = os.path.split(file_path)
    new_filename = filename.rsplit('.', 1)[0] + '_k_converted.' + filename.rsplit('.', 1)[1]
    new_file_path = os.path.join(directory, new_filename)

    try:
        # Load the data from the file (two or three columns)
        data = pd.read_csv(file_path, delim_whitespace=True, header=None)

        # Handle case where the file has 3 columns, discard the third column
        if data.shape[1] >= 3:
            data = data[[0, 1]]  # Keep only the first two columns (2theta and intensity)

        # Extract 2theta and intensity
        tth_degrees = data[0].values  # First column (2theta)
        intensity = data[1].values    # Second column (intensity)

        # Convert 2theta to k
        th_degrees = tth_degrees / 2
        theta_radians = np.radians(th_degrees)
        k_values = (2 * np.sin(theta_radians)) / wavelength_nm

        # Write the new data (k and intensity) to the new file
        with open(new_file_path, 'w') as new_file:
            for k, inten in zip(k_values, intensity):
                new_file.write(f"{k} {inten}\n")

        print(f"File saved successfully as {new_file_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


# In[6]:


file_path = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\Winter 2024\\New LXRD Data\\Nov 1 Data MA754'


# In[7]:


process_directory(file_path)


# In[10]:


tth_values = [44, 51.24, 75.527, 91.812, 44.146, 51.411, 75.7504, 92.086]
tth_values = np.array(tth_values)


# In[11]:


def tth_to_k(tth_values, wavelength_nm=0.154056):
    """
    Convert a list of 2theta (tth) values to k values.
    
    Args:
    tth_values (list or array-like): A list or array of 2theta values (in degrees).
    wavelength_nm (float): The wavelength in nm (default is 0.154056 nm).
    
    Returns:
    k_values (numpy array): An array of converted k values.
    """
    try:
        # Convert 2theta to theta
        th_degrees = np.array(tth_values) / 2
        theta_radians = np.radians(th_degrees)

        # Convert theta to k
        k_values = (2 * np.sin(theta_radians)) / wavelength_nm

        return k_values

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# In[12]:


k_values = tth_to_k(tth_values)


# In[13]:


k_values


# In[ ]:




