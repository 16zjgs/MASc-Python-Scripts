#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import os
import re
from difflib import SequenceMatcher
import itertools


# In[1]:


data_dict = {
    '617_SA': [None, 0.0, 0.0, 37.77, 999.0, 5064.0, 7.73, 1.44, 0.4, None, 11.85, 2.424],
    '617_SA_2dpa': [None, 2.0, 96.23, 666.2, 999.98, 1356.0, 0.01, 0.1669, 0.05, 0.6054, 3590.0, 18350.0],
    '617_SA_HPT_c': [0.0, 0.0, 83.36, 0.5017, 12.02, 0.06056, 221.2, 1.801, 0.7463, 0.001635, 12.25, 0.1073],
    '617_SA_HPT_m': [0.5, 0.0, 89.42, 0.6798, 10.58, 0.1173, 254.8, 5.268, 0.3476, 0.007042, 18.31, 0.6193],
    '617_SA_HPT_e': [1, 0.0, 93.23, 0.4504, 10.6, 0.116, 301.5, 4.372, 0.5051, 0.008573, 17.35, 0.5606],
    '617_SA_HPT_1dpa': [None, 1.0, 0.0, 2.771, 16.93, 0.3424, 62.18, 1.489, 0.5081, 0.01558, 15.27, 0.4489],
    '617_SA_HPT_2dpa': [None, 2.0, 16.75, 2.256, 16.16, 0.2957, 94.87, 2.08, 0.3365, 0.007381, 30.47, 1.119],
    '617_32000': [None, 0.0, 100.0, 1.534, 999.0, 2830.0, 0.01, 0.045, 0.4, None, 87.56, 17.09],
    '617_32000_2dpa': [None, 2.0, 100.0, 8.909, 330.76, 35.81, 38.78, 6.92, 0.063, 0.0064, 121.8, 19.38],
    '617_32000_HPT_c': [0, 0.0, 94.33, 0.6565, 11.01, 0.1398, 276.6, 6.444, 0.4033, 0.01018, 9.367, 0.1735],
    '617_32000_HPT_m': [0.5, 0.0, 77.38, 0.2154, 12.3, 0.09735, 298.5, 5.262, 0.3146, 0.00538, 11.97, 0.2018],
    '617_32000_HPT_e': [1, 0.0, 87.96, 0.2046, 13.65, 0.1497, 200.0, 3.161, 0.671, 0.01647, 11.58, 0.222],
    '617_32000_HPT_1dpa': [None, 1.0, 34.0, 2.17, 15.1, 0.2402, 74.27, 0.3847, 0.3646, 0.005985, 40.61, 2.637],
    '617_32000_HPT_2dpa': [None, 2.0, 78.9, 1.247, 14.17, 0.129, 94.24, 1.187, 0.2216, 0.001103, 25.53, 0.6861],
    'MA957': [None, 0.0, 100.0, 4.927, 59.783, 1.627, 0.103, 0.0324, 0.4, None, None, None],
    'MA957_1dpa': [None, 1.0, 100.0, 3.466, 25.616, 0.7872, 30.88, 0.5084, 0.143, 0.006807, None, None],
    'MA957_2dpa': [None, 2.0, 50.94, 4.641, 24.63, 0.3992, 87.66, 4.94, 0.106, 0.005, None, None],
    'MA957_HPT_c': [0.0, 0.0, 51.48, 1.224, 11.05, 0.04918, 163.2, 2.261, 0.4356, 0.00679, None, None],
    'MA957_HPT_m': [0.5, 0.0, 49.5, 1.037, 12.11, 0.05383, 184.0, 2.006, 0.4925, 0.005836, None, None],
    'MA957_HPT_e': [1, 0.0, 51.97, 0.4864, 12.42, 0.06655, 174.9, 1.749, 0.5697, 0.007, None, None],
    'MA957_HPT_1dpa': [None, 1.0, 60.58, 1.102, 17.93, 0.04907, 75.96, 0.5478, 0.1198, 0.001154, None, None],
    'MA957_HPT_2dpa': [None, 2.0, 71.47, 0.7286, 22.46, 0.05442, 17.39, 0.2582, 0.4571, 0.01039, None, None],
    'MA754': [None, 0.0, 100.0, 9.623, 97.4, 4.784, 2.79, 0.1243, 0.4, None, None, None],
    'MA754_1dpa': [None, 1.0, 12.09, 2.083, 999.7, 1005, 52.42, 1.826, 0.1016, 0.002, None, None],
    'MA754_2dpa': [None, 2.0, 100.0, 2.304, 199.7, 57.22, 93.9, 16.03, 0.06042, 0.005, None, None],
    'MA754_HPT_c': [0, 0.0, 100.0, 0.5959, 14.61, 0.1936, 137.5, 1.659, 0.6456, 0.01049, 30.82, 1.542],
    'MA754_HPT_m': [0.5, 0.0, 100.0, 0.06544, 15.25, 0.1627, 121.0, 1.411, 0.9408, 0.02045, 16.88, 0.5453],
    'MA754_HPT_e': [1, 0.0, 100.0, 0.5856, 18.64, 0.3745, 105.7, 1.408, 1.536, 0.04944, 12.47, 0.2976],
    'MA754_HPT_2dpa': [None, 2.0, 0.0, 17.14, 16.32, 1.218, 37.32, 4.05, 0.7934, 0.14, 21.67, 3.308],
}


# In[2]:


data_cols = [
    "radii",
    "dose [dpa]",
    "% screw",
    "% screw_err",
    "CSDS [nm]",
    "CSDS_err [nm]",
    "rho (1e-14) [m^-2]",
    "rho_err (1e-14) [m^-2]",
    "M",
    "M_err",
    "twin spacing [nm]",
    "twin spacing_err [nm]"
]


# In[3]:


save_dir = 'C:\\Users\\Zachary\\Desktop\\MASc\\Thesis\\MASc Thesis\\Results\\Results for paper\\CMWP Plots Nov\\Nov 11'


# In[4]:


def plot_selected_samples(data_dict, column_names, xy_columns, selected_sample_groups, error_columns, title, save_dir=None, filename=None):
    """
    Plots specified samples from the data_dict using the column indices for x and y axes,
    with optional error bars, and connects grouped samples with lines along the x-axis.
    Saves the plot as a PNG file with a dynamic name based on the first sample, x-axis, and y-axis if filename is not provided.

    Parameters:
    - data_dict: Dictionary containing sample data.
    - column_names: List of column names corresponding to data_dict columns.
    - xy_columns: Tuple (x_col_index, y_col_index) indicating the indices of the columns for the x and y axes.
    - selected_sample_groups: List of lists, where each sublist contains specific keys (samples) to be connected by lines.
    - error_columns: Tuple (x_err_col_index, y_err_col_index) indicating the columns for x and y error bars, or None if no error bars.
    - title: Title for the plot.
    - save_dir: Directory where the plot will be saved as a PNG file (optional).
    - filename: Custom filename for saving the plot. If not provided, a dynamic name is generated.
    """
    
    # Set the figure size to 16:9 and the font to Arial
    plt.figure(figsize=(16, 9))
    plt.rc('font', family='Arial')

    # Get the column indices for the x and y axis and error bars from the tuples
    x_col_index, y_col_index = xy_columns
    x_err_col_index, y_err_col_index = error_columns

    # Get the column names for labeling
    x_col_name = column_names[x_col_index]
    y_col_name = column_names[y_col_index]

    all_x_values = []
    all_y_values = []

    # Extract the name of the first sample for use in the filename if filename is not provided
    first_sample = selected_sample_groups[0][0] if selected_sample_groups and selected_sample_groups[0] else "sample"

    # Function to find the longest common substring in a list of strings
    def common_substring(strings):
        common = strings[0]
        for string in strings[1:]:
            match = SequenceMatcher(None, common, string).find_longest_match(0, len(common), 0, len(string))
            common = common[match.a: match.a + match.size]
        return common.strip("_")  # Remove trailing/leading underscores

    # Color cycle for each dataset
    color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    # Loop through the selected sample groups
    for group in selected_sample_groups:
        x_values = []
        y_values = []
        x_errors = []
        y_errors = []
        
        # Get the common substring for the legend from group labels
        group_labels = [key for key in group if key in data_dict]
        common_label = common_substring(group_labels)
        
        # Assign a consistent color for this group
        color = next(color_cycle)
        
        for key in group:
            # Check if the selected key exists in the dictionary and has valid data points for both x and y axes
            if key in data_dict:
                values = data_dict[key]
                # Ensure there are valid data points for x and y
                x_value = values[x_col_index] if len(values) > x_col_index else None
                y_value = values[y_col_index] if len(values) > y_col_index else None

                if x_value is not None and y_value is not None:
                    # Error bars for x and y (if specified)
                    x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else None
                    y_err = values[y_err_col_index] if y_err_col_index is not None and len(values) > y_err_col_index else None
                    
                    # Collect values for plotting
                    x_values.append(x_value)
                    y_values.append(y_value)
                    x_errors.append(x_err)
                    y_errors.append(y_err)
        
        # Collect all values to adjust axis limits later
        all_x_values.extend(x_values)
        all_y_values.extend(y_values)

        # Ensure that the data points are connected in the x direction
        if len(x_values) > 0:  # Ensure there's something to plot
            if any(x_errors) or any(y_errors):
                # Plot with error bars, using the assigned color consistently
                plt.errorbar(x_values, y_values, 
                             xerr=[x if x is not None else 0 for x in x_errors], 
                             yerr=[y if y is not None else 0 for y in y_errors], 
                             fmt='o', color=color, capsize=5, markersize=5, elinewidth=2,  
                             markerfacecolor=color, linewidth=1.5, label="_nolegend_")
                plt.plot(x_values, y_values, '-', color=color, linewidth=1.5, label=common_label)
            else:
                # Plot without error bars, using the assigned color consistently
                plt.plot(x_values, y_values, 'o-', color=color, markersize=5, 
                         markerfacecolor=color, linewidth=1.5, label=common_label)

    # Set the plot title and axis labels using the column names with larger font sizes
    plt.title(title, fontsize=20)  # Increased font size for title
    plt.xlabel(x_col_name, fontsize=20)  # Increased x-axis title font size to 18
    plt.ylabel(y_col_name, fontsize=20)  # Increased y-axis title font size to 18
    
    # Set tick label font size to make axis values bigger
    plt.xticks(fontsize=18)  # Increased x-axis value font size to 16
    plt.yticks(fontsize=18)  # Increased y-axis value font size to 16
    
    # Add faint grid lines
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # Calculate the axis limits: x-axis starts at -(0.1 * max(x_values))
    x_max = max(all_x_values) if all_x_values else 0
    plt.xlim(left=-0.1 * x_max)  # Start x-axis at -(0.1 * max(x_value))

    # Keep the y-axis as default (no changes to y_min)
    plt.ylim(auto=True)  # Default y-axis limits

    # Add a legend with a larger font size
    plt.legend(loc='best', fontsize=16)  # Increased legend font size to 16
    
    # Adjust the layout to fit everything properly
    plt.tight_layout()

    # Generate a dynamic filename if filename is not provided
    if save_dir:
        if filename is None:
            filename = f"{first_sample}_{x_col_name}_{y_col_name}.png"
        else:
            if not filename.endswith('.png'):
                filename += '.png'  # Append .png if not present

        save_path = os.path.join(save_dir, filename)

        # Ensure the directory exists
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)  # Create the directory if it doesn't exist

        # Save the plot to the specified directory
        plt.savefig(save_path, format='png', dpi=300)  # Save at 300 dpi for good resolution
        print(f"Plot saved to {save_path}")

    # Show the plot
    plt.show()


# In[ ]:





# In[5]:


def plot_selected_samples(data_dict, column_names, xy_columns, selected_sample_groups, error_columns, title, save_dir=None, filename=None):
    plt.figure(figsize=(16, 9))
    plt.rc('font', family='Arial')

    x_col_index, y_col_index = xy_columns
    x_err_col_index, y_err_col_index = error_columns
    x_col_name = column_names[x_col_index]
    y_col_name = column_names[y_col_index]

    all_x_values = []
    all_y_values = []
    first_sample = selected_sample_groups[0][0] if selected_sample_groups and selected_sample_groups[0] else "sample"

    def common_substring(strings):
        common = strings[0]
        for string in strings[1:]:
            match = SequenceMatcher(None, common, string).find_longest_match(0, len(common), 0, len(string))
            common = common[match.a: match.a + match.size]
        return common.strip("_")

    color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for group in selected_sample_groups:
        x_values, y_values, x_errors, y_errors = [], [], [], []
        group_labels = [key for key in group if key in data_dict]
        common_label = common_substring(group_labels)
        color = next(color_cycle)
        
        for key in group:
            if key in data_dict:
                values = data_dict[key]
                x_value = values[x_col_index] if len(values) > x_col_index else None
                y_value = values[y_col_index] if len(values) > y_col_index else None

                if x_value is not None and y_value is not None:
                    x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else None
                    y_err = values[y_err_col_index] if y_err_col_index is not None and len(values) > y_err_col_index else None
                    x_values.append(x_value)
                    y_values.append(y_value)
                    x_errors.append(x_err)
                    y_errors.append(y_err)
        
        all_x_values.extend(x_values)
        all_y_values.extend(y_values)

        if len(x_values) > 0:
            if any(x_errors) or any(y_errors):
                plt.errorbar(x_values, y_values, 
                             xerr=[x if x is not None else 0 for x in x_errors], 
                             yerr=[y if y is not None else 0 for y in y_errors], 
                             fmt='o', color=color, capsize=5, markersize=5, elinewidth=2,  
                             markerfacecolor=color, linewidth=1.5, label="_nolegend_")
                plt.plot(x_values, y_values, '-', color=color, linewidth=1.5, label=common_label)
            else:
                plt.plot(x_values, y_values, 'o-', color=color, markersize=5, 
                         markerfacecolor=color, linewidth=1.5, label=common_label)

    plt.title(title, fontsize=20)
    plt.xlabel(x_col_name, fontsize=20)
    plt.ylabel(y_col_name, fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    
    x_max = max(all_x_values) if all_x_values else 0
    plt.xlim(left=-0.1 * x_max)
    plt.ylim(auto=True)
    
    # Limit legend locations to the edges
    plt.legend(loc='upper right', fontsize=16)

    plt.tight_layout()

    if save_dir:
        if filename is None:
            filename = f"{first_sample}_{x_col_name}_{y_col_name}.png"
        else:
            if not filename.endswith('.png'):
                filename += '.png'

        save_path = os.path.join(save_dir, filename)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_path, format='png', dpi=300)
        print(f"Plot saved to {save_path}")

    plt.show()


# In[6]:


def plot_selected_samples(data_dict, column_names, xy_columns, selected_sample_groups, error_columns, title, save_dir=None, filename=None):
    import matplotlib.pyplot as plt
    import itertools
    from difflib import SequenceMatcher
    import os

    plt.figure(figsize=(16, 9))
    plt.rc('font', family='Arial')

    x_col_index, y_col_index = xy_columns
    x_err_col_index, y_err_col_index = error_columns
    x_col_name = column_names[x_col_index]
    y_col_name = column_names[y_col_index]

    all_x_values = []
    all_y_values = []
    all_x_errors = []
    all_y_errors = []
    first_sample = selected_sample_groups[0][0] if selected_sample_groups and selected_sample_groups[0] else "sample"

    def common_substring(strings):
        common = strings[0]
        for string in strings[1:]:
            match = SequenceMatcher(None, common, string).find_longest_match(0, len(common), 0, len(string))
            common = common[match.a: match.a + match.size]
        return common.strip("_")

    color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for group in selected_sample_groups:
        x_values, y_values, x_errors, y_errors = [], [], [], []
        group_labels = [key for key in group if key in data_dict]
        common_label = common_substring(group_labels)
        color = next(color_cycle)
        
        for key in group:
            if key in data_dict:
                values = data_dict[key]
                x_value = values[x_col_index] if len(values) > x_col_index else None
                y_value = values[y_col_index] if len(values) > y_col_index else None

                if x_value is not None and y_value is not None:
                    x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else None
                    y_err = values[y_err_col_index] if y_err_col_index is not None and len(values) > y_err_col_index else None
                    x_values.append(x_value)
                    y_values.append(y_value)
                    x_errors.append(x_err)
                    y_errors.append(y_err)
                    all_x_values.append(x_value)
                    all_y_values.append(y_value)
                    all_x_errors.append(x_err)
                    all_y_errors.append(y_err)

        if len(x_values) > 0:
            if any(x_errors) or any(y_errors):
                plt.errorbar(x_values, y_values, 
                             xerr=[x if x is not None else 0 for x in x_errors], 
                             yerr=[y if y is not None else 0 for y in y_errors], 
                             fmt='o', color=color, capsize=5, markersize=5, elinewidth=2,  
                             markerfacecolor=color, linewidth=1.5, label="_nolegend_")
                plt.plot(x_values, y_values, '-', color=color, linewidth=1.5, label=common_label)
            else:
                plt.plot(x_values, y_values, 'o-', color=color, markersize=5, 
                         markerfacecolor=color, linewidth=1.5, label=common_label)

    plt.title(title, fontsize=20)
    plt.xlabel(x_col_name, fontsize=20)
    plt.ylabel(y_col_name, fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    
    x_max = max(all_x_values) if all_x_values else 0
    plt.xlim(left=-0.1 * x_max)
    plt.ylim(auto=True)
    
    # Adjust legend location to avoid overlapping with data points or error bars
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()

    # Define corner regions (excluding bottom left)
    corner_regions = {
        'upper right': [(x_min + 0.7*(x_max - x_min), x_max), (y_min + 0.7*(y_max - y_min), y_max)],
        'upper left': [(x_min, x_min + 0.3*(x_max - x_min)), (y_min + 0.7*(y_max - y_min), y_max)],
        'lower right': [(x_min + 0.7*(x_max - x_min), x_max), (y_min, y_min + 0.3*(y_max - y_min))],
    }

    # Function to check if any data points or error bars are in the corner region
    def is_corner_empty(x_values, y_values, x_errors, y_errors, x_range, y_range):
        x_left, x_right = x_range
        y_bottom, y_top = y_range
        for x, y, xerr, yerr in zip(x_values, y_values, x_errors, y_errors):
            xerr = xerr if xerr is not None else 0
            yerr = yerr if yerr is not None else 0
            x_low = x - xerr
            x_high = x + xerr
            y_low = y - yerr
            y_high = y + yerr
            # Check if error bar overlaps with corner region
            if (x_high >= x_left) and (x_low <= x_right) and (y_high >= y_bottom) and (y_low <= y_top):
                return False  # Data point or error bar is in the corner
        return True  # No data points or error bars in the corner

    # Check each corner for data points
    available_corners = []
    for corner_name, (x_range, y_range) in corner_regions.items():
        corner_empty = is_corner_empty(
            all_x_values,
            all_y_values,
            [0 if x is None else x for x in all_x_errors],
            [0 if y is None else y for y in all_y_errors],
            x_range,
            y_range
        )
        if corner_empty:
            available_corners.append(corner_name)

    # Choose a corner
    if available_corners:
        legend_location = available_corners[0]  # Select the first available corner
    else:
        legend_location = 'upper right'  # Default to 'upper right' if all corners are occupied

    plt.legend(loc=legend_location, fontsize=16)

    plt.tight_layout()

    if save_dir:
        if filename is None:
            filename = f"{first_sample}_{x_col_name}_{y_col_name}.png"
        else:
            if not filename.endswith('.png'):
                filename += '.png'

        save_path = os.path.join(save_dir, filename)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_path, format='png', dpi=300)
        print(f"Plot saved to {save_path}")

    plt.show()


# In[7]:


def plot_selected_samples(data_dict, column_names, xy_columns, selected_sample_groups, error_columns, title, save_dir=None, filename=None):
    import matplotlib.pyplot as plt
    import itertools
    from difflib import SequenceMatcher
    import os

    plt.figure(figsize=(16, 9))
    plt.rc('font', family='Arial')

    x_col_index, y_col_index = xy_columns
    x_err_col_index, y_err_col_index = error_columns
    x_col_name = column_names[x_col_index]
    y_col_name = column_names[y_col_index]

    all_x_values = []
    all_y_values = []
    all_x_errors = []
    all_y_errors = []
    first_sample = selected_sample_groups[0][0] if selected_sample_groups and selected_sample_groups[0] else "sample"

    def common_substring(strings):
        common = strings[0]
        for string in strings[1:]:
            match = SequenceMatcher(None, common, string).find_longest_match(0, len(common), 0, len(string))
            common = common[match.a: match.a + match.size]
        return common.strip("_")

    color_cycle = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for group in selected_sample_groups:
        x_values, y_values, x_errors, y_errors = [], [], [], []
        group_labels = [key for key in group if key in data_dict]
        common_label = common_substring(group_labels)
        color = next(color_cycle)
        
        for key in group:
            if key in data_dict:
                values = data_dict[key]
                x_value = values[x_col_index] if len(values) > x_col_index else None
                y_value = values[y_col_index] if len(values) > y_col_index else None

                if x_value is not None and y_value is not None:
                    x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else None
                    y_err = values[y_err_col_index] if y_err_col_index is not None and len(values) > y_err_col_index else None
                    x_values.append(x_value)
                    y_values.append(y_value)
                    x_errors.append(x_err)
                    y_errors.append(y_err)
                    all_x_values.append(x_value)
                    all_y_values.append(y_value)
                    all_x_errors.append(x_err)
                    all_y_errors.append(y_err)

        if len(x_values) > 0:
            if any(x_errors) or any(y_errors):
                plt.errorbar(x_values, y_values, 
                             xerr=[x if x is not None else 0 for x in x_errors], 
                             yerr=[y if y is not None else 0 for y in y_errors], 
                             fmt='o', color=color, capsize=5, markersize=5, elinewidth=2,  
                             markerfacecolor=color, linewidth=1.5, label="_nolegend_")
                plt.plot(x_values, y_values, '-', color=color, linewidth=1.5, label=common_label)
            else:
                plt.plot(x_values, y_values, 'o-', color=color, markersize=5, 
                         markerfacecolor=color, linewidth=1.5, label=common_label)

    plt.title(title, fontsize=20)
    plt.xlabel(x_col_name, fontsize=20)
    plt.ylabel(y_col_name, fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    
    x_max = max(all_x_values) if all_x_values else 0
    plt.xlim(left=-0.1 * x_max)
    plt.ylim(auto=True)
    
    # Adjust legend location to avoid overlapping with data points or error bars
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()

    # Define right-side regions for legend placement
    right_regions = {
        'upper right': (1.0, 0.75),
        'center right': (1.0, 0.5),
        'lower right': (1.0, 0.25),
    }

    def is_region_empty(x_values, y_values, x_errors, y_errors, x_range, y_range):
        x_left, x_right = x_range
        y_bottom, y_top = y_range
        for x, y, xerr, yerr in zip(x_values, y_values, x_errors, y_errors):
            xerr = xerr if xerr is not None else 0
            yerr = yerr if yerr is not None else 0
            x_low = x - xerr
            x_high = x + xerr
            y_low = y - yerr
            y_high = y + yerr
            if (x_high >= x_left) and (x_low <= x_right) and (y_high >= y_bottom) and (y_low <= y_top):
                return False
        return True

    # Check each right-side region for data points
    for region_name, (x_anchor, y_anchor) in right_regions.items():
        x_range = (x_max - 0.1 * (x_max - x_min), x_max)
        y_range = (y_min + (y_anchor - 0.1) * (y_max - y_min), y_min + (y_anchor + 0.1) * (y_max - y_min))
        if is_region_empty(
            all_x_values,
            all_y_values,
            [0 if x is None else x for x in all_x_errors],
            [0 if y is None else y for y in all_y_errors],
            x_range,
            y_range
        ):
            plt.legend(bbox_to_anchor=(x_anchor, y_anchor), loc='center left', fontsize=16)
            break
    else:
        # Default to upper right if no empty region is found
        plt.legend(loc='upper right', fontsize=16)

    plt.tight_layout()

    if save_dir:
        if filename is None:
            filename = f"{first_sample}_{x_col_name}_{y_col_name}.png"
        else:
            if not filename.endswith('.png'):
                filename += '.png'

        save_path = os.path.join(save_dir, filename)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_path, format='png', dpi=300)
        print(f"Plot saved to {save_path}")

    plt.show()


# In[8]:


data_cols


# In[9]:


selected_samples = [['617_32000', '617_32000_2dpa'],
                    ['MA754', 'MA754_1dpa', 'MA754_2dpa'],
                    ['MA957', 'MA957_1dpa', 'MA957_2dpa']]


# In[12]:


plot_selected_samples(data_dict, data_cols, (1,4) , selected_samples, (None, 5), "CSDS vs Dose", save_dir, 'CSDSirrall')


# In[50]:


def plot_dual_y_axes(data_dict, column_names, xy_columns, selected_sample_groups, error_columns, title, save_dir=None, filename=None, secondary_y_axis=None, show_legend=False):
    """
    Plots specified samples from the data_dict using the column indices for x and y axes,
    with optional error bars, and connects grouped samples with lines along the x-axis.
    Supports an optional secondary Y-axis for plotting an additional parameter.
    
    Parameters:
    - data_dict: Dictionary containing sample data.
    - column_names: List of column names corresponding to data_dict columns.
    - xy_columns: Tuple (x_col_index, y_col_index) indicating the indices of the columns for the primary x and y axes.
    - selected_sample_groups: List of lists, where each sublist contains specific keys (samples) to be connected by lines.
    - error_columns: Tuple (x_err_col_index, y_err_col_index) indicating the columns for x and y error bars, or None if no error bars.
    - title: Title for the plot.
    - save_dir: Directory where the plot will be saved as a PNG file (optional).
    - filename: Custom filename for saving the plot. If not provided, a dynamic name is generated.
    - secondary_y_axis: Optional. Tuple (y2_col_index, y2_err_col_index) to define a second Y-axis with its own data and error.
    - show_legend: Boolean to control display of legend. Default is True.
    """
    
    # Set the figure size to 16:9 and the font to Arial
    fig, ax1 = plt.subplots(figsize=(16, 9))
    plt.rc('font', family='Arial')

    # Primary axis settings
    x_col_index, y_col_index = xy_columns
    x_err_col_index, y_err_col_index = error_columns
    x_col_name = column_names[x_col_index]
    y_col_name = column_names[y_col_index]

    # Secondary Y-axis (if provided)
    ax2 = ax1.twinx() if secondary_y_axis else None
    if secondary_y_axis:
        y2_col_index, y2_err_col_index = secondary_y_axis
        y2_col_name = column_names[y2_col_index]

    # Color cycles for primary and secondary Y-axes
    primary_color = 'tab:blue'
    secondary_color = 'tab:orange'
    color_cycle = itertools.cycle([primary_color, secondary_color])

    all_x_values = []
    all_y_values = []

    # Helper function to find the longest common substring
    def common_substring(strings):
        common = strings[0]
        for string in strings[1:]:
            match = SequenceMatcher(None, common, string).find_longest_match(0, len(common), 0, len(string))
            common = common[match.a: match.a + match.size]
        return common.strip("_")

    # Plot data for primary Y-axis
    for group in selected_sample_groups:
        x_values, y_values, x_errors, y_errors = [], [], [], []
        group_labels = [key for key in group if key in data_dict]
        common_label = common_substring(group_labels)
        
        # Assign color based on Y-axis
        color = primary_color if not secondary_y_axis else next(color_cycle)
        
        for key in group:
            if key in data_dict:
                values = data_dict[key]
                x_value = values[x_col_index] if len(values) > x_col_index else None
                y_value = values[y_col_index] if len(values) > y_col_index else None

                if x_value is not None and y_value is not None:
                    x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else np.nan
                    y_err = values[y_err_col_index] if y_err_col_index is not None and len(values) > y_err_col_index else np.nan
                    x_values.append(x_value)
                    y_values.append(y_value)
                    x_errors.append(x_err)
                    y_errors.append(y_err)
        
        all_x_values.extend(x_values)
        all_y_values.extend(y_values)

        # Plot primary Y-axis data with consistent color
        ax1.errorbar(x_values, y_values, xerr=x_errors, yerr=y_errors, fmt='o', color=color, capsize=5, markersize=5, elinewidth=2, markerfacecolor=color, linewidth=1.5, label=common_label)
        ax1.plot(x_values, y_values, '-', color=color, linewidth=1.5)

    # Plot data for secondary Y-axis if provided
    if secondary_y_axis:
        for group in selected_sample_groups:
            x_values, y2_values, x_errors, y2_errors = [], [], [], []
            group_labels = [key for key in group if key in data_dict]
            common_label = common_substring(group_labels)
            
            color = secondary_color

            for key in group:
                if key in data_dict:
                    values = data_dict[key]
                    x_value = values[x_col_index] if len(values) > x_col_index else None
                    y2_value = values[y2_col_index] if len(values) > y2_col_index else None

                    if x_value is not None and y2_value is not None:
                        x_err = values[x_err_col_index] if x_err_col_index is not None and len(values) > x_err_col_index else np.nan
                        y2_err = values[y2_err_col_index] if y2_err_col_index is not None and len(values) > y2_err_col_index else np.nan
                        x_values.append(x_value)
                        y2_values.append(y2_value)
                        x_errors.append(x_err)
                        y2_errors.append(y2_err)
            
            # Plot secondary Y-axis data with its color
            ax2.errorbar(x_values, y2_values, xerr=x_errors, yerr=y2_errors, fmt='s', color=color, capsize=5, markersize=5, elinewidth=2, markerfacecolor=color, linewidth=1.5, label=common_label)
            ax2.plot(x_values, y2_values, '-', color=color, linewidth=1.5)

    # Labels and title
    ax1.set_xlabel(x_col_name, fontsize=20)
    ax1.set_ylabel(y_col_name, fontsize=20, color=primary_color)
    if secondary_y_axis:
        ax2.set_ylabel(y2_col_name, fontsize=20, color=secondary_color)
        ax2.tick_params(axis='y', labelcolor=secondary_color)
    plt.title(title, fontsize=20)

    # Legends and layout
    if show_legend:
        ax1.legend(loc='upper left', fontsize=16)
        if secondary_y_axis:
            ax2.legend(loc='upper right', fontsize=16)
    plt.tight_layout()

    # Save the plot
    if save_dir:
        if filename is None:
            filename = f"{first_sample}_{x_col_name}_{y_col_name}.png"
        else:
            if not filename.endswith('.png'):
                filename += '.png'
        save_path = os.path.join(save_dir, filename)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_path, format='png', dpi=300)
        print(f"Plot saved to {save_path}")

    plt.show()


# In[51]:


data_cols


# In[82]:


selected_samples = [['MA957', 'MA957_1dpa', 'MA957_2dpa']]


# In[85]:


plot_dual_y_axes(
    data_dict=data_dict,
    column_names=data_cols,
    xy_columns=(1, 6),                # Primary X and Y columns
    selected_sample_groups=selected_samples,
    error_columns=(None, 7),          # Error columns for primary Y-axis
    title="MA957",
    save_dir=save_dir,
    filename='rhoMMA957',
    secondary_y_axis=(8, 9)           # Secondary Y-axis: column 6 for data and column 7 for error
)


# In[ ]:




