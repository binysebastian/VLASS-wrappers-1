import os
import glob
import pandas as pd

# Directory containing the files
directory = 'catalogue_output_files'

# Pattern to match the files
pattern = os.path.join(directory, '*_components.csv')

# List of files matching the pattern, sorted by modification time
files = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)

# Select the most recently modified file
if files:
    latest_file = files[0]
    print(f"Processing file: {latest_file}")

    # Read the CSV file
    data = pd.read_csv(latest_file)

    # Filter the data
    data = data[data.Duplicate_flag < 2]

    # Save the filtered data to a new file
    new_filename = latest_file.replace('.csv', '') + '_duplicate_free.csv'
    data.to_csv(new_filename, index=False)
    print(f"Saved filtered data to {new_filename}")
else:
    print("No matching files found.")
