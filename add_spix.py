          




import os
import pandas as pd
import numpy as np
from uncertainties import unumpy as unp
from mymodules import computemedian_spix, crossmatch
from mymodules import nearestneighbour
import glob

# Folder paths
folder = 'data/tiles/'
fol2 = 'data/alpha/'
final_cat = glob.glob('/tmp/vlass_cat/catalogue_output_files/VLASS*CIR_components.csv')[0]

# Load the master CSV file
master_df = pd.read_csv(final_cat)

# Initialize an empty DataFrame to store all matched data
all_matched_data = pd.DataFrame()

# Iterate over tiles
tiles = os.listdir(folder)
for t in tiles:
    files = os.listdir(folder + t)
    
    # Filter and process the CSV files
    file_names = [i.replace('.csv', '') for i in files if i.endswith('.csv') and not i.endswith('.vlad.csv')]
    for j in file_names:
        
        image_name = folder + t + '/' + j + '.fits'
        alpha_name = fol2 + j.replace('image.pbcor.tt0', 'alpha') + '.fits'
        alpha_err = fol2 + j.replace('image.pbcor.tt0', 'alpha.error') + '.fits'
        
        # Extract subtile from the filename
        subtile = j.split('.')[4]
        
        # Load data from the tile CSV
        file_data = pd.read_csv(folder + t + '/' + j + '.csv')
        
        # Filter master_df based on tile and subtile before performing crossmatch
        filtered_master_df = master_df[(master_df['Tile'] == t) & (master_df['Subtile'] == subtile)]
        try:    
            # Compute spectral index if not already present
            columns_to_drop = [col for col in file_data.columns if col.startswith('Alpha')]
            file_data = file_data.drop(columns=columns_to_drop)
            if 'Alpha' not in file_data.columns:
                size = 3
                _, _, weighted_spix_and_err, _, _ = computemedian_spix(file_data, 'RA', 'DEC', size, image_name, [25, 70], 'off', alpha_name, alpha_err)
                
                # Update the file_data DataFrame with new columns
                file_data['Alpha'] = [np.nan if unp.isnan(i) else round(i.nominal_value, 2) for i in weighted_spix_and_err]
                file_data['Alpha_err'] = [np.nan if unp.isnan(i) else round(i.std_dev, 4) for i in weighted_spix_and_err]
                
                # Save the updated DataFrame
                file_data.to_csv(folder + t + '/' + j + '.csv', index=False)
            
            if not filtered_master_df.empty:
                # Perform crossmatch and update the master DataFrame
                matched_data, _ = crossmatch(filtered_master_df, file_data, 'RA', 'DEC', 'RA', 'DEC', searchrad=0.2)  # Adjust searchrad as needed
                columns_to_drop = [col for col in matched_data.columns if col.endswith('_y')]
                matched_data = matched_data.drop(columns=columns_to_drop)
                
                # Rename columns by removing '_x' suffix
                matched_data.columns = [col.replace('_x', '') if '_x' in col else col for col in matched_data.columns]

                # Append the matched data to the all_matched_data DataFrame
                all_matched_data = pd.concat([all_matched_data, matched_data], ignore_index=True)
            
        except Exception as e:
            # Add new columns with NaN values to the filtered_master_df
            filtered_master_df['Alpha'] = np.nan
            filtered_master_df['Alpha_err'] = np.nan
            # filtered_master_df['Alpha_weighted'] = np.nan
            # filtered_master_df['Alpha_weighted_err'] = np.nan
            all_matched_data = pd.concat([all_matched_data, filtered_master_df], ignore_index=True)
            print(f"Error in spectral index_est: {str(e)}")

# At this point, all_matched_data contains all the matched rows. You can save or further process it.
secat=all_matched_data
secat['Alpha_quality_flag']=[1]*len(secat)

index=secat[(secat['S_Code']!='S') | (secat['DC_Maj']>1)].index

secat.loc[index,'Alpha_quality_flag']=3

nn=nearestneighbour(secat,'RA','DEC',5)


index=secat[secat.Component_name.isin(nn[0].Component_name_x)].index
secat.loc[index,'Alpha_quality_flag']=secat.loc[index,'Alpha_quality_flag']+4

# secat[secat.Alpha_quality_flag==7]
secat=secat.drop(['dist2d','indexx'],axis=1)
secat=secat.rename(columns={"Alpha_err": "Alpha_err_(TT1/TT0)"})


secat['Alpha_err_corr']=(secat['Alpha_err_(TT1/TT0)']**2+(0.144+0.08/np.log10(secat['Alpha_err_(TT1/TT0)']))**2)**0.5

index=secat[secat.Total_flux<1.3].index
secat.loc[index,'Alpha_err_corr']=secat.loc[index,'Alpha_err_(TT1/TT0)']

index=secat[(secat.Alpha<-50)].index
secat.loc[index,'Alpha_err_corr']=np.nan


secat.loc[index,'Alpha_err_(TT1/TT0)']=np.nan
secat.loc[index,'Alpha_quality_flag']=np.nan
secat.loc[index,'Alpha']=np.nan

lis=list(secat.columns)[0:-3]
lis.append(list(secat.columns)[-1])
lis.append(list(secat.columns)[-2])
lis.append(list(secat.columns)[-3])

secat=secat[lis]
secat.to_csv(final_cat.split('/')[-1], index=False)


