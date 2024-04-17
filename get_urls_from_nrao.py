

import sys
import requests
from bs4 import BeautifulSoup as bs

def get_soup(URL, filetype,imagetype='I'):
    # Retrieve the HTML from the URL using the requests library
    html = bs(requests.get(URL).text, 'html.parser')
    
    # Create an empty list to store the file URLs
    file_urls = []
    
    # Loop through all links on the page
    for link in html.find_all('a'):
        # Get the URL of the file
        file_link = link.get('href')
        # Check if the file type matches the specified type
        if imagetype=='I':
            if all(i in file_link for i in filetype) and 'rms' not in file_link and 'alpha' not in file_link and '.tt1' not in file_link:
            # If it matches, append the URL to the list
                file_urls.append(file_link)
        if imagetype=='alpha':
            if all(i in file_link for i in filetype) and 'pbcor' not in file_link and 'alpha' in file_link and '.tt1' not in file_link:
            # If it matches, append the URL to the list
                file_urls.append(file_link)
    
    # Return the list of file URLs
    return file_urls#[1:]  # slice to exclude the first element

def main():
    # Get the URL and mode from the command line arguments
    if len(sys.argv) != 4:
        print("Usage: python get_urls_from_nrao.py <URL (e.g. https://archive-new.nrao.edu/vlass/se_continuum_imaging/VLASS2.1)> <mode> <imagetype>")
        print("  mode: 'w' for write or 'a' for append")
        print("  imagetype: 'I' for total intensity images or 'alpha' for spectral index and error images")
        return
    URL = sys.argv[1]
    mode = sys.argv[2]
    imagetype=sys.argv[3]
    # Get a list of tile URLs that match the specified file type
    tile_list = get_soup(URL, ['/', 'T'])
    if 'QA_REJECTED/' in tile_list:
        tile_list.remove('QA_REJECTED/')
    # Open the CSV file in the appropriate mode
    with open('manifest.csv', mode) as fout:
        # If we're overwriting the file, write the header
        if mode == 'w':
            fout.write('file\n')

        # Loop through each tile URL
        for tile_url in tile_list:
            print(f"Processing tile URL: {tile_url}")
            try:
                # Get a list of subtile URLs that match the specified file type
                url_t = URL + '/' + tile_url
                subtile_list = get_soup(url_t, ['/', 'VLASS'])[1:]
                # print(get_soup(url_t, ['/', 'VLASS']))
                # Loop through each subtile URL
                for subtile_url in subtile_list:
                    try:
                    # Get the URL of the FITS file and append it to a list
                        url_st = url_t + subtile_url
                        # print(url_st)
                        # print(get_soup(url_st, ['.fits'],imagetype))
                        fits_url = get_soup(url_st, ['.fits'],imagetype)
                        for i in  fits_url:
                            fout.write(url_st + i + '\n')
                    except Exception as e:
                        print(e)                        
                    
            except Exception as e:
                # If there was an error, skip this tile
                print(e)

if __name__ == '__main__':
    main()
