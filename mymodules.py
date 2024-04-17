#!/usr/bin/env python
# coding: utf-8

# In[131]:


from matplotlib import pyplot as plt
import numpy as np
from astropy.coordinates import ICRS
from astropy import units as u
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord as SC
from astropy.io import ascii
from astropy.table import Table
from astropy import table
from astroquery.cadc import Cadc
import pandas as pd
import pandas as pd
from scipy.stats import norm

from astropy.wcs import WCS
from astropy import wcs
import os
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord as SC
from matplotlib import cm
#import casatasks as cts
import math
import uncertainties as un
import seaborn as sns
from uncertainties import umath
from uncertainties import unumpy
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator
from matplotlib.colors import ListedColormap
from scipy import stats
cadc = Cadc()


# In[132]:


def skycoord(data,RA1,DEC1):
    #Extracting the coordinates from the catalogs
    RA=data[RA1]
    DEC=data[DEC1]
    coordinates=[ICRS(i*u.degree,j*u.degree) for i,j in zip(RA,DEC)]
#     coordinates=[str(i.ra)+' '+str(i.dec) for i in coo]
    return coordinates


# ## Function to find offset from VLASS image header

# In[133]:


def find_offset(coords,*subtiles):
    radius = '0.0014 degree'
    if subtiles:
        coords=zip(coords,subtiles[0])
#         print('here')
    offsets_RA,offsets_DEC=[],[]
    for i in coords:
        try:
            if subtiles:
#                 print(i[0],i[1])
                results = cadc.query_region(i[0], radius, collection='VLASS')
                #print('Point here')
                image_list = cadc.get_image_list(results, i[0], '0.0014 degree')
                # print(image_list)
                image_list=[k for k in image_list if not 'VLASS2.' in k]
                #image_list=[k for k in image_list if j in k]
                # print(image_list)
                # subtile_of_images=[fits.open(k)[0].header['FILNAM05'] for k in image_list]
                image_list=[k for k in image_list if i[1][0:7] in k and i[1][8:14] in k]   
                image_list=image_list[0]
                # print(j[0:7]+j[8:14])

            # print(header)
            else:
                image_list=i
            
            header=fits.open(image_list)[0].header
            RA_orig=float(header['HISTORY'][-3].split()[3])
            Dec_orig=float(header['HISTORY'][-2].split()[3])
            RA_correct=header['CRVAL1']
            Dec_correct=header['CRVAL2']

            
            

            offsetstr=header['HISTORY'][-1]
            ind0=offsetstr.find('(')+1
            ind1=offsetstr[::-1].find(')')+1
            offset_dir=[float(i)/3600.0 for i in offsetstr[ind0:-ind1].replace('/cos(CRVAL2),','').split()]
            offset_dir[0]=offset_dir[0]/np.cos(np.deg2rad(Dec_orig))
            # print(offset_dir)
            offset_RA=RA_correct-RA_orig
            offset_DEC=Dec_correct-Dec_orig
            # offsets_RA.append(offset_RA)
            # offsets_DEC.append(offset_DEC)
            offsets_RA.append(offset_dir[0])
            offsets_DEC.append(offset_dir[1])

            print('success')
        except Exception as e:
            print(e)
            offsets_RA.append(np.nan)
            offsets_DEC.append(np.nan)

    #print(header)
    return offsets_RA,offsets_DEC
        


# ## Reading the subtile and the data catalog from CIRADA website

# In[134]:

def shorten_cat(data,RA,DEC,image):
    hdu1=fits.open(image)
    w=wcs.WCS(hdu1[0].header)

    if hdu1[0].data.ndim==4:
        size=len(hdu1[0].data[0][0][0])
        sky = w.pixel_to_world(0,0,1, 1)
        RA_max,DEC_min=sky[0].ra.deg,sky[0].dec.deg
        sky = w.pixel_to_world(size,size,1,1)
    elif hdu1[0].data.ndim==3:
        size=len(hdu1[0].data[0][0])
        sky = w.pixel_to_world(0,0,1,1)
        RA_max,DEC_min=sky[0].ra.deg,sky[0].dec.deg
        sky = w.pixel_to_world(size,size,1,1)
    elif hdu1[0].data.ndim==2:
        size=len(hdu1[0].data[0])
        sky = w.pixel_to_world(0,0,1,1)
        RA_max,DEC_min=sky[0].ra.deg,sky[0].dec.deg
        sky = w.pixel_to_world(size,size,1,1)
        
    
    
    
    RA_min,DEC_max=sky[0].ra.deg,sky[0].dec.deg
    # print(RA_min,RA_max,DEC_min,DEC_max)
    # hdu1.close()
    data_short=data[((data[RA]<RA_max) & (data[RA]>RA_min) & (data[DEC]<DEC_max) & (data[DEC]>DEC_min))]
    if RA_max<RA_min:
        data_short=data[((((data[RA]<RA_max) & (data[RA]>0)) | ((data[RA]>RA_min) & (data[RA]<360))) & (data[DEC]<DEC_max) & (data[DEC]>DEC_min))]
    # print(RA_min,RA_max)
    return data_short


# 
# - Set a bounding box around every source given the size, data catalog with positions, 
# - display alpha images, 
# - find median and weighted median for every subimage. 

# In[229]:



def skycoord2(data,RA1,DEC1):
    #Extracting the coordinates from the catalogs
    RA=list(map(str,data[RA1]))
    DEC=list(map(str,data[DEC1]))
    coordinates=SC([i+' '+j for i,j in zip(RA,DEC)],unit=(u.deg,u.deg))
    return coordinates





def crossmatch(data_1,data_2,RA1,DEC1,RA2,DEC2,searchrad):  #the radius is in arcsec
    
    coordinates_1=skycoord2(data_1,RA1,DEC1)
    coordinates_2=skycoord2(data_2,RA2,DEC2)
    
    # identifying the nearest source in 2nd catalog to every source in 1st catalog
    index,dist2d,dist3d = coordinates_1.match_to_catalog_sky(coordinates_2)
    
    
    # identifying sources within a search radius of searchradius
    ifst=dist2d.arcsec<searchrad
    match,index1,dist2d1=np.where(ifst,True,False),np.where(ifst,index,np.nan),np.where(ifst,dist2d,np.nan)
    
    
    # adding the index and the distance column to the first catalog
    # data_1.add_column(1,name='index')
    # data_1.add_column(1,name='dist2d')
    data_1['dist2d']=dist2d1.arcsec #*u.arcsec
    data_1['indexx']=index1
    
    # adding the index column to the second catalog
    # data_2.add_column(1,name='index')
    index_ql=[i for i in range(len(data_2))]
    data_2['indexx']=index_ql
    
    # Joining the first and second tables based on the index 
    data_matched=pd.merge(data_1,data_2,on='indexx',how='inner')
    
    
    #keeping only the relevant columns
    # data_matched.keep_columns([RA1,DEC1,'dist2d','Total_flux_source','Total_flux','Isl_Total_flux'])
    # data_matched['Isl_Total_flux']=data_matched['Isl_Total_flux']*1e3 *u.mJy
    return data_matched,match





def converttopix(data,RA,DEC,box_size,image): 
    if type(image)==str:
        hdu1=fits.open(image)
    else:
        hdu1=image    
 # Converting the box_size from arcsec of pixels
    w=wcs.WCS(hdu1[0].header)
    image_data=hdu1[0].data
    
    
    pix_width_arcsec=Angle(hdu1[0].header['CDELT2'],u.deg).arcsec
    #box_size=5 # u.arcsec
    box_size_pixel=int(box_size/pix_width_arcsec)
#     print(pix_width_arcsec)
#     box_size_pixel=3  # forcing the pixel size to 5x5 pixel
    box_size_pixel=int(box_size_pixel/2)
    
    # Finding the array indices and pixel values from source world coordinates
    array_coord=[list(w.world_to_array_index_values(i,j,1,1)) for i,j in zip(data[RA],data[DEC])]
    pix_coord=[list(w.wcs_world2pix(i,j,1,1,1)) for i,j in zip(data[RA],data[DEC])]
    
    return array_coord,pix_coord,box_size_pixel


def findconf(alpha,*args):
    
    import scipy
  
    # Generate some data for this 
    # demonstration.
    alpha_val=[i.nominal_value for i in alpha]
    data = alpha[np.logical_not(np.isnan(alpha_val))] 
    median=np.median(data)
    if args:
        median=args[0][-1]
    CI_p=median.nominal_value+median.std_dev
    CI_n=median.nominal_value-median.std_dev
    #np.random.normal(170, 10, 250)
    data=[i.nominal_value for i in data]
    
    # Fit a normal distribution to
    # the data:
    # mean and standard deviation
    mu, std = norm.fit(data) 
      
    # Plot the histogram.
    # hist=plt.hist(data, bins=10)
      
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
      
    # plt.plot(x, p, 'k', linewidth=2)
    # title = "Fit Values: {:.2f} and {:.2f}".format(mu, std)
    # plt.axvline(x=np.median(data),linestyle='--',color='red')
    # plt.axvline(x=np.mean(data),linestyle='--',color='black')
    
    # plt.axvline(x=CI_p,linestyle='--',color='green')
    # plt.axvline(x=CI_n,linestyle='--',color='green')
    # plt.title(title)
      
    # plt.show()
    gauss_approx=scipy.stats.norm.cdf(CI_p, loc=mu, scale=std)-scipy.stats.norm.cdf(CI_n, loc=mu, scale=std)

    simple_calc=len([i for i in data if (i> CI_n and i< CI_p)])/len(data)
    return simple_calc,gauss_approx
    
def visualizationandmedian_spix(full_data,RA,DEC,box_size,image,figsize,centering='off',*alpha_im): 
    
    hdu1=fits.open(image) #opening Itt0 image
    
    if alpha_im:
        hdu2=fits.open(alpha_im[0]) #opening alpha image
        alpha_data=hdu2[0].data
        try:
            label1=hdu2[0].header['BUNIT']
        except:
            pass    
        if len(alpha_im)>1:
            hdu3=fits.open(alpha_im[1]) #opening alpha err image
            alpha_err_data=hdu3[0].data
            # uncomment for masking based on alpha values
            # for i in range(len(alpha_err_data[0][0])):
            #     for j in range(len(alpha_err_data[0][0][0])):
            #         if alpha_err_data[0][0][i][j]>0.2:
            #             alpha_err_data[0][0][i][j]=np.nan
            #             alpha_data[0][0][i][j]=np.nan

    # shortening the data catalog to that with the range of image RA and Dec
    data=shorten_cat(full_data,RA,DEC,image)
    image_data=hdu1[0].data
    
    # Converting the box_size from arcsec of pixels
    array_coord,pix_coord,box_size_pixel=converttopix(data,RA,DEC,box_size,image)

    matched_len=len(array_coord)
    rem=np.sign(matched_len%3)+int(matched_len/3)
    
    
    fig = plt.figure(1)
    median=[]
    weighted_median=[]
    conflev=[]
    for j in range(matched_len):
        # Getting the peak flux location and changing the center
        indexx=[int(array_coord[j][2])-3,int(array_coord[j][2])+3,int(array_coord[j][3])-3,int(array_coord[j][3])+3]
        
#         print()
        if centering =='on':
            try:
                cent=np.where(image_data[0][0]==np.max(image_data[0][0][indexx[0]:indexx[1]+1,indexx[2]:indexx[3]+1]))
            except:
                print(indexx)
                print(image_data[0][0][indexx[0]:indexx[1]+1,indexx[2]:indexx[3]+1])
    #         print(array_coord[j][2])
            try:
                array_coord[j][2],array_coord[j][3]=int(cent[0]),int(cent[1])
            except Exception as e:
                # print(e)
                print((cent[0]),(cent[1]))
    
        ind=[int(array_coord[j][2])-box_size_pixel,int(array_coord[j][2])+box_size_pixel,int(array_coord[j][3])-box_size_pixel,int(array_coord[j][3])+box_size_pixel]
        weight=image_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
        if alpha_im:
            alpha=alpha_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            if len(alpha_im)>1:
                alpha_err=alpha_err_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
        
        #find weighted median
    

        vars()['ax'+str(j+1)] = fig.add_subplot(rem, 3, (j + 1))
        if alpha_im:
            vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(alpha)
        else:
            vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(weight,cmap=cm.jet)
            try:
                label1=hdu1[0].header['BUNIT']
            except:
                pass      
        try:                  
            plt.colorbar(vars()["im" + str(j+1)],label=label1)
        except:
            pass            
        try:
            names=['Component_name','src_name']
            columnname=[i for i in names if i in full_data.columns][0]
            # column=
            vars()['ax'+str(j+1)].set_title(full_data[columnname][full_data.index[j]],fontsize=20)
        except Exception as e:
            # print(e)
            pass
        try:
        # if alpha_im:
            ind=[int(array_coord[j][2])-box_size_pixel,int(array_coord[j][2])+box_size_pixel,int(array_coord[j][3])-box_size_pixel,int(array_coord[j][3])+box_size_pixel]
            weight=image_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            alpha=alpha_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            if len(alpha_im)>1:
                alpha_err=alpha_err_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
                alpha=np.array([un.ufloat(i,j) for i,j in zip (np.ravel(alpha),np.ravel(alpha_err))])
                # print(alpha)
            
        
            #find weighted median
            
            
            # vars()['ax'+str(j+1)] = fig.add_subplot(rem, 3, (j + 1))
            # vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(alpha)
            # plt.colorbar(vars()["im" + str(j+1)],label=hdu2[0].header['BUNIT'])
        
            df=pd.DataFrame([np.ravel(weight),np.ravel(alpha)]).transpose() # making a dataframe from 1D arrays of the flux densities and alpha 
            if len(alpha_im)>1:
                df=pd.DataFrame([np.ravel(weight),alpha]).transpose() # making a dataframe from 1D arrays of the flux densities and alpha 
            df.columns=['Weight','Alpha']

            df.sort_values('Alpha', inplace=True)
            msk=pd.isna(df['Alpha'])
            df=df[~msk] # dropping columns with nan alpha values
            
            cumsum = df.Weight.cumsum()
            cutoff = df.Weight.sum() / 2.0
            # print(cumsum,cutoff)
            weighted_median.append(df.Alpha[cumsum >= cutoff].iloc[0])
            median.append(np.nanmedian(alpha))
            if len(alpha_im)>1:
                conflev.append(findconf(alpha,weighted_median))
            # print(df)
            
        except Exception as e:
            weighted_median.append(np.nan)
            median.append(np.nan)
            conflev.append(np.nan)
            print(e)
            
    fig.set_size_inches(figsize[0],figsize[1], forward=True)
    plt.show()
    
#     Finding the median of the region centred around array_coord and have a size of 2*box_size_pixel
#     median=[np.nanmedian(image_data[0][0]) for j in range(matched_len)]
#     In case the image passed is total intensity, we use casa to estimate the statistics.
    try:
        print(nn)
        stat=[ctks.imstat(image,box=str(int(pix_coord[j][0])-box_size_pixel)+','+str(int(pix_coord[j][1])-box_size_pixel)+','+str(int(pix_coord[j][0])+box_size_pixel)+','+str(int(pix_coord[j][1])+box_size_pixel)) for j in range(matched_len)]
        df=np.nan
    except:
        stat=median

    if len(alpha_im)>1:   
        return stat,data,weighted_median,alpha,conflev
    else:
        return stat,data,weighted_median

# - Compare the median alpha with the alpha in the catalog by over plotting the individual histograms and the difference 

# In[230]:
import seaborn as sns
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator
from matplotlib.colors import ListedColormap
from scipy import stats
# In[230]:
import seaborn as sns
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter, AutoMinorLocator
from matplotlib.colors import ListedColormap
from scipy import stats

def compare_flux(flux_1,flux_2,bins,label1,label2,xlabel1):

    a=min(min(flux_1),min(flux_2))
    b=max(max(flux_1),max(flux_2))
    bins=np.histogram(np.hstack((a,b)), bins=bins)[1] #get the bin edges


#plt.ioff()
    plt.ion()
    plt.rc('font',family='sans-serif')
    plt.rc('font', serif='Helvetica')
    plt.rcParams['ytick.major.pad']='5'
    plt.rcParams['xtick.major.pad']='5'
    
    plt.rcParams['mathtext.fontset']='custom'
    plt.rcParams['mathtext.rm']='Bitstream Vera Sans'
    plt.rcParams['mathtext.it']='Bitstream Vera Sans:italic'
    plt.rcParams['mathtext.bf']='Bitstream Vera Sans:bold'
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.set_style("darkgrid")
    #sns.set_palette('muted')
    pal = sns.set_palette("muted")
    #hexcolors=pal.as_hex()
    #sns.set_style("darkgrid", {'axes.linewidth': 1.5, 'axes.grid': True,'axes.edgecolor': '.15','grid.color': '.8', 'xtick.color': '.15', 'ytick.color': '.15','xtick.direction': u'in','ytick.direction': u'in'})
    cmap = ListedColormap(sns.color_palette(n_colors=10))
    colors = cmap.colors
    
    hist_kwsa = dict(histtype='stepfilled',linewidth=3,alpha=0.05)
    hist_kwsb = dict(histtype='step',linewidth=3,alpha=1.0)


    plt.clf()
    fig1 = plt.figure()
    sub1 = fig1.add_subplot(1,1,1)
    diff=flux_2-flux_1
    
    #mygraph=sns.distplot(flux_FIRST_CNSS_all_plot, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[0], label=r'CNSS ($\flux_{\rm med} = ' + flux_FIRST_CNSS_all_median_str+')$')
    mygraph=sns.distplot(flux_1, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[3], label=label1+r' ($median = ' + str(round(np.median(flux_1),2))+'$\n$'+',mad = ' + str(round(stats.median_absolute_deviation(flux_1),2))+')$',bins=bins)
    mygraph=sns.distplot(flux_2, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[4],label=label2+r' ($median = ' + str(round(np.median(flux_2),2))+'$\n$'+',mad = ' + str(round(stats.median_absolute_deviation(flux_2),2))+')$',bins=bins)
    mygraph=sns.distplot(diff, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[5],label=r'Difference ($median = ' + str(round(np.median(diff),2))+'$\n$'+',mad = ' + str(round(stats.median_absolute_deviation(diff),2))+')$',bins=bins)


    leg=plt.legend(loc='upper right', fontsize=10, frameon=True)
    #myxlabel=plt.xlabel(r'$\alpha^{3.0}_{1.4}$',fontsize=36, fontweight='normal',labelpad=5)
    myxlabel=plt.xlabel(xlabel1,fontsize=26, fontweight='normal',labelpad=5)
    myylabel=plt.ylabel(r'$N$',fontsize=26, fontweight='normal',labelpad=5)
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=9, prune='lower'))
    plt.tick_params(which='major', axis='both',length=0, labelsize=16)
    plt.tick_params(which='minor', axis='both',length=0, labelsize=16)
    #plt.title('Spectral Indices with FIRST', fontsize=26, fontweight='normal', y=1.015)
    plt.xlim(min(min(flux_1),min(flux_2))-1,max(max(flux_1),max(flux_2))+1)
    

    
    # plt.hist(diff,alpha=0.5,label='difference',bins=bins)
    print(np.count_nonzero(np.isnan(diff)),len(flux_1),len(flux_2))
    # except Exception as e:
    #     print(e)
    # plt.legend()
    # plt.savefig('image.eps')
    plt.show()
#     return diff,frac_err


def compare_spix(alpha_1,alpha_2,bins,label1,label2):
    a=min(min(alpha_1),min(alpha_2))
    b=2
    bins=np.histogram(np.hstack((a,b)), bins=bins)[1] #get the bin edges


#plt.ioff()
    plt.ion()
    plt.rc('font',family='sans-serif')
    plt.rc('font', serif='Helvetica')
    plt.rcParams['ytick.major.pad']='5'
    plt.rcParams['xtick.major.pad']='5'
    
    plt.rcParams['mathtext.fontset']='custom'
    plt.rcParams['mathtext.rm']='Bitstream Vera Sans'
    plt.rcParams['mathtext.it']='Bitstream Vera Sans:italic'
    plt.rcParams['mathtext.bf']='Bitstream Vera Sans:bold'
    
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    sns.set_style("darkgrid")
    #sns.set_palette('muted')
    pal = sns.set_palette("muted")
    #hexcolors=pal.as_hex()
    #sns.set_style("darkgrid", {'axes.linewidth': 1.5, 'axes.grid': True,'axes.edgecolor': '.15','grid.color': '.8', 'xtick.color': '.15', 'ytick.color': '.15','xtick.direction': u'in','ytick.direction': u'in'})
    cmap = ListedColormap(sns.color_palette(n_colors=10))
    colors = cmap.colors
    
    hist_kwsa = dict(histtype='stepfilled',linewidth=3,alpha=0.05)
    hist_kwsb = dict(histtype='step',linewidth=3,alpha=1.0)


    plt.clf()
    fig1 = plt.figure()
    sub1 = fig1.add_subplot(1,1,1)
    diff=alpha_2-alpha_1
    #mygraph=sns.distplot(alpha_FIRST_CNSS_all_plot, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[0], label=r'CNSS ($\alpha_{\rm med} = ' + alpha_FIRST_CNSS_all_median_str+')$')
    mygraph=sns.distplot(alpha_1, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[3], label=label1+r' ($\alpha = ' + str(round(np.median(alpha_1),2))+'$\n$'+',\delta = ' + str(round(stats.median_absolute_deviation(alpha_1),2))+')$',bins=bins)
    mygraph=sns.distplot(alpha_2, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[4],label=label2+r' ($\alpha = ' + str(round(np.median(alpha_2),2))+'$\n$'+',\delta = ' + str(round(stats.median_absolute_deviation(alpha_2),2))+')$',bins=bins)
    mygraph=sns.distplot(diff, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsb, color=colors[5],label=r'Difference ($\alpha = ' + str(round(np.median(diff),2))+'$\n$'+',\delta = ' + str(round(stats.median_absolute_deviation(diff),2))+')$',bins=bins)

#     mygraph=sns.distplot(alpha_1, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsa)
#     mygraph=sns.distplot(alpha_2, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsa)
#     mygraph=sns.distplot(diff, hist=True, kde=False, norm_hist=False, hist_kws=hist_kwsa)

    leg=plt.legend(loc='upper right', fontsize=10, frameon=True)
    #myxlabel=plt.xlabel(r'$\alpha^{3.0}_{1.4}$',fontsize=36, fontweight='normal',labelpad=5)
    myxlabel=plt.xlabel(r'$\alpha$',fontsize=36, fontweight='normal',labelpad=5)
    myylabel=plt.ylabel(r'$N$',fontsize=26, fontweight='normal',labelpad=5)
    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=9, prune='lower'))
    plt.tick_params(which='major', axis='both',length=0, labelsize=16)
    plt.tick_params(which='minor', axis='both',length=0, labelsize=16)
    #plt.title('Spectral Indices with FIRST', fontsize=26, fontweight='normal', y=1.015)
    plt.xlim(min(min(alpha_1),min(alpha_2))-1,7)
    print(np.count_nonzero(np.isnan(diff)),len(alpha_1),len(alpha_2))
    plt.show()





# 
# - Print the flux-weighted median alpha, median alpha, and median alpha from the cube

# In[231]:

def spectral_plot(data):
    for i in data.index:
        plt.scatter(data.freq[i],data.stokesI_x[i])
        farr=np.linspace(min(data.freq[i]),max(data.freq[i]),100)
        ref_freq=np.array(data.reffreq_pol[i])
        plt.errorbar(data.freq[i],data.stokesI_x[i],yerr=data.stokesI_error[i],linestyle="None")
        plt.plot(farr,data.I_modelcoeff0[i]*(farr/ref_freq)**data.I_modelcoeff1[i],label=r'$\alpha$='+str(round(data.I_modelcoeff1[i],2)))
        plt.legend(prop={'size': 13})
        plt.title(data.src_name[i])
        plt.yscale('log',base=10) 
        plt.xscale('log',base=10) 
        plt.xlabel('Frequency (Hz)',size=15)
        plt.ylabel('Flux density (Jy/beam)',size=15)
        plt.show()



def make_cutout_using_rmtools(data_ext,RA,DEC,size,fits_file,folder):
    array_coord,pix_coord,box_size_pixel=converttopix(data_ext,RA,DEC,size,fits_file)
    a=[(int(round(float(pix_coord[i][0]),0)),box_size_pixel,int(round(float(pix_coord[i][1]),0)),box_size_pixel) for i in range(len(pix_coord))]
    b=[('rmtools_extractregion -c -o  '+fits_file+' '+folder+'cutout_'+str(round(data_ext.RA[data_ext.index[i]],2))+'_'+str(round(data_ext['DEC'][data_ext.index[i]],2))+'.fits '+ str(a[i]).replace('(','').replace(')','')).replace(',','') for i in range(len(pix_coord))]
    for i in b:
        os.system(i)
    

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u

def create_cutout(fits_file, ra, dec, box_size, output_file):
    # Open the FITS file
    with fits.open(fits_file) as hdulist:
        header = hdulist[0].header
        data = hdulist[0].data

    # Create a WCS object from the header
    wcs = WCS(header)

    # Extract the 2D data from the input data and adjust WCS
    if data.ndim == 4:
        data = data[0, 0]
        wcs = wcs.dropaxis(2)
        wcs = wcs.dropaxis(2)
    elif data.ndim == 3:
        data = data[0]
        wcs = wcs.dropaxis(2)
    elif data.ndim != 2:
        raise ValueError(f"Unsupported number of dimensions: {data.ndim}")

    # Convert the given RA and DEC to a SkyCoord object
    target_position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

    # Define the cutout size
    size = u.Quantity((box_size, box_size), u.arcsec)

    # Create the cutout
    cutout = Cutout2D(data, target_position, size, wcs=wcs)

    # Update the original header with the cutout WCS
    cutout_header = cutout.wcs.to_header()
    for key in cutout_header:
        header[key] = cutout_header[key]

    # Save the cutout to a new FITS file
    cutout_hdu = fits.PrimaryHDU(data=cutout.data, header=header)
    cutout_hdu.writeto(output_file, overwrite=True)

    print(f"Cutout saved to {output_file}")




def visualization(files,fol): 
    length=len(files)
    nocol=2
    # cent_p=8192
    # box_s=50
    rem=np.sign(length%nocol)+int(length/nocol)
    fig = plt.figure()
    for i,j in zip(files,range(length)):
        try:
            hdu=fits.open(fol+i)[0]
            wcs = WCS(hdu.header)
            # level=[2**(i)*3e-5 for i in range(10)]
            vars()['ax'+str(j+1)] = fig.add_subplot(rem, nocol, (j + 1),projection=wcs,slices=('x', 'y',1,1))
            # vars()['ax'+str(j+1)].contour(hdu.data[0][0][cent_p-box_s-1:cent_p+box_s,cent_p-box_s-1:cent_p+box_s], levels=level,colors='white', alpha=0.5)
            if hdu.data.ndim==2:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data)#,cmap=cm.jet)
            elif  hdu.data.ndim==3:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0])#,cmap=cm.jet)               
            elif  hdu.data.ndim==4:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0][0])#,cmap=cm.jet)               
            try:
                plt.colorbar(vars()["im" + str(j+1)],label=hdu.header['BUNIT'])
            except:
                plt.colorbar(vars()["im" + str(j+1)])              
            vars()['ax'+str(j+1)].set_title(hdu.header['OBJECT'],fontsize=20)
            plt.rcParams.update({'font.size': 22})
            # hdu.close()
        except Exception as e:
            print(e)
            # err.append(i)
    # fig.set_size_inches(int(length/nocol)*10,10*nocol, forward=True)    
    
    # print(err)




def dat_to_csv(inputf,outputf):
    import csv
    
    with open(inputf) as dat_file, open(outputf, 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
    
        for line in dat_file:
            row = [field.strip() for field in line.split()]
    #         if len(row) == 6 and row[3] and row[4]:
            csv_writer.writerow(row)   

    return outputf             




def nearestneighbour(data_1,RA1,DEC1,searchrad):  #the radius is in arcsec
    
    coordinates_1=skycoord2(data_1,RA1,DEC1)
    index,dist2d,dist3d = coordinates_1.match_to_catalog_sky(coordinates_1,nthneighbor=2)
    
    
    ifst=dist2d.arcsec<searchrad
    match,index1,dist2d1=np.where(ifst,True,False),np.where(ifst,index,np.nan),np.where(ifst,dist2d,np.nan)
    
    
    # adding the index and the distance column to the first catalog
    data_2=data_1.copy()
    data_1['dist2d']=dist2d1.arcsec #*u.arcsec
    data_1['indexx']=index1
    
    # adding the index column to the second catalog
    index_ql=[i for i in range(len(data_2))]
    data_2['indexx']=index_ql
    
    # Joining the first and second tables based on the index 
    data_matched=pd.merge(data_1,data_2,on='indexx',how='inner')
    # 
    
    return data_matched,match




def nthnearestneighbour(data_1,RA1,DEC1,searchrad,n):  #the radius is in arcsec
    
    coordinates_1=skycoord2(data_1,RA1,DEC1)
    index,dist2d,dist3d = coordinates_1.match_to_catalog_sky(coordinates_1,nthneighbor=n)
    
    
    ifst=dist2d.arcsec<searchrad
    match,index1,dist2d1=np.where(ifst,True,False),np.where(ifst,index,np.nan),np.where(ifst,dist2d,np.nan)
    
    
    # adding the index and the distance column to the first catalog
    data_2=data_1.copy()
    data_1['dist2d']=dist2d1.arcsec #*u.arcsec
    data_1['indexx']=index1
    
    # adding the index column to the second catalog
    index_ql=[i for i in range(len(data_2))]
    data_2['indexx']=index_ql
    
    # Joining the first and second tables based on the index 
    data_matched=pd.merge(data_1,data_2,on='indexx',how='inner')
    # 
    
    return data_matched,match


def computemedian_spix(full_data,RA,DEC,box_size,image,figsize,centering='off',*alpha_im): 
    
    hdu1=fits.open(image) #opening Itt0 image
    
    if alpha_im:
        hdu2=fits.open(alpha_im[0]) #opening alpha image
        alpha_data=hdu2[0].data
        try:
            label1=hdu2[0].header['BUNIT']
        except:
            pass    
        if len(alpha_im)>1:
            hdu3=fits.open(alpha_im[1]) #opening alpha err image
            alpha_err_data=hdu3[0].data
            # uncomment for masking based on alpha values
            # for i in range(len(alpha_err_data[0][0])):
            #     for j in range(len(alpha_err_data[0][0][0])):
            #         if alpha_err_data[0][0][i][j]>0.2:
            #             alpha_err_data[0][0][i][j]=np.nan
            #             alpha_data[0][0][i][j]=np.nan

    # shortening the data catalog to that with the range of image RA and Dec
    data=shorten_cat(full_data,RA,DEC,image)
    image_data=hdu1[0].data
    
    # Converting the box_size from arcsec of pixels
    array_coord,pix_coord,box_size_pixel=converttopix(data,RA,DEC,box_size,image)

    matched_len=len(array_coord)
    rem=np.sign(matched_len%3)+int(matched_len/3)
    
    
    # fig = plt.figure(1)
    median=[]
    weighted_median=[]
    conflev=[]
    for j in range(matched_len):
        # Getting the peak flux location and changing the center
        indexx=[int(array_coord[j][2])-3,int(array_coord[j][2])+3,int(array_coord[j][3])-3,int(array_coord[j][3])+3]
        
#         print()
        if centering =='on':
            try:
                cent=np.where(image_data[0][0]==np.max(image_data[0][0][indexx[0]:indexx[1]+1,indexx[2]:indexx[3]+1]))
            except:
                print(indexx)
                print(image_data[0][0][indexx[0]:indexx[1]+1,indexx[2]:indexx[3]+1])
    #         print(array_coord[j][2])
            try:
                array_coord[j][2],array_coord[j][3]=int(cent[0]),int(cent[1])
            except Exception as e:
                print(e)
                # print((cent[0]),(cent[1]))
    
        ind=[int(array_coord[j][2])-box_size_pixel,int(array_coord[j][2])+box_size_pixel,int(array_coord[j][3])-box_size_pixel,int(array_coord[j][3])+box_size_pixel]
        weight=image_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
        if alpha_im:
            alpha=alpha_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            if len(alpha_im)>1:
                alpha_err=alpha_err_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
        
        try:
        # if alpha_im:
            ind=[int(array_coord[j][2])-box_size_pixel,int(array_coord[j][2])+box_size_pixel,int(array_coord[j][3])-box_size_pixel,int(array_coord[j][3])+box_size_pixel]
            weight=image_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            alpha=alpha_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
            if len(alpha_im)>1:
                alpha_err=alpha_err_data[0][0][ind[0]:ind[1]+1,ind[2]:ind[3]+1]
                alpha=np.array([un.ufloat(i,j) for i,j in zip (np.ravel(alpha),np.ravel(alpha_err))])
                # print(alpha)
            
            
            #find weighted median
            
            
            
            df=pd.DataFrame([np.ravel(weight),np.ravel(alpha)]).transpose() # making a dataframe from 1D arrays of the flux densities and alpha 
            if len(alpha_im)>1:
                df=pd.DataFrame([np.ravel(weight),alpha]).transpose() # making a dataframe from 1D arrays of the flux densities and alpha 
            df.columns=['Weight','Alpha']
    
            df.sort_values('Alpha', inplace=True)
            msk=pd.isna(df['Alpha'])
            df=df[~msk] # dropping columns with nan alpha values
            
            cumsum = df.Weight.cumsum()
            cutoff = df.Weight.sum() / 2.0
            # print(cumsum,cutoff)
            weighted_median.append(df.Alpha[cumsum >= cutoff].iloc[0])
            
            if len(alpha_im)>1:
                conflev.append(findconf(alpha,weighted_median))
                median.append(np.median(alpha[~unumpy.isnan(alpha)]))
            else:
                median.append(np.nanmedian(alpha))    
            # print(df)
            
        except Exception as e:
            if alpha_im:
                if len(alpha_im)>1:
                    nan=float('nan')
                    weighted_median.append(un.ufloat(nan,nan))
                    median.append(un.ufloat(nan,nan))
                    # conflev.append(np.nan)

            else:                    
                weighted_median.append(np.nan)
                median.append(np.nan)

            conflev.append(np.nan)
            print(e)
            
    # 
#     Finding the median of the region centred around array_coord and have a size of 2*box_size_pixel
    try:
        print(nn)
        stat=[ctks.imstat(image,box=str(int(pix_coord[j][0])-box_size_pixel)+','+str(int(pix_coord[j][1])-box_size_pixel)+','+str(int(pix_coord[j][0])+box_size_pixel)+','+str(int(pix_coord[j][1])+box_size_pixel)) for j in range(matched_len)]
        df=np.nan
    except:
        stat=median

    if len(alpha_im)>1:   
        return stat,data,weighted_median,alpha,conflev
    else:
        return stat,data,weighted_median

# - Compare the median alpha with the alpha in the catalog by over plotting the individual histograms and the difference 

# In[230]:


import pandas as pd
def plotaitoff(data,ra,dec,colr):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import numpy as np
    if isinstance(data, pd.DataFrame):
        datafile=data
    else:    
        datafile=pd.read_csv(data)
    ra_random = np.array(datafile[ra]) * u.degree
    dec_random = np.array(datafile[dec])* u.degree

    c = SkyCoord(ra=ra_random, dec=dec_random, frame='icrs')
    
    
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,4.2))
    plt.subplot(111, projection="aitoff")
    # plt.title("Aitoff projection of our random data")
    plt.grid(True)
    plt.plot(ra_rad, dec_rad, 'o', markersize=0.05, alpha=0.3,color=colr)
    plt.subplots_adjust(top=0.95,bottom=0.0)
    plt.show()





def compare_spix_with_sns(data,alpha1,alpha2,*fluxclass):
    difference=data[alpha1]-data[alpha2]
    diff_string='('+alpha1+')'+'-'+'('+alpha2+')'
    data[diff_string]=difference
    cmap = ListedColormap(sns.color_palette(n_colors=5))
    colors = cmap.colors
    if fluxclass:
        g=sns.pairplot(data[[alpha1,alpha2,diff_string,fluxclass[0]]],  height=3,hue='Flux_class',corner=True,plot_kws={"s": 6})
        g.map_lower(sns.kdeplot, levels=4, color=".2")
    else:
        g=sns.pairplot(data[[alpha1,alpha2,diff_string]],  height=3,corner=True,plot_kws={"s": 6}) 
        g.map_lower(sns.kdeplot, levels=4,color='black')
    # g.axes

    for i,j in enumerate(g.axes.ravel()):
#         if i in [1,2,5]:
#             j.set_visible(False)
        if  i==3:
            x=np.linspace(np.min(data[alpha1])-0.5,np.max(data[alpha1])+0.5,100)
            j.plot(x,x,color='black',linestyle='--')
    #         j.set_linestyle("--")
        if i in [6,7]:
            j.axhline(0,color='black')
        if not fluxclass:
          if i in [0,4,8]:  
            if i==0:
                median_val=round(np.median(data[alpha1]),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[alpha1]),2)
                se_val=round(sig_val/np.sqrt(len(data[alpha1])),2)
            if i==4:
                median_val=round(np.median(data[alpha2]),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[alpha2]),2)
                se_val=round(sig_val/np.sqrt(len(data[alpha1])),2)
            if i==8:
                median_val=round(np.median(data[diff_string]),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[diff_string]),2)
                se_val=round(sig_val/np.sqrt(len(data[alpha1])),2)
                
            j.annotate(r'$\alpha = ' + str(median_val)+'\pm'+str(se_val)+'$,\n$'+'\sigma = ' + str(sig_val)+'$', xy = (0, 0), xytext = (0.9, 0.9), size=10, textcoords = 'axes    fraction', ha = 'left', va = 'center', color='k', weight='normal')
            j.axvline(median_val,color='black')
            j.axvline(median_val+sig_val,color='black',linestyle='--')
            j.axvline(median_val-sig_val,color='black',linestyle='--')
        else:
          if i in [0,4,8]:  
            if i==0:
                median_val=round(np.median(data[data[fluxclass[0]]=='Flux > 10mJy']),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha1]),2)
                se_val=round(sig_val/np.sqrt(len(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha1])),2)
            if i==4:
                median_val=round(np.median(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha2]),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha2]),2)
                se_val=round(sig_val/np.sqrt(len(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha1])),2)
            if i==8:
                median_val=round(np.median(data[data[fluxclass[0]]=='Flux > 10mJy'][diff_string]),2)
                sig_val=round( 1.48* stats.median_absolute_deviation(data[data[fluxclass[0]]=='Flux > 10mJy'][diff_string]),2)
                se_val=round(sig_val/np.sqrt(len(data[data[fluxclass[0]]=='Flux > 10mJy'][alpha1])),2)
                
            j.annotate(r'$\alpha = ' + str(median_val)+'\pm'+str(se_val)+'$,\n$'+'\sigma = ' + str(sig_val)+'$', xy = (0, 0), xytext = (0.9, 0.9), size=10, textcoords = 'axes    fraction', ha = 'left', va = 'center', color='k', weight='normal')
            j.axvline(median_val,color='black')
            j.axvline(median_val+sig_val,color='black',linestyle='--')
            j.axvline(median_val-sig_val,color='black',linestyle='--')

    g.savefig(r'Comparison_alpha.eps')
#     g.title('Comparison between alpha values')        
            
            
    #         j.plot(x,x)
        
#     if not 
    

from matplotlib.colors import LogNorm


def visualization(files,fol): 
    length=len(files)
    nocol=9
    # cent_p=8192
    # box_s=50
    rem=np.sign(length%nocol)+int(length/nocol)
    fig = plt.figure()
    for i,j in zip(files,range(length)):
        try:
            hdu=fits.open(fol+i)[0]
            wcs = WCS(hdu.header)
            # level=[2**(i)*3e-5 for i in range(10)]
            # 

            # vars()['ax'+str(j+1)].contour(hdu.data[0][0][cent_p-box_s-1:cent_p+box_s,cent_p-box_s-1:cent_p+box_s], levels=level,colors='white', alpha=0.5)
            if hdu.data.ndim==2:
                vars()['ax'+str(j+1)] = fig.add_subplot(rem, nocol, (j + 1),projection=wcs)
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data,cmap=cm.jet)
            elif  hdu.data.ndim==3:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0],cmap=cm.jet)#,cmap=cm.jet)               
            elif  hdu.data.ndim==4:
                vars()['ax'+str(j+1)] = fig.add_subplot(rem, nocol, (j + 1),projection=wcs,slices=('x', 'y',1,1))
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0][0],cmap=cm.jet)#,cmap=cm.jet)               
            # try:
                # plt.colorbar(vars()["im" + str(j+1)],label=hdu.header['BUNIT'])
            # except:
                # plt.colorbar(vars()["im" + str(j+1)])              
            # vars()['ax'+str(j+1)].set_title(hdu.header['OBJECT'],fontsize=20)
            plt.rcParams.update({'font.size': 10})
            # hdu.close()
        except Exception as e:
            print(e)
            # err.append(i)
    fig.set_size_inches(int(length/nocol)*10,10*nocol, forward=True)    
    
    # print(err)



def VLASS_cutout(data,RA,DEC,vlassver,pixsize,outputdir):
    link1='https://www.legacysurvey.org/viewer/cutout.fits?ra='
    link2='&dec='
    link3='&layer=vlass'+vlassver+'&size='+str(pixsize)
    cutouts=[link1+str(data.loc[j][RA])+link2+str(data.loc[j][DEC])+link3 for i,j in zip(data,data.index)]
    os.system('rm -rf '+outputdir)
    os.system('mkdir '+outputdir)
    cwd=os.getcwd()
    os.chdir(''+outputdir)
    [os.system('wget "'+i+'"') for i in cutouts]
    os.chdir(cwd)









def visualization_image(images): 
    length=len(images)
    nocol=3
    # cent_p=8192
    # box_s=50
    rem=np.sign(length%nocol)+int(length/nocol)
    fig = plt.figure()
    for i,j in zip(images,range(length)):
        try:
            hdu=i[0]
            wcs = WCS(hdu.header)
            # level=[2**(i)*3e-5 for i in range(10)]
            vars()['ax'+str(j+1)] = fig.add_subplot(rem, nocol, (j + 1),projection=wcs,slices=('x', 'y',1,1))
            # vars()['ax'+str(j+1)].contour(hdu.data[0][0][cent_p-box_s-1:cent_p+box_s,cent_p-box_s-1:cent_p+box_s], levels=level,colors='white', alpha=0.5)
            if hdu.data.ndim==2:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data,cmap=cm.jet)
            elif  hdu.data.ndim==3:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0],cmap=cm.jet)               
            elif  hdu.data.ndim==4:
                vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(hdu.data[0][0],cmap=cm.jet)               
            try:
                plt.colorbar(vars()["im" + str(j+1)],label=hdu.header['BUNIT'])
            except:
                plt.colorbar(vars()["im" + str(j+1)])              
            vars()['ax'+str(j+1)].set_title(hdu.header['OBJECT'],fontsize=4*nocol)
            plt.rcParams.update({'font.size': 4*nocol})
            # hdu.close()
        except Exception as e:
            print(e)
            # err.append(i)
    # if 
    fig.set_size_inches(np.ceil(length/nocol)*10,5*nocol, forward=True)    
    
    # print(err)

def VLASS_cutout_cadc(data,RA,DEC,epoch,size,name):
    coords=skycoord(data,RA,DEC)
    
    # epoch='VLASS2'
    # qlorse='se'
    from astropy import units as u
    from astroquery.cadc import Cadc
    cadc = Cadc()
    # coords = '01h45m07.5s +23d18m00s'
    radius = size*u.arcsec
    # results = cadc.query_region(coords, radius, collection='VLASS')
    
    # image_list = cadc.get_image_list(results, coords, radius)
    # print(image_list)   
    
    for j,jj,kk in zip(coords,data[name],data['Subtile_x']):
        try:
            imageneeded=[]
            images = cadc.get_images(j, radius, collection='VLASS')
            
            for i in range(len(images)):
                if images[i][0].header['FILNAM05']==kk and  not images[i][0].header['FILNAM01']=='s13_0':
                    # if epoch==images[i][0].header['FILNAM01'] and qlorse==images[i][0].header['FILNAM03']:
                    # if epoch==images[i][0].header['FILNAM01']:
                    imageneeded.append(images[i])
                    versions=images[i][0].header['FILNAM02']
                    images[i][0].header['OBJECT']='QL_'+images[i][0].header['FILNAM01']+'_'+jj[12:]
                    
                elif images[i][0].header['FILNAM06']==kk:
                    # if epoch==images[i][0].header['FILNAM02'] and qlorse==images[i][0].header['FILNAM04'] and 'alpha' not in images[i][0].header['FILNAM12']:
                    if 'alpha' not in images[i][0].header['FILNAM12'] and 'tt1' not in images[i][0].header['FILNAM13']:  
                        imageneeded.append(images[i])
                        versions=images[i][0].header['FILNAM03']
                        images[i][0].header['OBJECT']='SE_'+images[i][0].header['FILNAM02']+'_'+jj[12:]
            visualization_image(imageneeded)
        except Exception as e:
            print(e)        
    return imageneeded

from matplotlib import colors
def visualization_array_10by10(images,names,namesave): 
    length=len(images)
    nocol=10
    rem=np.sign(length%nocol)+int(length/nocol)
    fig = plt.figure()
    for i,j in zip(images,range(length)):
        try:
            cent_p=int(i.shape[0]/2)
            box_s=50
            vars()['ax'+str(j+1)] = fig.add_subplot(rem, nocol, (j + 1))
            vars()["im" + str(j+1)] = vars()['ax'+str(j+1)].imshow(i[cent_p-box_s-1:cent_p+box_s,cent_p-box_s-1:cent_p+box_s],cmap=cm.bwr,norm=colors.PowerNorm(1/2))
            vars()['ax'+str(j+1)].text(5, 5, names[j], fontsize=100, bbox={'facecolor': 'white', 'pad': 20})

            vars()["ax" + str(j+1)].axis('off')
        except Exception as e:
            print(e)
    fig.set_size_inches(nocol*10,np.ceil(length/nocol)*10, forward=True) 
    
    plt.tight_layout()
    fig.savefig(namesave)
    plt.show()

    
import glob
def make2darray_fromfits(folder):
    if folder[-1]=='/':
        names=glob.glob(folder+'*.fits')
    else:
        names=glob.glob(folder+'/*.fits')
    names.sort()
    im_array=[]
    
    # print(names)
    for ii in names:
        hdu=fits.open(ii)[0]
        if hdu.data.ndim==2:
            out=hdu.data
        elif  hdu.data.ndim==3:
            out=hdu.data[0]
        elif  hdu.data.ndim==4:
            out=hdu.data[0][0]
        im_array.append(out)

    return im_array,names
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
import pandas as pd

def crossmatch_using_vizier(cat, ra_col, dec_col, vizcat, searchrad):
    catalogs = {
        "tgss": "J/A+A/598/A78/table3",
        "nvss": "VIII/65/nvss",
        "first": "VIII/92/first14",
        "vlassql1": "J/ApJS/255/30/comp"
    }

    # Check if the requested Vizier catalog exists
    if vizcat not in catalogs:
        return "Catalog not found"

    # Create a SkyCoord object from the catalog
    coords = SkyCoord(cat[ra_col], cat[dec_col], unit=(u.deg, u.deg), frame='icrs')

    # Query Vizier catalog
    vizier_catalog = catalogs[vizcat]
    v = Vizier(columns=['**'], row_limit=-1)
    matched_results = []

    for index, coord in enumerate(coords):
        result = v.query_region(coord, radius=searchrad*u.arcsec, catalog=vizier_catalog)
        if result:
            for table in result:
                # Convert to Pandas DataFrame
                df = table.to_pandas()

                # Add original catalog data to the DataFrame
                for col in cat.columns:
                    df[col] = cat.iloc[index][col]

                matched_results.append(df)

    # Concatenate all matched results into a single DataFrame
    if matched_results:
        return pd.concat(matched_results, ignore_index=True)
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no matches were found

