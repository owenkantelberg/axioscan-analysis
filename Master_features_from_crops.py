#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:00:15 2021
@author: Mathew & Expanded by Owen
"""

import numpy as np
from skimage.io import imread
import matplotlib.pyplot as plt
from skimage import filters,measure
from PIL import Image
import pandas as pd
from scipy.spatial import distance
import csv
from pathlib import Path
#import tifffile
#import czifile


# Function to load images:
    
def load_image(toload):
    
    image=imread(toload)

    return image

def FileName(string):
    if string == "/Nuc.tif":
        return "DAPI"
    elif string == "/Ab.tif":
        return "Antibody"
    elif string == "/Apt.tif":
        return "Aptamer"
    else:
        return filename

# Function to make area histogram from region_props_table

def area_hist(title,areas):
    plt.hist(areas, bins = 50,range=[0,100], rwidth=0.7,color='blue', edgecolor='black', linewidth=0.5)
    plt.xlabel('Area (pixels)')
    plt.ylabel('Number of Features')
    plt.title(title + " - Pixels/Feature")
    plt.savefig(parent_directory + "/" + title + " - PixelsPerFeature.pdf")
  
    
def area_hist_from_csv(title,csv_file,feature):
    case_no = parent_directory[35:39]
    df = np.loadtxt(parent_directory + "/" + csv_file, delimiter="\t", skiprows=1, usecols=1)
    plt.hist(df, bins = 100,range=[10,2000], rwidth=0.7,color='blue', edgecolor='black', linewidth=0.5)
    plt.ylim((0, 1000))
    if title == 'Aptamer':
        plt.ylim((0, 200))
    plt.xlabel('Area (pixels)')
    plt.ylabel('Number of Features')
    plt.title('Case: ' + case_no + '     ' + title + " - Pixels/Feature      Features: " + str(feature))
    plt.savefig(parent_directory + "/SD" + case_no + "_All crops " + title + " - PixelsPerFeature.pdf")
    plt.show()



# Threshold image using otsu method and output the filtered image along with the threshold value applied:
    
def threshold_image_otsu(input_image):
    threshold_value=filters.threshold_otsu(input_image)    
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Define set-threshold for apt

def threshold_image_apt(input_image):
    threshold_value=30000
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Define set-threshold for ab

def threshold_image_ab(input_image):
    threshold_value=10000  
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Label and count the features in the thresholded image:
def label_image(input_image):
    labelled_image=measure.label(input_image)
    number_of_features=labelled_image.max()
 
    return number_of_features,labelled_image
    
# Function to show the particular image:
def show(input_image):
    plt.imshow(input_image,cmap="Reds")
    plt.show()

# Take a labelled image and the original image and measure intensities, sizes etc.
def analyse_labelled_image(labelled_image,original_image):
    measure_image=measure.regionprops_table(labelled_image,intensity_image=original_image,properties=('area','perimeter','centroid','orientation','major_axis_length','minor_axis_length','mean_intensity','max_intensity'))
    measure_dataframe=pd.DataFrame.from_dict(measure_image)
    #area_hist_(measure_image['area'])
    return measure_dataframe

# This is to look at coincidence purely in terms of pixels

def coincidence_analysis_pixels(binary_image1,binary_image2):
    pixel_overlap_image=binary_image1&binary_image2         
    pixel_overlap_count=pixel_overlap_image.sum()
    pixel_fraction=pixel_overlap_image.sum()/binary_image1.sum()
    
    return pixel_overlap_image,pixel_overlap_count,pixel_fraction

# Look at coincidence in terms of features. Needs binary image input 

def feature_coincidence(binary_image1,binary_image2):
    number_of_features,labelled_image1=label_image(binary_image1)          # Labelled image is required for this analysis
    coincident_image=binary_image1 & binary_image2        # Find pixel overlap between the two images
    coincident_labels=labelled_image1*coincident_image   # This gives a coincident image with the pixels being equal to label
    coinc_list, coinc_pixels = np.unique(coincident_labels, return_counts=True)     # This counts number of unique occureences in the image
    # Now for some statistics
    total_labels=labelled_image1.max()
    total_labels_coinc=len(coinc_list)
    fraction_coinc=total_labels_coinc/total_labels
    
    # Now look at the fraction of overlap in each feature
    # First of all, count the number of unique occurances in original image
    label_list, label_pixels = np.unique(labelled_image1, return_counts=True)
    fract_pixels_overlap=[]
    for i in range(len(coinc_list)):
        overlap_pixels=coinc_pixels[i]
        label=coinc_list[i]
        total_pixels=label_pixels[label]
        fract=1.0*overlap_pixels/total_pixels
        fract_pixels_overlap.append(fract)
    
    
    # Generate the images
    coinc_list[0]=1000000   # First value is zero- don't want to count these. 
    coincident_features_image=np.isin(labelled_image1,coinc_list)   # Generates binary image only from labels in coinc list
    coinc_list[0]=0
    non_coincident_features_image=~np.isin(labelled_image1,coinc_list)  # Generates image only from numbers not in coinc list.    
    return coinc_list,coinc_pixels,fraction_coinc,coincident_features_image,non_coincident_features_image,fract_pixels_overlap

# Function to measure minimum distances between two sets of data
def minimum_distance(measurements1,measurements2):
    s1 = measurements1[["centroid-0","centroid-1"]].to_numpy()
    s2 = measurements2[["centroid-0","centroid-1"]].to_numpy()
    minimum_lengths=distance.cdist(s1,s2).min(axis=1)
    return minimum_lengths


# Function to create master .csv files for each channel's region_props data
def master_file_region_props(name):
    with open(parent_directory+'/'+name+'.csv', 'w', newline='') as file:
        csv.writer(file, delimiter='\t')
        csv.writer(file, delimiter='\t').writerow(['area','perimeter','centroid','orientation','major_axis_length','minor_axis_length','mean_intensity','max_intensity'])
# Function to create master .csv file compiling coincidence data for all crops    
def master_file_coincidence():
    with open(parent_directory+'/master_coincidence.csv', 'w', newline='') as file:
        csv.writer(file, delimiter='\t')
        csv.writer(file, delimiter='\t').writerow(['Apt pixel fraction coincident with Ab pixels','Ab pixel fraction coincident with Apt pixels','Fraction of Apt features coincident with Ab features','Average Apt Ab feature overlap'])#,'Aptamer_min_dist_nuc','Antibody_min_dist_nuc','Nuc_min_dist_to_nuc'])
    with open(parent_directory+'/master_features.csv', 'w', newline='') as file:
        csv.writer(file, delimiter='\t')
        csv.writer(file, delimiter='\t').writerow(['apt_features','ab_features','nuc_features'])

def append_coincidence_master():
    avg_feature_overlap_apt=(sum(ab_fraction_pixels_overlap)/len(ab_fraction_pixels_overlap))
    with open(parent_directory+'//master_coincidence.csv', 'a', newline='') as file:
        csv.writer(file, delimiter='\t').writerow([apt_pixel_fraction,ab_pixel_fraction,apt_fraction_coinc,avg_feature_overlap_apt])#,aptamer_distance_to_nuc,antibody_distance_to_nuc,nuc_distance_to_nuc])

def append_feature_master(feature_1,feature_2,feature_3):
    with open(parent_directory+'//master_features.csv', 'a', newline='') as file:
        csv.writer(file, delimiter='\t').writerow([feature_1,feature_2,feature_3])

pathList=[]


pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD008_15_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD015_16_crops")
#pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD003_17_crops")
'''pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD016_14_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD016_16_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD022_17_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD024_18_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD025_17_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD026_18_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD027_18_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD030_18_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD036_17_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD041_19_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD050_19_crops")
pathList.append(r"C:/Data/ECAS_cohort/2021_01_19\SD053_16_crops")'''



for i in range(len(pathList)):
    parent_directory = pathList[i]
    
    #Run mast_file (Only needs to be run once per folder)
    master_file_region_props('master_ab_region_props')
    master_file_region_props('master_apt_region_props')
    master_file_region_props('master_nuc_region_props')
    master_file_coincidence()



    sub_pathList=[]
    
    sub_pathList.append(parent_directory + "/1")
    sub_pathList.append(parent_directory + "/2")
    sub_pathList.append(parent_directory + "/3")
    sub_pathList.append(parent_directory + "/4")
    sub_pathList.append(parent_directory + "/5")
    sub_pathList.append(parent_directory + "/6")
    sub_pathList.append(parent_directory + "/7")
    sub_pathList.append(parent_directory + "/8")
    sub_pathList.append(parent_directory + "/9")
    sub_pathList.append(parent_directory + "/10")
    
    total_apt_number = 0
    total_ab_number = 0
    total_nuc_number = 0
    
    
    for i in range(len(sub_pathList)):
        
        directory=sub_pathList[i]
        parent_dir = str(Path(directory).parent)
        
        
    
        
        filename="/Apt.tif"
        apt_image=load_image(directory+filename)
        apt_threshold,apt_binary=threshold_image_apt(apt_image)
        apt_number,apt_labelled=label_image(apt_binary)
        print("%d feautres were detected in the aptamer image."%apt_number)
        apt_measurements=analyse_labelled_image(apt_labelled,apt_image)
        apt_measurements.to_csv(directory + '/' + 'all_aptamer_metrics.csv', sep = '\t')
        apt_measurements.to_csv(parent_dir + '/master_apt_region_props.csv', sep = '\t', mode = 'a', header = False)
        
        
        
        filename="/Ab.tif"
        ab_image=load_image(directory+filename)
        ab_threshold,ab_binary=threshold_image_ab(ab_image)
        ab_number,ab_labelled=label_image(ab_binary)
        print("%d feautres were detected in the antibody image."%ab_number)
        ab_measurements=analyse_labelled_image(ab_labelled,ab_image)
        ab_measurements.to_csv(directory + '/' + 'all_ab_metrics.csv', sep = '\t')
        ab_measurements.to_csv(parent_dir + '/master_ab_region_props.csv', sep = '\t', mode = 'a', header = False)
        
        
        filename="/Nuc.tif"
        nuc_image=load_image(directory+filename)
        nuc_threshold,nuc_binary=threshold_image_otsu(nuc_image)
        nuc_number,nuc_labelled=label_image(nuc_binary)
        print("%d nuclei were detected in the image."%nuc_number)
        nuc_measurements=analyse_labelled_image(nuc_labelled,nuc_image)
        nuc_measurements.to_csv(directory + '/' + 'all_nuc_metrics.csv', sep = '\t')
        nuc_measurements.to_csv(parent_dir + '/master_nuc_region_props.csv', sep = '\t', mode = 'a', header = False)
        append_feature_master(apt_number,ab_number,nuc_number)
        
        
        # Coincidence functions
        
        apt_pixel_coincident_image,apt_pixel_overal_count,apt_pixel_fraction=coincidence_analysis_pixels(apt_binary,ab_binary)
        print("%.2f of aptamer pixels had coincidence with the antibody image."%apt_pixel_fraction)
        
        ab_pixel_coincident_image,ab_pixel_overal_count,ab_pixel_fraction=coincidence_analysis_pixels(ab_binary,apt_binary)
        print("%.2f of antibody pixels had coincidence with the aptamer image."%ab_pixel_fraction)
    
        apt_coinc_list,apt_coinc_pixels,apt_fraction_coinc,apt_coincident_features_image,apt_noncoincident_features_image,apt_fraction_pixels_overlap=feature_coincidence(apt_binary,ab_binary)
        print("%.2f of aptamer features had coincidence with features in antibody image. Average overlap was %2f."%(apt_fraction_coinc,sum(apt_fraction_pixels_overlap)/len(apt_fraction_pixels_overlap)))
        
        aptamer_coincident_tosave=apt_coincident_features_image*apt_image
        im = Image.fromarray(aptamer_coincident_tosave)
        im.save(directory+'/Aptamer_features_coincident.tif')
        
        aptamer_noncoincident_tosave=apt_noncoincident_features_image*apt_image
        im = Image.fromarray(aptamer_noncoincident_tosave)
        im.save(directory+'/Aptamer_features_noncoincident.tif')
        
        
        ab_coinc_list,ab_coinc_pixels,ab_fraction_coinc,ab_coincident_features_image,ab_noncoincident_features_image,ab_fraction_pixels_overlap=feature_coincidence(ab_binary,apt_binary)
        print("%.2f of antibody features had coincidence with features in aptamer image. Average overlap was %2f."%(ab_fraction_coinc,sum(ab_fraction_pixels_overlap)/len(ab_fraction_pixels_overlap)))
        
        antibody_coincident_tosave=ab_coincident_features_image*ab_image
        im = Image.fromarray(antibody_coincident_tosave)
        im.save(directory+'/Antibody_features_coincident.tif')
        
        antibody_noncoincident_tosave=ab_noncoincident_features_image*ab_image
        im = Image.fromarray(antibody_noncoincident_tosave)
        im.save(directory+'/Antibody_features_noncoincident.tif')
    
        
        aptamer_distance_to_nuc=minimum_distance(apt_measurements,nuc_measurements)
        antibody_distance_to_nuc=minimum_distance(ab_measurements,nuc_measurements)  
        nuc_distance_to_nuc=minimum_distance(nuc_measurements,nuc_measurements)
        
        append_coincidence_master()
        
        total_apt_number += apt_number
        total_ab_number += ab_number
        total_nuc_number += nuc_number
    
    
    #Create master histograms of all crop areas
    area_hist_from_csv('Aptamer','master_apt_region_props.csv',total_apt_number)
    area_hist_from_csv('Antibody','master_ab_region_props.csv',total_ab_number)
    area_hist_from_csv('Nucleus','master_nuc_region_props.csv',total_nuc_number)
        
        
