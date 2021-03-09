#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:52:02 2021
@author: Owen Kantelberg
"""


import numpy as np
from PIL import Image
import os
import slideio
import napari

def make_crops(directory,file_to_load,folder_name):
    # Make a new folder to store the output images:
    directory=path+'/'+folder_name+"_crops/"
    if not os.path.exists(directory):
        os.mkdir(directory)
            
            
    # Open the slide
    slide = slideio.open_slide(file_to_load, "CZI")
    scene = slide.get_scene(0)
    
    
    # The next section is to make a scaled version and then open it using napari. 
    # Need to get the size of the image:
    im_width,im_height=scene.size
    
    # Scale image by factor:
    scale_factor=50
    
    new_width=int(im_width/scale_factor)
    new_height=int(im_height/scale_factor)
    
    # Size of images to generate
    image_size=2000
    scaled_size=int(image_size/(2*scale_factor))
    
    nuc_zoom = scene.read_block(size=(new_width,0), channel_indices=[0])
    ab_zoom = scene.read_block(size=(new_width,0), channel_indices=[1])
    apt_zoom = scene.read_block(size=(new_width,0), channel_indices=[2])
    
    
    # Save the scaled images:
    
    im = Image.fromarray(nuc_zoom)
    im.save(directory+'Nuc_scaled.tif')
    
    im = Image.fromarray(apt_zoom)
    im.save(directory+'Apt_scaled.tif')
    
    im = Image.fromarray(ab_zoom)
    im.save(directory+'Ab_scaled.tif')
    
    
    
    # Now open napari to click the points:
        
    with napari.gui_qt():
            viewer = napari.Viewer()
            viewer.add_image(nuc_zoom,colormap="blue",opacity=0.5)
            viewer.add_image(ab_zoom,colormap="green",opacity=0.5)
            viewer.add_image(apt_zoom,colormap="red",opacity=0.5)
    
    
    # This array contains all of the points (x and y coords)
    positions = viewer.layers["Points"].data
    
    j=0
    
    # Show the ROIs
    
    rois=np.zeros((new_height,new_width))
    
    for i in range(len(positions)):
    
        xpos=int(positions[i,0]*scale_factor-(image_size/2))
        ypos=int(positions[i,1]*scale_factor-(image_size/2))
    
        xposroi=int(xpos/scale_factor)
        yposroi=int(ypos/scale_factor)
        
        rois[(xposroi):(xposroi+2*scaled_size),(yposroi):(yposroi+2*scaled_size)]=i+1
        
        
        roi_directory=path+"/"+folder_name+"_crops/"+str(i+1)+"/"
        if not os.path.exists(roi_directory):
            os.mkdir(roi_directory)
            
        
        image_nuc = scene.read_block((ypos,xpos,image_size,image_size),channel_indices=[0])
        image_ab = scene.read_block((ypos,xpos,image_size,image_size),channel_indices=[1])
        image_apt = scene.read_block((ypos,xpos,image_size,image_size),channel_indices=[2])
        
        im = Image.fromarray(image_nuc)
        im.save(roi_directory+'Nuc.tif')
          
        im = Image.fromarray(image_ab)
        im.save(roi_directory+'Ab.tif')
        
        im = Image.fromarray(image_apt)
        im.save(roi_directory+'Apt.tif')
    
        
    im = Image.fromarray(rois)
    im.save(directory+'Rois.tif')


# Path and file to load
#path=r"C:/Data/ECAS_cohort/2021_01_26/"
#filename="2021_01_26_SD058_13.czi"
#file_to_load=path+filename

pathList=[]

pathList.append(r"C:/Data/ECAS_cohort/2021_01_29")


for i in range(len(pathList)):
    path=pathList[i]
    for root, dirs, files in os.walk(path):
        for name in files:
            if '.czi' in name:
                #for files in path:
                working_file = path + '/' + name
                folder_name = name[11:19]
                make_crops(path, working_file,folder_name)
    

