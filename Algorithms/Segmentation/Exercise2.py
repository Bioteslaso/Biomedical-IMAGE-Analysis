from skimage import io
import numpy as np
#%pylab inline 
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
import matplotlib.image as mpimg
import cv2
from pydicom import dcmread
import math
from skimage.morphology import reconstruction, square, disk, cube
import imimposemin as im
from skimage.segmentation import watershed

from scipy.ndimage import gaussian_filter
import scipy.ndimage.filters as filters
from skimage.filters import sobel, prewitt
import imageio

###########################################################################################################################################

def imimposemin(I, BW, conn=None, max_value=255):
    if not I.ndim in (2, 3):
        raise Exception("'I' must be a 2-D or 3D array.")

    if BW.shape != I.shape:
        raise Exception("'I' and 'BW' must have the same shape.")

    if BW.dtype is not bool:
        BW = BW != 0

    # set default connectivity depending on whether the image is 2-D or 3-D
    if conn == None:
        if I.ndim == 3:
            conn = 26
        else:
            conn = 8
    else:
        if conn in (4, 8) and I.ndim == 3:
            raise Exception("'conn' is invalid for a 3-D image.")
        elif conn in (6, 18, 26) and I.ndim == 2:
            raise Exception("'conn' is invalid for a 2-D image.")

    # create structuring element depending on connectivity
    if conn == 4:
        selem = disk(1)
    elif conn == 8:
        selem = square(3)
    elif conn == 6:
        selem = ball(1)
    elif conn == 18:
        selem = ball(1)
        selem[:, 1, :] = 1
        selem[:, :, 1] = 1
        selem[1] = 1
    elif conn == 26:
        selem = cube(3)

    fm = I.astype(float)

    try:
        fm[BW]                 = -math.inf
        fm[np.logical_not(BW)] = math.inf
    except:
        fm[BW]                 = -float("inf")
        fm[np.logical_not(BW)] = float("inf")

    if I.dtype == float:
        I_range = np.amax(I) - np.amin(I)

        if I_range == 0:
            h = 0.1
        else:
            h = I_range*0.001
    else:
        h = 1

    fp1 = I + h

    g = np.minimum(fp1, fm)

    # perform reconstruction and get the image complement of the result
    if I.dtype == float:
        J = reconstruction(1 - fm, 1 - g, selem=selem)
        J = 1 - J
    else:
        J = reconstruction(255 - fm, 255 - g, method='dilation', selem=selem)
        J = 255 - J

    try:
        J[BW] = -math.inf
    except:
        J[BW] = -float("inf")

    return J

###########################################################################################################################################

def RegionGrowingP2(tipo, image_in, threshold_sup, threshold_inf, seed):
    '''
    The function takes 5 different arguments. 
    Tipo is either "dicom" or "png", depending on what type of image you want to segment.
    Image_in is just the name of the image.
    Threshold_sup, and threshold_inf are 2 values (int or float) to select what pixels we want to examine.
    Seed is the list with the coordinates of the initial pixel of the area of interest. It has to be of the form [(x,y)]
    '''
    if tipo == "dicom":
        ds = dcmread(image_in)
        img_in = ds.pixel_array
        img = img_in
        img = img/np.max(img)
        img_in = img_in/np.max(img_in)
        
    if tipo == "png":
        img=mpimg.imread(image_in)
        img=img[:,:,0] 
        img_in = img
        img = img/np.max(img)
        img_in = img_in/np.max(img_in)
        
    pad_width=1
    padded_image = np.pad(img, pad_width = pad_width, mode='constant',constant_values=100000) # We add a very high constant value to the   padding of the image so that the algorithm will not use the borders of the image for the connectivity
    
    row, column = padded_image.shape # remember to add thrid argument for 3D structures!
    row_original, column_original = img.shape # remember to add thrid argument for 3D structures!
    img=padded_image
    initial_matrix = img
    
    threshold_inf = np.abs(threshold_inf) # We get the absolute value of the thresholds chosen by the user to avoid mistakes.
    threshold_sup = np.abs(threshold_sup)

    #seed = coord
    a, b = seed[0][0]+pad_width, seed[0][1]+pad_width # We must consider the padding to correctly select the initial seed

    zeros_matrix = np.zeros(img_in.shape,np.int8) # Zeros_matrix to generate the mask that will multiply the original image afterwards

    zeros_matrix[a-1,b-1] = 1 

    iteration_matrix = submatrix(img,a,b,1) # We create the initial patch to for which the first connectivity of 8 is going to be done 
    row_patch, column_patch = iteration_matrix.shape # 3x3

    seeds_list = [] # create new list with all the values that fullfil the condition of the threshold, to afterwards perform the connectivity of 8 to them
    black_list= [] # new list to append all the "seeds" that have already been studied

    seeds_list.append([a,b]) # add initial pixel (the seed)

    while bool(seeds_list): # iterate until the list of seeds is empty:
        black_list.append([seeds_list[0][0],seeds_list[0][1]]) # append to the black_list those pixels already studied
        for t in range(0,row_patch): #iterate over all the pixels of the patch (connectivity of 8)
            for w in range(0,column_patch):
                iteration_matrix = submatrix(initial_matrix,seeds_list[0][0],seeds_list[0][1],1) # We create a patch centered in the new pixel to be studied
                central_pixel=initial_matrix[seeds_list[0][0],seeds_list[0][1]]

                if (initial_matrix[a,b] - threshold_inf) <=  iteration_matrix[t,w] <= (initial_matrix[a,b] + threshold_sup) :  # check if the current pixel in 8 connectivity patch is in the desired range 
                                    #seed = [k,l]

                    if [t+seeds_list[0][0]-1,w+seeds_list[0][1]-1] not in black_list:
                        if [t+seeds_list[0][0]-1,w+seeds_list[0][1]-1] not in seeds_list:
                            seeds_list.append([t+seeds_list[0][0]-1,(w+seeds_list[0][1]-1)]) # add  (current iterating) homogeneous pixel to list
                            zeros_matrix[t+seeds_list[0][0]-1-pad_width,(w+seeds_list[0][1]-1-pad_width)] = 1 #iteration_matrix[t,w] # assign value of (current iterating) homogeneous pixel to zeros_matrix (new matrix for initial_matrix)

                else:
                    pass #zeros_matrix[t+a-1,w+b-1] = 0  
        
        seeds_list.remove(seeds_list[0]) # done with this pixel being used as central in 8-connectivity patch. We keep iterating th rest
        
    filtered_image = np.multiply(zeros_matrix, img_in) #Multiply the obtained binary mask (zeros_matrix) with the original image to obtain the desired region.

    return filtered_image, zeros_matrix
    
###########################################################################################################################################   
def WatershedP2(image_in, seed, tipo):
    '''
    The function takes 3 different arguments. 
    Tipo is either "dicom" or "png", depending on what type of image you want to segment.
    Image_in is just the name of the image.
    Seed is the list with the coordinates of the initial pixel of the area of interest. It has to be of the form [(y, x)]
    '''
    
    if tipo == "dicom":
        ds = dcmread(image_in)
        img_in = ds.pixel_array
        img = img_in
        img = img/np.max(img)
        img_in = img_in/np.max(img_in)
        
    if tipo == "png":
        img=mpimg.imread(image_in)
        img=img[:,:,0] 
        img_in = img
        img = img/np.max(img)
        img_in = img_in/np.max(img_in)
        
    sobel_image = sobel(img) # We perform the derivative filter to our original image to obtain the edges.
    
    zeros_matrix = np.zeros(img.shape, np.int8) 
    zeros_matrix = zeros_matrix*255
    for i in range(0, len(seed)):
        zeros_matrix[seed[i][0],seed[i][1]] = 255 # This is the mask with the positions we want to select as the local minima, to afterwards perform the watershed algorithm 
    
    J = im.imimposemin(I=sobel_image, BW=zeros_matrix) 

    watershed_image_mask =watershed(J)
    watershed_sobel = watershed(sobel_image)  ########################### PREGUNTAR
    
    fig = plt.figure(figsize=(15, 15))
    plt.subplot(141)
    plt.imshow(watershed_image_mask)
    plt.title("Watershed Segmentation")
    plt.subplot(142)
    plt.imshow(watershed_image_mask)
    plt.imshow(img_in,alpha = 0.5, cmap = "gray")
    plt.title("Original + Segmentation")
    plt.subplot(143)
    plt.imshow(watershed_sobel, cmap = "gray")
    plt.title("Segmentation w/o imimposemin")
    
    return watershed_image_mask, watershed_sobel, img_in

###########################################################################################################################################
def submatrix( matrix, centerRow, centerCol, pad_width):
    
    # Function to obtain the patches required. Depending on the parameter pad_width, patches will
    # be of different dimensions.
    
    return matrix[centerRow-pad_width:centerRow+1+pad_width,centerCol-pad_width:centerCol+1+pad_width]


