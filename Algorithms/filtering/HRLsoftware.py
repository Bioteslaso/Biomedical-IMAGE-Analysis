###################
#   HRLSoftware   #
###################

# h.calero.2017@alumnos.urjc.es
# r.lopezgon.2017@alumnos.urjc.es             
# p.laso.2017@alumnos.urjc.es #

# 100:150, 150:200 # [same as pyhton in-built fucntion]   (acceptable size >> meaningful results!)

# Modules that we use in this practice

from skimage import io
import numpy as np
#%pylab inline 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import cv2
from tqdm import tqdm
import winsound
from skimage import data, img_as_float
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.metrics import peak_signal_noise_ratio
from skimage.util import random_noise
from scipy.ndimage import gaussian_filter
import scipy.ndimage.filters as filters
from skimage.filters import sobel
import imageio
import acoustics
import random


##########################################################################################
##########################################################################################
###########################  NL Means algorithm  #########################################
##########################################################################################
##########################################################################################

def NL_Means(img, pad_width,h):
    
    image=img #unless debugging

    #padding
    padded_image = np.pad(image, pad_width = pad_width, mode='edge') 

    row, column = padded_image.shape # remember to add thrid argument for 3D structures!
    row_original, column_original = image.shape # remember to add thrid argument for 3D structures!
    denoised_image= np.zeros((row_original, column_original))

    # go over every pixel in our image (2D --> i,j)
    pbar= tqdm(range(0,row_original))
    for i in pbar:
        for j in range(0,column_original):
            pbar.set_description('Re-defining pixel # %s in image. Please, wait' % str(i))
            # create patch for pixel to be filtered (re-defined):
            patch_initial = submatrix(padded_image,i+pad_width,j+pad_width,pad_width) # centered on current image pixel p(i,j)

            # compare current pixel with all other pixels in our image (2D --> k,l)
            distance_matrix_peq=np.zeros((row_original,column_original)) # stores eculidean distance (between pathces) for each central                                                                                  pixel 
            weight_matrix= np.zeros((row_original,column_original))        
            for k in range(0,row_original):
                for l in range (0,column_original):
                
                    patch = submatrix(padded_image,k+pad_width,l+pad_width,pad_width) # centered on image pixel used for comparison p(k,l)
                    eud = euclidean_distance(patch_initial, patch) # compute euclidean distance between pacthes (how gray-similar both                                                                        patches are)
                    distance_matrix_peq[k,l] = eud 
                    weight_matrix[k,l] = np.exp(-(eud)/(h**2)) # compute weight (higher for pixels with a more similar patch)
                    
            # weights correction:
            weight_matrix[weight_matrix==1 ] = 0 # get rid of same patches' distance (eud=0 --> w=1)
            weight_matrix[weight_matrix==0 ] = weight_matrix.max() # re-define as the next maximum weight (other than 1)
            
            # normalization:
            Z = sum(sum(weight_matrix)) # normalize weights (all toghether must sum up to 1)
            weight_matrix_normalized = weight_matrix/Z  
            
            new_pixel_value=0
            for k in range(0,row_original):
                for l in range (0,column_original): 
                    new_pixel_value += image[k,l]*weight_matrix_normalized[k,l] # re-assign pixel value by a weighted average of all other                                                                                 pixels, favouring similarity of patches
            denoised_image[i,j]= new_pixel_value
            
    print('Process ended succesfully.')
    show_results(image,denoised_image)
    
    return denoised_image
	
##########################################################################################
##########################################################################################
###########################  NL Means CPP algorithm  #####################################
##########################################################################################
##########################################################################################

def NL_Means_CPP(image, pad_width,h,D_0,alpha):
    

    #image = mpimg.imread(img)

    #padding
    padded_image = np.pad(image, pad_width = pad_width, mode='edge') 

    row, column = padded_image.shape # remember to add thrid argument for 3D structures!
    row_original, column_original = image.shape # remember to add thrid argument for 3D structures!
    denoised_image= np.zeros((row_original, column_original))


    # go over every pixel in our image (2D --> i,j)
    pbar= tqdm(range(0,row_original))
    for i in pbar:
        for j in range(0,column_original):
            pbar.set_description('Re-defining pixel # %s in image. Please, wait' % str(i))
            # create patch for pixel to be filtered (re-defined):
            patch_initial = submatrix(padded_image,i+pad_width,j+pad_width,pad_width) # centered on current image pixel p(i,j)

            # compare current pixel with all other pixels in our image (2D --> k,l)
            distance_matrix_peq=np.zeros((row_original,column_original)) # stores eculidean distance (between pathces) for each central                                                                            pixel 
            weight_matrix= np.zeros((row_original,column_original))        
            for k in range(0,row_original):
                for l in range (0,column_original):
                                            
                    patch = submatrix(padded_image,k+pad_width,l+pad_width,pad_width) # centered on image pixel used for comparison p(k,l)
                    eud = euclidean_distance(patch_initial, patch) # compute euclidean distance between pacthes (how gray-similar both                                                                            patches are)
                    
                    
                    mu= (1/(1+((image[i,j]-image[k,l])/D_0)**(2*alpha))) # compute weight similarity coeficient from CPP theory
                    distance_matrix_peq[k,l] = eud 
                    weight_matrix[k,l] = np.exp(-(eud)/(h**2))*mu # compute new weight matrix (higher for pixels with a more similar patch)
                    
            # weights correction:
            weight_matrix[weight_matrix==1 ] = 0 # get rid of same patches' distance (eud=0 --> w=1)
            weight_matrix[weight_matrix==0 ] = weight_matrix.max() # re-define as the next maximum weight (other than 1)
            
            # normalization:
            Z = sum(sum(weight_matrix)) # normalize weights (all toghether must sum up to 1)
            weight_matrix_normalized = weight_matrix/Z   
            
            new_pixel_value=0
            for k in range(0,row_original):
                for l in range (0,column_original): 
                    new_pixel_value += image[k,l]*weight_matrix_normalized[k,l] # re-assign pixel value by a weighted average of all other                                                                                 pixels, favouring similarity of patches
            denoised_image[i,j]= new_pixel_value
            
    print('Process ended succesfully.')
    show_results(image,denoised_image)
    
    return denoised_image

##########################################################################################
##########################################################################################
###########################  Anisotropic  ################################################
##########################################################################################
##########################################################################################

def Anisotropic(img, pad_width, treshold, smoothing):
    
    image=img #unless debugging
    
    image_sobel = sobel(image)

    if smoothing == "gaussian":   # In the case of choosing gaussian filter, request pop up to the user in order to know which standard 
                                  # deviation is wanted
            
        std = int((input((print(("enter std"))))))

    
    #padding
    padded_image = np.pad(image, pad_width = pad_width, mode='edge') 

    row, column = padded_image.shape # remember to add thrid argument for 3D structures!
    row_original, column_original = image.shape # remember to add thrid argument for 3D structures!
    
    denoised_image= np.zeros((row_original, column_original))
    
    # go over every pixel in our image (2D --> i,j)
    pbar= tqdm(range(0,row_original))
    for i in pbar:
        for j in range(0,column_original):
            pbar.set_description('Re-defining pixel # %s in image. Please, wait' % str(i))
            # create patch for pixel to be filtered (re-defined):
            patch_no_sobel = submatrix(padded_image,i+pad_width,j+pad_width, pad_width) # centered on current image pixel p(i,j)
            
            # If statements depending on the filter chosen, so that we can apply the desired one
            
            if smoothing == "mean": 
                if image_sobel[i,j] <= treshold: # Comparing from each pixel of the derivatived image its value with our treshold
                                                 # and then decide if filter is  going to be applied or not in the original image
                        
                    denoised_image[i,j] = np.mean(patch_no_sobel)
                else: 
                    denoised_image[i,j] = image[i,j]
            elif smoothing == "median":
                if image_sobel[i,j] <= treshold:
                    denoised_image[i,j] = np.median(patch_no_sobel)
                else: 
                    denoised_image[i,j] = image[i,j]
            elif smoothing == "gaussian":
                if image_sobel[i,j] <= treshold:
                    blur = cv2.GaussianBlur(patch_no_sobel,((3,3)), std)
                    denoised_image[i,j] = blur[1,1]
                else: 
                    denoised_image[i,j] = image[i,j]
                    
                    
    print('Process ended succesfully.')
    show_results(image,denoised_image)
    
    
    return denoised_image



def Anisotropic_2(img, pad_width, treshold, smoothing):  
    
    # The main difference from this anisotropic_2 filter and the previous one is the treshold. In this code, the treshold is going to be         compared with the sum of values of the patch used to iterate over the derived image instead of single pixels.
    
    
    image=img #unless debugging
    
    image_sobel = sobel(image)

    if smoothing == "gaussian":    # In the case of choosing gaussian filter, request pop up to the user in order to know which standard 
                                   # deviation is wanted
            
        std = int((input((print(("enter std"))))))

    
    #padding
    padded_image = np.pad(image, pad_width = pad_width, mode='edge') 
    padded_image_sobel = np.pad(image_sobel, pad_width = pad_width, mode='edge') 

    row, column = padded_image.shape # remember to add thrid argument for 3D structures!
    row_original, column_original = image.shape # remember to add thrid argument for 3D structures!
    
    denoised_image= np.zeros((row_original, column_original))
    
    # go over every pixel in our image (2D --> i,j)
    pbar= tqdm(range(0,row_original))
    for i in pbar:
        for j in range(0,column_original):
            pbar.set_description('Re-defining pixel # %s in image. Please, wait' % str(i))
            # create patch for pixel to be filtered (re-defined):
            patch_no_sobel = submatrix(padded_image,i+pad_width,j+pad_width, pad_width) # centered on current image pixel p(i,j)
            
            patch_sobel = submatrix(padded_image_sobel,i+pad_width,j+pad_width, pad_width)
            
             # If statements depending on the filter chosen, so that we can apply the desired one
            
            if smoothing == "mean":
                if sum(sum(patch_sobel)) <= treshold:   # Comparing from each pixel of the derivatived image its value with our treshold
                                                        # and then decide if filter is  going to be applied or not in the original image
                        
                    denoised_image[i,j] = np.mean(patch_no_sobel)
                else: 
                    denoised_image[i,j] = image[i,j]
            elif smoothing == "median":
                if sum(sum(patch_sobel)) <= treshold:
                    denoised_image[i,j] = np.median(patch_no_sobel)
                else: 
                    denoised_image[i,j] = image[i,j]
            elif smoothing == "gaussian":
                if sum(sum(patch_sobel)) <= treshold:
                    blur = cv2.GaussianBlur(patch_no_sobel,((3,3)), std)
                    denoised_image[i,j] = blur[1,1]
                else: 
                    denoised_image[i,j] = image[i,j]
                    
    print('Process ended succesfully.')
    show_results(image,denoised_image)
    
    
    return denoised_image



def Anisotropic_3(img, iterations, pad_width, treshold, smoothing):  
    
    # The main difference from this anisotropic_2 filter and the previous one is the treshold. In this code, the treshold is going to be         compared with the sum of values of the patch used to iterate over the derived image instead of single pixels.
    
    image=img #unless debugging
    
    if smoothing == "gaussian":    # In the case of choosing gaussian filter, request pop up to the user in order to know which standard 
                                       # deviation is wanted

            std = int((input((print(("enter std"))))))

    
    for i in range(0, iterations):
        
    
        

        image_sobel = sobel(image)

        

        #padding
        padded_image = np.pad(image, pad_width = pad_width, mode='edge') 
        padded_image_sobel = np.pad(image_sobel, pad_width = pad_width, mode='edge') 

        row, column = padded_image.shape # remember to add thrid argument for 3D structures!
        row_original, column_original = image.shape # remember to add thrid argument for 3D structures!

        denoised_image= np.zeros((row_original, column_original))



        # go over every pixel in our image (2D --> i,j)
        pbar= tqdm(range(0,row_original))
        for i in pbar:
            for j in range(0,column_original):
                pbar.set_description('Re-defining pixel # %s in image. Please, wait' % str(i))
                # create patch for pixel to be filtered (re-defined):
                patch_no_sobel = submatrix(padded_image,i+pad_width,j+pad_width, pad_width) # centered on current image pixel p(i,j)

                patch_sobel = submatrix(padded_image_sobel,i+pad_width,j+pad_width, pad_width)

                 # If statements depending on the filter chosen, so that we can apply the desired one

                if smoothing == "mean":
                    if sum(sum(patch_sobel)) <= treshold:   # Comparing from each pixel of the derivatived image its value with our                     treshold
                                                            # and then decide if filter is  going to be applied or not in the original       image

                        denoised_image[i,j] = np.mean(patch_no_sobel)
                    else: 
                        denoised_image[i,j] = image[i,j]
                elif smoothing == "median":
                    if sum(sum(patch_sobel)) <= treshold:
                        denoised_image[i,j] = np.median(patch_no_sobel)
                    else: 
                        denoised_image[i,j] = image[i,j]
                elif smoothing == "gaussian":
                    if sum(sum(patch_sobel)) <= treshold:
                        blur = cv2.GaussianBlur(patch_no_sobel,((3,3)), std)
                        denoised_image[i,j] = blur[1,1]
                    else: 
                        denoised_image[i,j] = image[i,j]
                

        print('Process ended succesfully.')
        show_results(image,denoised_image)
        
        image = denoised_image


    return denoised_image


##########################################################################################
##########################################################################################
###########################  Other functions  ############################################
##########################################################################################
##########################################################################################

def euclidean_distance(patch, patch_initial): 
    a = patch-patch_initial
    b = a**2
    c = sum(sum(b))
    d = np.sqrt(c)
    return d

def submatrix( matrix, centerRow, centerCol, pad_width):
    
    # Function to obtain the patches required. Depending on the parameter pad_width, patches will be of different dimensions.
    
    return matrix[centerRow-pad_width:centerRow+1+pad_width,centerCol-pad_width:centerCol+1+pad_width]

def show_results(image,output_image):
    
 
    fig, axes = plt.subplots(nrows=1,ncols=2, figsize= (12, 5))
    axes[0].imshow(image,cmap='gray')
    axes[1].imshow(output_image,cmap='gray')
    axes[0].axis('off')
    axes[1].axis('off')
    axes[0].set_title('Original')
    axes[1].set_title('Denoised')
    plt.show()
   
    
def notify_me():          # Function that warns when the code has been applied.
    frequency = 1500  
    duration = 2000  
    winsound.Beep(frequency, duration)
    print('>> Notification here')


def show_images(im_out, noise_img):

    import cv2
    import imageio
    import matplotlib.pyplot as plt

    cv2.imwrite(im_out, noise_img)
    im1 = imageio.imread(im_out)
    plt.axis('off')
    plt.imshow(im1, cmap='gray')
    
    return


##########################################################################################
##########################################################################################
###########################  Noises  #####################################################
##########################################################################################
##########################################################################################

def white_noise(image, intensity):

    image = cv2.imread(image)
    image = cv2.resize(image,(240,240)) # resizing for visualization
    
    # Because of the difficulty of implementing a uniform noise in any kind of image, python helps from gaussian function in order to
    # apply the tool acoustics.generator.white
    
    noise = intensity*acoustics.generator.white(image.size).reshape(*image.shape)
    white_image=image+noise
    return white_image


#################################################
    
def gnoise(image, mean, var, intensity):
    
    # Parameters of the gaussian 
    image = cv2.imread(image)
    sigma = var**0.5
    gaussian = intensity* np.random.normal(mean, sigma, image.shape) # Creates gaussian

    noisy_img = image + gaussian # Unify image with the noise

    plt.imshow((noisy_img * 255).astype(np.uint8), cmap=plt.cm.gray)
    plt.title("Gaussian")

    return noisy_img


#################################################

def impulsive_noise(image_in,noise_intensity): # salt-and-pepper

    input_image = cv2.imread(image_in,0)
    untouch_image = input_image
    output = np.zeros(input_image.shape,np.uint8)
    threshold = 1 - noise_intensity 
    
    for i in range(input_image.shape[0]):         # Iteration of whole image to add 0/255 (black/white) pixels values at random positions
        for j in range(input_image.shape[1]):
            rdn = random.random()
            if rdn < noise_intensity:
                output[i][j] = 0
            elif rdn > threshold:
                output[i][j] = 255
            else:
                output[i][j] = input_image[i][j]
                
    untouch_image = input_image
    return output
