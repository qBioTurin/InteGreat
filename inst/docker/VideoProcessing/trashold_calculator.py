import cv2
import skimage as ski
import numpy as np
import matplotlib.pyplot as plt
import sys
# Definisci una funzione per scrivere sulla console di RStudio



def calculate_trashold(im,trash_val):
    

    blurred = cv2.GaussianBlur(im, (3, 3), 0) 
    kernel_sharpening = np.array([[-1, -1, -1],
                    [-1, 9, -1],
                    [-1, -1, -1]]) 
    sharpened = cv2.filter2D(blurred, -1, kernel_sharpening)
    _, thresholded = cv2.threshold(sharpened, trash_val, 255, cv2.THRESH_BINARY)
    return thresholded

if __name__ == "__main__":
    trash_val = sys.argv[1] 
    trash_val = int(trash_val)
    image_phat = sys.argv[2]
    image =cv2.imread(image_phat)
   
    tresh_list = calculate_trashold(im=image,trash_val=trash_val)
    tresh_phat = "../home/tresh.png"
    cv2.imwrite(tresh_phat,tresh_list)
    

    
