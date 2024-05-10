import sys
import skimage as ski
import numpy as np
import matplotlib.pyplot as plt
import json
from readlif.reader import LifFile


# fileImportato = LifFile('./movies.lif')

def first_frame(path_fileImportato):
    fileImportato = LifFile(path_fileImportato)
    images = []
    for i in range(0, len(fileImportato.image_list)):
        video = fileImportato.get_image(i)
        frame = video.get_frame(z=0, t=0, c=0)
        images.append((np.array(frame)).tolist())
    return images


# BORDER DETECTION

def border_detection(path_fileImportato, videos):
    fileImportato = LifFile(path_fileImportato)
    # Morphological Image Processing
    images_processed = []
    borders = []
    for i in range(0, len(fileImportato.image_list)):
        if videos[i] == 1:
            video = fileImportato.get_image(i)
            frame_mean = video.get_frame(z=0, t=0, c=0)
            frame_mean = np.array(frame_mean)
            for j in range(1, video.dims_n[4]):
                frame_j = video.get_frame(z=0, t=j, c=0)
                frame_j = np.array(frame_j)
                frame_mean = (frame_j + frame_mean) / 2
            frame_mean = ski.filters.unsharp_mask(frame_mean, radius=30, amount=1)
            image_processed = ski.morphology.remove_small_holes(frame_mean.astype(bool),
                                                                area_threshold=3000)  # remove small holes and first detection of nucleus
            image_processed = ski.morphology.area_opening(image_processed,
                                                          area_threshold=15000)  # opening alto per eliminare aree di sfondo che si riempivano
            images_processed.append(image_processed)
    # Segmentation
    for i in range(0, len(images_processed)):
        edge_sobel = ski.filters.sobel(images_processed[i])
        borders.append(edge_sobel.tolist())
    return borders


# NUCLEUS DETECTION


def nucleus_detection(path_fileImportato, videos):
    fileImportato = LifFile(path_fileImportato)
    nucleus = []
    for i in range(0, len(fileImportato.image_list)):
        if videos[i] == 1:
            image = fileImportato.get_image(i).get_frame(z=0, t=0, c=0)
            image = np.array(image)
            image = ski.filters.unsharp_mask(image, radius=50, amount=1)
            image=ski.morphology.area_opening(image, area_threshold=1000) #opening
            image = ski.morphology.area_closing(image, area_threshold=1500)  # closing
            
            canny = ski.feature.canny(image, sigma=4)
            canny = ski.morphology.remove_small_holes(canny, area_threshold=2000)
            canny = ski.filters.gaussian(canny, sigma=0.3)  # rendo blurrato per chiudere eventuali buchi per il flood
            mask = ski.morphology.flood(canny, (0, 0))
            mask = np.array(mask)
            mask = np.bitwise_not(mask)  # preparo per opening
            mask = ski.morphology.area_opening(mask, area_threshold=2000)  # tolgo bordini
            border = ski.filters.sobel(mask)  # bordo
            nucleus.append(border.tolist())
    return nucleus


def read_checkbox_states(file_path):
    with open(file_path, "r") as file:
        data = file.read().strip()  # Read the contents of the file
        # Split the string into individual boolean values
        checkbox_states_str = data.split()
        # Convert the boolean strings to actual boolean values
        checkbox_states = [True if state == "TRUE" else False for state in checkbox_states_str]
    return checkbox_states


if __name__ == "__main__":
    if sys.argv[1] == "first_frame":
        video_path = sys.argv[2]  # Path to the uploaded video
        ris = first_frame(video_path)
        with open("../home/resultFirstFrames.json", "w") as f:
            json.dump(ris, f)
    elif sys.argv[1] == "border_detection":
        video_path = sys.argv[2]  # Path to the uploaded video
        checkbox_states_path = sys.argv[3]  # Path to the checkbox states file
        checkbox_states = read_checkbox_states(checkbox_states_path)
        borders = border_detection(video_path, checkbox_states)
        nucleus = nucleus_detection(video_path, checkbox_states)
        with open("../home/resultNucleus.json", "w") as f:
            json.dump(nucleus, f)
        with open("../home/resultBorder.json", "w") as f:
            json.dump(borders, f)
