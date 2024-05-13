import numpy as np
import matplotlib.pyplot as plt
import skimage as ski
from skimage.segmentation import flood_fill
import cv2
import imageio
from math import sqrt
from skimage.feature import blob_log
from skimage import img_as_ubyte
import sys
from readlif.reader import LifFile
import json
# Definisci le variabili globali per le coordinate del frame precedente
previous_x = 0
previous_y = 0
Vescicole_erratiche = 0
Vescicole_Spola_locale = 0
Vescicole_Spola_Nucleo_Membrana = 0
Vescicole_totali = 0
distance_trashold = 50
raggio_di_trashold = 2
previous_x_arrow = [None] 
previous_y_arrow = [None] 
PrePrivius_x = [None]
PrePrivius_y = [None]
color_arrow = [None]
distances_vescicles = [0]
processed_ids = []
cheked_id_pattern = []
check_terget = []
VescicleType = [None]
Vescicle_target_color = [None]
last_arrow = []
arrow_points = []
def read_checkbox_states(file_path):
    with open(file_path, "r") as file:
        data = file.read().strip()  # Read the contents of the file
        # Split the string into individual boolean values
        checkbox_states_str = data.split()
        # Convert the boolean strings to actual boolean values
        checkbox_states = [True if state == "TRUE" else False for state in checkbox_states_str]
    return checkbox_states
def nucleus_target_detection(Nucleus):
    
    
    
    kernel = np.ones((5, 5), np.uint8)
    Nucleus_target = cv2.dilate(Nucleus, kernel, iterations=2)
    Nucleus_target = (Nucleus_target>0)
    Nucleus_target = Nucleus_target * 6#identificativo del bordo
    
            
    return Nucleus_target

  



    
    
    
    
    

def check_distance_to_higth(x_precedente, y_precedente, x_corrente, y_corrente, soglia_base):
    distanza = euclidean_distance(x_precedente, y_precedente, x_corrente, y_corrente)
    soglia_attuale = soglia_base * distanza
    return distanza > soglia_attuale

# Definisci la funzione per calcolare la distanza euclidea tra due punti
def euclidean_distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

# Definisci la funzione per calcolare la direzione tra due punti
def calculate_direction(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    direction = np.arctan2(dy, dx)  # Calcola l'angolo della direzione in radianti
    
    return direction

# Definisci la funzione per verificare se delle coordinate cadono in una direzione specifica
def check_direction(x, y, x1, y1, x2, y2, tolerance=0.1):
    direction = calculate_direction(x1, y1, x2, y2)
    angle_tolerance = tolerance * np.pi  # Converte la tolleranza da gradi a radianti
    target_direction = np.arctan2(y - y1, x - x1)  # Calcola l'angolo rispetto al punto iniziale (x1, y1)
    
    # Calcola la differenza tra l'angolo della direzione e l'angolo del punto rispetto al punto iniziale
    angle_difference = np.abs(direction - target_direction)
    
    # Assicurati che l'angolo della differenza sia compreso tra -pi e pi
    if angle_difference > np.pi:
        angle_difference = 2 * np.pi - angle_difference
    
    # Verifica se l'angolo della differenza è minore della tolleranza
    return angle_difference < angle_tolerance

# Definisci la funzione per predire la posizione futura
def predict_origin_position(image,current_x, current_y, possible_positions,frame,show_ids,spola_target):
    global previous_x, previous_y,processed_ids,arrow_points,previous_x_arrow,previous_y_arrow,color_arrow, distances_vescicles,check_terget,Vescicle_target_color,VescicleType# Utilizza le variabili globali per le coordinate del frame precedente
    global PrePrivius_x,PrePrivius_y,Vescicole_erratiche,Vescicole_Spola_locale,Vescicole_Spola_Nucleo_Membrana,cheked_id_pattern,raggio_di_trashold,Vescicole_totali
    
    # Calcola la distanza euclidea tra la posizione corrente e tutte le possibili posizioni future
    distances = [euclidean_distance(current_x, current_y, blob[0], blob[1]) for blob in possible_positions]
    
    # Calcola le probabilità delle posizioni basate sul reciproco della distanza euclidea
    # Più vicino è il punto, più alta è la probabilità
    total_distance = sum(distances)
    probabilities = [1 - dist / total_distance for dist in distances]
    
    # Verifica se le coordinate cadono nella direzione desiderata
    direction_check_results = [check_direction(blob[0], blob[1],previous_x, previous_y, current_x, current_y) for blob in possible_positions]
    
    # Considera solo le posizioni che soddisfano il check di direzione
    possible_positions_filtered = [pos for pos, check_result in zip(possible_positions, direction_check_results) if check_result]
    probabilities_filtered = [prob for prob, check_result in zip(probabilities, direction_check_results) if check_result]
    
     # Se ci sono posizioni valide, seleziona la posizione con la probabilità più alta come posizione futura
    if possible_positions_filtered:
        max_probability_index = np.argmax(probabilities_filtered)
        origin_x, origin_y, origin_radius, origin_id = possible_positions_filtered[max_probability_index]
    else:
        # Se non ci sono posizioni valide, restituisci None per indicare che non è possibile predire una posizione futura
        origin_x, origin_y, origin_radius, origin_id = None, None, None, None
    
    
    # Disegnamento Frecce
    
    if origin_id not in processed_ids and origin_id:

        
        if  frame == 1:
            
            previous_x_arrow[origin_id] = current_x
            previous_y_arrow[origin_id] = current_y
            distance = euclidean_distance(origin_x,origin_y,current_x,current_y)
            if distance<distance_trashold and color_arrow[origin_id] is not None :
                distances_vescicles[origin_id] = distance
                arrow_points.append(((origin_x,origin_y),(current_x,current_y),origin_id,color_arrow[origin_id],VescicleType[origin_id]))

           
                
        else: #dal terzo in poi 
            
           
                
            if previous_x_arrow[origin_id] is not None and previous_y_arrow[origin_id] is not None and color_arrow[origin_id] is not None:
                distance = euclidean_distance(previous_x_arrow[origin_id],previous_y_arrow[origin_id],current_x,current_y)
                if distance<distance_trashold:
                    distances_vescicles[origin_id] += distance
                    arrow_points.append(((previous_x_arrow[origin_id],previous_y_arrow[origin_id]),(current_x,current_y),origin_id,color_arrow[origin_id],VescicleType[origin_id]))
                    # Definisci il font, la scala, il colore e lo spessore del testo
                    PrePrivius_x[origin_id] = previous_x_arrow[origin_id]
                    PrePrivius_y[origin_id] = previous_y_arrow[origin_id]
                    previous_x_arrow[origin_id] = current_x
                    previous_y_arrow[origin_id] = current_y
       
            # dal terzo frame in poi perchè altrimenti non ci sono abbastanza informazioni
            if PrePrivius_x[origin_id] is not None and PrePrivius_y[origin_id] is not None:
                
                if origin_id not in check_terget:
                    
                    id_pos = spola_target[current_y,current_x]
                   
                    if id_pos == 11:## zona di overlapp conta come se facesse spola 
                        Vescicole_Spola_Nucleo_Membrana +=1
                        VescicleType[origin_id] = "Target"#Tipo Target 
                        Vescicole_totali += 1
                        check_terget.append(origin_id)
                    
                    if id_pos == 6 and Vescicle_target_color[origin_id] is None:#non ha toccato la membrana
                        Vescicle_target_color[origin_id] = 1#RED
                        
                        
                    if id_pos == 6 and Vescicle_target_color[origin_id] == 2 :#ha toccato nucleo e membrana
                        Vescicole_Spola_Nucleo_Membrana += 1 
                        VescicleType[origin_id] = "Target"#Tipo Target 
                        Vescicole_totali += 1
                        check_terget.append(origin_id)
                        
                    if id_pos == 5 and Vescicle_target_color[origin_id] is None:#non ha toccato il nucleo 
                        Vescicle_target_color[origin_id] = 2#BLUE
                        
                    if id_pos == 5 and Vescicle_target_color[origin_id] == 1:#ha toccato membrana e nucleo
                        Vescicole_Spola_Nucleo_Membrana +=1
                        VescicleType[origin_id] = "Target"#Tipo Target
                        Vescicole_totali += 1
                        check_terget.append(origin_id)
                        
                if origin_id not in cheked_id_pattern and origin_id not in check_terget and euclidean_distance(current_x,current_y,PrePrivius_x[origin_id],PrePrivius_y[origin_id]) < raggio_di_trashold:
                    Vescicole_totali += 1
                    Vescicole_Spola_locale += 1
                    VescicleType[origin_id] = "Local"
                    cheked_id_pattern.append(origin_id)
                   
                if origin_id not in cheked_id_pattern:
                    Vescicole_totali += 1    
                    Vescicole_erratiche += 1
                    VescicleType[origin_id] = "Erratic"
                    cheked_id_pattern.append(origin_id)
                
               
                
               

        if show_ids:    
                cv2.putText(image, str(origin_id), (current_x,current_y),cv2.FONT_HERSHEY_SIMPLEX,0.3,(0,0,255))
        
    
   
        processed_ids.append(origin_id)
    
    previous_x = current_x
    previous_y = current_y

def boreder_cleaner_and_target(border):
    border = border.astype(np.uint8)*255
    
    ##parte di eliminazione membrane di cellule inutili###############
    # Applica un'operazione di erosione seguita dalla dilatazione per ridurre i dettagli
    kernel = np.ones((3, 3), np.uint8)
    image_processed = cv2.morphologyEx(border, cv2.MORPH_CLOSE, kernel)

    # Trova i contorni nell'immagine
    contours, _ = cv2.findContours(image_processed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    
    # Trova il contorno con l'area più grande# si presuppone sia la cellula in esame
    largest_contour = max(contours, key=cv2.contourArea)

    # Crea un'immagine vuota con lo stesso shape dell'originale
    erasing_mask = np.zeros_like(border)

    # Trova i contorni che toccano il bordo dell'immagine in almeno due punti, escludendo il contorno più grande
    for contour in contours:
        if contour is not largest_contour:
            for point in contour:
                x, y = point[0]
                if x == 0 or x == border.shape[1]-1 or y == 0 or y == border.shape[0]-1:
                    cv2.drawContours(erasing_mask, [contour], -1, 255, thickness=cv2.FILLED)
                    break
    kernel = np.ones((8, 8), np.uint8)
    erasing_mask = 255-cv2.morphologyEx(erasing_mask, cv2.MORPH_DILATE, kernel)               
    border_clean = (border * erasing_mask)*255
    
    #####################################
    ####Parte di creazione del bordo target########
    border_target = ski.filters.gaussian(border_clean, sigma=0.5)
    border_target = (border_target > 0)
    
    border_target = border_target * 5#identificativo del bordo 
    border_clean = border_clean *255
    
    return border_clean,border_target
        
def Remove_usless_vescicles(im,trash_val,edge_canny):
    ###prendere in input il bordo pulito 
    
    blurred = cv2.GaussianBlur(im, (3, 3), 0) 
    kernel_sharpening = np.array([[-1, -1, -1],
                              [-1, 9, -1],
                              [-1, -1, -1]]) 
    sharpened = cv2.filter2D(blurred, -1, kernel_sharpening)
    _, thresholded = cv2.threshold(sharpened, trash_val, 255, cv2.THRESH_BINARY)
   
    ###da mantenere pulizia delle vesciole inutili 
    kernel = np.ones((10, 10), np.uint8)
    edge_canny = cv2.morphologyEx(edge_canny, cv2.MORPH_DILATE, kernel)    
    Vescicle_cleaner = 255-flood_fill(edge_canny, (0,0), 255)
    denoise_img = cv2.bitwise_and(Vescicle_cleaner, thresholded)
    denoise_img = img_as_ubyte(denoise_img)
    
    return denoise_img
def read_input_values_tracking(path_file):
   # Apre il file txt e legge i dati
    with open(path_file, "r") as file:
        data = file.readlines()

    # Inizializza le liste per i valori numerici e booleani
    valori_numerici = []
    valori_booleani = []
    
    i = 0
    j = 0
    while i < len(data):
        img_index = int(data[i].strip().split("_")[1])
        valori_numerici.insert(img_index,int(data[i+1].strip()))
        valori_booleani.insert(img_index,[x.strip()=="TRUE" for x in data[i+2:i+8]])
        
        i+=8
    return valori_numerici,valori_booleani


def Run_detector_tracker(show_Ids,tresh_val,Fast_mode,Visibility,file_tif,output_file,Show_erratic,Show_traget,Show_local,border,Nucleus):
    global arrow_points,previous_x_arrow,previous_y_arrow, color_arrow, distances_vescicles,PrePrivius_x,PrePrivius_y,Vescicole_totali,Vescicle_target_color
    global VescicleType
    ##border

    
   
    # Definisci un dizionario per memorizzare i tag associati a ciascun blob
    blob_IDs = []
    Video_time = 0
    frame_rate = 2 
    frame_index = 0
    # Crea un writer per il video TIFF
    Nucleus = np.array(Nucleus)
    Nucleus[Nucleus>0] = 255
    Nucleus = Nucleus.astype(np.uint8)
    border = np.array(border)
    border[border>0] = 255
    
    border_clean, border_target = boreder_cleaner_and_target(border)
    Nucleus_target = nucleus_target_detection(Nucleus)
    
    fourcc = cv2.VideoWriter_fourcc(*'VP90')
    output_video = cv2.VideoWriter(output_file, fourcc, frame_rate, (512, 512))
    num_frames =file_tif.dims.t
    
    for i in range(num_frames):

        image = file_tif.get_frame(t=i)
        image = np.array(image)
       
    
        denoise_img  = Remove_usless_vescicles(image,tresh_val,border_clean)
       
        if Fast_mode:
            blobs= blob_log(denoise_img,min_sigma=3, max_sigma=7, threshold=0.05)
        else:
            blobs= blob_log(denoise_img,min_sigma=2, max_sigma=7, threshold=0.05)

        blobs[:, 2] = blobs[:, 2] * sqrt(2)
        
        
            
        if Visibility:
            image = denoise_img + Nucleus + border_clean
        else:
            
            image = Nucleus + border_clean
       
        

        # Combina le immagini risultanti
        Spola_target = Nucleus_target + border_target
      
        image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)
        
        ##inizializzazione delle grafiche "freccie"
        
        #liste ids
        previous_x_arrow  += [None] * (len(blobs) - len(previous_x_arrow))
        previous_y_arrow  += [None] * (len(blobs) - len(previous_y_arrow))
        PrePrivius_x += [None] * (len(blobs)-len(PrePrivius_x))
        PrePrivius_y += [None] * (len(blobs)-len(PrePrivius_y))
        Vescicle_target_color += [None] * (len(blobs)-len(Vescicle_target_color))
        color_arrow += [None] * (len(blobs) - len(color_arrow))
        distances_vescicles += [0] * (len(blobs) - len(distances_vescicles)) 
        VescicleType += [None] * (len(blobs)-len(VescicleType)) 

        blob_id_counter = 0
        if frame_index % frame_rate == 0:
            Video_time +=1
        for blob_id_counter, blob in enumerate(blobs):
            y, x, r = blob
            
            
            # Verifica se il blob è già stato assegnato a un tag
            if frame_index == 0:
                
                color1 = (list(np.random.choice(range(256), size=3)))  
                color =int(color1[0]), int(color1[1]), int(color1[2])
                
                    
                         
                blob_IDs.append((int(x), int(y), int(r), blob_id_counter))
                previous_x_arrow[blob_id_counter] = int(x)
                previous_y_arrow[blob_id_counter] = int(y)
                color_arrow[blob_id_counter] = color
                if blob_IDs:#escludo l'inizializzazione 
                    
                    if show_Ids:
                        cv2.putText(image, str(blob_IDs[blob_id_counter][3]), (int(x),int(y)),cv2.FONT_HERSHEY_SIMPLEX,0.3,(0,0,255))
                    
                    

            else:
            
                if blob_IDs:
                        
                    predict_origin_position(image,int(x),int(y),blob_IDs,frame_index,show_Ids,Spola_target)
                     
                blob_id_counter +=1
            
        
        if arrow_points :
            arrow_points= sorted(arrow_points, key=lambda x: x[2])
            
            for i in range (len(arrow_points)):
                if Show_traget and Show_local and Show_erratic:
                    cv2.arrowedLine(image,arrow_points[i][0],arrow_points[i][1],arrow_points[i][3],1)
                if Show_traget and arrow_points[i][4] == "Target":
                    cv2.arrowedLine(image,arrow_points[i][0],arrow_points[i][1],arrow_points[i][3],1)
                
                if Show_local and arrow_points[i][4] == "Local": 
                    cv2.arrowedLine(image,arrow_points[i][0],arrow_points[i][1],arrow_points[i][3],1)
                
                if Show_erratic and arrow_points[i][4] == "Erraic": 
                    cv2.arrowedLine(image,arrow_points[i][0],arrow_points[i][1],arrow_points[i][3],1)

        
        
        frame_index = frame_index + 1      
        processed_ids.clear()
        
        
        output_video.write(image)
    
    output_video.release()
                        ##Metriche
                        #velocità media 
    distances_vescicles = list(filter(lambda x: x != 0, distances_vescicles))
    Mean_velocity = np.average(distances_vescicles)/Video_time#nanometri/secondo
              
    return Mean_velocity,Vescicole_Spola_locale,Vescicole_Spola_Nucleo_Membrana,Vescicole_erratiche,Vescicole_totali
                    
def save_metrics(file_phat,Mean_vel,Erratic,Local,Target,Total):
      # Apri il file in modalità scrittura
    with open(file_phat, "w") as file:
        # Scrivi le metriche nel file
        file.write("Mean velocity: {} nanometers/seconds\n".format(Mean_vel))
        file.write("N° Erratic Vescicles: {}\n".format(Erratic))
        file.write("N° 'back and forth' Vescicles: {}\n".format(Local))
        file.write("N° 'back and forth between Nucleus -> Membrane' Vescicles: {}\n".format(Target))
        file.write("N° Total Vescicles: {}\n".format(Total))



def CalcCount(chek_v,index):
    i =0
    count = 0
    while index >= 0:
      if chek_v[i] != True:
        count +=1
      else:
        count+=1
        index -=1
      i+=1
    return count-1


if __name__ == "__main__":
    phat_lif = sys.argv[1]
    File_lif = LifFile(phat_lif)
    img_index = int(sys.argv[2]) - 1 
    checkbox_slider_phat = sys.argv[3]
    thres_val_list,input_boolean_values = read_input_values_tracking(checkbox_slider_phat)
    phat_border = sys.argv[4]
    with open(phat_border, 'r') as file:
         border = json.load(file)
    phat_nucleus = sys.argv[5]
    with open(phat_nucleus, 'r') as file:
         Nucleus = json.load(file)     
    checkbox_states_path = sys.argv[6]  # Path to the checkbox states file
    checkbox_states = read_checkbox_states(checkbox_states_path)
    j = CalcCount(checkbox_states,img_index)
    file_tif = File_lif.get_image(j)
    output_file_path =f"../home/Tracked_{img_index+1}.webm"
    output_Metrics_path ="../home/Metrics.txt"
    Mean_velocity,Vescicole_Spola_locale,Vescicole_Spola_Nucleo_Membrana,Vescicole_erratiche,Vescicole_totali = Run_detector_tracker(show_Ids=input_boolean_values[img_index][0],tresh_val=thres_val_list[img_index],Fast_mode=input_boolean_values[img_index][1],Visibility = input_boolean_values[img_index][2],file_tif=file_tif,output_file=output_file_path,Show_erratic = input_boolean_values[img_index][3],Show_traget=input_boolean_values[img_index][4],Show_local=input_boolean_values[img_index][5],border=border[img_index],Nucleus=Nucleus[img_index])
    save_metrics(output_Metrics_path,Mean_velocity,Vescicole_erratiche,Vescicole_Spola_locale,Vescicole_Spola_Nucleo_Membrana,Vescicole_totali)



