print('running')

from cellpose import models
import cellpose
import skimage
import numpy as np
import os
import matplotlib.pyplot as plt
import cv2
import math
from tqdm import tqdm
import yaml
import pandas as pd
from aicsimageio.writers import OmeTiffWriter
from roifile import ImagejRoi
import easygui
from skan import Skeleton, summarize



def SegmentWithCellpose(im, model='nuclei'):
    model = models.Cellpose(model_type=model)
    channels = [[0,0]]  
    masks, flows, styles, diams = model.eval(im, diameter=None, flow_threshold=None, channels=channels)
    return masks

def SegmentWithCellposeCustom(im, model_path, diam):
    model = models.CellposeModel(pretrained_model=model_path)
    #channels = [[0,0]]
    masks, flows, styles = model.eval(im, diameter=diam, flow_threshold=0.4)
    return masks

def RemoveSmallObjecs(label_im, threshold):
    props = skimage.measure.regionprops(label_im)
    for prop in props:
            if prop['area'] < threshold:
                label_im[label_im == prop['label']] = 0
    return label_im


def GetPerinuclearRegion(nuc_im, lab, radius=15):
    Im = np.zeros(nuc_im.shape, dtype=np.uint8)
    Im[nuc_im == lab] = 1
    rad = radius
    Im1_d = cv2.dilate(Im, skimage.morphology.disk(rad), iterations=1)
    Im1_r = np.logical_and(Im1_d, np.logical_not(Im))
    Im1_r[nuc_im > 0] = 0
    return Im1_r

def GetOutlineImage(image, struct):
    
    if image.dtype == 'bool':
         image = image.astype(np.uint8)
    im = cv2.dilate(image, struct, iterations=1)
    im_sub = np.logical_and(im, np.logical_not(image))
    return im_sub


def AddMeasurementToImage(baseim, measureim, mean, struct):
    outline_im = GetOutlineImage(measureim, struct)
    baseim[outline_im > 0] = 255
    labs = skimage.measure.label(measureim)
    tempprops = skimage.measure.regionprops(labs)
    text_string = str(round(mean,0))
    location = (int(tempprops[0]['centroid'][1]-20),int(tempprops[0]['centroid'][0]))
    baseim = cv2.putText(baseim, text_string, location, font, .5,(255,255,255),1)
    return baseim

def swap_in(lst, fro, to):
    lst[fro], lst[to] = lst[to], lst[fro]


def get_roi_from_im(image):
    con = skimage.measure.find_contours(image, level=0.9999)
    if len(con) > 1:
        concat_list = np.concatenate([con[0], con[1]], axis=0).tolist()
    else:
        concat_list = con[0]
    for sublist in concat_list:
        swap_in(sublist, 0, 1)
    roi = ImagejRoi.frompoints(concat_list)
    return roi

def box_with_skimage_bbox(image, bbox, val=1):
    xmin = bbox[0] - 5
    ymin = bbox[1] -5
    xmax = bbox[2] + 5
    ymax = bbox[3] + 5

    xx,yy = skimage.draw.polygon_perimeter([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin])

    try:
        image[xx, yy] = val
    except: 
        print('caught exception')
    return image

def AddTextToImage(im, length, mean, cilia_prop):
    text_string = 'length = ' + str(round(length,2))
    location = (int(cilia_prop['bbox'][1]-20),int(cilia_prop['bbox'][0]-40))
    im = cv2.putText(im, text_string, location, font, .5,(255,255,255),1)
    text_string = 'kif int = ' + str(round(mean, 0))
    location = (int(cilia_prop['bbox'][1]-20),int(cilia_prop['bbox'][0]-20))
    im = cv2.putText(im, text_string, location, font, .5,(255,255,255),1)
    return im

font = cv2.FONT_HERSHEY_SIMPLEX
struct = np.ones((3,3))


CONFIG_NAME = 'config.yaml'
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

cwd = os.getcwd()
print(cwd)

with open(CONFIG_NAME, "r") as f:
        config = yaml.safe_load(f)


def MakeDividingCellMask(labs, raw):
    props = skimage.measure.regionprops(labs)
    mean_vals = []
    for prop in props:
        meanval = np.mean(raw[labs == prop['label']])
        mean_vals.append(meanval)
    
        plt.hist(mean_vals)

raw_data_location = easygui.diropenbox('Select raw data location')
output_data_location = easygui.diropenbox('Select output location')
cilia_location = easygui.diropenbox('Select classified cilia location')

###for testing purposes
#raw_data_location = 'test_data'
#output_data_location = 'output_5chan'
#cilia_location = 'classified_cilia'

cellpose_model_path_cell = config['cellpose_model_path_cell']
cellpose_model_path_nuc = config['cellpose_model_path_nuc']


min_cilia_area = config['min_cilia_area']
min_cell_area = config['min_cell_area']
pernuclear_expansion = config['perinuc_expansion']

nuc_channel = config['nuc_channel'] - 1
kif_channel = config['kif_channel'] - 1
mem_channel = config['mem_channel'] - 1
arl_channel = config['arl_channel'] - 1
second_measure_channel = config['second_measure_channel'] - 1


input_data_list = os.listdir(raw_data_location)
print(input_data_list)

for file in input_data_list:
    if file == "Thumbs.db":
        input_data_list.remove(file)
    
print(input_data_list)

for i, file in enumerate(input_data_list):

    if i > -1:
        print('===================================================================')
        print('running file ' + file)
        df = pd.DataFrame()

        output_data_location_file = output_data_location + os.path.sep + file
        os.makedirs(output_data_location_file, exist_ok=True)
        os.makedirs(output_data_location_file + os.path.sep + 'cilia_images', exist_ok=True)
        os.makedirs(output_data_location_file + os.path.sep + 'measurement_rois', exist_ok=True)
        mean_array = []
        second_mean_array = []
        cilia_array = []
        cilia_length_array = []
        nuclear_label_array = []
        cilia_label_array = []

        positive_nuclei = []  #keep track of the nuclei with associated cilia already
    

        im = skimage.io.imread(raw_data_location + os.path.sep + file)
        
        print('data loaded, starting segmentation')
        nuc_masks = SegmentWithCellposeCustom(im[:,:,nuc_channel], model_path=cellpose_model_path_nuc, diam=78.24)
        nuc_masks_noedge = np.copy(nuc_masks)
        nuc_masks_noedge = cellpose.utils.remove_edge_masks(nuc_masks_noedge, change_index=False)
        nuc_props = skimage.measure.regionprops(nuc_masks)
 
        cell_masks = SegmentWithCellposeCustom(im[:,:,kif_channel], model_path=cellpose_model_path_cell, diam=175)
        #cell_masks = cellpose.utils.remove_edge_masks(cell_masks, change_index=True)
        cell_masks_filtered = RemoveSmallObjecs(cell_masks, min_cell_area)
        print('segmentation finished, beginning processing')

        MakeDividingCellMask(cell_masks, im[:,:,nuc_channel])


        measurement_image = np.zeros(cell_masks.shape, dtype = np.uint8)
        
        
        skel_all = np.zeros(cell_masks.shape, dtype=np.uint16)
        skel_overlay = np.copy(im[:,:,arl_channel])
        if  np.amax(skel_all) > 65500:
            max_skel_value = np.amax(skel_all) * 3
        else:
            max_skel_value = 65500
        cilia = skimage.io.imread(cilia_location + os.path.sep + 'classified_'+ file)

        cilia_labeled = skimage.measure.label(cilia)
        cilia_labeled = RemoveSmallObjecs(cilia_labeled, min_cilia_area)
        cilia_props = skimage.measure.regionprops(cilia_labeled)


        print('total number of cilia detected: ' + str(len(cilia_props)))
        print('processing detected cilia')
        for cilia_prop in tqdm(cilia_props):
            #print('processing cilia ' + str(cilia_prop['label']))
            cilia_copy = np.zeros(cilia_labeled.shape)
            cilia_copy[cilia_labeled == cilia_prop['label']] = 1


            #get cilia length using major axis
            #cilia_length = cilia_prop['axis_major_length']
            
            
            cilia_copy = np.zeros(cilia_labeled.shape)
            cilia_copy[cilia_labeled == cilia_prop['label']] = 1
            skel = skimage.morphology.skeletonize(cilia_copy)
            skel_all[skel > 0] = max_skel_value
            skel_overlay[skel > 0] = max_skel_value
            branch_data = summarize(Skeleton(skel))
            cilia_length = branch_data['branch-distance']
            #print('cilia length: ' + str(cilia_length))
        


            containing_cell_label = np.argmax(np.bincount(cell_masks[cilia_copy == 1]))
            #print('matching cell label is: ' + str(containing_cell_label))
      
            if containing_cell_label > 0:
                #find the matching nucleus

                all_nuclei_in_cell = np.unique([nuc_masks[cell_masks == containing_cell_label]])
                all_nuclei_in_cell = all_nuclei_in_cell[all_nuclei_in_cell>0]

                #find the nucleus with the maximum percentage overlap
                max_overlap = 0
                max_nucleus = -1
                for nuclab in all_nuclei_in_cell:
        
                    nucleus_size = np.count_nonzero(nuc_masks[nuc_masks == nuclab])
                    masked_nuclear_image = np.copy(nuc_masks)
                    masked_nuclear_image[cell_masks != containing_cell_label] = 0
                    overlap_size = np.count_nonzero(masked_nuclear_image[masked_nuclear_image == nuclab])
                    overlap_percent = overlap_size / nucleus_size
                    if overlap_percent > max_overlap:
                        max_overlap = overlap_percent
                        max_nucleus = nuclab
                #print('matching nucleus label is: ' + str(max_nucleus))

                #check if the nucleus is on an edge, ignore if it is:
                
                if np.sum(nuc_masks_noedge[nuc_masks == max_nucleus]) == 0:
                    print('edge mask at label ' + str(max_nucleus))

                nuc_to_dilate = np.zeros(nuc_masks.shape)
                nuc_to_dilate[nuc_masks == max_nucleus] = 1
  
                measure_im = GetPerinuclearRegion(nuc_masks, max_nucleus, radius = pernuclear_expansion)

          
                roi = get_roi_from_im(measure_im)
                roi.tofile(output_data_location_file + os.path.sep + 'measurement_rois' + os.path.sep + str(max_nucleus) + '.roi')


                mean = round(np.mean(im[:,:,kif_channel][measure_im > 0]), 2)
                mean_array.append(mean)
                if second_measure_channel > 0:
                    mean_ch2 = round(np.mean(im[:,:,second_measure_channel][measure_im > 0]), 2)
                    second_mean_array.append(mean_ch2)

                cilia_array.append("Y")
                cilia_length_array.append(cilia_length)
                nuclear_label_array.append(max_nucleus) 
                cilia_label_array.append(str(cilia_prop['label']))

                positive_nuclei.append(max_nucleus)

                single_nuc_im = np.zeros(nuc_masks.shape, dtype = np.uint16)
                single_cell_im = np.zeros(nuc_masks.shape, dtype = np.uint16)
                single_measure_im = np.zeros(nuc_masks.shape, dtype = np.uint16)
                single_nuc_im[nuc_masks == max_nucleus] = max_nucleus
                single_cell_im[cell_masks == containing_cell_label] = containing_cell_label
                single_measure_im = AddMeasurementToImage(single_measure_im, measure_im, mean, struct)

                memtemp = np.copy(im[:,:,mem_channel])
                mem_annotated = box_with_skimage_bbox(memtemp, cilia_prop['bbox'], val=np.amax(memtemp))

                kiftemp = np.copy(im[:,:,kif_channel])

                con = skimage.measure.find_contours(single_nuc_im, level=0.9999)
                for c in con[0]:
                    kiftemp[int(c[0]), int(c[1])] = np.amax(kiftemp)

                con = skimage.measure.find_contours(single_cell_im, level=0.9999)
                for c in con[0]:
                    kiftemp[int(c[0]), int(c[1])] = np.amax(kiftemp)

                b = cv2.normalize(kiftemp, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                g = cv2.normalize(im[:,:,arl_channel], None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                r = cv2.normalize(mem_annotated, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                bgr = cv2.merge((b,g,r))
                bgr = AddTextToImage(bgr, cilia_length, mean, cilia_prop)
                cv2.imwrite(output_data_location_file + os.path.sep + 'cilia_images' + os.path.sep + 'cilia_' + str(cilia_prop['label']) + '_' + file + '.png', bgr)  

                measurement_image = AddMeasurementToImage(measurement_image, measure_im, mean, struct)
             
                   
            if containing_cell_label == 0:
                #find the closest nucleus
                min_dist = 9999
                for nucprop in nuc_props:
                    nuc_centroid = nucprop['centroid']
                    dist = math.dist(nuc_centroid, cilia_prop['centroid'])
                    if dist < min_dist:
                        min_dist = dist
                        closest_nucleus = nucprop['label']

                measure_im = GetPerinuclearRegion(nuc_masks, closest_nucleus, radius = pernuclear_expansion)
                roi = get_roi_from_im(measure_im)
                roi.tofile(output_data_location_file + os.path.sep + 'measurement_rois' + os.path.sep + str(closest_nucleus) + '.roi')
                mean = round(np.mean(im[:,:,kif_channel][measure_im > 0]), 2)

                if second_measure_channel > 0:
                    mean_ch2 = round(np.mean(im[:,:,second_measure_channel][measure_im > 0]), 2)
                    second_mean_array.append(mean_ch2)

                #print('mean value: ' + str(mean))
                mean_array.append(mean)

                cilia_array.append("Y")
                cilia_length_array.append(cilia_length)
                nuclear_label_array.append(closest_nucleus) 
                cilia_label_array.append(str(cilia_prop['label']))

                positive_nuclei.append(closest_nucleus)

                measurement_image = AddMeasurementToImage(measurement_image, measure_im, mean, struct)

                single_nuc_im = np.zeros(nuc_masks.shape, dtype = np.uint16)
                single_measure_im = np.zeros(nuc_masks.shape, dtype = np.uint16)
                single_nuc_im[nuc_masks == closest_nucleus] = closest_nucleus
                single_measure_im = AddMeasurementToImage(single_measure_im, measure_im, mean, struct)

                memtemp = np.copy(im[:,:,mem_channel])
                kiftemp = np.copy(im[:,:,kif_channel])

                mem_annotated = box_with_skimage_bbox(memtemp, cilia_prop['bbox'], val=np.amax(memtemp))
            

                con = skimage.measure.find_contours(single_nuc_im, level=0.9999)
                for c in con[0]:
                    kiftemp[int(c[0]), int(c[1])] = np.amax(kiftemp)

                b = cv2.normalize(kiftemp, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                g = cv2.normalize(im[:,:,arl_channel], None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                r = cv2.normalize(mem_annotated, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
                bgr = cv2.merge((b,g,r))
                bgr = AddTextToImage(bgr, cilia_length, mean, cilia_prop)
                cv2.imwrite(output_data_location_file + os.path.sep + 'cilia_images' + os.path.sep + 'cilia_' + str(cilia_prop['label']) + '_' + file + '.png', bgr)  


        #got through the rest of the images
        print('running nonciliated nuclei')
        for nucprop in tqdm(nuc_props):
            if nucprop['label'] not in positive_nuclei:
                measure_im = GetPerinuclearRegion(nuc_masks, nucprop['label'], radius = pernuclear_expansion)
                mean = round(np.mean(im[:,:,kif_channel][measure_im > 0]), 2)
                mean_array.append(mean)
                if second_measure_channel > 0:
                    mean_ch2 = round(np.mean(im[:,:,second_measure_channel][measure_im > 0]), 2)
                    second_mean_array.append(mean_ch2)
                cilia_array.append("N")
                cilia_length_array.append("NA")
                nuclear_label_array.append(nucprop['label']) 
                cilia_label_array.append('NA')
                measurement_image = AddMeasurementToImage(measurement_image, measure_im, mean, struct)
                roi = get_roi_from_im(measure_im)
                roi.tofile(output_data_location_file + os.path.sep + 'measurement_rois' + os.path.sep + str(nucprop['label']) + '.roi')
    
        stacked = np.stack((im[:,:,nuc_channel], im[:,:,kif_channel], im[:,:,mem_channel], im[:,:,arl_channel], nuc_masks.astype(np.uint16), nuc_masks_noedge.astype(np.uint16), cell_masks.astype(np.uint16), cilia_labeled.astype(np.uint16), skel_all.astype(np.uint16),  skel_overlay.astype(np.uint16), measurement_image.astype(np.uint16)))
        OmeTiffWriter.save(stacked, output_data_location_file + os.path.sep + 'summary_' + file, dim_order="CYX")

        df['Mean Intensity'] = mean_array
        if second_measure_channel > 0:
            df['Mean Intensity Second Channel'] = second_mean_array
        df['Cilia Present'] = cilia_array
        df['Cilia Length'] = cilia_length_array
        df['Nucleus Label'] = nuclear_label_array
        df['Cilia Label'] = cilia_label_array

        df.to_csv(output_data_location_file + os.path.sep + 'output.csv')
        


print("Pipeline Complete")