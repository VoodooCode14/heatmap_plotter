'''
Created on Jul 27, 2022

@author: voodoocode
'''

import numpy as np

import finnpy.file_io.csv_data
import finnpy.file_io.nifti_data
import finnpy.visualization.volumetric_plots

import mayavi.mlab

#===============================================================================
# # Default
# # CSV_FILE = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_all_shifted.csv"
# # DATA_LEVEL = -2
#===============================================================================

CSV_FILE = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_all_shifted.csv"
#CSV_FILE = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_all.csv"

#DATA_LEVEL = -2
DATA_LEVEL = -1

def read_csv_data(path = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_all_shifted.csv", level = -2):
    data = finnpy.file_io.csv_data.run(path)
    data = np.asarray(data)
    
    return data[:, [3, 4, 5, level]]    

def read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/STN.nii", threshold = 0.05):
    data = finnpy.file_io.nifti_data.run(path)
    values = list(); coords = list()
    for x in range(0, data.shape[0]):
        for y in range(0, data.shape[1]):
            for z in range(0, data.shape[2]):
                if (data.get_fdata()[x, y, z] > threshold):
                    coords.append([x, y, z]); values.append(data.get_fdata()[x, y, z])
    coords = np.asarray(coords); values = np.asarray(values)
    trans = data._affine
    coords = np.dot(coords, trans[:3, :3].T)
    coords += trans[:3, 3]
    
    return np.concatenate((coords, np.expand_dims(values, axis = 1)), axis = 1)
    
def plot_erna_pts():
    for (row_idx, values) in enumerate(read_csv_data(path = CSV_FILE, level = DATA_LEVEL)):
        if (row_idx == 0):
            continue
        finnpy.visualization.volumetric_plots.plot_circle(x = float(values[0]), y = float(values[1]), z = float(values[2]), r = .2, c_val = float(values[3]))

# Need to be updated to new csv file
#===============================================================================
# def plot_updrs_pts():
#     x_pts = list(); y_pts = list(); z_pts = list(); c_pts = list()
#     for (row_idx, values) in enumerate(read_csv_data(path = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_max.csv", level = -3)):
#         if (row_idx == 0):
#             continue
#         x_pts.append(float(values[0])); y_pts.append(float(values[1])); z_pts.append(float(values[2])); c_pts.append(-float(values[3]))
#     finnpy.visualization.volumetric_plots.plot_circle(x = x_pts, y = y_pts, z = z_pts, r = .2, c_val = c_pts, colormap = "RdBu")
#     
# def plot_rigid_pts():
#     for (row_idx, values) in enumerate(read_csv_data(path = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_max.csv", level = -1)):
#         if (row_idx == 0):
#             continue
#         finnpy.visualization.volumetric_plots.plot_circle(x = float(values[0]), y = float(values[1]), z = float(values[2]), r = .2, c_val = float(values[3]))
#===============================================================================

def perspective_tester():
    
    import mayavi.api
     
    engine = mayavi.api.Engine()
    engine.start()
    scene = engine.new_scene()
     
    scene.scene.camera.position = [-10.853422308877313, -29.88176360985174, 16.32684811852981]
    scene.scene.camera.focal_point = [11.0964053, -13.89933787, -9.06410754]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.412711194846655, 0.5664308609453176, 0.7133200890322611]
    scene.scene.camera.clipping_range = [18.019153692902677, 61.376498749008256]
    
    print(mayavi.mlab.view())

def set_perspective():
    mayavi.mlab.view(-143.94034704238993, 46.91965685864011, 37.1743660852843, np.asarray([ 11.0964053 , -13.89933787,  -9.06410754]))

#def main(magnification = 2, filter_sz = 1.5):
#def main(magnification = 5, filter_sz = 5, gaussian_filter_sz = .75, 
def main(magnification = 5, filter_sz = 5, gaussian_filter_sz = 1.25, 
         otf_pts = [[1, 1], [.9, 1], [.1, .4], [.09, 0], [0, 0]]):
    
    finnpy.visualization.volumetric_plots.create_figure()
    
    #otf_pts = [[1, .66], [.9, .66], [.5, .25], [.49, 0], [0, 0]]
    #otf_pts = [[1, .66], [.9, .66], [.1, .25], [.09, 0], [0, 0]]
    otf_pts = [[1, 1], [.75, .66], [.70, .33], [.1, .125], [.09, 0], [0, 0]]
    #otf_pts = [[1, 1], [.75, .66], [.70, .4], [.1, .125], [.09, 0], [0, 0]]
    #otf_pts = [[1, 1], [.2, .4], [0, 0]]
    otf_pts = [[1, 1], [.75, .66], [.70, .5], [.1, .25], [.09, 0], [0, 0]]
    otf_pts = [[1, 1], [.75, .8], [.70, .25], [.1, .125], [.09, 0], [0, 0]]
    otf_pts = [[1, 1], [.75, .9], [.70, .6], [.45, .5], [.1, .25], [.09, 0], [0, 0]]
    
    #ctf_pts = [[1, 255/255, 0/255, 0/255], [.75, 255/255, 255/255, 0/255], [.65, 0/255, 255/255, 127/255], [0, 0/255, 0/255, 255/255]]
    
    #log config
    ctf_pts = [[1, 255/255, 255/255, 0/255], [0.7, 0/255, 255/255, 255/255], [0, 0/255, 0/255, 255/255]]
    if (DATA_LEVEL == -2):
        c_min_val = 0.25
        c_values = [1, -1, (c_min_val + 1)/2, -1, c_min_val]; c_values[1] = (c_values[0] + c_values[2])/2; c_values[3] = (c_values[2] + c_values[4])/2
        ctf_pts = [[c_values[0], 255/255, 0/255, 0/255], [c_values[1], 255/255, 255/255, 0/255], [c_values[2], 0/255, 255/255, 0/255], [c_values[3], 0/255, 255/255, 255/255], [c_values[4], 255/255, 0/255, 0/255]]
        otf_pts = [[1, 1], [c_min_val, 0], [0, 0]]
    if (DATA_LEVEL == -1):
        c_min_val = 0.3
        c_values = [1, -1, (c_min_val + 1)/2, -1, c_min_val]; c_values[1] = (c_values[0] + c_values[2])/2; c_values[3] = (c_values[2] + c_values[4])/2
        ctf_pts = [[c_values[0], 255/255, 0/255, 0/255], [c_values[1], 255/255, 255/255, 0/255], [c_values[2], 0/255, 255/255, 0/255], [c_values[3], 0/255, 255/255, 255/255], [c_values[4], 255/255, 0/255, 0/255]]
        otf_pts = [[1, 1], [c_min_val, 0], [0, 0]]
    
    
    
    input_data = read_csv_data(CSV_FILE, DATA_LEVEL)[1:, :]; input_data = np.asarray(input_data, dtype = float);
    if (DATA_LEVEL == -1):
        non_zero_idx = np.argwhere(input_data[:, -1] != 0).squeeze(1)
        input_data[non_zero_idx, -1] = np.log(input_data[non_zero_idx, -1])# / np.log(10000)
        input_data[:, -1] += 1
    volume = finnpy.visualization.volumetric_plots.plot_volumetric_data(input_data,
                                                                        magnification, nanmean_filter_sz = filter_sz,
                                                                        gaussian_filter_sz = gaussian_filter_sz,
                                                                        ctf_pts = ctf_pts, otf_pts = otf_pts, mode = "mean")
 
    #plot_erna_pts()
    #plot_updrs_pts()
    #plot_rigid_pts()
    
    finnpy.visualization.volumetric_plots.plot_contour3d_data(read_nii_data(), magnification = magnification, color = (.7, .7, .7))
    
    set_perspective()
    #perspective_tester()
    
    #finnpy.visualization.volumetric_plots.plot_circle(x = 12.5, y = -12.995, z = -5.895, r = 1, c_val = 1.0)
    finnpy.visualization.volumetric_plots.plot_circle(x = 12.5, y = -12.995, z = -5.895, r = .2, c_val = 1.0, color = (0, 0, 0))
    
    mayavi.mlab.colorbar(volume)
    
    finnpy.visualization.volumetric_plots.show_figure()

main()
print("terminated successfully")

