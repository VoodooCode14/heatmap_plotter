'''
Created on Jul 27, 2022

@author: voodoocode
'''

import pyexcel
import numpy as np

import scipy.ndimage
import plotly.graph_objects

import nibabel

import plotly.io as pio
pio.renderers.default = "chromium"

import dash.html
import dash.dcc
import dash.dependencies

def read_xlsx_data(path = "/mnt/data/Professional/UHN/projects/data/minor/snr_artur.xlsx", feature_idx = 8):
    data = pyexcel.get_sheet(file_name = path)
    data = np.asarray(data)
    
    data = data[:, [2, 3, 4, feature_idx]]
    
    return data

def read_csv_data(path = "/mnt/data/Professional/UHN/projects/data/minor/POINT_CLOUD_all_flipped_to_right_all.csv"):
    data = pyexcel.get_sheet(file_name = path)
    data = np.asarray(data)
    
    return data[:, [3, 4, 5, -1]]
#    return data[:, [3, 4, 5, -2]]

def read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/STN.nii", threshold = 0.05):
    data = nibabel.load(path)
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

def print_data(data, res = 60, scale = None, threshold = None, surface_count = 20, colorscale = None, 
#def print_data(data, res = 50, scale = None, threshold = None, surface_count = 20, colorscale = None, 
               magnification = 3, filter_sz = 2, 
               isomin = None):
    res = res * magnification
    res_offset = res // 2
    (x_pos, y_pos, z_pos) = np.mgrid[:(res+1), :(res+1), :(res+1)]
    x_pos -= res_offset; y_pos -= res_offset; z_pos -= res_offset; 
    c_pos = np.empty(x_pos.shape)
    
    pre_x = np.asarray(data[1:, 0], dtype = float); pre_y = np.asarray(data[1:, 1], dtype = float)
    pre_z = np.asarray(data[1:, 2], dtype = float); pre_c = np.asarray(data[1:, 3], dtype = float)
    
    for idx in range(1, data.shape[0] - 1, 1):
        c_pos[int(pre_x[idx]*magnification + res_offset), int(pre_y[idx]*magnification + res_offset), int(pre_z[idx]*magnification + res_offset)] = pre_c[idx]
    
    if (filter_sz != 0):
        c_pos = scipy.ndimage.gaussian_filter(c_pos, filter_sz)
    
    if (scale is not None):
        c_pos *= scale
    
    if (isomin is None):
        if (threshold is not None):
            c_pos[c_pos < threshold] = 0
        
        c_pos -= np.min(c_pos)
        c_pos /= np.max(c_pos)
    
        if (threshold is None):
            threshold = 0.05
            
        if (colorscale is None):
            volume = plotly.graph_objects.Volume(x=x_pos.flatten(), y=y_pos.flatten(), z=z_pos.flatten(), value=c_pos.flatten(), 
                                                  isomin = threshold, isomax = 1, surface_count=surface_count, opacity = .5)
        else:
            volume = plotly.graph_objects.Volume(x=x_pos.flatten(), y=y_pos.flatten(), z=z_pos.flatten(), value=c_pos.flatten(), 
                                                  isomin = threshold, isomax = 1, surface_count=surface_count, opacity = .5, colorscale = colorscale)
            
    else:
        if (colorscale is None):
            volume = plotly.graph_objects.Volume(x=x_pos.flatten(), y=y_pos.flatten(), z=z_pos.flatten(), value=c_pos.flatten(), 
                                                  isomin = 4.5, surface_count=surface_count, opacity = .5)
        else:
            volume = plotly.graph_objects.Volume(x=x_pos.flatten(), y=y_pos.flatten(), z=z_pos.flatten(), value=c_pos.flatten(), 
                                                  isomin = 4.5, surface_count=surface_count, opacity = .5, colorscale = colorscale)
    
    return volume

def ms(x, y, z, radius, resolution=20):
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)

def plot_circle(x = 0, y = 0, z = 0, r = .2, opacity = 0.5, colorscale = "Blues", magnification = 1):
    (x_pns_surface, y_pns_surface, z_pns_suraface) = ms(x * magnification, y * magnification, z * magnification, r)
    return plotly.graph_objects.Surface(x=x_pns_surface, y=y_pns_surface, z=z_pns_suraface, opacity=opacity, colorscale = colorscale)
    
def run_in_dash_server(fig, port):
    
    app = dash.Dash()
    app.layout = dash.html.Div([
        dash.html.Div(id="output"),
        dash.dcc.Graph(id="fig", figure=fig)
    ])
    
    
    @app.callback(
        dash.dependencies.Output("output", "children"),
        dash.dependencies.Input("fig", "relayoutData")
    )
    def show_data(data):
        return [str(data)]
    
    
    app.run_server(debug=False, use_reloader=False, port = port)

X = [15.43039994, 15.25782405, 15.08524815, 14.91267225, 14.74009635, 14.56752045, 14.39494455, 14.22236865, 14.04979276, 13.87721686, 13.70464096, 13.53206506, 13.35948916, 13.18691326, 13.01433736, 12.84176147, 12.66918557, 13.6946996, 13.49337908, 13.29205856, 13.09073803, 12.88941751, 12.68809699, 12.48677647, 12.28545595, 12.08413542, 11.8828149, 11.68149438, 11.48017386, 11.27885334, 11.07753281, 10.87621229, 10.60518129, 10.41622208, 10.22726286, 10.03830365, 9.849344431, 9.741955334, 9.510326611, 9.278697889, 9.047069166, 8.815440444, 8.583811721, 8.352182999, 8.120554276, 7.888925554, 7.657296831, 7.425668109, 7.194039386, 6.962410664, 8.457567789, 8.283084758, 8.108601727, 7.934118696, 12.25835368, 12.21902678, 12.17969987, 12.14037297, 12.10104606, 12.06171916, 12.02239225, 11.98306535, 11.94373844, 11.90441154, 11.86508463, 11.82575773, 11.78643083, 11.74710392, 11.70777702, 11.66845011, 11.62912321, 13.88457881, 13.66007691, 13.43557502, 13.21107312, 12.98657122, 12.76206933, 12.53756743, 12.31306553, 12.08856364, 11.86406174, 11.63955985, 11.41505795, 12.30524332, 12.17717552, 12.04910772, 11.92103992, 11.79297212, 11.66490432, 11.53683652, 11.40876872, 11.28070092, 11.15263312, 11.02456532, 12.30023483, 12.1248106, 11.94938637, 11.77396215, 11.59853792, 11.42311369, 11.24768947, 11.07226524, 11.10152366, 10.90002883, 10.698534, 10.49703917, 10.29554434, 10.0940495, 9.892554671, 9.691059839, 9.489565008, 9.288070176, 9.086575344, 8.885080512, 11.93882864, 11.85422783, 11.76962702, 11.68502621, 11.60042539, 11.51582458, 11.43122377, 11.34662296, 11.26202215, 11.17742134, 14.93108495, 14.73549476, 14.53990457, 14.34431438, 14.14872418, 13.95313399, 13.7575438, 13.56195361, 13.36636342, 13.17077323, 12.97518303, 12.77959284, 11.22880306, 11.01095859, 10.79311413, 10.57526967, 10.3574252, 10.13958074, 9.921736276, 9.703891812, 9.486047349, 9.268202886, 9.050358423, 8.832513959, 12.9051, 12.73841667, 12.57173333, 12.40505, 12.23836667, 12.07168333, 11.905, 11.73831667, 11.57163333, 11.40495, 11.23826667, 10.6071, 10.43667778, 10.26625556, 10.09583333, 9.925411111, 9.754988889, 9.584566667, 9.414144444, 9.243722222, 9.0733, 8.902877778, 8.732455556, 13.25772811, 12.99683469, 12.73594127, 12.47504785, 12.21415443, 11.95326101, 11.69236758, 11.43147416, 11.17058074, 10.90968732, 13.01813273, 12.79979231, 12.58145189, 12.36311147, 12.14477105, 11.92643063, 11.7080902, 11.48974978, 11.27140936, 11.05306894, 10.83472852, 10.6163881, 10.39804768, 10.17970726, 9.961366835, 9.743026414, 11.2022, 11.04197414, 10.88174828, 10.72152243, 10.56129657, 10.40107071, 10.24084485, 13.64850622, 13.47079483, 13.29308344, 13.11537205, 12.93766066, 12.75994927, 12.58223788, 12.40452649, 12.2268151, 12.04910371, 11.87139232, 11.69368093]
Y = [-11.89599637, -12.29592816, -12.69585996, -13.09579175, -13.49572354, -13.89565534, -14.29558713, -14.69551892, -15.09545072, -15.49538251, -15.8953143, -16.29524609, -16.69517789, -17.09510968, -17.49504147, -17.89497327, -18.29490506, -10.91840855, -11.3258153, -11.73322205, -12.1406288, -12.54803555, -12.9554423, -13.36284904, -13.77025579, -14.17766254, -14.58506929, -14.99247604, -15.39988279, -15.80728954, -16.21469629, -16.62210303, -12.15036356, -12.467792, -12.78522044, -13.10264887, -13.42007731, -10.70164543, -10.92597428, -11.15030313, -11.37463198, -11.59896083, -11.82328968, -12.04761853, -12.27194738, -12.49627622, -12.72060507, -12.94493392, -13.16926277, -13.39359162, -11.24850162, -11.49353539, -11.73856917, -11.98360295, -11.28347575, -11.59060434, -11.89773293, -12.20486152, -12.51199011, -12.8191187, -13.1262473, -13.43337589, -13.74050448, -14.04763307, -14.35476166, -14.66189025, -14.96901884, -15.27614743, -15.58327602, -15.89040461, -16.1975332, -9.672345655, -10.06986318, -10.4673807, -10.86489822, -11.26241574, -11.65993326, -12.05745079, -12.45496831, -12.85248583, -13.25000335, -13.64752087, -14.04503839, -11.71541843, -12.10116479, -12.48691115, -12.87265751, -13.25840387, -13.64415023, -14.02989659, -14.41564295, -14.80138931, -15.18713567, -15.57288203, -14.19477469, -14.4697737, -14.74477271, -15.01977172, -15.29477073, -15.56976974, -15.84476875, -16.11976776, -14.21807991, -14.48508805, -14.75209619, -15.01910434, -15.28611248, -15.55312062, -15.82012876, -16.0871369, -16.35414504, -16.62115319, -16.88816133, -17.15516947, -12.93932316, -13.24350189, -13.54768062, -13.85185935, -14.15603809, -14.46021682, -14.76439555, -15.06857428, -15.37275301, -15.67693175, -11.70451005, -12.01712485, -12.32973965, -12.64235445, -12.95496924, -13.26758404, -13.58019884, -13.89281364, -14.20542844, -14.51804324, -14.83065804, -15.14327284, -13.25977241, -13.56610848, -13.87244454, -14.17878061, -14.48511667, -14.79145274, -15.09778881, -15.40412487, -15.71046094, -16.016797, -16.32313307, -16.62946914, -12.285, -12.65310556, -13.02121111, -13.38931667, -13.75742222, -14.12552778, -14.49363333, -14.86173889, -15.22984444, -15.59795, -15.96605556, -12.3434, -12.67767778, -13.01195556, -13.34623333, -13.68051111, -14.01478889, -14.34906667, -14.68334444, -15.01762222, -15.3519, -15.68617778, -16.02045556, -9.630198349, -9.979760385, -10.32932242, -10.67888446, -11.02844649, -11.37800853, -11.72757057, -12.0771326, -12.42669464, -12.77625667, -10.44883616, -10.77354756, -11.09825895, -11.42297035, -11.74768174, -12.07239314, -12.39710453, -12.72181593, -13.04652732, -13.37123872, -13.69595011, -14.02066151, -14.34537291, -14.6700843, -14.9947957, -15.31950709, -11.6864, -12.0063966, -12.32639321, -12.64638981, -12.96638641, -13.28638302, -13.60637962, -9.962584067, -10.35647674, -10.75036942, -11.1442621, -11.53815478, -11.93204745, -12.32594013, -12.71983281, -13.11372549, -13.50761816, -13.90151084, -14.29540352]
Z = [-3.93683456, -4.348831902, -4.760829244, -5.172826586, -5.584823928, -5.99682127, -6.408818613, -6.820815955, -7.232813297, -7.644810639, -8.056807981, -8.468805323, -8.880802665, -9.292800007, -9.704797349, -10.11679469, -10.52879203, -4.675839286, -5.089901536, -5.503963787, -5.918026037, -6.332088288, -6.746150538, -7.160212789, -7.574275039, -7.98833729, -8.40239954, -8.816461791, -9.230524041, -9.644586292, -10.05864854, -10.47271079, -7.291704669, -7.71511097, -8.138517272, -8.561923574, -8.985329876, -7.969865572, -8.472430054, -8.974994535, -9.477559016, -9.980123498, -10.48268798, -10.98525246, -11.48781694, -11.99038142, -12.49294591, -12.99551039, -13.49807487, -14.00063935, -8.970555348, -9.4601128, -9.949670252, -10.4392277, -5.82204472, -6.214632829, -6.607220938, -6.999809048, -7.392397157, -7.784985266, -8.177573375, -8.570161484, -8.962749594, -9.355337703, -9.747925812, -10.14051392, -10.53310203, -10.92569014, -11.31827825, -11.71086636, -12.10345447, -5.056881919, -5.460178434, -5.863474949, -6.266771464, -6.670067979, -7.073364493, -7.476661008, -7.879957523, -8.283254038, -8.686550552, -9.089847067, -9.493143582, -6.144744851, -6.559441323, -6.974137795, -7.388834266, -7.803530738, -8.21822721, -8.632923682, -9.047620154, -9.462316625, -9.877013097, -10.29170957, -6.571604228, -7.030078869, -7.488553509, -7.94702815, -8.40550279, -8.863977431, -9.322452071, -9.780926712, -7.915549519, -8.374194475, -8.832839431, -9.291484386, -9.750129342, -10.2087743, -10.66741925, -11.12606421, -11.58470916, -12.04335412, -12.50199908, -12.96064403, -6.453833134, -6.892786566, -7.331739999, -7.770693431, -8.209646864, -8.648600296, -9.087553729, -9.526507161, -9.965460593, -10.40441403, -3.927575729, -4.372350616, -4.817125502, -5.261900389, -5.706675276, -6.151450162, -6.596225049, -7.040999936, -7.485774823, -7.930549709, -8.375324596, -8.820099483, -7.04407659, -7.491774771, -7.939472951, -8.387171132, -8.834869313, -9.282567494, -9.730265674, -10.17796386, -10.62566204, -11.07336022, -11.5210584, -11.96875658, -5.48466, -5.96351, -6.44236, -6.92121, -7.40006, -7.87891, -8.35776, -8.83661, -9.31546, -9.79431, -10.27316, -7.03489, -7.506149778, -7.977409556, -8.448669333, -8.919929111, -9.391188889, -9.862448667, -10.33370844, -10.80496822, -11.276228, -11.74748778, -12.21874756, -4.898813341, -5.253327839, -5.607842338, -5.962356837, -6.316871336, -6.671385835, -7.025900334, -7.380414832, -7.734929331, -8.08944383, -5.136006905, -5.55589233, -5.975777755, -6.39566318, -6.815548605, -7.235434031, -7.655319456, -8.075204881, -8.495090306, -8.914975731, -9.334861157, -9.754746582, -10.17463201, -10.59451743, -11.01440286, -11.43428828, -6.95158, -7.397947751, -7.844315502, -8.290683253, -8.737051004, -9.183418755, -9.629786505, -4.670876243, -5.070992601, -5.471108959, -5.871225318, -6.271341676, -6.671458034, -7.071574393, -7.471690751, -7.87180711, -8.271923468, -8.672039826, -9.072156185]

X = np.asarray(X)
Y = np.asarray(Y)
Z = np.asarray(Z)

def main(magnification = 2, filter_sz = 1.5):
    volume1 = print_data(read_xlsx_data(feature_idx = 8), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume2 = print_data(read_xlsx_data(feature_idx = 9), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume3 = print_data(read_xlsx_data(feature_idx = 10), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume4 = print_data(read_xlsx_data(feature_idx = 11), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume5 = print_data(read_xlsx_data(feature_idx = 12), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume6 = print_data(read_xlsx_data(feature_idx = 13), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume7 = print_data(read_xlsx_data(feature_idx = 14), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    volume8 = print_data(read_xlsx_data(feature_idx = 15), surface_count = 20, colorscale = None, threshold = 0.0, magnification = magnification, filter_sz = filter_sz, isomin = 4.5)
    
    #snr1 = print_data(read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/snr/rh_SN_cluster3_1_ncut.nii"), surface_count = 1, colorscale = "Reds", magnification = magnification, filter_sz = filter_sz)
    #snr2 = print_data(read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/snr/rh_SN_cluster3_1_ncut.nii"), surface_count = 1, colorscale = "Reds", magnification = magnification, filter_sz = filter_sz)
    #snr3 = print_data(read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/snr/rh_SN_cluster3_1_ncut.nii"), surface_count = 1, colorscale = "Reds", magnification = magnification, filter_sz = filter_sz)
    
    snr = print_data(read_nii_data(path = "/mnt/data/Professional/UHN/projects/data/minor/snr/substantia_nigra_rh.nii"), surface_count = 1, colorscale = "Reds", magnification = magnification, filter_sz = filter_sz)
    
    scene = 0
    cam = 1
    
    if (scene == 0):
        fig = plotly.graph_objects.Figure(data = [snr, volume1])
    if (scene == 1):
        fig = plotly.graph_objects.Figure(data = [snr, volume2])
    if (scene == 2):
        fig = plotly.graph_objects.Figure(data = [snr, volume3])
    if (scene == 3):
        fig = plotly.graph_objects.Figure(data = [snr, volume4])
    if (scene == 4):
        fig = plotly.graph_objects.Figure(data = [snr, volume5])
    if (scene == 5):
        fig = plotly.graph_objects.Figure(data = [snr, volume6])
    if (scene == 6):
        fig = plotly.graph_objects.Figure(data = [snr, volume7])
    if (scene == 7):
        fig = plotly.graph_objects.Figure(data = [snr, volume8])
    
    if (cam == 1):
        camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0.21810523329138587, y=-0.3223437443149643, z=-0.20916416611249108),
            eye=dict(x=-0.042510031300014406, y=0.13708987387861754, z=-0.11215032051561241)
        )
    if (cam == 2):
        camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0.21810523329138587, y=-0.3223437443149643, z=-0.20916416611249108),
            eye=dict(x=-0.062414808939955047, y=-0.6009302179431355, z=0.15430431191543018)
        )
    
    fig.update_layout(scene_camera=camera)
    
    fig.update_layout(autosize=False, width=1920, height=1000, margin=dict(l=0, r=0, b=0, t=0, pad=0),)
    
    run_in_dash_server(fig, port = 9050 + cam)

main()
print("terminated successfully")


#===============================================================================
# Data:
#  - Norm 
#  - Norm2
# 
# Angle [1, 2, 3, 4]:
#  - {'scene.camera': {'up': {'x': 0, 'y': 0, 'z': 1}, 'center': {'x': 0.16760867111505937, 'y': -0.2842659394019215, 'z': -0.14406079269444813}, 'eye': {'x': -0.009349092920459656, 'y': -0.4571475498663773, 'z': 0.02890372582033912}, 'projection': {'type': 'perspective'}}}
#  - {'scene.camera': {'up': {'x': 0, 'y': 0, 'z': 1}, 'center': {'x': 0.1581529296852149, 'y': -0.27680225021166854, 'z': -0.18098527774133144}, 'eye': {'x': -0.16716140549851857, 'y': -0.3060092562947876, 'z': -0.10358145700536278}, 'projection': {'type': 'perspective'}}}
#  - {'scene.camera': {'up': {'x': 0, 'y': 0, 'z': 1}, 'center': {'x': 0.2217874339950377, 'y': -0.27431829834283367, 'z': -0.18199247994440507}, 'eye': {'x': 0.22190798952539384, 'y': -0.6017262122608962, 'z': -0.10797990861299861}, 'projection': {'type': 'perspective'}}}
#  - {'scene.camera': {'up': {'x': 0, 'y': 0, 'z': 1}, 'center': {'x': 0.22205540256695705, 'y': -0.2537726899329517, 'z': -0.17960391762565964}, 'eye': {'x': 0.22241370927615767, 'y': -0.29253996422096545, 'z': 0.1538189181332006}, 'projection': {'type': 'perspective'}}}
# 
# Scenes [1, 2, 3, 4]:
#  - All
#  - STN center only (1 green)
#  - Other STN center only (1 green)
#  - No dots
# 
# Norm_[1, 2, 3, 4]_[A, B, C, D]
#===============================================================================

