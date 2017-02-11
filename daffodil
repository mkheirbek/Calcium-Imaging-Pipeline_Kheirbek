
# coding: utf-8

# # Daffodil
# 
# Daffodil extracts df/f signals from the motion corrected videos. It uses 3 steps:
# 
#  - automatic detection using ICA/PCA from Dombeck's scripts
#  - roi extraction using custom algorithm
#  - manual selection of final ROIs with ImageJ

# In[91]:

import os
import sima
import glob
import pylab as pl
import numpy as np
import time
import tifffile
from scipy.io import loadmat
from ipywidgets import interact

import skimage
from skimage.measure import find_contours
from skimage.filters import gaussian
from skimage.filters import gaussian
from skimage.segmentation import active_contour, clear_border
from skimage.measure import label, regionprops, approximate_polygon
from skimage.morphology import closing, square, opening, erosion, disk
from skimage.color import label2rgb

get_ipython().magic(u'matplotlib inline')
# pl.rcParams["figure.facecolor"] = 'white'
pl.rcParams['savefig.dpi'] = 120


# In[92]:

def calc_dFoF_oldschool(traces):
    return (traces-traces.mean(0))/traces.mean(0)

def calc_dFoF_jia(traces, tau0=0.2, tau1=0.75, tau2=3):
    ntau1 = int(tau1 * fps)
    fbar = traces.copy()
    for i in xrange(traces.shape[1]):
        fbar[:, i] = np.convolve(traces[:, i], np.ones(ntau1).astype(float)/ntau1, mode='same')
    ntau2 = int(tau2 * fps)
    fzero = fbar.copy()
    for i in xrange(traces.shape[1]):
        for t in xrange(1, len(fzero[:, i])):
#             print fzero[t, i]
#             print np.min(fbar[max(t-ntau2, 0):t][:, i])
            fzero[t, i] = np.min(fbar[np.clip(t-ntau2, 0):t][:, i])
    # fix the first ntau2 time bins
    fzero[:ntau2] = fbar[:ntau2]
    R = (traces - fzero)/fzero
    ntau0 = tau0 * fps
    w_tau = np.exp(-np.arange(-20, 20)/fps/tau0)
    w_tau[:20] = 0
    for i in xrange(traces.shape[1]):
        R[:, i] = np.convolve(R[:, i], w_tau/np.sum(w_tau), mode='same')
    return R

def calc_dFoF_sima(traces, t1=0.1, t2=60, convzero=True):
    traces_smooth = traces.copy()
    uni_win = np.ones(int(t1*fps))
    uni_win /= uni_win.sum()
    for i in xrange(traces.shape[1]):
        traces_smooth[:, i] = np.convolve(traces[:, i], uni_win, mode='same')
    fzero = traces_smooth.copy()
    nt2 = int(t2 * fps)
    for i in xrange(traces.shape[1]):
        for t in xrange(1, len(traces[:, i])):
#             print min(traces_smooth[:, i][max((t-nt2, 0)):t])
            if t-nt2 < 0:
                minzero = np.min(traces_smooth[:, i][nt2:2*nt2])
                fzero[t, i] = minzero
            else:
                fzero[t, i] = np.min(traces_smooth[:, i][t-nt2:t])
#             print np.min(traces_smooth[:, i][max((t-nt2, 0)):t])
    if convzero:
        for i in xrange(traces.shape[1]):
            fzero[:, i] = np.convolve(fzero[:, i], uni_win, mode='same')
    return (traces_smooth - fzero)/fzero, fzero, traces_smooth

def calc_dFoF_dombeck(traces, time_ax, percentile=8, time_window=15):
    fzero = []
    traces_ = traces.copy()
    fzeros = np.zeros_like(traces)
    for cell in xrange(traces.shape[1]):
        fzero = []
        for t0 in time_ax:
            t_start = t0 - time_window/2.
            t_stop = t_start + time_window
            time_bool = (time_ax>=t_start) * (time_ax<t_stop)
            fzero.append(np.percentile(traces[:, cell][time_bool], percentile))
        fzeros[:, cell] = fzero
    return traces_-fzeros, fzeros, traces

def extract_candidate_regions(img, opening_square_size=3, closing_square_size=3, min_area=150, thr=0.3):
    bw = opening(img > thr, square(opening_square_size))
    bw = closing(bw, square(closing_square_size))
    label_image = label(bw)
    regions = []
    for region in regionprops(label_image):
        if region.area > min_area:
            regions.append(region)
    return regions

def create_imgroi(img, region):
    imgroi = np.zeros(img.shape)
    for r in region.coords:
        imgroi[r[0], r[1]] = 1
    return imgroi

def erode_region(img, region, erosion_disk_size):
    imgroi = create_imgroi(img, region)
    eds = max(1, erosion_disk_size)
    regions = regionprops(label(erosion(imgroi, disk(eds))))
    return regions[0] if len(regions)>0 else []

def extract_roi(img, region):
    imgroi = create_imgroi(img, region)
    return np.fliplr(find_contours(imgroi, 0.5)[0])

def gen_candidates(rois_matlab, default_parameters, parameters={}, auto_scale_regions=False,
                   min_area=200, max_area=500):
    candidates = []
    regions = []
    ids = []
    erosion_sizes = []
    if auto_scale_regions:
        print "Scaling regions to max_area..."
    for i, img in enumerate(rois_matlab):
        if i in parameters.keys():
            pars = parameters[i]
        else:
            pars = default_parameters
        regs = extract_candidate_regions(img, thr=pars[0])
        for r in regs:
            erosion_size = pars[1]
            r_e = erode_region(img, r, erosion_size)
            if not hasattr(r_e, 'coords') or r_e.area < min_area:
                continue
            if auto_scale_regions:
                while r_e.area > max_area:
                    r_e = erode_region(img, r, erosion_size)
                    erosion_size += .1
                if np.random.random() >= 0.5:
                    print "boom",
                else:
                    print "baboooom",
                print r_e.area,
        #         erosion_sizes.append(erosion_size)
        
            if hasattr(r_e, 'coords'):
                candidates.append(extract_roi(img, r_e))
                regions.append(r)
                ids.append(i)
    return regions, candidates, ids


# In[93]:

## THIS IS THE MAIN FOLDER
#You should only need to change root_folder
root_folder = "/home/sebi/data/sebi/experiment_name/mouse_number/date/session"
xmlfile = os.path.join(root_folder, 'tseries/raw/tseries.xml')


# In[94]:

filenames = 'tseries/normcorre/mc.tif'


# In[95]:

#rois_matlab_all = loadmat(os.path.join(root_folder, 'tseries/turboreg/ica_filters.mat'))['ica_filters']


# # From here on, Shift+Enter only

# In[96]:

# these are loaded as different sequences but it should be 1 big sequence of N frames
# sequences = [sima.Sequence.create('TIFF', filename) for filename in glob.glob(root_folder+'/TS*.tif')]
sequences = [sima.Sequence.create('TIFF', filename) for filename in glob.glob(os.path.join(root_folder, filenames))]


# In[97]:

# preparing folders

traces_output_folder = os.path.join(root_folder, "traces/")
if not os.path.isdir(traces_output_folder):
    os.makedirs(traces_output_folder)

sima_folder = os.path.join(root_folder, 'traces/sima_files')
try:
    dataset = sima.ImagingDataset(sequences, sima_folder)
except OSError, e:
    print e
    sima_folder += time.strftime("__%h%d_%H%M%S")
    print "Writing on %s instead" % sima_folder
    dataset = sima.ImagingDataset(sequences, sima_folder)


# In[98]:

# def plotme(cell=0, thr=0.5, threshold=False):
#     if not threshold:
#         pl.imshow(rois_matlab_all[cell], cmap=pl.cm.gray_r)
#     else:
#         pl.imshow(rois_matlab_all[cell]>thr, cmap=pl.cm.gray_r)
# interact(plotme, cell=(0, rois_matlab_all.shape[0]-1, 1), thr=(0, 2, 0.1))


# # Write the list of keep and delete

# In[99]:

# # write in only one of these two
# exclude_these = []
# keep_these = [13, 17, 19, 22, 23, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 45, 46, 47, 50, 52, 54,
#               57, 58, 61, 64, 67, 74]


# In[100]:

# if len(exclude_these) == 0:
#     exclude_these = np.delete(np.arange(rois_matlab_all.shape[0]), keep_these)
# else:
#     keep_these = np.delete(np.arange(rois_matlab_all.shape[0]), exclude_these)


# In[101]:

# exclude_these


# In[102]:

# rois_matlab = rois_matlab_all[keep_these]


# # Play with dffault parameters here

# In[103]:

# # these are applied to all cells
# default_parameters = [0.1, 1]  # [threshold, erosion]


# In[104]:

average = dataset.time_averages[0, :, :, 0]


# In[105]:

# def plot_me(threshold=default_parameters[0], erosion=default_parameters[1], mode='averageroi', smooth=True):
#     regions, candidates, ids = gen_candidates(rois_matlab, default_parameters=[threshold, erosion])
#     fig, ax = pl.subplots(1, 1, figsize=(5, 5))
#     colors = pl.cm.rainbow(np.linspace(0, 1, len(candidates)))
#     [pl.plot(c[:, 0], c[:, 1], color=col) for c, col in zip(candidates, colors)]
#     [pl.text(r.centroid[1], r.centroid[0], i, color=col) for i, r, col in zip(ids, regions, colors)]
#     ax.set_xticks(())
#     ax.set_yticks(())
#     if mode == 'averageroi':
#         if smooth:
#             ax.imshow(gaussian(rois_matlab.mean(0), 5), alpha=1, cmap=pl.cm.gray)
#         else:
#             ax.imshow(rois_matlab.mean(0), alpha=1, cmap=pl.cm.gray)
#     if mode == 'stdroi':
#         if smooth:
#             ax.imshow(gaussian(np.std(rois_matlab, 0), 5), alpha=1, cmap=pl.cm.gray)
#         else:
#             ax.imshow(np.std(rois_matlab, 0), alpha=1, cmap=pl.cm.gray)
#     elif mode == 'average':
#         ax.imshow(average, alpha=1, cmap=pl.cm.gray)
# interact(plot_me, threshold=(0.01, 3, 0.01), erosion=(0., 5, 0.1), mode=['average', 'averageroi', 'stdroi'])


# # Write the new dffault parameters here

# In[106]:

# # these are applied to all cells
# default_parameters = [0.07, 3.6]  # [threshold, erosion]


# In[107]:

# def plot_me(image=0, threshold=default_parameters[0], erosion=default_parameters[1], auto_erode=False):

#     fig, ax = pl.subplots(1, 1)
    
#     regions, candidates, ids = gen_candidates(rois_matlab, parameters=parameters, default_parameters=default_parameters,
#                                      auto_scale_regions=True, min_area=min_area, max_area=max_area)
    
#     img = rois_matlab[image]
#     pl.imshow(img, cmap=pl.cm.gray_r, alpha=0.5)
#     pl.ylim(365, 0)
#     pl.xlim(0, 456)
#     regs = extract_candidate_regions(img, thr=threshold)
#     for region in regs:
# #         roi = extract_roi(img, region)
# #         pl.plot(roi[:, 0], roi[:, 1], 'k--') 
#         if erosion >= 1:
#             region_e = erode_region(img, region, erosion)
#         else:
#             region_e = region
#         try:
#             roi = extract_roi(img, region_e)
#             pl.plot(roi[:, 0], roi[:, 1]) 
#         except AttributeError:
#             print "destroyed..."
#         pl.text(region.centroid[1]+20, region.centroid[0]-20, region_e.area, color='blue')
#     #             poly_coords = approximate_polygon(region.coords, tolerance)
#     #             ax.add_patch(Polygon(np.fliplr(poly_coords), color='r', alpha=0.2))
#     #             pl.plot(region.perimeter[:, 0], region.perimeter[:, 1])
    
# interact(plot_me, image=(0, rois_matlab.shape[0]-1, 1), threshold=(0, 2, 0.01), erosion=(0, 20, .1))


# In[108]:

# # but you can still modify the parameters to some cells
# parameters = {16: [0.1, 3.6],  # image : [threshold, erosion]
#               19: [0.15, 3.1],
#               28: [0.03, 2.1]
#              }


# # Save the ROIs

# In[109]:

# max_area = 200
# min_area = 10
# mode = 'averageroi'
# smooth = True
# auto_scale_regions = True  # set to false to avoid automatic erosion to max area

# regions, candidates, ids = gen_candidates(rois_matlab, parameters=parameters, default_parameters=default_parameters,
#                                         auto_scale_regions=auto_scale_regions, min_area=min_area, max_area=max_area)

# fig, ax = pl.subplots(1, 1, figsize=(5, 5))
# colors = pl.cm.rainbow(np.linspace(0, 1, len(candidates)))
# [pl.plot(c[:, 0], c[:, 1], color=col) for c, col in zip(candidates, colors)]
# [pl.text(r.centroid[1], r.centroid[0], i, color=col) for i, r, col in zip(ids, regions, colors)]
# ax.set_xticks(())
# ax.set_yticks(())
# if mode == 'averageroi':
#     if smooth:
#         ax.imshow(gaussian(rois_matlab.mean(0), 5), alpha=1, cmap=pl.cm.gray)
#     else:
#         ax.imshow(rois_matlab.mean(0), alpha=1, cmap=pl.cm.gray)
# if mode == 'stdroi':
#     if smooth:
#         ax.imshow(gaussian(np.std(rois_matlab, 0), 5), alpha=1, cmap=pl.cm.gray)
#     else:
#         ax.imshow(np.std(rois_matlab, 0), alpha=1, cmap=pl.cm.gray)
# elif mode == 'average':
#     ax.imshow(average, alpha=1, cmap=pl.cm.gray)
    
# myrois_folder = os.path.join(root_folder, 'traces/myrois')
# if not os.path.isdir(myrois_folder):
#     os.makedirs(myrois_folder)
# os.system('rm %s/*'%myrois_folder)
# [np.savetxt(os.path.join(myrois_folder, '%04d-%04d.roi') %
#            (region.centroid[0], region.centroid[1]), c, fmt='%d')
#  for region, c in zip(regions, candidates)];


# # WAIT! Go to ImageJ and save your ROIs from there
# 
# Open ImageJ and open the video. Then `Plugins->Macros->ROI_select`.
# 
# Edit the ROIs and once satisfied select all ROIs then, in the ROI manager, `Properties...` and set `Position` to `none` to all of them.
# Then also in the ROI manager window hit `More...` then `Save` and save the ROIs as `traces/RoiSet.zip` (case sensitive).

# In[110]:

# reading ROIs

from sima.ROI import ROIList
rois_filename = os.path.join(root_folder, 'traces/RoiSet.zip')

rois = ROIList.load(rois_filename, fmt='ImageJ')
# rois = ROIList.load('/home/fabios/data/testingmanualroi/jbiane_test/sess01/tseries/RoiSet_29_122116.zip', fmt='ImageJ')

# for some reason, the ROI im_shape attributes may be screwed up so here we override them
for r in rois:
    r.im_shape = dataset.frame_shape[1:-1]
dataset.add_ROIs(rois, 'from_ImageJ')


# In[111]:

fig, ax = pl.subplots(1, 1, figsize=(5, 5))
colors = pl.cm.rainbow(np.linspace(0, 1, len(rois)))
for i, r in enumerate(rois):
    mask = r.mask[0]
    ax.imshow(np.where(mask.todense()==True, i, np.nan),
              cmap=pl.cm.rainbow, alpha=0.5, zorder=19, vmin=0, vmax=len(rois))
    x, y = max(r.coords[0][:, 0]), max(r.coords[0][:, 1])
    ax.text(x, y, i+1, color='y')
    ax.set_xticks(())
    ax.set_yticks(())
ax.imshow(average, alpha=1, cmap=pl.cm.gray)


# In[112]:

# signals_auto = dataset.extract(rois_auto)


# In[113]:

# np.r_[signals_auto['raw']].shape


# In[114]:

# grab time axis from the xml file

import xml.etree.ElementTree as ET

print "I infer the time axis from:\n", xmlfile
tree = ET.parse(xmlfile)
root = tree.getroot()

# unfortunately we miss the first frame
time_ax = np.r_[[child.attrib['absoluteTime']
                 for child in root.iter('Frame')]].astype(float)


# In[115]:

dataset.num_frames


# In[116]:

fps = 1./np.diff(time_ax)[0]


# In[117]:

signals = dataset.extract(rois)


# In[118]:

nseq, nrois, nframes = np.r_[signals['raw']].shape
print "We have %d sequences with %d ROIs and %d time bins." % (nseq, nrois, nframes)


# In[119]:

# these are df/f signals, we take only the first sequence!

sequence_id = 0
traces = np.r_[signals['raw']][sequence_id, :, :].T
df_f = (traces-traces.mean(0))/traces.mean(0)


# In[120]:

# trim the time_ax
time_ax = time_ax[:traces.shape[0]]


# In[121]:

time_ax.shape


# In[122]:

traces[:, 0]


# In[123]:

df_f, fzero, traces = calc_dFoF_dombeck(np.nan_to_num(traces), time_ax, time_window=30)


# In[124]:

def plot_me(cell=0, tstart=0, deltat=120):

    pl.plot(time_ax[:traces.shape[0]], traces[:, cell], 'k-', lw=2, alpha=0.3)
    pl.plot(time_ax[:traces.shape[0]], df_f[:, cell], 'r-', lw=2)
    pl.plot(time_ax[:traces.shape[0]], fzero[:, cell], 'b-', lw=2)

    pl.xlim(tstart, tstart+deltat)
    pl.ylim(-0.5, 3)
    pl.xlabel("time (s)")
    pl.ylabel("df/f")
interact(plot_me, cell=(0, len(rois)-1, 1), tstart=(0, 3000, 1), deltat=(0, 500, 0.1))


# In[125]:

colors = pl.cm.rainbow(np.linspace(0, 1, len(rois)))
# print "********* WE ARE SKIPPING LAST FRAME FOR NOW, NICK KNOWS WHY! *********"
# [pl.plot(time_ax[:traces.shape[0]]/60., (traces[:, i]-traces[:, i].mean())+i, color=c, alpha=0.3) for i, c in enumerate(colors)];
[pl.plot(time_ax[:traces.shape[0]]/60., df_f[:, i]+i, color=c) for i, c in enumerate(colors)];
pl.xlabel('Time (m)')
pl.ylabel("df/f from SIMA")
# pl.xlim(0, 1)


# # Save ROIs

# In[126]:

if not os.path.isdir(traces_output_folder):
    os.mkdir(traces_output_folder)
    print "Created folder", traces_output_folder
np.savetxt(os.path.join(traces_output_folder, "sima_manual_ROIs_traces.txt"), traces)
np.savetxt(os.path.join(traces_output_folder, "sima_manual_ROIs_df_f.txt"), df_f)
np.savetxt(os.path.join(traces_output_folder, "sima_manual_ROIs_time_ax.txt"), time_ax)


# In[ ]:



