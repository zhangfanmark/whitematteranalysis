#!/usr/bin/env python
import argparse
import os
import glob
import nibabel
import vtk
import numpy
import time
from nibabel.affines import apply_affine
from joblib import Parallel, delayed

try:
    import whitematteranalysis as wma
except:
    print "<wm_harden_transform_with_slicer> Error importing white matter analysis package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Endpoint analysis using the labels already stored in the fiber cluster vtk/vtp files.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'inputDirectory',
    help='Directory of the subject-specific fiber clustering results.')
parser.add_argument(
    'regionList', type=int, nargs='+',
    help='A list of regions to be analyzed, such as, 8,1024')
parser.add_argument(
    '-hemi', type=str, dest="hemi", default="whole",
    help='Hemispheres or whole brain fiber clusters to be analyzed.')
parser.add_argument(
    '-fn', type=str, dest="folder_name", default="",
    help='String in the folder name of the subject-specific fiber clustering')


args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

hemi = ''
if args.hemi == "comm":
    hemi = 'tracts_commissural'
elif args.hemi == "left":
    hemi = 'tracts_left_hemisphere'
elif args.hemi == 'right':
    hemi = 'tracts_right_hemisphere'
elif args.hemi == 'whole':
    hemi = ''
else:
    print "Error: Hemisphere or whole brain: ", args.hemi, "does not exist."
    exit()

folder_name = args.folder_name

def list_fc_folders(input_dir, hemi, folder_name):
    # Find input files
    input_mask = ("{0}/*"+folder_name+"*/"+hemi).format(input_dir)
    input_fc_folders = glob.glob(input_mask)
    input_fc_folders = sorted(input_fc_folders)
    return (input_fc_folders)

fc_folders = list_fc_folders(args.inputDirectory, hemi, folder_name)

num_subjects = len(fc_folders)
num_clusters = 0
flag_first = 0
for fc_folder in fc_folders:
    pd_fc_paths = wma.io.list_vtk_files(fc_folder)
    #print "Subject", fc_folder, 'has', len(pd_fc_paths), 'clusters.'
    if flag_first == 0:
        num_clusters = len(pd_fc_paths)
    else:
        if num_clusters != len(pd_fc_paths):
            print "Error: All subjects should have the same number of clusters"
            exit()

def get_endpoint_labels_per_cluster(pd_fc):

    pd_fc.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()
    inpointdata = pd_fc.GetPointData()

    point_data_array_indices = range(inpointdata.GetNumberOfArrays())
    label_array = None
    for idx in point_data_array_indices:
        label_array = inpointdata.GetArray(idx)
        if label_array.GetName().lower() == 'region_label':
            break

    eps_1 = []
    eps_2 = []
    eps_all = []
    num_fibers = pd_fc.GetNumberOfLines()

    mask = numpy.zeros(num_fibers)
    flag_tmp = 0

    for lidx in range(0, num_fibers):
        pd_fc.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        ep_1 = [0]
        ep_2 = [0]
        label_array.GetTupleValue(line_ptids.GetId(0), ep_1)
        label_array.GetTupleValue(line_ptids.GetId(line_length-1), ep_2)

        ep_1 = ep_1[0]
        ep_2 = ep_2[0]

        if ep_1 > 3000:
            ep_1 = ep_1 - 2000
        if ep_2 > 3000:
            ep_2 = ep_2 - 2000
        #
        # if ep_1 > ep_2:
        #     tmp = ep_2
        #     ep_2 = ep_1
        #     ep_1 = tmp


        # if ep_1 == 1033 or ep_2 == 1033:# and flag_tmp ==0:
        #     mask[lidx] = 1

        # eps_1.append(ep_1)
        # eps_2.append(ep_2)

        eps_all.append(ep_1)
        eps_all.append(ep_2)

    tmp = numpy.argsort(numpy.bincount(eps_all))[-2:]

    try:
        ep_1_label = min(tmp)
        ep_1_label_percent = numpy.sum(numpy.sum(eps_all == ep_1_label)) / float(len(eps_all))
    except:
        ep_1_label = -1
        ep_1_label_percent = -1

    try:
        ep_2_label = max(tmp)
        ep_2_label_percent = numpy.sum(numpy.sum(eps_all == ep_2_label)) / float(len(eps_all))
    except:
        ep_2_label = -1
        ep_2_label_percent = -1

    # print ep_1_label, ep_1_label_percent, ep_2_label, ep_2_label_percent

    # tmp = wma.filter.mask(pd_fc, mask, preserve_point_data=True, preserve_cell_data=True, verbose=False)
    # wma.io.write_polydata(tmp, '/data/lmi/projects/HCP/Tractography/UKF_2T_fs_label_mask/clusters_l40_k2000_f5000/cluster_atlas_01_00002_remove_outliers/tmp.vtk')

    return (ep_1_label, ep_2_label, ep_1_label_percent, ep_2_label_percent)

def get_endpoint_labels_per_subject(fc_folder):
    print fc_folder

    pd_fc_paths = wma.io.list_vtk_files(fc_folder)

    ep_1_labels = []
    ep_2_labels = []
    ep_1_labels_percent = []
    ep_2_labels_percent = []

    count  = 0
    for pd_fc_path in pd_fc_paths:
        # print pd_fc_path
        count = count + 1
        if count == 10:
            break

        pd_fc = wma.io.read_polydata(pd_fc_path)

        ep_1_label, ep_2_label, ep_1_label_percent, ep_2_label_percent = get_endpoint_labels_per_cluster(pd_fc)

        ep_1_labels.append(ep_1_label)
        ep_2_labels.append(ep_2_label)
        ep_1_labels_percent.append(ep_1_label_percent)
        ep_2_labels_percent.append(ep_2_label_percent)

    return numpy.concatenate((ep_1_labels, ep_2_labels, ep_1_labels_percent, ep_2_labels_percent))

print 'A total of', num_subjects, 'subjects'
print 'Eech subject has', num_clusters, 'clusters'

mat_ep_labels = \
    Parallel(n_jobs=10, verbose=1)(
        delayed(get_endpoint_labels_per_subject)(fc_folder)
        for fc_folder in fc_folders)

mat_ep_labels = numpy.array(mat_ep_labels)

len_mat = mat_ep_labels.shape[1] / 4

num_clusters = len_mat

mat_ep_1_labels = mat_ep_labels[:, 0:num_clusters]
mat_ep_2_labels = mat_ep_labels[:, num_clusters:2*num_clusters]
mat_ep_1_labels_percent = mat_ep_labels[:, 2*num_clusters:3*num_clusters]
mat_ep_2_labels_percent = mat_ep_labels[:, 3*num_clusters:4*num_clusters]

for c_idx in range(num_clusters):

    print 'Cluster', c_idx+1
    ep_1_labels_per_cluster = mat_ep_1_labels[:, c_idx]
    ep_2_labels_per_cluster = mat_ep_2_labels[:, c_idx]

    ep_1_labels_per_cluster = ep_1_labels_per_cluster.astype(numpy.int64)
    ep_2_labels_per_cluster = ep_2_labels_per_cluster.astype(numpy.int64)

    ep_1_label_per_cluster = numpy.argmax(numpy.bincount(ep_1_labels_per_cluster))
    ep_2_label_per_cluster = numpy.argmax(numpy.bincount(ep_2_labels_per_cluster))

    ep_1_label_per_cluster_percent = numpy.sum(numpy.sum(ep_1_labels_per_cluster == ep_1_label_per_cluster)) / float(len(ep_1_labels_per_cluster))
    ep_2_label_per_cluster_percent = numpy.sum(numpy.sum(ep_2_labels_per_cluster == ep_2_label_per_cluster)) / float(len(ep_2_labels_per_cluster))

    print ep_1_label_per_cluster, ep_2_label_per_cluster, ep_1_label_per_cluster_percent, ep_2_label_per_cluster_percent