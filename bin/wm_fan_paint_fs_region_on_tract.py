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
    description="Extract region (e.g. freesurfer) of each fiber point and paint the region label on the fiber.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'input',
    help='Directory or path of input VTK/VTP file(s) that are going to be painted.')
parser.add_argument(
    'label_map_file',
    help='Label map file in nifti (default for freesurfer result).')
parser.add_argument(
    'outputDirectory',
    help='Directory of output transformed results.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.exists(args.input):
    print "Error: Input directory", args.input, "does not exist."
    exit()


pd_tract_list = []
if os.path.isfile(args.input):
    flag_single = True
    pd_tract_list.append(args.input)
elif os.path.isdir(args.input):
    flag_single = False
    def list_cluster_files(input_dir):
        # Find input files
        input_mask = "{0}/cluster_*.vtk".format(input_dir)
        input_mask2 = "{0}/cluster_*.vtp".format(input_dir)
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_pd_fnames = sorted(input_pd_fnames)
        return (input_pd_fnames)

    pd_tract_list = list_cluster_files(args.input)
    if len(pd_tract_list) == 0:
        print "Error: No cluster files found."
        exit()

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

if not os.path.exists(args.label_map_file):
    print "Label map", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

lb = nibabel.load(args.label_map_file)
voxel_data = lb.get_data()
print args.label_map_file, ', label volume shape: ', lb.get_data().shape

number_of_endpoints = 5

for pd_tract_path in pd_tract_list:
    pd_tract = wma.io.read_polydata(pd_tract_path)

    num_fibers = pd_tract.GetNumberOfLines()
    print pd_tract_path, ', fiber number:', num_fibers

    pd_tract.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()
    inpoints = pd_tract.GetPoints()
    num_points = inpoints.GetNumberOfPoints()

    label_list = numpy.zeros(num_points)
    if flag_single:
        t_start = time.time()

    mask_touch_0_fibers = numpy.zeros(num_fibers)
    for lidx in range(0, num_fibers):

        pd_tract.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        ptidx_list = numpy.zeros(line_length)
        line_labels = numpy.zeros(line_length)
        for pidx in range(0, line_length):
            ptidx = line_ptids.GetId(pidx)
            ptidx_list[pidx] = ptidx
            if pidx < number_of_endpoints or pidx > line_length - number_of_endpoints - 1:
                point_ras = inpoints.GetPoint(ptidx)
                point_ijk = apply_affine(numpy.linalg.inv(lb.affine), point_ras)
                point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
                label = voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])]
                line_labels[pidx] = label

        # if endpoint touches 0 region
        if line_labels[0] < 1:
            for t_idx in range(0,number_of_endpoints):
                if line_labels[t_idx] > 0:
                    line_labels[0:t_idx] = line_labels[t_idx]
                    break

        if line_labels[-1] < 1:
            for t_idx in range(0,number_of_endpoints):
                if line_labels[-t_idx-1] > 0:
                    line_labels[-t_idx:] = line_labels[-t_idx-1]
                    break

        if line_labels[0] < 1 or line_labels[-1] < 1:
            mask_touch_0_fibers[lidx] = 1

        for pidx in range(0, line_length):
            label = line_labels[pidx]
            ptidx = ptidx_list[pidx]
            label_list[ptidx] = label

        step = 1000
        if lidx % step == 0:
            print 'Fiber %8d / %d, length, %4d, endpoint1: %5d, endpoint2: %5d' % (lidx, num_fibers, line_length, line_labels[0], line_labels[-1])
            if flag_single:
                t_end = time.time()
                print '   - time elapsed:', t_end - t_start, ', estimated total time:', num_fibers / step * (t_end - t_start)
                t_start = t_end

    label_array = vtk.vtkIntArray()
    label_array.SetName('region_label')
    label_mask = vtk.vtkIntArray()
    label_mask.SetName('region_mask')
    for val in label_list:
        label_array.InsertNextValue(val)
        if val > 0:
            label_mask.InsertNextValue(1)
        else:
            label_mask.InsertNextValue(0)

    inpointsdata = pd_tract.GetPointData()
    inpointsdata.AddArray(label_array)
    inpointsdata.AddArray(label_mask)
    inpointsdata.Update()
    pd_tract.Update()

    vtk_name = os.path.split(pd_tract_path)[1]
    vtk_name = pd_tract_path[:-4] + '_with_region' + pd_tract_path[-4:]
    wma.io.write_polydata(pd_tract, os.path.join(args.outputDirectory, vtk_name))

    pd_tract_not_touch = wma.filter.mask(pd_tract, mask_touch_0_fibers, color=None, preserve_point_data=True, preserve_cell_data=True, verbose=True)
    vtk_name = pd_tract_path[:-4] + '_with_region_not_touch' + pd_tract_path[-4:]
    wma.io.write_polydata(pd_tract_not_touch, os.path.join(args.outputDirectory, vtk_name))
