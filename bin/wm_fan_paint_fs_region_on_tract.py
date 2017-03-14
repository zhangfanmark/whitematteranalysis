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
    'outputDirectory',
    help='Directory of output transformed results.')
parser.add_argument(
    '-lm', dest="label_map_file",
    help='Label map file in nifti (default for freesurfer result).')
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

for pd_tract_path in pd_tract_list:
    pd_tract = wma.io.read_polydata(pd_tract_path)

    num_fibers = pd_tract.GetNumberOfLines()
    print pd_tract_path, ', fiber number:', num_fibers

    pd_tract.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()
    inpoints = pd_tract.GetPoints()
    inpointsdata = pd_tract.GetPointData()

    label_array = vtk.vtkIntArray()
    label_array.SetNumberOfComponents(inpoints.GetNumberOfPoints())
    label_array.SetName('region_label')

    if flag_single:
        t_start = time.time()
    for lidx in range(0, num_fibers):

        pd_tract.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        fiber_array_r = numpy.zeros(line_length)
        fiber_array_a = numpy.zeros(line_length)
        fiber_array_s = numpy.zeros(line_length)
        ptidx_list = []
        for pidx in range(0, line_length):
            ptidx = line_ptids.GetId(pidx)
            ptidx_list.append(ptidx)
            point = inpoints.GetPoint(ptidx)

            fiber_array_r[pidx] = point[0]
            fiber_array_a[pidx] = point[1]
            fiber_array_s[pidx] = point[2]

        array_ras = numpy.column_stack([fiber_array_r, fiber_array_a, fiber_array_s])
        array_ijk = apply_affine(numpy.linalg.inv(lb.affine), array_ras)
        array_ijk = numpy.rint(array_ijk).astype(numpy.int32)

        line_labels = []
        for pidx in range(0, line_length):
            label = voxel_data[(array_ijk[pidx, 0], array_ijk[pidx, 1], array_ijk[pidx, 2])]
            line_labels.append(label)
        if line_labels[0] < 1:
            line_labels[0] = line_labels[1]
        if line_labels[-1] < 1:
            line_labels[-1] = line_labels[-2]

        for pidx in range(0, line_length):
            label = line_labels[pidx]
            ptidx = ptidx_list[pidx]
            label_array.InsertComponent(ptidx, 0, label)

        if lidx % 20 == 0:
            print 'Fiber %8d, length, %4d, endpoint1: %5d, endpoint2: %5d' % (lidx, line_length, line_labels[0], line_labels[-1])
            if flag_single:
                t_end = time.time()
                print ' - time elapsed:', t_end - t_start
                t_start = t_end


    inpointsdata.AddArray(label_array)
    inpointsdata.Update()
    pd_tract.Update()

    vtk_name = os.path.split(pd_tract_path)[1]
    vtk_name = pd_tract_path[:-4] + '_with_region' + pd_tract_path[-4:]
    wma.io.write_polydata(pd_tract, os.path.join(args.outputDirectory, vtk_name))