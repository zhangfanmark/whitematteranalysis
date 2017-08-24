#!/usr/bin/env python
# If the cluster was created in bilateral clustering and is not commissural, separate into separate directories
# output should be a right_hem and a left_hem directory
# copy MRML files from toplevel directory into this one
# perhaps make a commissural directory also. why not.
# note that it could be cleaner to have a cell data array in the polydata and the option to 
# measure left and right hem parts of a cluster separately.
# Note: if this is performed at the atlas clustering stage, it can be used to separate clusters into groups,
# and this can be learned. At that point all data are midsagitally aligned, which this requires.
# For running this per-subject, the alignment should be performed to handle the tracts near the midline better.
# That should be added as an option.
import numpy
import argparse
import os
import vtk
import glob

try:
    import whitematteranalysis as wma
except:
    print "<wm_access_cluster_location.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Access the location of each atlas fiber cluster, which can be used to guide the cluster separation in individual subjects.",
    epilog="Written by Fan Zhang")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='A directory of clustered whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='output temp files')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "<wm_access_cluster_location.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory."
    exit()


print "<wm_access_cluster_location.py> Starting computation."
print ""
print "=====input directory ======\n", args.inputDirectory
print "=========================="
print ""

# relatively high number of points for accuracy
points_per_fiber = 40

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/cluster_*.vtk".format(input_dir)
    input_mask2 = "{0}/cluster_*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

input_polydatas = list_cluster_files(args.inputDirectory)

number_of_subjects = len(input_polydatas)

print "<wm_access_cluster_location.py> Input number of vtk/vtp files: ", number_of_subjects


def get_label_array(input_data, verbose=True):

    input_data.GetLines().InitTraversal()
    num_fibers = input_data.GetNumberOfLines()
    inpoints = input_data.GetPoints()
    inpointsdata = input_data.GetPointData()
    num_points = inpoints.GetNumberOfPoints()

    if verbose:
        print 'Total fiber number in the cluster:', num_fibers, ', point number in the cluster:', num_points

    if inpointsdata.GetNumberOfArrays() > 0:
        point_data_array_indices = range(inpointsdata.GetNumberOfArrays())
        for idx in point_data_array_indices:
            array = inpointsdata.GetArray(idx)
            if array.GetName().lower() == 'region_label':
                label_array = array

    return label_array

def get_along_tract_region(input_data, verbose=True):

    label_array = get_label_array(input_data, verbose=False)

    num_fibers = input_data.GetNumberOfLines()

    # loop over lines
    input_data.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()

    along_tract_regions = [] # this is unique region: set
    along_tract_regions_all = [] # all regions all points.
    for lidx in range(0, num_fibers):

        input_data.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        # loop over the indices that we want and get those points
        region_one_fiber = []
        for line_index in range(line_length):
            # do nearest neighbor interpolation: round index
            ptidx = line_ptids.GetId(int(round(line_index)))
            array_label = label_array.GetTuple(ptidx)[0]

            region_one_fiber.append(array_label)

        region_one_fiber = numpy.array(region_one_fiber)
        along_tract_regions_all.append(region_one_fiber[numpy.where(region_one_fiber>0)])

        region_one_fiber = numpy.unique(region_one_fiber)
        region_one_fiber = region_one_fiber[numpy.where(region_one_fiber>0)]

        if verbose:
            print ' --Fiber', lidx, 'regions:', region_one_fiber

        along_tract_regions.append(region_one_fiber)

    return along_tract_regions, along_tract_regions_all

def remove_BS_points(inpd):

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()

    array = inpointdata.GetArray(0)

    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outpoints = vtk.vtkPoints()

    inpd.GetLines().InitTraversal()
    outlines.InitTraversal()

    for lidx in range(0, inpd.GetNumberOfLines()):
        ptids = vtk.vtkIdList()
        inpd.GetLines().GetNextCell(ptids)

        cellptids = vtk.vtkIdList()

        num_points_before = ptids.GetNumberOfIds()
        num_points_kept = 0
        for pidx in range(0, ptids.GetNumberOfIds()):
            point = inpoints.GetPoint(ptids.GetId(pidx))

            fs_label = array.GetTuple(ptids.GetId(pidx))

            if fs_label[0] != 16:
                idx_ = outpoints.InsertNextPoint(point)
                cellptids.InsertNextId(idx_)
                num_points_kept += 1

        if num_points_kept > 0:
            outlines.InsertNextCell(cellptids)

        # print 'Line:', lidx, ', before removal:', num_points_before, ', after removal:', num_points_kept,
        # print ', point removed:', num_points_before - num_points_kept

    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)

    return outpd

def decide_location_OLD(pd, verbose=False):

    hemisphere_percent_threshold = 60

    fibers = wma.fibers.FiberArray()
    fibers.points_per_fiber = points_per_fiber
    fibers.hemisphere_percent_threshold = hemisphere_percent_threshold / 100.0
    fibers.hemispheres = True
    fibers.convert_from_polydata(pd)

    num_fibers_h_60 = len(fibers.index_left_hem) + len(fibers.index_right_hem)
    num_fibers_c_60 = len(fibers.index_commissure)

    hemisphere_percent_threshold = 90

    fibers = wma.fibers.FiberArray()
    fibers.points_per_fiber = points_per_fiber
    fibers.hemisphere_percent_threshold = hemisphere_percent_threshold / 100.0
    fibers.hemispheres = True
    fibers.convert_from_polydata(pd)

    num_fibers_h_90 = len(fibers.index_left_hem) + len(fibers.index_right_hem)
    num_fibers_c_90 = len(fibers.index_commissure)

    # if a cluster was too small, we need to check.
    if num_fibers_h_90 + num_fibers_c_90 < 10:
        decision_60 = 'ng'
        decision_90 = 'ng'
    else:
        if num_fibers_h_60 > 4 * num_fibers_c_60 and (num_fibers_c_60 < 20 or num_fibers_h_60 > 8 * num_fibers_c_60):
            decision_60 = 'h'
        elif num_fibers_c_60 > 4 * num_fibers_h_60 and (num_fibers_h_60 < 20 or num_fibers_c_60 > 8 * num_fibers_h_60):
            decision_60 = 'c'
        else:
            decision_60 = 'ng'

        if num_fibers_h_90 > 4 * num_fibers_c_90 and (num_fibers_c_90 < 20 or num_fibers_h_90 > 8 * num_fibers_c_90):
            decision_90 = 'h'
        elif num_fibers_c_90 > 4 * num_fibers_h_90 and (num_fibers_h_90 < 20 or num_fibers_c_90 > 8 * num_fibers_h_90):
            decision_90 = 'c'
        else:
            decision_90 = 'ng'

    if verbose:
        print ' - th = 60: %5d:%5d - %3s' % (num_fibers_h_60, num_fibers_c_60, decision_60)
        print ' - th = 90: %5d:%5d - %3s' % (num_fibers_h_90, num_fibers_c_90, decision_90)


    if decision_60 == 'c' and decision_90 == 'c':
        # both settings have comm
        decision_final = 'c'
    elif decision_60 == 'h' and decision_90 == 'h':
        # both settings have hemi
        decision_final = 'h'
    elif decision_90 == 'ng' and decision_60 == 'h':
        # If with th=90, we are not sure what the cluster is;
        # Then, we see the lower th=60, if it is hemi, we consider it is a hemi.
        decision_final = 'ng_h'
    elif decision_60 == 'ng' and decision_90 == 'c':
        # If with th=60, we are not sure what the cluster is;
        # Then, we see the higher th=90, if it is comm, we consider it is a comm.
        decision_final = 'ng_c'
    else:
        # Otherwise, it is a Not Given
        decision_final = 'ng'

    return decision_final


def decide_location(pd, verbose=False):

    hemisphere_percent_threshold = 80

    fibers = wma.fibers.FiberArray()
    fibers.points_per_fiber = points_per_fiber
    fibers.hemisphere_percent_threshold = hemisphere_percent_threshold / 100.0
    fibers.hemispheres = True
    fibers.convert_from_polydata(pd)

    num_fibers_h = len(fibers.index_left_hem) + len(fibers.index_right_hem)
    num_fibers_c = len(fibers.index_commissure)


    if num_fibers_h > num_fibers_c:
        location_decision = 'h'
    elif num_fibers_c > num_fibers_h:
        location_decision = 'c'
    else:
        location_decision = 'ng'

    if verbose:
        print ' - hemisphere percent threshold = %3d: hemispheric %5d - commissural %5d - %3s' % (hemisphere_percent_threshold, num_fibers_h, num_fibers_c, location_decision)

    return location_decision

def output_mrml(cluster_list, mrml_filename):
    number_of_files = len(cluster_list)
    if number_of_files > 0:
        if number_of_files > 1:
            step = int(100 * 255.0 / (number_of_files - 1))
        elif number_of_files == 1:
            step = int(100 * 255.0)

        R = numpy.array(range(0, 100 * 255 + 1, step)) / 100.0
        G = numpy.abs(range(100 * -127, 100 * 128 + 1, step)) * 2.0 / 100.0
        B = numpy.array(range(100 * 255 + 1, 0, -step)) / 100.0

        colors = list()
        for idx in range(len(cluster_list)):
            colors.append([R[idx], G[idx], B[idx]])

        colors = numpy.array(colors)
        wma.mrml.write(cluster_list, colors, mrml_filename, ratio=1.0)

cluster_location_list = []
BS_related_list = []

outstr = 'Cluster Index' + '\t'+ 'Location Label' + '\n'

for fname in input_polydatas:

    # figure out filename and extension
    fname_base = os.path.basename(fname)
    print "\n***", fname_base

    # read data
    pd = wma.io.read_polydata(fname)

    if pd.GetNumberOfLines() == 0:
        cluster_location_list.append('ng')
        BS_related_list.append('No')
        outstr = outstr + fname_base + '\t' + 'ng' + '\n'

        continue

    # get region label
    along_tract_regions, along_tract_regions_all = get_along_tract_region(pd, verbose=False)

    # decide if a cluster is just from BS
    along_tract_regions_all_all = numpy.concatenate(along_tract_regions_all).astype(int)
    num_points_in_BS = numpy.where(along_tract_regions_all_all == 16)[0].shape[0]

    if num_points_in_BS > along_tract_regions_all_all.shape[0]*3/4:
        print ' - ', num_points_in_BS, 'out of ', along_tract_regions_all_all.shape[0], 'points are in the brain stem.'
        decision_final = 'c'
        BS_related_list.append('Yes')
    else:
        # decide if a cluster passing the brain stem (BS)
        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)
        unique_regions = numpy.unique(all_regions_all_fibers)
        region_occurrence = numpy.bincount(all_regions_all_fibers)

        if numpy.where(unique_regions == 16)[0].shape[0] > 0 and region_occurrence[16] > pd.GetNumberOfLines() / 4:
            print ' - ', region_occurrence[16], ' out of ',  pd.GetNumberOfLines(), 'fibers passing through the brain stem.'
            # remove the BS points
            pd = remove_BS_points(pd)
            BS_related_list.append('Yes')
        else:
            BS_related_list.append('No')

        decision_final = decide_location(pd, verbose=True)

    cluster_location_list.append(decision_final)

    outstr = outstr + fname_base + '\t' + decision_final + '\n'

    print '\n == hemisphere location:', decision_final

output_file = open(os.path.join(args.inputDirectory, 'cluster_hemisphere_location.txt'), 'w')
output_file.write(outstr)
output_file.close()


if args.flag_verbose:

    yes_h_list = []
    yes_c_list = []
    yes_ng_list = []

    no_h_list = []
    no_c_list = []
    no_ng_list = []

    outdir = os.path.join(args.inputDirectory, 'tmp_cluster_location')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for fname, BS_yes_no, decision_label in zip(input_polydatas, BS_related_list, cluster_location_list) :

        # figure out filename and extension
        fname_base = os.path.basename(fname)
        print "\n***", fname_base

        if BS_yes_no == 'Yes':
            if decision_label == 'h':
                yes_h_list.append(fname_base)
            elif decision_label == 'c':
                yes_c_list.append(fname_base)
            elif decision_label == 'ng':
                yes_ng_list.append(fname_base)

        elif BS_yes_no == 'No':
            if decision_label == 'h':
                no_h_list.append(fname_base)
            elif decision_label == 'c':
                no_c_list.append(fname_base)
            elif decision_label == 'ng':
                no_ng_list.append(fname_base)


    output_mrml(yes_h_list, os.path.join(outdir, "m_yes_h_list"+str(len(yes_h_list))+".mrml"))
    output_mrml(yes_c_list, os.path.join(outdir, "m_yes_c_list"+str(len(yes_c_list))+".mrml"))
    output_mrml(yes_ng_list, os.path.join(outdir, "m_yes_ng_list"+str(len(yes_ng_list))+".mrml"))

    output_mrml(no_h_list, os.path.join(outdir, "m_no_h_list"+str(len(no_h_list))+".mrml"))
    output_mrml(no_c_list, os.path.join(outdir, "m_no_c_list"+str(len(no_c_list))+".mrml"))
    output_mrml(no_ng_list, os.path.join(outdir, "m_no_ng_list"+str(len(no_ng_list))+".mrml"))

