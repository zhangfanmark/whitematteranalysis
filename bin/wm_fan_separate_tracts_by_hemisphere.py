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
import shutil
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_separate_hemispheres.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Separate each cluster into left/right/commissural tracts according to the percentage of each fiber. The output is three directories of fiber bundles according to left hemisphere, right hemisphere, and commissural tracts. Also copies any Slicer MRML scene file that is given as input into the three separate directories.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputDirectory',
    help='A directory of clustered whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-pthresh', action="store", dest="hemispherePercentThreshold", type=float,
    help='The percent of a fiber that has to be in one hemisphere to consider the fiber as part of that hemisphere (rather than a commissural fiber). Default number is 0.5, where a higher number will tend to label fewer fibers as hemispheric and more fibers as commmissural (not strictly in one hemisphere or the other), while a lower number will be stricter about what is classified as commissural.')

args = parser.parse_args()


if not os.path.isdir(args.inputDirectory):
    print "<wm_separate_hemispheres.py> Error: Input directory", args.inputDirectory, "does not exist or is not a directory."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_separate_hemispheres.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

        
# default to be changed if user input is there
hemisphere_percent_threshold = 0.500001
if args.hemispherePercentThreshold is not None:
    if (args.hemispherePercentThreshold > 0.5) & (args.hemispherePercentThreshold <= 1.0):
        hemisphere_percent_threshold = args.hemispherePercentThreshold
    else:
        print "<wm_separate_hemispheres.py> Hemisphere fiber percent threshold", args.hemispherePercentThreshold, "must be between 0.5 and 1. (0.75 to 0.95 recommmended range)."
        exit

print "<wm_separate_hemispheres.py> Hemisphere fiber percent threshold", hemisphere_percent_threshold
        
print "<wm_separate_hemispheres.py> Starting computation."
print ""
print "=====input directory ======\n", args.inputDirectory
print "=====output directory =====\n", args.outputDirectory
print "=========================="
print ""

# relatively high number of points for accuracy
points_per_fiber = 40

input_polydatas = wma.io.list_vtk_files(args.inputDirectory)
number_of_subjects = len(input_polydatas)

print "<wm_separate_hemispheres.py> Input number of vtk/vtp files: ", number_of_subjects


# read in data
input_pds = list()
for fname in input_polydatas:

    # figure out filename and extension
    fname_base = os.path.basename(fname)

    if fname_base == 'clustered_whole_brain.vtp':
        continue

    if fname_base.startswith('T_CC'):
        fname_output = os.path.join(outdir, fname_base)
        shutil.copyfile(fname, fname_output)
        continue

    # read data
    print "<wm_separate_hemispheres.py> Reading input file:", fname
    pd = wma.io.read_polydata(fname)

    # internal representation for fast similarity computation
    # this also detects which hemisphere fibers are in
    fibers = wma.fibers.FiberArray()
    fibers.points_per_fiber = points_per_fiber
    fibers.hemisphere_percent_threshold = hemisphere_percent_threshold
    # must request hemisphere computation from object
    fibers.hemispheres = True
    # Now convert to array with points and hemispheres as above
    fibers.convert_from_polydata(pd)
        
    # separate into right and left hemispheres
    # note: this assumes RAS coordinates as far as right and left labels
    # LPS would be switched
    # -------------------------
    mask_right = numpy.zeros(pd.GetNumberOfLines())
    mask_right[fibers.index_right_hem] = 1
    pd_right = wma.filter.mask(pd, mask_right, preserve_point_data=True, preserve_cell_data=True,verbose=False)

    mask_left = numpy.zeros(pd.GetNumberOfLines())
    mask_left[fibers.index_left_hem] = 1
    pd_left = wma.filter.mask(pd, mask_left, preserve_point_data=True, preserve_cell_data=True,verbose=False)

    mask_commissure = numpy.zeros(pd.GetNumberOfLines())
    mask_commissure[fibers.index_commissure] = 1
    pd_commissure = wma.filter.mask(pd, mask_commissure, preserve_point_data=True, preserve_cell_data=True,verbose=False)

    print ' - fiber number: left %5d - right %5d - comm %5d' % (pd_left.GetNumberOfLines(), pd_right.GetNumberOfLines(), pd_commissure.GetNumberOfLines())

    # output polydatas into the appropriate subdirectories
    fname_output = os.path.join(outdir, fname_base)
    fname_output = fname_output.replace('.vtp', '_right.vtp')
    wma.io.write_polydata(pd_right, fname_output)

    fname_output = os.path.join(outdir, fname_base)
    fname_output = fname_output.replace('.vtp', '_left.vtp')
    wma.io.write_polydata(pd_left, fname_output)

print 'Done! output to', outdir
