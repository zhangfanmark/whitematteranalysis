#!/usr/bin/env python
import numpy
import argparse
import os
import multiprocessing
import time

import vtk

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise

try:
    from joblib import Parallel, delayed
    USE_PARALLEL = 1
except ImportError:
    USE_PARALLEL = 0
    print "<cluster.py> Failed to import joblib, cannot multiprocess."
    print "<cluster.py> Please install joblib for this functionality."

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Pre compute the pairwise fiber distance B in cluster.spectral_atlas_label(). "
                "This script divide the whole brain tractography into multiple tracts"
                "It takes the same parameters as in wm_cluster_from_atlas.py",
    epilog="Written by Fan Zhang. fzhang@bwh.harvard.edu")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputFile',
    help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
parser.add_argument(
    'atlasDirectory',
    help='The directory where the atlas is stored. Must contain atlas.p and atlas.vtp')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    '-f', action="store", dest="numberOfFibers", type=int,
    help='Number of fibers to analyze from each subject. If this parameter is not used, all fibers will be analyzed by default.')
parser.add_argument(
    '-l', action="store", dest="fiberLength", type=int, default=60,
    help='Minimum length (in mm) of fibers to analyze. 60mm is default.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int, default=1,
    help='Number of processors to use.')
parser.add_argument(
    '-verbose', action='store_true', dest="flag_verbose",
    help='Verbose. Run with -verbose for more text output.')
parser.add_argument(
    '-mrml_fibers', action="store", dest="showNFibersInSlicer", type=float,
    help='Approximate upper limit on number of fibers to show when MRML scene of clusters is loaded into slicer.Default is 10000 fibers; increase for computers with more memory. Note this can be edited later in the MRML file by searching for SubsamplingRatio and editing that number throughout the file. Be sure to use a text editor program (save as plain text format). An extra MRML file will be saved for visualizing 100%% of fibers.')
parser.add_argument(
    '-reg', action='store_true', dest="registerAtlasToSubjectSpace",
    help='To cluster in individual subject space, register atlas polydata to subject. Otherwise, by default this code assumes the subject has already been registered to the atlas.')
parser.add_argument(
    '-norender', action='store_true', dest="flag_norender",
    help='No Render. Prevents rendering of images that would require an X connection.')

args = parser.parse_args()

if not os.path.exists(args.inputFile):
    print "<wm_cluster_from_atlas.py> Error: Input file", args.inputFile, "does not exist."
    exit()

if not os.path.isdir(args.atlasDirectory):
    print "<wm_cluster_from_atlas.py> Error: Atlas directory", args.atlasDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "<wm_cluster_from_atlas.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

fname = args.inputFile
subject_id = os.path.splitext(os.path.basename(fname))[0]
outdir = os.path.join(outdir, subject_id)
if not os.path.exists(outdir):
    print "<wm_cluster_from_atlas.py> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "\n=========================="
print "input file:", args.inputFile
print "atlas directory:", args.atlasDirectory
print "output directory:", args.outputDirectory

if args.numberOfFibers is not None:
    print "fibers to analyze per subject: ", args.numberOfFibers
else:
    print "fibers to analyze per subject: ALL"
number_of_fibers = args.numberOfFibers

fiber_length = args.fiberLength
print "minimum length of fibers to analyze (in mm): ", fiber_length

number_of_jobs = args.numberOfJobs
print 'Using N jobs:', number_of_jobs

if args.flag_verbose:
    print "Verbose ON."
else:
    print "Verbose OFF."
verbose = args.flag_verbose

if args.showNFibersInSlicer is not None:
    show_fibers = args.showNFibersInSlicer
else:
    show_fibers = 10000.0
print "Maximum total number of fibers to display in MRML/Slicer: ", show_fibers

if args.registerAtlasToSubjectSpace:
    print "Registration of atlas fibers to subject fibers is ON."
    print "Warning: the registration pipeline is under improvements--the registration implemented here is not currently supported. Please register to atlas first, then call this script without the -reg option."
    exit()
else:
    print "Registration of atlas fibers to subject fibers is OFF. Subject must be in atlas space before calling this script."

if args.flag_norender:
    print "No rendering (for compute servers without X connection)."
else:
    print "Rendering. After clustering, will create colorful jpg images."
render = not args.flag_norender

print "==========================\n"

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

# read atlas
print "<wm_cluster_from_atlas.py> Loading input atlas:", args.atlasDirectory
atlas = wma.cluster.load_atlas(args.atlasDirectory, 'atlas')

# read data
print "<wm_cluster_from_atlas.py> Reading input file:", args.inputFile
input_data = wma.io.read_polydata(args.inputFile)

#-----------------
# Cluster the data using clusters from the atlas
#-----------------

# num_lines = input_data.GetNumberOfLines()
# fiber_mask = numpy.ones(num_lines)
# input_data_line_only = wma.filter.mask(input_data, fiber_mask, preserve_point_data=False, preserve_cell_data=False, verbose=False)

input_polydata = input_data

number_fibers = input_polydata.GetNumberOfLines()
sz = atlas.nystrom_polydata.GetNumberOfLines()

if not os.path.exists(os.path.join(outdir, 'tmp')):
    os.makedirs(os.path.join(outdir, 'tmp'))

masks_n = numpy.divide(range(number_fibers), number_fibers / 99)
pd_sub_n_list = wma.cluster.mask_all_clusters(input_polydata, masks_n, numpy.max(masks_n) + 1, color=None,
                                              preserve_point_data=False, preserve_cell_data=False, verbose=False)

for c in range(len(pd_sub_n_list)):
    pd_c = pd_sub_n_list[c]

    fname_c = os.path.join(outdir, 'tmp') + '/tract_subdivision_{0:05d}.vtp'.format(c)
    print fname_c
    wma.io.write_polydata(pd_c, fname_c)


