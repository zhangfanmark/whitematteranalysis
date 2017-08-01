#!/usr/bin/env python
import numpy
import argparse
import os
import multiprocessing
import time

import vtk
import glob

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
                "This script divide the whole computation into multiple small jobs, os the whole process could be finished soon. "
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
input_data = None # wma.io.read_polydata(args.inputFile)

def _rectangular_similarity_matrix(input_polydata_n, input_polydata_m, threshold, sigma,
                                number_of_jobs=3, landmarks_n=None, landmarks_m=None, distance_method='Hausdorff',
                                bilateral=False, tmp_folder=None):

    """ Internal convenience function available to clustering
    routines.

    Computes similarity matrix (nxn) for all n fibers in input
    polydata.  Calls function _pairwise_distance_matrix first.

    For distance computation and conversion to similarity uses
    parameters threshold and sigma.

    """

    distances = _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                                             number_of_jobs, landmarks_n, landmarks_m, distance_method, bilateral=bilateral, tmp_folder=tmp_folder)

    if distance_method == 'StrictSimilarity':
        similarity_matrix = distances
    else:
        # similarity matrix
        sigmasq = sigma * sigma
        # similarity_matrix_ = wma.similarity.distance_to_similarity(distances, sigmasq)

        similarity_matrix = numpy.zeros((distances.shape[0], distances.shape[1]))
        for r in range(distances.shape[0]):
            similarity_matrix[r, :] = numpy.exp(-distances[r, :] / (sigmasq))
        # print similarity_matrix == similarity_matrix_

    return similarity_matrix


def _rectangular_distance_matrix(input_polydata_n, input_polydata_m, threshold,
                                 number_of_jobs=3, landmarks_n=None, landmarks_m=None,
                                 distance_method='Hausdorff', bilateral=False, tmp_folder=None):
    """ Internal convenience function available to clustering
    routines.

    Computes distance matrix (nxm) for all n+m fibers in input
    polydata. each fiber in input_polydata_n is compared to each fiber
    in input_polydata_m.


    """

    if distance_method == 'Frechet':

        distances = wma.similarity.rectangular_frechet_distances(input_polydata_n, input_polydata_m)
        distances = numpy.array(distances)

    else:
        if False:
            fiber_array_n = wma.fibers.FiberArray()
            fiber_array_n.convert_from_polydata(input_polydata_n, points_per_fiber=15)
            fiber_array_m = wma.fibers.FiberArray()
            fiber_array_m.convert_from_polydata(input_polydata_m, points_per_fiber=15)

            if landmarks_n is None:
                landmarks_n = numpy.zeros((fiber_array_n.number_of_fibers, 3))

            # pairwise distance matrix
            all_fibers_n = range(0, fiber_array_n.number_of_fibers)

            distances = Parallel(n_jobs=number_of_jobs,
                                 verbose=0)(
                delayed(wma.similarity.fiber_distance)(
                    fiber_array_n.get_fiber(lidx),
                    fiber_array_m,
                    threshold, distance_method=distance_method,
                    fiber_landmarks=landmarks_n[lidx, :],
                    landmarks=landmarks_m, bilateral=bilateral)
                for lidx in all_fibers_n)

            distances = numpy.array(distances).T

        else:

            num_fibers_m = input_polydata_m.GetNumberOfLines() # atlas
            #num_fibers_n = input_polydata_n.GetNumberOfLines() # new subject

            # masks_n = numpy.divide(range(num_fibers_n), num_fibers_n / 5)
            # pd_sub_n_list = wma.cluster.mask_all_clusters(input_polydata_n, masks_n, numpy.max(masks_n) + 1, color=None,
            #                                   preserve_point_data=False, preserve_cell_data=False, verbose=False)

            input_mask = "{0}/tract_subdivision_*.vtp".format(tmp_folder)
            input_pd_fnames = glob.glob(input_mask)
            input_pd_fnames = sorted(input_pd_fnames)
            n_tract_subdivision = len(input_pd_fnames)

            fiber_array_m = wma.fibers.FiberArray()
            fiber_array_m.convert_from_polydata(input_polydata_m, points_per_fiber=15)

            # dirpath = tempfile.mkdtemp()
            # print 'Create a temp dir:', dirpath

            # try:
            for m_n in range(n_tract_subdivision):
                tmp_distance_npy = os.path.join(tmp_folder, 'distances_sub_n_' + str(m_n) + '.npy')

                if not os.path.exists(tmp_distance_npy):
                    numpy.save(tmp_distance_npy, [])
                    print ' -Divide ploydata m and n:', m_n, 'in', n_tract_subdivision
                    fname_c = os.path.join(outdir, 'tmp') + '/tract_subdivision_{0:05d}.vtp'.format(m_n)
                    print '  Read', fname_c
                    #input_polydata_n_sub = pd_sub_n_list[m_n]
                    input_polydata_n_sub = wma.io.read_polydata(fname_c)

                    fiber_array_n_sub = wma.fibers.FiberArray()
                    fiber_array_n_sub.convert_from_polydata(input_polydata_n_sub, points_per_fiber=15)

                    all_fibers_n_sub = range(0, fiber_array_n_sub.number_of_fibers)

                    distances_sub_n = Parallel(n_jobs=number_of_jobs, verbose=0)(
                        delayed(wma.similarity.fiber_distance)(
                            fiber_array_n_sub.get_fiber(lidx),
                            fiber_array_m,
                            threshold, distance_method=distance_method,
                            fiber_landmarks=None,
                            landmarks=None, bilateral=bilateral)
                        for lidx in all_fibers_n_sub)

                    distances_sub_n = numpy.array(distances_sub_n).T
                    numpy.save(tmp_distance_npy, distances_sub_n)

            distances = []
            for m_n in range(n_tract_subdivision):
                distances_sub_n = numpy.load(os.path.join(tmp_folder, 'distances_sub_n_' + str(m_n) + '.npy'))
                if distances_sub_n.shape[0] == 0:
                    print 'There are still unfinished jobs!'
                    exit()
                distances.append(distances_sub_n)
                # distances[:, mask_n] = distances_sub_n
            distances = numpy.concatenate(distances, axis=1)
            print 'Distances computed', distances.shape
            # except:
            #     print 'Error: exit!'
            #     exit()

    return distances

#-----------------
# Cluster the data using clusters from the atlas
#-----------------

# num_lines = input_data.GetNumberOfLines()
# fiber_mask = numpy.ones(num_lines)
# input_data_line_only = wma.filter.mask(input_data, fiber_mask, preserve_point_data=False, preserve_cell_data=False, verbose=False)

input_polydata = input_data

# number_fibers = input_polydata.GetNumberOfLines()
# sz = atlas.nystrom_polydata.GetNumberOfLines()

if not os.path.exists(os.path.join(outdir, 'tmp')):
    os.makedirs(os.path.join(outdir, 'tmp'))

if os.path.exists(os.path.join(outdir, 'tmp/B.npy')):
    print 'Matrix B has been computed!'
    exit()

B = _rectangular_similarity_matrix(input_polydata, atlas.nystrom_polydata,
                                   atlas.threshold, atlas.sigma, number_of_jobs=1, distance_method=atlas.distance_method,
                                   bilateral=atlas.bilateral, tmp_folder=os.path.join(outdir, 'tmp'))

if os.path.exists(os.path.join(outdir, 'tmp/B.npy')):
    print 'Matrix B has been computed!'
else:
    numpy.save(os.path.join(outdir, 'tmp/B.npy'), B)


