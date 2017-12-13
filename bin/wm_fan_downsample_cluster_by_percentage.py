#!/usr/bin/env python


import glob
import os
import argparse
import multiprocessing
import math

try:
    import whitematteranalysis as wma
except:
    print "<wm_laterality.py> Error importing white matter analysis package\n"
    raise

try:
    from joblib import Parallel, delayed
except:
    print "<wm_laterality.py> Error importing joblib package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Applies preprocessing to input directory. Downsamples, removes short fibers. Preserves tensors and scalar point data along retained fibers.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'input',
    help='Contains whole-brain tractography as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    '-p', action="store", dest="percentageOfFibers", type=float,
    help='Percentage of fibers to keep from each dataset.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')
parser.add_argument(
    '-save_data', action='store_true', dest="save_data",
    help='Save along tract data or not')

args = parser.parse_args()


if not os.path.exists(args.input):
    print "Error: Input directory", args.input, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

print "wm_laterality. Starting white matter laterality computation."
print ""
print "=====input directory======\n", args.input
print "=====output directory=====\n", args.outputDirectory
print "=========================="

if args.percentageOfFibers is not None:
    print "fibers to retain per subject: ", args.percentageOfFibers
else:
    print "fibers to retain per subject: ALL"
    args.percentageOfFibers = 1

print 'CPUs detected:', multiprocessing.cpu_count()
if args.numberOfJobs is not None:
    parallel_jobs = args.numberOfJobs
else:
    parallel_jobs = 1
print 'Using N jobs:', parallel_jobs


print "=========================="

# =======================================================================
# Above this line is argument parsing. Below this line is the pipeline.
# =======================================================================

# Loop over input DWIs
if os.path.isdir(args.input):
    inputPolyDatas = wma.io.list_vtk_files(args.input)
else:
    inputPolyDatas = []
    inputPolyDatas.append(args.input)


print "<wm_preprocess.py> Input number of files: ", len(inputPolyDatas)

# for testing
#inputPolyDatas = inputPolyDatas[0:2]

def pipeline(inputPolyDatas, sidx, args):
    # get subject identifier from unique input filename
    # -------------------
    subjectID = os.path.splitext(os.path.basename(inputPolyDatas[sidx]))[0]
    id_msg = "<wm_preprocess.py> ", sidx + 1, "/", len(inputPolyDatas)  
    msg = "**Starting subject:", subjectID
    print(id_msg + msg)

    fname = os.path.join(args.outputDirectory, subjectID + '.vtp')
    if os.path.exists(fname):
        print 'Already processed!'
        return

    # read input vtk data
    # -------------------
    msg = "**Reading input:", subjectID
    print(id_msg + msg)

    wm = wma.io.read_polydata(inputPolyDatas[sidx])

    num_lines = wm.GetNumberOfLines()
    print "Input number of fibers", num_lines

    num_fibers_of_cluster = math.ceil(wm.GetNumberOfLines() * args.percentageOfFibers)
    pd_fcluster_ds, pd_indices = wma.filter.downsample(wm, num_fibers_of_cluster, \
                                                       return_indices=True,
                                                       preserve_point_data=False, preserve_cell_data=False,
                                                       verbose=False)

    # outputs
    # -------------------
    msg = "**Writing output data for subject:", subjectID
    print id_msg, msg

    try:
        print "Writing output polydata", fname, "..."
        wma.io.write_polydata(pd_fcluster_ds, fname)
        print "Wrote output", fname, "."
    except:
        print "Unknown exception in IO"
        raise
    del pd_fcluster_ds

# loop over all inputs
Parallel(n_jobs=parallel_jobs, verbose=0)(
        delayed(pipeline)(inputPolyDatas, sidx, args)
        for sidx in range(0, len(inputPolyDatas)))

exit()
