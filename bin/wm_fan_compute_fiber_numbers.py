import connectivity_constraints
import numpy

import os
import argparse
import glob

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Calculate number of fibers per cluster.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

print inputdir

stat_file = os.path.join(inputdir, 'Stat_FiberNumber.txt')

if os.path.exists(stat_file):
    print "\n Already Computed!!! \n"
    exit()

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/cluster_*.vtk".format(input_dir)
    input_mask2 = "{0}/cluster_*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

def list_tract_files(input_dir):
    # Find input files
    input_mask = "{0}/*.vtk".format(input_dir)
    input_mask2 = "{0}/*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

cluster_paths = list_cluster_files(inputdir)

if len(cluster_paths) == 0:
    cluster_paths = list_tract_files(inputdir)

outstr_num_fibers = 'Cluster Index' + '\t'+ 'Fiber Number' + '\n'
for cluster_path in cluster_paths:

    cluster_file_name = os.path.split(cluster_path)[1]
    print '\n', cluster_file_name

    pd_cluster = wma.io.read_polydata(cluster_path)

    num_fibers = pd_cluster.GetNumberOfLines()
    print ' - Number of fibers:', num_fibers
    outstr_num_fibers = outstr_num_fibers + cluster_file_name + '\t' + str(num_fibers) + '\n'

output_file = open(stat_file, 'w')
output_file.write(outstr_num_fibers)
output_file.close()

print 'Done! Result is in', stat_file


