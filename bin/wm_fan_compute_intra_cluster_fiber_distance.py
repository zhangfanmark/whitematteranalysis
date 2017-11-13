import connectivity_constraints
import numpy

import os
import argparse
import glob

from joblib import Parallel, delayed

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Calculate within-cluster fiber distance.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

def calculate_mean_pairwise_distance(pd_cluster):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:

        fiber_array = wma.fibers.FiberArray()
        fiber_array.convert_from_polydata(pd_cluster, points_per_fiber=15)

        # pairwise distance matrix
        all_fibers = range(0, fiber_array.number_of_fibers, 5)

        threshold = 0.0
        landmarks = None
        distance_method = 'Hausdorff'
        bilateral = True
        sigmasq = 6400

        if landmarks is None:
            landmarks2 = numpy.zeros((fiber_array.number_of_fibers, 3))
        else:
            landmarks2 = landmarks

        distances = Parallel(n_jobs=1, verbose=0)(
            delayed(wma.similarity.fiber_distance)(
                fiber_array.get_fiber(lidx),
                fiber_array,
                threshold, distance_method=distance_method,
                fiber_landmarks=landmarks2[lidx, :],
                landmarks=landmarks, bilateral=bilateral, sigmasq=sigmasq)
            for lidx in all_fibers)

        distances = numpy.array(distances)

        dis_cluster = numpy.sqrt(distances)

        mean_dis = numpy.mean(dis_cluster)
        std_dis = numpy.std(dis_cluster)
    else:
        mean_dis = 0
        std_dis = 0

    print ' - Mean pairwise distance: %.2f +/- ' % mean_dis, std_dis

    return mean_dis, std_dis

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

print inputdir

stat_file = os.path.join(inputdir, 'Stat_IntraClusterDistance.txt')

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

cluster_paths = list_cluster_files(inputdir)

outstr_distances = 'Cluster Index' + '\t'+ 'Mean Intra-cluster Distance' + '\t' + 'Std' + '\n'

for cluster_path in cluster_paths:

    cluster_file_name = os.path.split(cluster_path)[1]
    print '\n', cluster_file_name

    pd_cluster = wma.io.read_polydata(cluster_path)

    num_fibers = pd_cluster.GetNumberOfLines()
    print ' - Number of fibers:', num_fibers

    mean_pairwise_distance, std_pairwise_distance = calculate_mean_pairwise_distance(pd_cluster)
    outstr_distances = outstr_distances + cluster_file_name + '\t' + str(mean_pairwise_distance) + '\t' + str(std_pairwise_distance) + '\n'

output_file = open(stat_file, 'w')
output_file.write(outstr_distances)
output_file.close()

print 'Done! Result is in', stat_file

