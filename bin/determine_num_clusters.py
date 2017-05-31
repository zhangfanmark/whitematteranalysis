
import whitematteranalysis as wma
import os
import connectivity_constraints
import numpy

import os
import argparse
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_append_cluster.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append multiple fiber clusters into one cluster.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

def calculate_endpoint(pd_cluster):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        endpoint_regions_per_cluster = connectivity_constraints.get_endpiont_region(pd_cluster)

        all_endpoint_regions = []
        for ep in endpoint_regions_per_cluster:
            all_endpoint_regions.append(ep[0])
            all_endpoint_regions.append(ep[1])

        ep_label_occurrence = numpy.bincount(all_endpoint_regions)
        top_two_labels = numpy.argsort(ep_label_occurrence)[-2:]

        mean_percent = sum(ep_label_occurrence[top_two_labels] / float(num_fibers)) / 2
    else:
        top_two_labels = []
        mean_percent = 0

    print ' - Mean endpoint percentage:', mean_percent, ', top two endpoint regions:', top_two_labels

    return mean_percent

def calculate_along_tract(pd_cluster, threshold_percentage=0.5, verbose=False):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, verbose=False)

        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)

        unique_regions =  numpy.unique(all_regions_all_fibers)

        region_occurrence = numpy.bincount(all_regions_all_fibers)

        tract_regions = unique_regions[region_occurrence[unique_regions] / float(num_fibers) > threshold_percentage]

        if verbose:
            print 'All possible regions', unique_regions
            print 'Region occurrence:', region_occurrence[unique_regions]

        dice_all = []
        for region_one_fiber in along_tract_regions:
            dice_one_fiber =  len(numpy.intersect1d(region_one_fiber, tract_regions)) * 2.0  / (len(region_one_fiber) + len(tract_regions))

            dice_all.append(dice_one_fiber)
            if verbose:
                print 'Fiber region', region_one_fiber, tract_regions, ', Dice Score:', dice_one_fiber

        mean_percent = numpy.mean(dice_all)
        std_percent = numpy.std(dice_all)

    else:
        tract_regions = []
        mean_percent = 0
        std_percent = 0

    print ' - Mean along tract FS percent:', mean_percent, '+/-', std_percent, ', tract regions:', tract_regions

    return mean_percent, std_percent

def calculate_mean_pairwise_distance(pd_cluster):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        dis_cluster = wma.cluster._pairwise_distance_matrix(pd_cluster, threshold=0)

        dis_cluster = numpy.sqrt(dis_cluster)

        mean_dis = numpy.mean(dis_cluster)
        std_dis = numpy.std(dis_cluster)
    else:
        mean_dis = 0
        std_dis = 0

    print ' - Mean pairwise distance:', mean_dis, '+/-', std_dis

    return mean_dis, std_dis

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

cluster_paths = wma.io.list_vtk_files(inputdir)

outstr = 'cluster index' + '\t' + 'fiber number' + '\t' + 'dice along tract mean' + '\t' + 'dice along tract std' + '\t' + 'pairwise distance mean' + '\t' + 'pairwise distance std' + '\t' + 'endpoint mean' + '\n'

for cluster_path in cluster_paths:
    print 'Analyzing', os.path.split(cluster_path)[1]

    if cluster_path.find('clustered_whole_brain.vtp') != -1 or cluster_path.find('atlas.vtp') != -1:
        continue

    pd_cluster = wma.io.read_polydata(cluster_path)

    num_fibers = pd_cluster.GetNumberOfLines()
    print ' - Number of fibers:', num_fibers

    mean_endpoint_percent = calculate_endpoint(pd_cluster)

    mean_along_tract_percent, std_along_tract_percent = calculate_along_tract(pd_cluster)

    mean_pairwise_distance, std_pairwise_distance = calculate_mean_pairwise_distance(pd_cluster)

    outstr = outstr + os.path.split(cluster_path)[1] + '\t' + str(num_fibers) + '\t' + str(mean_along_tract_percent) + '\t' + str(std_along_tract_percent)  \
             + '\t' + str(mean_pairwise_distance) + '\t' + str(std_pairwise_distance) + '\t' + str(mean_endpoint_percent) + '\n'

output_file = open(os.path.join(inputdir, 'statistics_cluster.txt'), 'w')
output_file.write(outstr)
output_file.close()

