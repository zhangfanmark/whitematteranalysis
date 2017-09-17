import connectivity_constraints
import numpy

import os
import argparse
import glob

try:
    import whitematteranalysis as wma
except:
    print "<wm_append_cluster.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Calculate statistics for each fiber cluster.",
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
        top_two_labels = numpy.array([])
        mean_percent = 0

    print ' - Mean endpoint percentage: %.2f, top two endpoint regions:' % mean_percent, top_two_labels

    return mean_percent, top_two_labels

def calculate_along_tract(pd_cluster, threshold_percentage=0.5, verbose=False):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, verbose=False)

        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)

        for r_idx in range(all_regions_all_fibers.shape[0]):
            all_regions_all_fibers[r_idx] = connectivity_constraints.region_label(all_regions_all_fibers[r_idx])

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
        tract_regions = numpy.array([])
        mean_percent = 0
        std_percent = 0

    print ' - Mean along-tract percentage: %.2f, with tract regions' % (mean_percent), tract_regions

    return mean_percent, std_percent, tract_regions

def calculate_mean_pairwise_distance(pd_cluster):

    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        dis_cluster = wma.cluster._pairwise_distance_matrix(pd_cluster,
                        threshold=0.0, number_of_jobs=1, landmarks=None, distance_method='Hausdorff', bilateral=True)

        dis_cluster = numpy.sqrt(dis_cluster)

        mean_dis = numpy.mean(dis_cluster)
        std_dis = numpy.std(dis_cluster)
    else:
        mean_dis = 0
        std_dis = 0

    print ' - Mean pairwise distance: %.2f' % mean_dis

    return mean_dis, std_dis

def calculate_number_of_subjects(pd_cluster):

    if pd_cluster.GetNumberOfLines() > 0:
        number_of_subjects = 0

        subjects_IDs = connectivity_constraints.get_subject_ID(pd_cluster)

        if subjects_IDs is not None:
            data_array = []
            l = [0]
            for lidx in range(subjects_IDs.GetNumberOfTuples()):
                subjects_IDs.GetTupleValue(lidx, l)
                data_array.append(l[0])

            data_array = numpy.array(data_array)

            number_of_subjects = numpy.unique(data_array).shape[0]
    else:
        number_of_subjects = 0

    print ' - Number of subjects:', number_of_subjects

    return number_of_subjects

def calculate_all_regions_on_cluster(pd_cluster):
    num_fibers = pd_cluster.GetNumberOfLines()

    if num_fibers > 0:
        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, verbose=False)

        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)

        for r_idx in range(all_regions_all_fibers.shape[0]):
            all_regions_all_fibers[r_idx] = connectivity_constraints.region_label(all_regions_all_fibers[r_idx])

        unique_regions = numpy.unique(all_regions_all_fibers)

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

print inputdir
if os.path.exists(os.path.join(inputdir, 'Stat_IntraClusterDistance.txt')):
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

outstr_num_fibers = 'Cluster Index' + '\t'+ 'Fiber Number' + '\n'
outstr_num_subjects = 'Cluster Index' + '\t'+ 'Subject Number' + '\n'
outstr_distances = 'Cluster Index' + '\t'+ 'Mean Intra-cluster Distance' + '\t' + 'Std' + '\n'
outstr_endpoints = 'Cluster Index' + '\t'+ 'Endpoint FS Dice Score' + '\t' + 'Regions' + '\n'
outstr_along_tract = 'Cluster Index' + '\t'+ 'Along-tract FS Dice Score' + '\t' + 'Regions' + '\n'

for cluster_path in cluster_paths:

    cluster_file_name = os.path.split(cluster_path)[1]
    print '\n', cluster_file_name

    pd_cluster = wma.io.read_polydata(cluster_path)

    num_fibers = pd_cluster.GetNumberOfLines()
    print ' - Number of fibers:', num_fibers
    outstr_num_fibers = outstr_num_fibers + cluster_file_name + '\t' + str(num_fibers) + '\n'

    number_of_subjects = calculate_number_of_subjects(pd_cluster)
    outstr_num_subjects = outstr_num_subjects + cluster_file_name + '\t' + str(number_of_subjects) + '\n'

    mean_pairwise_distance, std_pairwise_distance = calculate_mean_pairwise_distance(pd_cluster)
    outstr_distances = outstr_distances + cluster_file_name + '\t' + str(mean_pairwise_distance) + '\t' + str(std_pairwise_distance) + '\n'

    mean_endpoint_percent, top_two_labels = calculate_endpoint(pd_cluster)
    outstr_endpoints = outstr_endpoints + cluster_file_name + '\t' + str(mean_endpoint_percent) + '\t' + numpy.array_str(top_two_labels) + '\n'

    mean_along_tract_percent, std_along_tract_percent, tract_regions = calculate_along_tract(pd_cluster)
    outstr_along_tract = outstr_along_tract + cluster_file_name + '\t' + str(mean_along_tract_percent) + '\t' + numpy.array_str(tract_regions) + '\n'


output_file = open(os.path.join(inputdir, 'Stat_FiberNumber.txt'), 'w')
output_file.write(outstr_num_fibers)
output_file.close()

output_file = open(os.path.join(inputdir, 'Stat_SubjectNumber.txt'), 'w')
output_file.write(outstr_num_subjects)
output_file.close()

output_file = open(os.path.join(inputdir, 'Stat_IntraClusterDistance.txt'), 'w')
output_file.write(outstr_distances)
output_file.close()

output_file = open(os.path.join(inputdir, 'Stat_Endpoints.txt'), 'w')
output_file.write(outstr_endpoints)
output_file.close()

output_file = open(os.path.join(inputdir, 'Stat_AlongTractRegions.txt'), 'w')
output_file.write(outstr_along_tract)
output_file.close()

