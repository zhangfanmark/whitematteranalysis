import argparse
import os
import numpy
import glob

try:
    import whitematteranalysis as wma
except:
    print "<wm_cluster_wmql_region> Error importing white matter analysis package\n"
    raise

parser = argparse.ArgumentParser(
    description="Find the WMQL regions per cluster across multiple subjects.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'inputDir',
    help='Contains all txt files obtained from wm_cluster_wmql_region_per_subject.py.')

args = parser.parse_args()

input_mask_comm = "{0}/*commissural.txt".format(args.inputDir)
region_per_subject_comm_paths = glob.glob(input_mask_comm)

input_mask_left = "{0}/*left.txt".format(args.inputDir)
region_per_subject_left_paths = glob.glob(input_mask_left)

input_mask_right = "{0}/*right.txt".format(args.inputDir)
region_per_subject_right_paths = glob.glob(input_mask_right)

num_subjects = len(region_per_subject_comm_paths)

region_per_subject_comm_list = []
region_per_subject_left_list = []
region_per_subject_right_list = []
for comm_file, left_file, right_file in zip(region_per_subject_comm_paths, region_per_subject_left_paths, region_per_subject_right_paths):
    print '<wm_cluster_wmql_region> Loading', comm_file.replace('_commissural.txt', '*')
    comm = numpy.genfromtxt(comm_file, delimiter='\t', dtype="|S")
    left = numpy.genfromtxt(left_file, delimiter='\t', dtype="|S")
    right = numpy.genfromtxt(right_file, delimiter='\t', dtype="|S")

    region_per_subject_comm_list.append(comm)
    region_per_subject_left_list.append(left)
    region_per_subject_right_list.append(right)

region_list = comm[0, 2:]
num_regions = len(region_list);

comm_cluster_list = comm[1:, 0]
hemi_cluster_list = left[1:, 0]

num_clusters_comm = len(comm_cluster_list)
num_clusters_hemi = len(hemi_cluster_list)

num_clusters = num_clusters_comm + num_clusters_hemi

hemi_or_comm_list = ['hemi'] * num_clusters
for c_comm_idx in range(num_clusters_comm):
    c_idx = int(comm_cluster_list[c_comm_idx].replace('cluster_', '')) - 1
    hemi_or_comm_list[c_idx] = 'comm'


def region_of_one_cluster(region_per_subject_list, c_sub_idx, num_subjects, num_regions):
    subject_region_matrix_of_one_cluster = numpy.zeros((num_subjects, num_regions))

    for s_idx in range(num_subjects):
        region_one_subject = region_per_subject_list[s_idx][c_sub_idx + 1, 1:]
        region_one_subject = [float(r) for r in region_one_subject] # convert into float

        if region_one_subject[0] == 0: # if cluster is empty
            region_belong = []
        else:
            region_one_subject_percentage = numpy.divide(region_one_subject, region_one_subject[0])

            if numpy.sum(region_one_subject_percentage[1:]) == 0: # if no WMQL result
                region_belong = []
            else:
                region_belong = numpy.argsort(region_one_subject_percentage[1:])[-1] # select the top region the cluster belongs to

        subject_region_matrix_of_one_cluster[s_idx, region_belong] = 1
    #     print 's'+str(s_idx), region_belong,
    # print ''

    num_subjects_per_region = numpy.sum(subject_region_matrix_of_one_cluster, axis=0)
    result_region_index = numpy.argsort(num_subjects_per_region)[-1]

    if num_subjects_per_region[result_region_index] == 0:
        result_region_index = None

    return result_region_index

print '\n<wm_cluster_wmql_region> Extracting WMQL regions: \n'
for c_idx in range(num_clusters):
    cluster_name = 'cluster_{0:05d}'.format(c_idx+1)

    if hemi_or_comm_list[c_idx] == 'comm':
        continue
        # need to update as multiple parcellation are given
        c_comm_idx = [i for i, s in enumerate(comm_cluster_list) if cluster_name in s][0]

        result_region_index_comm = region_of_one_cluster(region_per_subject_comm_list, c_comm_idx, num_subjects, num_regions)

        if result_region_index_left is not None:
            print region_list[result_region_index_left],
        else:
            print None,

    elif hemi_or_comm_list[c_idx] == 'hemi':
        print '<wm_cluster_wmql_region>', cluster_name, 'from', hemi_or_comm_list[c_idx],

        c_hemi_idx = [i for i, s in enumerate(hemi_cluster_list) if cluster_name in s][0]

        result_region_index_left = region_of_one_cluster(region_per_subject_left_list, c_hemi_idx, num_subjects, num_regions)
        result_region_index_right = region_of_one_cluster(region_per_subject_right_list, c_hemi_idx, num_subjects, num_regions)

        if result_region_index_left is not None:
            print region_list[result_region_index_left],
        else:
            print None,

        if result_region_index_right is not None:
            print region_list[result_region_index_right]
        else:
            print None