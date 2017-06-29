
import whitematteranalysis as wma
import numpy
import connectivity_constraints
import os
import argparse


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="wm_fan_subject_tract_regions",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

def extract_along_tract_regions(pd_cluster):

    along_tract_regions_per_fiber = connectivity_constraints.get_along_tract_region(pd_cluster, verbose=False)
    along_tract_regions_all_fibers = numpy.concatenate(along_tract_regions_per_fiber).astype(int)

    along_tract_regions_unique = numpy.unique(along_tract_regions_all_fibers)
    occurrence_per_region = numpy.bincount(along_tract_regions_all_fibers)[along_tract_regions_unique]

    idx_sorted = numpy.argsort(occurrence_per_region)[::-1][:]

    along_tract_regions = along_tract_regions_unique[idx_sorted]
    region_occurrence = occurrence_per_region[idx_sorted]

    return (along_tract_regions, region_occurrence)


def extract_endpoint_regions(pd_cluster):
    endpoint_regions_per_fiber = connectivity_constraints.get_endpiont_region(pd_cluster)

    endpoint_regions_all_fibers = []
    for ep in endpoint_regions_per_fiber:
        endpoint_regions_all_fibers.append(ep[0])
        endpoint_regions_all_fibers.append(ep[1])

    endpoint_regions_all_fibers = numpy.array(endpoint_regions_all_fibers).astype(int)

    endpoint_regions_all_fibers = endpoint_regions_all_fibers[numpy.where(endpoint_regions_all_fibers>0)]

    endpoint_regions_unique = numpy.unique(endpoint_regions_all_fibers)
    occurrence_per_region = numpy.bincount(endpoint_regions_all_fibers)[endpoint_regions_unique]

    idx_sorted = numpy.argsort(occurrence_per_region)[::-1][:]

    endpoint_regions = endpoint_regions_unique[idx_sorted]
    region_occurrence = occurrence_per_region[idx_sorted]

    return (endpoint_regions, region_occurrence)

def output_regions(region_data_list, output_file, max_region):

    outstr = 'Cluster Idx' + '\t' + 'Fiber Number'
    for r in range(max_region):
        outstr += '\t' + 'Region' + str(r)
    outstr += '\n'

    for file_name, num_fibers, regions, occurrence in region_data_list:
        substr1 = file_name + '\t' + str(num_fibers)
        substr2 = file_name + '\t' + str(num_fibers)
        tmp_c = 0
        for r, o in zip(regions, occurrence):
            substr1 += '\t' + str(r)
            substr2 += '\t' + str(o)
            tmp_c += 1

        for r in range(tmp_c, max_region):
            substr1 += '\t' + ' '
            substr2 += '\t' + ' '

        outstr += substr1 + '\n'
        outstr += substr2 + '\n'

    output_file = open(output_file, 'w')
    output_file.write(outstr)
    output_file.close()

pd_paths = wma.io.list_vtk_files(args.inputDirectory)
along_tract_list = []
endpoint_list = []
max_along_tract_region = 0
max_endpoint_region = 0
for pd_path in pd_paths:

    fiber_name = os.path.split(pd_path)[1]
    pd_cluster = wma.io.read_polydata(pd_path)
    num_fibers = pd_cluster.GetNumberOfLines()

    print fiber_name, ', num fiber:', num_fibers

    if num_fibers > 0:

        along_tract_regions, at_region_occurrence = extract_along_tract_regions(pd_cluster)
        print ' - Along-tract:', along_tract_regions, numpy.divide(at_region_occurrence.astype(float), num_fibers)
        if len(along_tract_regions) > max_along_tract_region:
            max_along_tract_region = len(along_tract_regions)

        endpoint_regions, ep_region_occurrence = extract_endpoint_regions(pd_cluster)
        print ' - Endpoint:', endpoint_regions, numpy.divide(ep_region_occurrence.astype(float), num_fibers)
        if len(endpoint_regions) > max_endpoint_region:
            max_endpoint_region = len(endpoint_regions)

    else:
        print ' - Empty'
        along_tract_regions = []
        at_region_occurrence = []
        endpoint_regions = []
        ep_region_occurrence = []

    along_tract_list.append((fiber_name, num_fibers, along_tract_regions, at_region_occurrence))
    endpoint_list.append((fiber_name, num_fibers, endpoint_regions, ep_region_occurrence))

output_regions(along_tract_list, os.path.join(args.inputDirectory, 'region_along_tract.txt'), max_along_tract_region)
output_regions(endpoint_list, os.path.join(args.inputDirectory, 'region_endpoint.txt'), max_endpoint_region)

exit()

if True:
    dice_cluster_matrix = numpy.zeros([len(tract_regions_list), len(tract_regions_list)])

    for pd_1_idx in range(len(tract_regions_list)):
        print pd_1_idx
        tract_1_regions = tract_regions_list[pd_1_idx]
        for pd_2_idx in range(pd_1_idx, len(tract_regions_list)):
            tract_2_regions = tract_regions_list[pd_2_idx]

            try:
                if len(tract_1_regions) > 3 or len(tract_2_regions) > 3:
                    dice_cluster_matrix[pd_1_idx, pd_2_idx] = len(numpy.intersect1d(tract_1_regions, tract_2_regions)) * 2.0 / (len(tract_1_regions) + len(tract_2_regions))

                    dice_cluster_matrix[pd_2_idx, pd_1_idx] = dice_cluster_matrix[pd_1_idx, pd_2_idx]
            except:
                None

    print dice_cluster_matrix

    numpy.save(os.path.join(args.inputDirectory, 'dice_cluster_matrix'), dice_cluster_matrix)

if True:
    dice_fiber_matrix = numpy.zeros([len(tract_regions_list), len(tract_regions_list)])
    for pd_idx in range(len(tract_regions_list)):
        tract_regions = tract_regions_list[pd_idx]
        print pd_idx
        for pd_sub_idx in range(len(tract_regions_list)):

            if dice_cluster_matrix[pd_idx, pd_sub_idx] < 0.5:
                along_tract_regions = along_tract_regions_list[pd_sub_idx]

                dice_all = []
                for region_one_fiber in along_tract_regions:
                    dice_one_fiber = len(numpy.intersect1d(region_one_fiber, tract_regions)) * 2.0 / (len(region_one_fiber) + len(tract_regions))

                    dice_all.append(dice_one_fiber)

                mean_percent = numpy.mean(dice_all)

                dice_fiber_matrix[pd_idx, pd_sub_idx] = mean_percent

    print dice_fiber_matrix
    numpy.save(os.path.join(args.inputDirectory, 'dice_fiber_matrix'), dice_fiber_matrix)


