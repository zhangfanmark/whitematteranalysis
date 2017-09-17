import connectivity_constraints
import numpy

import os
import argparse
import glob
import csv

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
parser.add_argument(
    'atlasDirectory',
    help='Contains Stat_AllRegionProbability.csv')

args = parser.parse_args()

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/cluster_*.vtk".format(input_dir)
    input_mask2 = "{0}/cluster_*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

def calculate_region_probability(pd_cluster, region_prob_cluster_atlas):
    num_fibers = pd_cluster.GetNumberOfLines()

    mean_prob_list_atlas = 0
    mean_prob_list_self = 0
    if num_fibers > 0:
        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, verbose=False)

        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)
        unique_regions =  range(NUMBER_OF_REGION + 1)
        region_occurrence = numpy.bincount(all_regions_all_fibers)
        region_prob_cluster_self = region_occurrence / float(num_fibers)
        # print numpy.where(region_prob_cluster_self > 0)[0]
        # print region_prob_cluster_self[numpy.where(region_prob_cluster_self > 0)[0]]
        # print numpy.where(region_prob_cluster_atlas > 0)[0]
        # print region_prob_cluster_atlas[numpy.where(region_prob_cluster_atlas > 0)[0]]

        prob_list_atlas = []
        prob_list_self = []
        for fiber_regions in along_tract_regions:
            # print fiber_regions.astype(int)
            # print region_prob_cluster_atlas[fiber_regions.astype(int)]
            # print numpy.sum(region_prob_cluster_atlas[fiber_regions.astype(int)])
            # print numpy.sum(region_prob_cluster_atlas)
            # exit()

            fiber_prob_atlas = numpy.sum(region_prob_cluster_atlas[fiber_regions.astype(int)]) / numpy.sum(region_prob_cluster_atlas)
            fiber_prob_self = numpy.sum(region_prob_cluster_self[fiber_regions.astype(int)]) / numpy.sum(region_prob_cluster_self)
            prob_list_atlas.append(fiber_prob_atlas)
            prob_list_self.append(fiber_prob_self)

        mean_prob_list_atlas = numpy.mean(prob_list_atlas)
        mean_prob_list_self = numpy.mean(prob_list_self)

    print ' - Mean region probability: atlas %.2f, self %.2f' % (mean_prob_list_atlas, mean_prob_list_self)
    
    return mean_prob_list_atlas, mean_prob_list_self

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

NUMBER_OF_REGION = 5002
inputcsv = os.path.join(args.atlasDirectory, 'Stat_AllRegionProbability.csv')
if not os.path.exists(inputcsv):
    print "Error: Stat_AllRegionProbability.csv does not exist."
    exit()
else:
    region_prob_atlas_list = []
    with open(inputcsv, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='\'')
        first_row = True
        for row in spamreader:
            if first_row:
                first_row = False
                continue

            region_prob_cluster_atlas = numpy.zeros(NUMBER_OF_REGION+1) # region_prob_cluster_atlas[0] will never be used as region start with 1
            for ele in row[1:-1]:
                region = ele.split(':')[0]
                prob = ele.split(':')[1]
                region_prob_cluster_atlas[int(region)] = prob

            region_prob_atlas_list.append(region_prob_cluster_atlas)

cluster_paths = list_cluster_files(inputdir)

outstr_region_prob = 'Cluster Index' + '\t'+ 'Region Probability Atlas' + '\t'+ 'Region Probability Self' + '\n'

for cluster_path, region_prob_cluster_atlas in zip(cluster_paths, region_prob_atlas_list):

    cluster_file_name = os.path.split(cluster_path)[1]
    print '\n', cluster_file_name

    pd_cluster = wma.io.read_polydata(cluster_path)

    prob_atlas, prob_self = calculate_region_probability(pd_cluster, region_prob_cluster_atlas)
    outstr_region_prob = outstr_region_prob + cluster_file_name + '\t' + str(prob_atlas) + '\t' + str(prob_self) + '\n'


output_file = open(os.path.join(inputdir, 'Stat_RegionProbability.txt'), 'w')
output_file.write(outstr_region_prob)
output_file.close()



