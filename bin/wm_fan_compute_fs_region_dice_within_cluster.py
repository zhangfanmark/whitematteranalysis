import connectivity_constraints
import numpy

import os
import argparse
import glob
import csv

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Calculate the consistency of the FS regions of the fibers in each cluster.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/*.vtk".format(input_dir)
    input_mask2 = "{0}/*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

def calculate_region_probability(pd_cluster, region_prob_cluster_atlas, occ_profile='Occ'):
    num_fibers = pd_cluster.GetNumberOfLines()

    mean_prob_list_atlas = 0
    mean_prob_list_self = 0
    if num_fibers > 0:

        if occ_profile == 'Occ':
            along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions='LR', verbose=False)
        elif occ_profile == 'Profile':
            along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions='GW', verbose=False)

        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)
        region_occurrence = numpy.bincount(all_regions_all_fibers)
        region_prob_cluster_self = region_occurrence / float(num_fibers)

        prob_list_atlas = []
        prob_list_self = []
        for fiber_regions in along_tract_regions:

            fiber_prob_atlas = numpy.sum(region_prob_cluster_atlas[fiber_regions.astype(int)]) / numpy.sum(region_prob_cluster_atlas)
            fiber_prob_self = numpy.sum(region_prob_cluster_self[fiber_regions.astype(int)]) / numpy.sum(region_prob_cluster_self)
            prob_list_atlas.append(fiber_prob_atlas)
            prob_list_self.append(fiber_prob_self)

        mean_prob_list_atlas = numpy.mean(prob_list_atlas)
        mean_prob_list_self = numpy.mean(prob_list_self)

    print ' - Mean region probability: atlas %.2f, self %.2f' % (mean_prob_list_atlas, mean_prob_list_self)
    
    return mean_prob_list_atlas, mean_prob_list_self


def calculate_self_region_dice(pd_cluster, combine_regions='LR'):
    num_fibers = pd_cluster.GetNumberOfLines()

    step = num_fibers / 10

    self_dice = numpy.nan
    if num_fibers != 0:
        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions=combine_regions,  verbose=False)

        mean_dice = []
        for f_idx_1 in range(0, len(along_tract_regions), step):
            fiber_regions_1 = along_tract_regions[f_idx_1]
            for f_idx_2 in range(f_idx_1, len(along_tract_regions), step):
                fiber_regions_2 = along_tract_regions[f_idx_2]
                mean_dice.append(len(numpy.intersect1d(fiber_regions_1, fiber_regions_2)) * 2.0 / (len(fiber_regions_1) + len(fiber_regions_2)))

        self_dice = numpy.mean(mean_dice)

    return self_dice

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

cluster_paths = list_cluster_files(inputdir)
cluster_paths = cluster_paths[0:20]
prob_file = os.path.join(inputdir, 'Stat_SelfRegionDice.txt')

if not os.path.exists(prob_file):
    outstr_region_prob = 'Cluster Index' + '\t'+ 'Region Dice Self' + '\n'

    for cluster_path in cluster_paths:

        cluster_file_name = os.path.split(cluster_path)[1]
        print '\n', cluster_file_name

        pd_cluster = wma.io.read_polydata(cluster_path)

        self_dice = calculate_self_region_dice(pd_cluster)

        print 'Self region dice', self_dice

        outstr_region_prob = outstr_region_prob + cluster_file_name  + '\t' + str(self_dice) + '\n'


    output_file = open(prob_file, 'w')
    output_file.write(outstr_region_prob)
    output_file.close()

print 'Done! Result is in', output_file