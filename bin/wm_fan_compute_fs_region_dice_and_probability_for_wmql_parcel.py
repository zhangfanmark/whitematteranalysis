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
parser.add_argument(
    'inputAtlasFile',
    help='Stat_AllRegionProbability.csv or Stat_AtlasTractProfile.csv')

args = parser.parse_args()

def list_cluster_file(input_dir, reg):
    # Find input files
    str_reg = "{0}/*_"+reg+".vtp"
    input_mask = str_reg.format(input_dir)
    input_pd_fnames = glob.glob(input_mask)
    return (input_pd_fnames)

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


def calculate_region_dice(pd_cluster, region_prob_cluster_atlas, thres_ranges, occ_profile='Occ'):
    num_fibers = pd_cluster.GetNumberOfLines()

    atlas_regions_different_thres = list()
    for thres in thres_ranges:
        atlas_regions_different_thres.append(numpy.where(region_prob_cluster_atlas >= thres / 100.0)[0])

    dice_list = []
    if num_fibers > 0:
        if occ_profile == 'Occ':
            along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions='LR', verbose=False)
        elif occ_profile == 'Profile':
            along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions='GW', verbose=False)

        for atlas_regions in atlas_regions_different_thres:
            mean_dice = []
            for fiber_regions in along_tract_regions:
                if (len(fiber_regions) + len(atlas_regions)) > 0:

                    mean_dice.append(len(numpy.intersect1d(fiber_regions, atlas_regions)) * 2.0 / (len(fiber_regions) + len(atlas_regions)))
                else:
                    mean_dice.append(0)

            dice_list.append(numpy.mean(mean_dice))

    print ' - Mean region dice at different thres: ', dice_list
    return dice_list

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

NUMBER_OF_REGION = 5002
inputcsv = args.inputAtlasFile
if not os.path.exists(inputcsv):
    print "Error: ", inputcsv, " does not exist."
    exit()
else:
    region_prob_atlas_list = []
    parcels_to_be_analyzed = []
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
            parcels_to_be_analyzed.append(row[0])

# cluster_paths = list_cluster_files(inputdir)

input_file_name = os.path.split(args.inputAtlasFile)[1]

occ_profile = ''
if 'AllRegionProbability' in input_file_name:
    occ_profile = 'Occ'
elif 'Profile' in input_file_name:
    occ_profile = 'Profile'

# probability scores
if occ_profile == 'Occ':
    prob_file = os.path.join(inputdir, 'Stat_RegionProbability.txt')
elif occ_profile == 'Profile':
    prob_file = os.path.join(inputdir, 'Stat_ProfileProbability.txt')

if not os.path.exists(prob_file):
    outstr_region_prob = 'Cluster Index' + '\t'+ 'Region Probability Atlas' + '\t'+ 'Region Probability Self' + '\n'

    for parce_pre, region_prob_cluster_atlas in zip(parcels_to_be_analyzed, region_prob_atlas_list):

        cluster_path = list_cluster_file(inputdir, parce_pre)

        if len(cluster_path) > 0:
            cluster_file_name = os.path.split(cluster_path[0])[1]
            print '\n', cluster_file_name

            pd_cluster = wma.io.read_polydata(cluster_path[0])

            prob_atlas, prob_self = calculate_region_probability(pd_cluster, region_prob_cluster_atlas, occ_profile)
        else:
            prob_atlas = 0
            prob_self = 0

        outstr_region_prob = outstr_region_prob + cluster_file_name + '\t' + str(prob_atlas) + '\t' + str(prob_self) + '\n'

    output_file = open(prob_file, 'w')
    output_file.write(outstr_region_prob)
    output_file.close()

# dice score with threshold
if occ_profile == 'Occ':
    dice_file = os.path.join(inputdir, 'Stat_RegionDice.txt')
elif occ_profile == 'Profile':
    dice_file = os.path.join(inputdir, 'Stat_ProfileDice.txt')

thres_ranges = range(30, 90, 10)
if not os.path.exists(dice_file):
    outstr_region_dice = 'Cluster Index'
    for thres in thres_ranges:
        outstr_region_dice += '\t' + 'Region Dice : Threshold ' + str(thres) + '%'
    outstr_region_dice += '\n'

    for parce_pre, region_prob_cluster_atlas in zip(parcels_to_be_analyzed, region_prob_atlas_list):

        cluster_path = list_cluster_file(inputdir, parce_pre)

        if len(cluster_path) > 0:
            cluster_file_name = os.path.split(cluster_path[0])[1]
            print '\n', cluster_file_name

            pd_cluster = wma.io.read_polydata(cluster_path[0])

            dice_list = calculate_region_dice(pd_cluster, region_prob_cluster_atlas, thres_ranges, occ_profile)
        else:
            dice_list = []

        outstr_region_dice += cluster_file_name
        for dice in dice_list:
            outstr_region_dice += '\t' + str(dice)
        outstr_region_dice += '\n'

    output_file = open(dice_file, 'w')
    output_file.write(outstr_region_dice)
    output_file.close()

print ''
print 'Done! Result is in', prob_file
print 'Done! Result is in', dice_file