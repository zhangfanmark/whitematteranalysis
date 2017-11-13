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
    'atlasDirectory',
    help='Contains Stat_AllRegionProbability.csv')

parser.add_argument(
    'outDirectory',
    help='Contains Stat_AllRegionProbability.csv')

args = parser.parse_args()

outdir = os.path.abspath(args.outDirectory)
if not os.path.isdir(args.outDirectory):
    print "Error: Output directory", args.outDirectory, "does not exist."
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


outstr = 'Cluster Index' + '\t'+ 'Region' + '\n'
thres = 40
for c_idx, region_prob_cluster_atlas in zip(range(len(region_prob_atlas_list)), region_prob_atlas_list):
    regions = numpy.where(region_prob_cluster_atlas >= thres / 100.0)[0]

    outstr = outstr + 'cluster_{0:05d}'.format(c_idx+1) + '\t' + numpy.array_str(regions) + '\n'

output_file = open(os.path.join(outdir, 'FiberClusterAnnotation.txt'), 'w')
output_file.write(outstr)
output_file.close()

print ''
print 'Done! Result is in', os.path.join(outdir, 'FiberClusterAnnotation.txt')
