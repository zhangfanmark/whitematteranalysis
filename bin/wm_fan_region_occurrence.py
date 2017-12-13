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
    description="calculate tract profile",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')

args = parser.parse_args()

## THIS IS THE SAME CODE TO THE TRACT PROFILE. ONLY get_along_tract_region WAS MODIFIED TO combine_regions='LR'

def calculate_all_regions(pd_cluster):
    num_fibers = pd_cluster.GetNumberOfLines()

    str_occurrence = ''
    if num_fibers > 0:

        along_tract_regions = connectivity_constraints.get_along_tract_region(pd_cluster, combine_regions='LR', verbose=False)
        all_regions_all_fibers = numpy.concatenate(along_tract_regions).astype(int)

        unique_regions = numpy.unique(all_regions_all_fibers)
        region_occurrence = numpy.bincount(all_regions_all_fibers)[unique_regions]

        idx = numpy.argsort(-region_occurrence)
        region_occurrence = numpy.divide(region_occurrence[idx].astype(float), num_fibers)
        unique_regions = unique_regions[idx]

        for r_l, r_o in zip(unique_regions, region_occurrence):
            str_occurrence += str(r_l) + ':' + str(r_o) + ','

    print ' - Region occurrence:', str_occurrence

    return str_occurrence

###################################################################

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

print inputdir

stat_file = os.path.join(inputdir, 'Stat_AllRegionOccurrence.txt')

if os.path.exists(stat_file):
    print "\n Already Computed!!! \n"
    exit()

def list_cluster_files(input_dir):
    # Find input files
    input_mask = "{0}/*.vtk".format(input_dir)
    input_mask2 = "{0}/*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

cluster_paths = list_cluster_files(inputdir)

outstr_all_region_occurrence = 'Cluster Index' + '\t'+ 'Region Occurrence' + '\n'

for cluster_path in cluster_paths:

    cluster_file_name = os.path.split(cluster_path)[1]
    print '\n', cluster_file_name

    pd_cluster = wma.io.read_polydata(cluster_path)

    region_str = calculate_all_regions(pd_cluster)
    outstr_all_region_occurrence = outstr_all_region_occurrence + cluster_file_name + '\t' + region_str + '\n'

output_file = open(stat_file, 'w')
output_file.write(outstr_all_region_occurrence)
output_file.close()

print 'Done! Result is in', stat_file
