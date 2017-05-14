import os
import argparse
import numpy

try:
    import whitematteranalysis as wma
except:
    print "<wm_cluster_wmql_region_per_subject.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Find the WMQL regions per cluster of an individual subject.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'inputClusterDir',
    help='Contains the vtk files per cluster.')
parser.add_argument(
    'inputWMQLDir',
    help='Contains the WMQL results.')
parser.add_argument(
    'regionListFile',
    help='Contains the WMQL ROI list')
parser.add_argument(
    '-subID',
    type=str,
    help='To be included in the output file name. Normally, this is the subject ID')
parser.add_argument(
    '-hemi',
    type=str, choices=['right', 'left', 'commissural'],
    help='Specify if the clusters are the hemisphere-based separated clusters or the bilateral clusters.')

args = parser.parse_args()

if not os.path.isdir(args.inputClusterDir):
    print "Error: Input cluster directory", args.inputClusterDir, "does not exist."
    exit()

if not os.path.isdir(args.inputWMQLDir):
    print "Error: Input WMQL directory", args.inputWMQLDir, "does not exist."
    exit()

if not os.path.exists(args.regionListFile):
    print "Error: Input regions list", args.regionListFile, "does not exist."
    exit()

pd_cluster_paths = wma.io.list_vtk_files(args.inputClusterDir)
print '<wm_cluster_wmql_region_per_subject.py> Find a total of', len(pd_cluster_paths), 'clusters.'

def list_subdirs(input_dir, with_base=False):

    if with_base is False:
        subdirs = [subdir_name for subdir_name in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, subdir_name)) and subdir_name.find('cluster') != -1]
    elif with_base is True:
        subdirs = [os.path.join(input_dir, subdir_name) for subdir_name in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, subdir_name))]

    subdirs = sorted(subdirs)
    return subdirs

wmql_clusters = list_subdirs(args.inputWMQLDir)
print '<wm_cluster_wmql_region_per_subject.py> Find a total of', len(wmql_clusters), 'WMQL results.'

region_list = numpy.genfromtxt(args.regionListFile, dtype="|S")
print '<wm_cluster_wmql_region_per_subject.py> Find a total of', len(region_list), 'WMQL regions of interest.'
print '  ', region_list

file_name = 'cluster_wmql_region'
if args.subID is not None:
    file_name += '_' + args.subID
if args.hemi is not None:
    file_name += '_' + args.hemi
file_name += '.txt'

output_file = open(os.path.join(args.inputWMQLDir, file_name), 'w')

output_str = 'cluster\ttotal fiber number'
for region in region_list:
    output_str += '\t'+region
output_str += '\n'

for pd_cluster_path in pd_cluster_paths:
    cluster_name = os.path.split(pd_cluster_path)[1][:13]

    cluster_str_list = [''] * (len(region_list) + 1)
    if any(cluster_name in s for s in wmql_clusters):
        print '-', cluster_name, ':',

        pd_cluster = wma.io.read_polydata(pd_cluster_path)
        num_pd_cluster = pd_cluster.GetNumberOfLines()
        cluster_str_list[0] = str(num_pd_cluster)

        wmql_cluster_path = os.path.join(args.inputWMQLDir, cluster_name)
        wmql_tract_paths = wma.io.list_vtk_files(wmql_cluster_path)

        for wmql_tract_path in wmql_tract_paths:
            pd_wmql_tract = wma.io.read_polydata(wmql_tract_path)
            num_pd_wmql_tract = pd_wmql_tract.GetNumberOfLines()

            wmql_region_name = os.path.split(wmql_tract_path)[1].replace(cluster_name+'_', '').replace('.vtp', '')
            wmql_region_index = [i for i, s in enumerate(region_list) if wmql_region_name in s][0]
            print wmql_region_index, '-', wmql_region_name, ',',

            cluster_str_list[wmql_region_index+1] = str(num_pd_wmql_tract)
        print ' '

    output_str += cluster_name
    for cluster_str in cluster_str_list:
        if len(cluster_str) == 0:
            output_str += '\t0'
        else:
            output_str += '\t' + cluster_str
    output_str += '\n'

output_file.write(output_str)
output_file.close()

print '<wm_cluster_wmql_region_per_subject> Find the result at', os.path.join(args.inputWMQLDir, file_name)
print '\n<wm_cluster_wmql_region_per_subject> Done!'
