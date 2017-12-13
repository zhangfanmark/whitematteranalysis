import os
import argparse
import vtk

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise

import glob
#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append multiple fiber clusters defined in a MRML file into one fiber tract.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    '-tractMRML', dest="tractMRML", 
    help='A MRML file that contains the cluster id')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

if not os.path.exists(args.tractMRML):
    print "Error: tract MRML", args.tractMRML, "does not exist, creating it."
    exit()

outdir = os.path.join(args.outputDirectory)
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

cluster_vtp_list = []

f = open(args.tractMRML, 'r') 
for line in f:
    idx = line.find('.vtp')
    if idx > 0:
        cluster_vtp_list.append(line[idx-13:idx+4])

print 'Cluster to be combined:', cluster_vtp_list

input_mask = "{0}/*".format(args.inputDirectory)
input_directories = sorted(glob.glob(input_mask))

for c_idx in range(len(cluster_vtp_list)):
    cluster_vtp = cluster_vtp_list[c_idx]

    print cluster_vtp
    # Loop over inputs and try to find clusters
    appender = vtk.vtkAppendPolyData()
    fname_list = list()
    for dir in input_directories:
        if os.path.isdir(dir):
            fname = os.path.join(dir, cluster_vtp)
            pd_cluster = wma.io.read_polydata(fname)

            if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
                appender.AddInputData(pd_cluster)
            else:
                appender.AddInput(pd_cluster)

    appender.Update()

    pd_appended_cluster = appender.GetOutput()
    output_file = os.path.join(outdir, cluster_vtp)
    wma.io.write_polydata(pd_appended_cluster, output_file)

    print 'Save result to', output_file
