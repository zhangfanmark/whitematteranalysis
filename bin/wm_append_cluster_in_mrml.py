import os
import argparse
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_append_cluster.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append multiple fiber clusters into one cluster.",
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
    help='mrml file that contains the cluster id')

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

appender = vtk.vtkAppendPolyData()
for c_idx in range(len(cluster_vtp_list)):
    cluster_vtp = cluster_vtp_list[c_idx]
    pd_cluster = wma.io.read_polydata(os.path.join(inputdir, cluster_vtp))

    vtk_array = vtk.vtkIntArray()
    vtk_array.SetName('cluster_idx')
    for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
        vtk_array.InsertNextTuple1(int(cluster_vtp_list[c_idx][8:12]))

    pd_cluster.GetPointData().AddArray(vtk_array)
    pd_cluster.Update()

    print '<wm_append_cluster>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_cluster)
    else:
        appender.AddInput(pd_cluster)

appender.Update()
pd_appended_cluster = appender.GetOutput()

output_file = os.path.split(args.tractMRML)[1][:-5] + '.vtp'
wma.io.write_polydata(pd_appended_cluster, os.path.join(outdir, output_file))

print '<wm_append_cluster> Appended clusters , number of fibers', pd_appended_cluster.GetNumberOfLines()
print ''
print '<wm_append_cluster> Save result to', os.path.join(outdir, output_file)
