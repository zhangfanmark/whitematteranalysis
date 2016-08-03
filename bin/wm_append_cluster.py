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
    'clusterList',
    help='A list of clusters to be appended, separated by comma,  e.g., 1,2,4,5')


args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "<wm_append_cluster> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "<wm_append_cluster> Error: Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

cluster_list = args.clusterList.split(",")
if len(cluster_list) == 0:
    print "<wm_append_cluster> Error: input list", args.clusterList, "is not correct."
    exit()

cluster_vtp_list = list()
output_suffix = ''
for c_idx in range(len(cluster_list)):

    if len(cluster_list[c_idx]) > 5:
        print "<wm_append_cluster> Error: input cluster index", cluster_list[c_idx], "should less that 99999."
        exit()

    try:
        int(cluster_list[c_idx])
    except:
        print "<wm_append_cluster> Error: input cluster index ", cluster_list[c_idx], "should be an integer"
        exit()

    cluster_vtp_list.append('cluster_'+cluster_list[c_idx].zfill(5)+'.vtp')
    output_suffix = output_suffix+cluster_list[c_idx]+'_'

    if not os.path.exists(os.path.join(inputdir, cluster_vtp_list[c_idx])):
        print "<wm_append_cluster> Error: cluster", cluster_vtp_list[c_idx], "does not exist."
        exit()

print ""
print "<wm_append_cluster> Starting appending cluster."
print ""
print "=====input directory======\n", inputdir
print "=====output directory=====\n", outdir
print "=====clusters to be appended====\n", cluster_vtp_list
print ""

appender = vtk.vtkAppendPolyData()
for c_idx in range(len(cluster_vtp_list)):
    cluster_vtp = cluster_vtp_list[c_idx]
    pd_cluster = wma.io.read_polydata(os.path.join(inputdir, cluster_vtp))

    vtk_array = vtk.vtkIntArray()
    vtk_array.SetName('cluster_idx')
    for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
        vtk_array.InsertNextTuple1(int(cluster_list[c_idx]))

    pd_cluster_point = pd_cluster.GetPointData().AddArray(vtk_array)
    pd_cluster.GetPointData().AddArray(vtk_array)
    pd_cluster.Update()

    print '<wm_append_cluster>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_cluster)
    else:
        appender.AddInput(pd_cluster)

appender.Update()
pd_appended_cluster = appender.GetOutput()

output_file = 'append_cluster_' + output_suffix[:-1] + '.vtp'
wma.io.write_polydata(pd_appended_cluster, os.path.join(outdir, output_file))

print '<wm_append_cluster> Appended clusters , number of fibers', pd_appended_cluster.GetNumberOfLines()
print ''
print '<wm_append_cluster> Save result to', os.path.join(outdir, output_file)
print '<wm_append_cluster> Use scalar \'cluster_idx\' to see the fibers from different subjects in 3D Slicer.'
