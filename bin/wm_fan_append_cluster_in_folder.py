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
    'outputFile',
    help='The output directory should be a new empty directory. It will be created if needed.')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "<wm_append_cluster> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

cluster_vtp_list = wma.io.list_vtk_files(inputdir)


print ""
print "<wm_append_cluster> Starting appending cluster."
print ""
print "=====input directory======\n", inputdir
print "=====clusters to be appended====\n", cluster_vtp_list
print ""

appender = vtk.vtkAppendPolyData()
for c_idx in range(len(cluster_vtp_list)):
    cluster_vtp = cluster_vtp_list[c_idx]
    pd_cluster = wma.io.read_polydata(os.path.join(inputdir, cluster_vtp))

    vtk_array = vtk.vtkIntArray()
    vtk_array.SetName('cluster_idx')
    for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
        vtk_array.InsertNextTuple1(int(c_idx))

    pd_cluster.GetPointData().AddArray(vtk_array)
    pd_cluster.Update()

    print '<wm_append_cluster>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_cluster)
    else:
        appender.AddInput(pd_cluster)

appender.Update()
pd_appended_cluster = appender.GetOutput()

output_file = args.outputFile
wma.io.write_polydata(pd_appended_cluster, output_file)

print '<wm_append_cluster> Appended clusters , number of fibers', pd_appended_cluster.GetNumberOfLines()
print ''
print '<wm_append_cluster> Save result to', output_file
