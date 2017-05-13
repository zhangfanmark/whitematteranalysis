import os
import argparse
import vtk

try:
    import whitematteranalysis as wma
except:
    print "<wm_combine_all_clusters.py> Error importing white matter analysis package\n"
    raise

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Combine all clusters in a folder into one vtk file. This is used for transforming the fiber clusters per subject into tje DWI space. "
                "Performing the transform on individual clusters would take a long time. So, all clusters first will be combined as one vtk file using this script. "
                "Then, the wm_harden_transform.py script will be used to do the transform from the atlas space to the DWI space."
                "Finally, the wm_separate_all_clusters.py will be used to separated the transformed vtk into multiple fiber clusters. ",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')
parser.add_argument(
    '-subID',
    type=str, default='NG',
    help='File name of the output vtk file. Normally, this is the subject ID')
parser.add_argument(
    '-hemi',
    type=str, default='bilateral', choices=['left', 'right', 'commissural', 'bilateral'],
    help='Specify if the clusters are the hemisphere-based separated clusters or the bilateral clusters.')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "<wm_combine_all_clusters> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

cluster_vtp_list = wma.io.list_vtk_files(inputdir)

if len(cluster_vtp_list) == 0:
    print "<wm_combine_all_clusters> Error: no cluster files found!"
    exit()

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "<wm_combine_all_clusters> Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

print ""
print "<wm_combine_all_clusters> Start combining all clusters."
print ""
print "=====input directory======\n", inputdir
print "=====output directory=====\n", outdir
print "=====number of clusters====\n", len(cluster_vtp_list)
print ""

appender = vtk.vtkAppendPolyData()
for c_idx in range(len(cluster_vtp_list)):
    cluster_vtp = cluster_vtp_list[c_idx]
    pd_cluster = wma.io.read_polydata(os.path.join(inputdir, cluster_vtp))

    vtk_array = vtk.vtkIntArray()
    vtk_array.SetName('cluster_idx')
    for p_idx in range(0, pd_cluster.GetNumberOfPoints()):
        vtk_array.InsertNextTuple1(int(c_idx + 1))

    pd_cluster.GetPointData().AddArray(vtk_array)
    pd_cluster.Update()

    print '<wm_combine_all_clusters>', cluster_vtp, ', number of fibers', pd_cluster.GetNumberOfLines()

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_cluster)
    else:
        appender.AddInput(pd_cluster)

appender.Update()
pd_appended_cluster = appender.GetOutput()

output_file = 'combined_clusters_'+args.subID+'_'+args.hemi+'.vtp'
wma.io.write_polydata(pd_appended_cluster, os.path.join(outdir, output_file))

print '<wm_combine_all_clusters> Combined clusters , number of fibers', pd_appended_cluster.GetNumberOfLines()
print ''
print '<wm_combine_all_clusters> Save result to', os.path.join(outdir, output_file)





