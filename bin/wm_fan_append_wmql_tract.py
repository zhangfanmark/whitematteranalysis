import os
import argparse
import vtk
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
    description="Append multiple fiber clusters into one cluster.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDirectory',
    help='Contains fiber clusters as vtkPolyData file(s).')
parser.add_argument(
    'outputDirectory',
    help='The output directory should be a new empty directory. It will be created if needed.')

args = parser.parse_args()

inputdir = os.path.abspath(args.inputDirectory)
if not os.path.isdir(args.inputDirectory):
    print "<wm_append_cluster> Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "<wm_append_cluster> Error: Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

def list_tract_files(input_dir, regstr):
    # Find input files
    input_mask = ("{0}/*"+regstr+"*.vtk").format(input_dir)
    input_mask2 = ("{0}/*"+regstr+"*.vtp").format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

# left_tract_paths = list_tract_files(args.inputDirectory, 'left')
# righ_tract_paths = list_tract_files(args.inputDirectory, 'right')
#
# for left_tract_path, righ_tract_path in zip(left_tract_paths, righ_tract_paths):
#
#     print 'Working on:'
#     print left_tract_path
#     print righ_tract_path
#     pd_left_tract = wma.io.read_polydata(left_tract_path)
#     pd_righ_tract = wma.io.read_polydata(righ_tract_path)
#
#     appender = vtk.vtkAppendPolyData()
#
#     if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
#         appender.AddInputData(pd_left_tract)
#         appender.AddInputData(pd_righ_tract)
#     else:
#         appender.AddInput(pd_left_tract)
#         appender.AddInput(pd_righ_tract)
#
#     appender.Update()
#     pd_appended_tract = appender.GetOutput()
#
#     output_file = os.path.split(left_tract_path)[1].replace('.left', '')
#
#     wma.io.write_polydata(pd_appended_tract, os.path.join(outdir, output_file))

cst_tract_paths = list_tract_files(args.inputDirectory, 'cst')

appender = vtk.vtkAppendPolyData()
for cst_tract_path in cst_tract_paths:

    pd_cst_tract = wma.io.read_polydata(cst_tract_path)

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_cst_tract)
    else:
        appender.AddInput(pd_cst_tract)

appender.Update()
pd_appended_tract = appender.GetOutput()

output_file = os.path.split(cst_tract_path)[1].replace('_pre', '')
output_file = output_file.replace('.right', '')

wma.io.write_polydata(pd_appended_tract, os.path.join(outdir, output_file))


print '<wm_append_cluster> Save result to', outdir
