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
    'tract',
    help='tract name')
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


tract_paths = list_tract_files(args.inputDirectory, args.tract)

appender = vtk.vtkAppendPolyData()
for tract_path in tract_paths:

    pd_tract = wma.io.read_polydata(tract_path)

    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        appender.AddInputData(pd_tract)
    else:
        appender.AddInput(pd_tract)

appender.Update()
pd_appended_tract = appender.GetOutput()

output_file = os.path.join(args.outputDirectory, args.tract+'.vtp')

wma.io.write_polydata(pd_appended_tract, output_file)


print '<wm_append_cluster> Save result to', outdir
