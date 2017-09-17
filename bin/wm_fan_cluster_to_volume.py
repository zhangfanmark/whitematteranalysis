import argparse
import os
import nibabel
import vtk
import numpy
from nibabel.affines import apply_affine
from joblib import Parallel, delayed

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Convert a fiber cluster (vtk) to a volume image.",
    epilog="Written by Fan Zhang")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'input',
    help='Directory or path of input VTK/VTP file(s) that are going to be converted.')
parser.add_argument(
    'volume',
    help='A volume image that the cluster will be converted on. Note this an image no matter it is.')
parser.add_argument(
    'outputDirectory',
    help='Directory of converted volumes.')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int, default=1,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.exists(args.input):
    print "Error: Input directory", args.input, "does not exist."
    exit()

if not os.path.exists(args.volume):
    print "Error: Input volume", args.volume, "does not exist, creating it."
    exit()

if os.path.isdir(args.input):
    pd_tract_list = wma.io.list_vtk_files(args.input)
    if len(pd_tract_list) == 0:
        print "Error: No cluster files found."
        exit()
else:
    pd_tract_list = []
    pd_tract_list.append(args.input)

outdir = os.path.abspath(args.outputDirectory)
if not os.path.exists(args.outputDirectory):
    print "Output directory", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

volume = nibabel.load(args.volume)
print args.volume, ', input volume shape: ', volume.get_data().shape

def convert_cluster_to_volume(inpd, volume_shape):

    new_voxel_data = numpy.zeros(volume_shape)

    inpoints = inpd.GetPoints()

    inpd.GetLines().InitTraversal()
    for lidx in range(0, inpd.GetNumberOfLines()):

        ptids = vtk.vtkIdList()
        inpd.GetLines().GetNextCell(ptids)

        for pidx in range(0, ptids.GetNumberOfIds()):
            point = inpoints.GetPoint(ptids.GetId(pidx))

            point_ijk = apply_affine(numpy.linalg.inv(volume.affine), point)
            point_ijk = numpy.rint(point_ijk).astype(numpy.int32)

            new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1

    return new_voxel_data

def run_one_cluster(pd_tract_path, ourdir, volume):

    cluster_file_name = os.path.split(pd_tract_path)[1]

    pd_tract = wma.io.read_polydata(pd_tract_path)

    num_fibers = pd_tract.GetNumberOfLines()
    print cluster_file_name, ", fiber number:", '{0:05d}'.format(num_fibers)

    new_voxel_data = convert_cluster_to_volume(pd_tract, volume.get_data().shape)

    volume_cluster = nibabel.Nifti1Image(new_voxel_data, volume.affine, volume.header)

    nibabel.save(volume_cluster, os.path.join(outdir, cluster_file_name.replace('.vtp', '.nii.gz')))


Parallel(n_jobs=args.numberOfJobs, verbose=0)(
            delayed(run_one_cluster)(pd_tract_path, args.outputDirectory, volume)
            for pd_tract_path in pd_tract_list)

print 'Done!'