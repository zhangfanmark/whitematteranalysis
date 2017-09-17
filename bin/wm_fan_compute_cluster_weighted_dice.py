import argparse
import glob
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
    'input1',
    help='Directory or path of input VTK/VTP file(s) that are going to be converted.')

parser.add_argument(
    'input2',
    help='Directory or path of input VTK/VTP file(s) that are going to be converted.')

parser.add_argument(
    'volume',
    help='A volume image that the cluster will be converted on. Note this an image no matter it is.')

args = parser.parse_args()


if not os.path.exists(args.input1):
    print "Error: Input directory", args.input1, "does not exist."
    exit()

if not os.path.exists(args.input1):
    print "Error: Input directory", args.input2, "does not exist."
    exit()

volume = nibabel.load(args.volume)
print args.volume, ', input volume shape: ', volume.get_data().shape

input1_vtk_list = wma.io.list_vtk_files(args.input1)
input2_vtk_list = wma.io.list_vtk_files(args.input2)

print 'Number of clusters', len(input1_vtk_list)


def convert_cluster_to_volume(inpd, volume):

    volume_shape  =volume.get_data().shape
    new_voxel_data = numpy.zeros(volume_shape)

    resampler = vtk.vtkPolyDataPointSampler()
    resampler.GenerateEdgePointsOn()
    resampler.GenerateVertexPointsOff()
    resampler.GenerateInteriorPointsOff()
    resampler.GenerateVerticesOff()
    resampler.SetDistance(0.5)

    inpoints = inpd.GetPoints()

    inpd.GetLines().InitTraversal()
    for lidx in range(0, inpd.GetNumberOfLines()):

        ptids = vtk.vtkIdList()
        inpd.GetLines().GetNextCell(ptids)

        tmpPd = vtk.vtkPolyData()
        tmpPoints = vtk.vtkPoints()
        tmpCellPtIds = vtk.vtkIdList()
        tmpLines =  vtk.vtkCellArray()
        
        for pidx in range(0, ptids.GetNumberOfIds()):
            point = inpoints.GetPoint(ptids.GetId(pidx))
            idx_ = tmpPoints.InsertNextPoint(point)
            tmpCellPtIds.InsertNextId(idx_)

        tmpLines.InsertNextCell(tmpCellPtIds)

        tmpPd.SetLines(tmpLines)
        tmpPd.SetPoints(tmpPoints)

        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            resampler.SetInputData(tmpPd)
        else:
            resampler.SetInput(tmpPd)

        resampler.Update()

        sampledCellPts = resampler.GetOutput().GetPoints()
        sampledNpts = resampler.GetOutput().GetNumberOfPoints()

        for pidx in range(0, sampledNpts):
            point = sampledCellPts.GetPoint(pidx)

            point_ijk = apply_affine(numpy.linalg.inv(volume.affine), point)
            point_ijk = numpy.rint(point_ijk).astype(numpy.int32)

            new_voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])] += 1

    return new_voxel_data

outstr = 'Cluster Index' + '\t'+ 'wDice' + '\t'+ 'sDice' + '\n'

for input1_vtk_path, input2_vtk_path in zip(input1_vtk_list, input2_vtk_list):

    print os.path.split(input1_vtk_path)[1]

    input1_vtk = wma.io.read_polydata(input1_vtk_path)
    input2_vtk = wma.io.read_polydata(input2_vtk_path)

    n_1_fibers = input1_vtk.GetNumberOfLines()
    n_2_fibers = input2_vtk.GetNumberOfLines()

    print " - Number of fibers: p1 %5d, p2 %5d" %(n_1_fibers, n_2_fibers)

    vox_1 = convert_cluster_to_volume(input1_vtk, volume)
    vox_2 = convert_cluster_to_volume(input2_vtk, volume)

    if n_1_fibers > 0 and n_2_fibers > 0:

        voxel_data_1 = vox_1 / n_1_fibers
        voxel_data_2 = vox_2 / n_2_fibers

        mask_data_1 = numpy.sign(voxel_data_1)
        mask_data_2 = numpy.sign(voxel_data_2)

        n_1 = numpy.sum(mask_data_1)
        n_2 = numpy.sum(mask_data_2)

        w_1 = numpy.sum(voxel_data_1)
        w_2 = numpy.sum(voxel_data_2)

        mask_intersection = numpy.logical_and(mask_data_1, mask_data_2)

        n_ = numpy.sum(mask_intersection)

        w_1_ = numpy.sum(voxel_data_1[mask_intersection])
        w_2_ = numpy.sum(voxel_data_2[mask_intersection])

        print ' - cluster 1: %5d voxels, %5d fibers, %5d intersected voxels' % (n_1, n_1_fibers, n_)
        print ' - cluster 2: %5d voxels, %5d fibers, %5d intersected voxels' % (n_2, n_2_fibers, n_)

        dice_standard = 2 * n_ / (n_1 + n_2)
        dice_weighted = (w_1_ + w_2_) / (w_1 + w_2)
    else:
        dice_standard = 0
        dice_weighted = 0

    print ' ** wDice = %2f, sDice = %2f' % (dice_weighted, dice_standard)

    outstr = outstr + os.path.split(input1_vtk_path)[1] + '\t' + str(dice_weighted) + '\t' + str(dice_standard) + '\n'

output_file = open(os.path.join(args.input2, 'Stat_RetestDice_RS.txt'), 'w')
output_file.write(outstr)
output_file.close()














