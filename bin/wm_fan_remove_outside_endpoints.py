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
    print "<wm_remove_outlier_endpoints> Error importing white matter analysis package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Remove endpoints that are outside the provided brain mask or freesurfer label",
    epilog="Written by Ye Wu and Fan Zhang")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'input',
    help='Directory or path of input VTK/VTP file(s) that are going to be painted.')
parser.add_argument(
    'label_map_file',
    help='Label map file in nifti (default for freesurfer result).')
parser.add_argument(
    'outputDirectory',
    help='Directory of EP removed results.')
parser.add_argument(
    '-r', type=int, dest="regionList", nargs='+', default=[0],
    help='regions to remove, e.g. 0')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int, default=1,
    help='Number of processors to use.')
parser.add_argument(
    '-removeData', action='store_true',
    help='Remove all data if given')
args = parser.parse_args()

if not os.path.exists(args.input):
    print "Error: Input directory", args.input, "does not exist."
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

if not os.path.exists(args.label_map_file):
    print "Label map", args.outputDirectory, "does not exist, creating it."
    os.makedirs(outdir)

lb = nibabel.load(args.label_map_file)
voxel_data = lb.get_data()
print args.label_map_file, ', label volume shape: ', lb.get_data().shape\


def CorrectStreamlines_From_YE(oldstreamlinesPolyData, mask_img):
    """Remove start and end point not in mask

    Pararmeters
    -----------
    oldstreamlinesPolyData : vtkPolyData,
        vtk fomat about fiber information
    fmask : str,
        path of mask file

    Returns
    -------
    newstreamlinesPolyData : vtkPolyData,
        vtk fomat about fiber information, which is corrected streamlines
    """
    #mask_img = nib.load(fmask)
    mask = mask_img.get_data()
    affine = mask_img.get_affine()
    affine_inv = numpy.linalg.inv(affine)

    newstreamlinesPolyData = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    tractsId = vtk.vtkIdList()
    newstreamlinesPolyData.DeepCopy(oldstreamlinesPolyData)
    newstreamlinesPolyData.BuildCells()
    tractsnum = oldstreamlinesPolyData.GetNumberOfCells()
    cellsnum = newstreamlinesPolyData.GetNumberOfCells()
    old_streamlines = []
    new_streamlines = []

    #Remove all cells in newstreamlines
    for i in range(cellsnum):
        idlist = vtk.vtkIdList()
        newstreamlinesPolyData.GetCellPoints(i, idlist)
        newstreamlinesPolyData.DeleteCell(i)
        del idlist
    newstreamlinesPolyData.RemoveDeletedCells()

    #all cells in oldstreamlines
    for i in range(tractsnum):
        k = 0
        oldstreamlinesPolyData.GetCellPoints(i, tractsId)
        pointsnum = tractsId.GetNumberOfIds()
        fibermrtrix = numpy.zeros((pointsnum, 3))
        for j in range(pointsnum):
            old_point = oldstreamlinesPolyData.GetPoint(tractsId.GetId(j))
            point = numpy.array(list(old_point))
            #point = numpy.dot(affine_inv[:-1, :-1], point.T).T + affine_inv[:-1, -1]
            point = numpy.round(numpy.dot(affine_inv[:-1, :-1], point.T).T + affine_inv[:-1, -1]).astype(int)

            if False:
                point_ras = point
                point_ijk = apply_affine(numpy.linalg.inv(lb.affine), point_ras)
                point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
                print point_ijk, point
                #label = voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])]

            if mask[point[0],point[1],point[2]] > 0:
                fibermrtrix[k] = old_point
                k += 1
            else:
                continue
        #one of newstreamlines cells
        fibermrtrix = numpy.delete(fibermrtrix, numpy.array([i for i in range(k, pointsnum)]), axis=0)
        old_streamlines.append(fibermrtrix)

    for i in range(len(old_streamlines)):
        fibermrtrix = old_streamlines[i]
        for j in range(old_streamlines[i].shape[0]):
            #point = numpy.dot(affine_inv[:-1, :-1], old_streamlines[i][j].T).T + affine_inv[:-1, -1]
            point = numpy.round(numpy.dot(affine_inv[:-1, :-1], old_streamlines[i][j].T).T + affine_inv[:-1, -1]).astype(int)
            if mask[point[0],point[1],point[2]] < 1:
                fibermrtrix = numpy.delete(fibermrtrix, j, 0)
            else:
                break
        for j in range(old_streamlines[i].shape[0])[::-1]:
            #point = numpy.dot(affine_inv[:-1, :-1], old_streamlines[i][j].T).T + affine_inv[:-1, -1]
            point = numpy.round(numpy.dot(affine_inv[:-1, :-1], old_streamlines[i][j].T).T + affine_inv[:-1, -1]).astype(int)
            if mask[point[0],point[1],point[2]] < 1:
                fibermrtrix = numpy.delete(fibermrtrix, j, 0)
            else:
                break
        fibermrtrix = numpy.delete(fibermrtrix, numpy.array([i for i in range(k, pointsnum)]), axis=0)
        new_streamlines.append(fibermrtrix)

    for j in range(len(new_streamlines)):
        #add points
        coord = new_streamlines[j]
        for i in range(coord.shape[0]):
            points.InsertNextPoint(coord[i][0], coord[i][1], coord[i][2])

    a = 0
    for j in range(len(new_streamlines)):
        #line points
        polyline = vtk.vtkPolyLine()
        coord = new_streamlines[j]
        polyline.GetPointIds().SetNumberOfIds(coord.shape[0])

        for i in range(coord.shape[0]):
            polyline.GetPointIds().SetId(i, a)
            a += 1

        cells.InsertNextCell(polyline)
        del polyline

    #set newstreamlinesPolyData
    newstreamlinesPolyData.SetPoints(points)
    newstreamlinesPolyData.SetLines(cells)

    return newstreamlinesPolyData


def remove_endpoints(inpd, mask, region_to_remove, preserve_cell_data=True, preserve_point_data=True, verbose=False):

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()
    incelldata = inpd.GetCellData()

    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outpoints = vtk.vtkPoints()

    if preserve_cell_data:
        if incelldata.GetNumberOfArrays() > 0:
            cell_data_array_indices = range(incelldata.GetNumberOfArrays())
            for idx in cell_data_array_indices:
                array = incelldata.GetArray(idx)
                dtype = array.GetDataType()

                if dtype == 10:
                    out_array = vtk.vtkFloatArray()
                elif dtype == 6:
                    out_array = vtk.vtkIntArray()
                elif dtype == 3:
                    out_array = vtk.vtkUnsignedCharArray()
                else:
                    out_array = vtk.vtkFloatArray()
                out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                out_array.SetName(array.GetName())

                # if verbose:
                #     print "Cell data array found:", array.GetName(), array.GetNumberOfComponents()

                outpd.GetCellData().AddArray(out_array)

    if preserve_point_data:
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = range(inpointdata.GetNumberOfArrays())
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)
                out_array = vtk.vtkFloatArray()
                out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                out_array.SetName(array.GetName())

                # if verbose:
                #     print "Point data array found:", array.GetName(), array.GetNumberOfComponents()

                outpd.GetPointData().AddArray(out_array)

    inpd.GetLines().InitTraversal()
    outlines.InitTraversal()

    for lidx in range(0, inpd.GetNumberOfLines()):
        ptids = vtk.vtkIdList()
        inpd.GetLines().GetNextCell(ptids)

        cellptids = vtk.vtkIdList()

        num_points_before = ptids.GetNumberOfIds()
        num_points_kept = 0
        for pidx in range(0, ptids.GetNumberOfIds()):
            point = inpoints.GetPoint(ptids.GetId(pidx))

            point_ijk = apply_affine(numpy.linalg.inv(lb.affine), point)
            point_ijk = numpy.rint(point_ijk).astype(numpy.int32)
            label = voxel_data[(point_ijk[0], point_ijk[1], point_ijk[2])]

            if not any(label == region_to_remove):
                idx_ = outpoints.InsertNextPoint(point)
                cellptids.InsertNextId(idx_)
                num_points_kept += 1
                if preserve_point_data:
                    for idx in point_data_array_indices:
                        array = inpointdata.GetArray(idx)
                        outpd.GetPointData().GetArray(idx).InsertNextTuple(array.GetTuple(ptids.GetId(pidx)))

        if verbose:
            print 'Line:', lidx, ', before removal:', num_points_before, ', after removal:', num_points_kept,
            print ', point removed:', num_points_before - num_points_kept
        outlines.InsertNextCell(cellptids)

        if preserve_cell_data:
            if incelldata.GetNumberOfArrays() > 0:
                for idx in cell_data_array_indices:
                    array = incelldata.GetArray(idx)
                    out_array = outpd.GetCellData().GetArray(idx)
                    out_array.InsertNextTuple(array.GetTuple(lidx))

    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)

    return outpd

def run_one_cluster(pd_tract_path, ourdir, lb, regionList, removeData):

    cluster_file_name = os.path.split(pd_tract_path)[1]

    pd_tract = wma.io.read_polydata(pd_tract_path)

    num_fibers = pd_tract.GetNumberOfLines()
    print cluster_file_name, ", fiber number:", '{0:05d}'.format(num_fibers)

    print ' - point number before:', pd_tract.GetNumberOfPoints(),
    pd_tract_removed = remove_endpoints(pd_tract, lb, regionList, preserve_point_data = not removeData, preserve_cell_data = not removeData, verbose=True)
    print ', point number after:', pd_tract_removed.GetNumberOfPoints()

    wma.io.write_polydata(pd_tract_removed, os.path.join(outdir, cluster_file_name))

# for pd_tract_path in pd_tract_list:
#     cluster_file_name = os.path.split(pd_tract_path)[1]
#
#     pd_tract = wma.io.read_polydata(pd_tract_path)
#
#     num_fibers = pd_tract.GetNumberOfLines()
#     print cluster_file_name, ", fiber number:", '{0:05d}'.format(num_fibers)
#
#     print ' - point number before:', pd_tract.GetNumberOfPoints(),
#     pd_tract_removed = remove_endpoints(pd_tract, lb, args.regionList)
#     #pd_tract_removed = CorrectStreamlines(pd_tract, lb)
#     print ', point number after:', pd_tract_removed.GetNumberOfPoints()
#
#     wma.io.write_polydata(pd_tract_removed, os.path.join(args.outputDirectory, cluster_file_name))

Parallel(n_jobs=args.numberOfJobs, verbose=1)(
            delayed(run_one_cluster)(pd_tract_path, args.outputDirectory, lb, args.regionList, args.removeData)
            for pd_tract_path in pd_tract_list)

print 'Done!'