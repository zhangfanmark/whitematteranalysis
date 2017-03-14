import nibabel
import os
import numpy
import vtk
import whitematteranalysis as wma
from nibabel.affines import apply_affine

def paint_scalars_on_point(pd, scalars, scalar_name):
    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetName(scalar_name)

    for p_idx in range(0, pd.GetNumberOfPoints()):
        vtk_array.InsertNextTuple1(float(scalars[p_idx]))

    pd.GetPointData().AddArray(vtk_array)
    pd.Update()

    return pd


label_map = os.path.join('/Users/fan/Desktop', '100307_wmparc.nii.gz')
tract = os.path.join('/Users/fan/Desktop', 'TP_SMG_100307-ukftrack_b3000_default7_minGA006_minFA008_seedFALimit010_Ql0_Qm0_Rs0.vtk')

lb = nibabel.load(label_map)
pd = wma.io.read_polydata(tract)

print lb
print pd

fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(pd, points_per_fiber=2) # Endpoints only

print fiber_array.fiber_array_r[0, 0], fiber_array.fiber_array_a[0, 0], fiber_array.fiber_array_s[0, 0]

print lb.affine
print lb.get_data().shape

voxel_data = lb.get_data()

pd.GetLines().InitTraversal()
line_ptids = vtk.vtkIdList()
inpoints = pd.GetPoints()

num_fibers = pd.GetNumberOfLines()

for lidx in range(0, num_fibers):

    pd.GetLines().GetNextCell(line_ptids)
    line_length = line_ptids.GetNumberOfIds()

    pidx = 0
    fiber_array_r_1 = numpy.zeros(line_length)
    fiber_array_a_1 = numpy.zeros(line_length)
    fiber_array_s_1 = numpy.zeros(line_length)
    for line_index in range(0, line_length):
        # do nearest neighbor interpolation: round index
        ptidx = line_ptids.GetId(line_index)

        point = inpoints.GetPoint(ptidx)

        fiber_array_r_1[pidx] = point[0]
        fiber_array_a_1[pidx] = point[1]
        fiber_array_s_1[pidx] = point[2]
        pidx = pidx + 1


    merged = numpy.column_stack([fiber_array_r_1, fiber_array_a_1, fiber_array_s_1])

    ijk = apply_affine(numpy.linalg.inv(lb.affine), merged)
    ijk = numpy.rint(ijk).astype(numpy.int32)

    line_labels = []
    for line_index in range(0, line_length):
        label = voxel_data[(ijk[line_index, 0], ijk[line_index, 1], ijk[line_index, 2])]
        line_labels.append(label)


    print 'output', lidx
    fiber_mask = numpy.zeros(num_fibers)
    fiber_mask[lidx] = 1
    pd_fiber = wma.filter.mask(pd, fiber_mask, preserve_point_data=True, preserve_cell_data=True)
    paint_scalars_on_point(pd_fiber, line_labels, 'fs_label')
    wma.io.write_polydata(pd_fiber, '/Users/fan/Desktop/fibers/' + str(lidx) + '.vtp')