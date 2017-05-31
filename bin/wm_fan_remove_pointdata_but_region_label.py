#!/usr/bin/env python
import numpy
import argparse
import os
import multiprocessing
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print "Error importing white matter analysis package\n"
    raise


#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Append multiple fiber clusters into one cluster.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputDir',
    help='input')
parser.add_argument(
    'outDir',
    help='output')

args = parser.parse_args()

def mask_remove(inpd, fiber_mask, color=None, preserve_point_data=True, preserve_cell_data=False, verbose=True):
    """ Keep lines and their points where fiber_mask == 1.

     Unlike vtkMaskPolyData that samples every nth cell, this function
     uses an actual mask, and also gets rid of points that
     are not used by any cell, reducing the size of the polydata file.

     This code also sets scalar cell data to input color data.  This
     input data is expected to be 1 or 3 components.

     If there is no input cell scalar color data, existing cell
     scalars that we have created (EmbeddingColor, ClusterNumber,
     EmbeddingCoordinate) are looked for and masked as well.

     """

    inpoints = inpd.GetPoints()
    inpointdata = inpd.GetPointData()
    incelldata = inpd.GetCellData()

    # output and temporary objects
    ptids = vtk.vtkIdList()
    outpd = vtk.vtkPolyData()
    outlines = vtk.vtkCellArray()
    outpoints = vtk.vtkPoints()
    outcolors = None
    outpointdata = outpd.GetPointData()
    outcelldata = outpd.GetCellData()
    tensor_names = []

    if color is not None:
        # if input is RGB
        if len(color.shape) == 2:
            if color.shape[1] == 3:
                outcolors = vtk.vtkUnsignedCharArray()
                outcolors.SetNumberOfComponents(3)

        # otherwise output floats as colors
        if outcolors == None:
            outcolors = vtk.vtkFloatArray()

    # check for cell data arrays to keep
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
                if verbose:
                    print "Cell data array found:", array.GetName(), array.GetNumberOfComponents()
                outcelldata.AddArray(out_array)
                # make sure some scalars are active so rendering works
                # outpd.GetCellData().SetActiveScalars(array.GetName())

                # if inpd.GetCellData().GetArray('ClusterNumber'):
                #    # this will be active unless we have embedding colors
                #    outpd.GetCellData().SetActiveScalars('ClusterNumber')
                # if inpd.GetCellData().GetArray('EmbeddingColor'):
                #    # Note Slicer can't display these cell scalars (all is red)
                #    outpd.GetCellData().SetActiveScalars('EmbeddingColor')

        else:
            preserve_cell_data = False

    # check for point data arrays to keep
    if preserve_point_data:
        if inpointdata.GetNumberOfArrays() > 0:
            point_data_array_indices = range(inpointdata.GetNumberOfArrays())
            for idx in point_data_array_indices:
                array = inpointdata.GetArray(idx)

                if array.GetName() == 'region_label':
                    out_array = vtk.vtkFloatArray()
                    out_array.SetNumberOfComponents(array.GetNumberOfComponents())
                    out_array.SetName(array.GetName())
                    if verbose:
                        print "Point data array found:", array.GetName(), array.GetNumberOfComponents()
                    outpointdata.AddArray(out_array)
                    # make sure some scalars are active so rendering works
                    # outpd.GetPointData().SetActiveScalars(array.GetName())
                    # keep track of tensors to choose which is active
                    if array.GetNumberOfComponents() == 9:
                        tensor_names.append(array.GetName())
        else:
            preserve_point_data = False

    # Set up scalars and tensors attributes for correct visualization in Slicer.
    # Slicer works with point data and does not handle cell data well.
    # This set of changes breaks old internal wma default visualization of cell scalars.
    # Changes must be propagated through wma so that render is called with the name of the field to visualize.
    # the new way in wma is like this line, below.
    # ren = wma.render.render(output_polydata_s, 1000, data_mode="Cell", data_name='EmbeddingColor')

    # For Slicer: First set one of the expected tensor arrays as default for vis
    tensors_labeled = False
    for name in tensor_names:
        if name == "tensors":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("tensors"))
            tensors_labeled = True
        if name == "Tensors":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("Tensors"))
            tensors_labeled = True
        if name == "tensor1":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("tensor1"))
            tensors_labeled = True
        if name == "Tensor1":
            outpd.GetPointData().SetTensors(outpd.GetPointData().GetArray("Tensor1"))
            tensors_labeled = True
    if not tensors_labeled:
        if len(tensor_names) > 0:
            print "Data has unexpected tensor name(s). Unable to set active for visualization:", tensor_names
    # now set cell data visualization inactive.
    outpd.GetCellData().SetActiveScalars(None)

    # loop over lines
    inpd.GetLines().InitTraversal()
    outlines.InitTraversal()

    for lidx in range(0, inpd.GetNumberOfLines()):
        inpd.GetLines().GetNextCell(ptids)

        if fiber_mask[lidx]:

            if verbose:
                if lidx % 100 == 0:
                    print "<filter.py> Line:", lidx, "/", inpd.GetNumberOfLines()

            # get points for each ptid and add to output polydata
            cellptids = vtk.vtkIdList()

            for pidx in range(0, ptids.GetNumberOfIds()):
                point = inpoints.GetPoint(ptids.GetId(pidx))
                idx = outpoints.InsertNextPoint(point)
                cellptids.InsertNextId(idx)
                if preserve_point_data:
                    for idx in point_data_array_indices:
                        array = inpointdata.GetArray(idx)
                        if array.GetName() == 'region_label':
                            out_array = outpointdata.GetArray(0)
                            out_array.InsertNextTuple(array.GetTuple(ptids.GetId(pidx)))

            outlines.InsertNextCell(cellptids)

            if color is not None:
                # this code works with either 3 or 1 component only
                if outcolors.GetNumberOfComponents() == 3:
                    outcolors.InsertNextTuple3(color[lidx, 0], color[lidx, 1], color[lidx, 2])
                else:
                    outcolors.InsertNextTuple1(color[lidx])

            if preserve_cell_data:
                for idx in cell_data_array_indices:
                    array = incelldata.GetArray(idx)
                    out_array = outcelldata.GetArray(idx)
                    out_array.InsertNextTuple(array.GetTuple(lidx))

    # put data into output polydata
    outpd.SetLines(outlines)
    outpd.SetPoints(outpoints)

    # if we had an input color requested during masking, set that to be the default scalar for vis
    if color is not None:
        outpd.GetCellData().SetScalars(outcolors)

    if verbose:
        print "<filter.py> Fibers sampled:", outpd.GetNumberOfLines(), "/", inpd.GetNumberOfLines()

    return outpd

if not os.path.isdir(args.outDir):
    os.mkdir(args.outDir)

if os.path.isdir(args.inputDir):
    pd_paths = wma.io.list_vtk_files(args.inputDir)
    for pd_path in pd_paths:

        pd_input = wma.io.read_polydata(pd_path)
        vtk_name = os.path.split(pd_path)[1]
        print vtk_name

        num_lines = pd_input.GetNumberOfLines()
        fiber_mask = numpy.ones(num_lines)
        outpd = mask_remove(pd_input, fiber_mask, preserve_point_data=True, preserve_cell_data=False, verbose=False)

        wma.io.write_polydata(outpd, os.path.join(args.outDir, vtk_name))
else:

    vtk_name = os.path.split(args.inputDir)[1]
    if not os.path.exists(os.path.join(args.outDir, vtk_name)):
        pd_input = wma.io.read_polydata(args.inputDir)

        print vtk_name

        num_lines = pd_input.GetNumberOfLines()
        fiber_mask = numpy.ones(num_lines)
        outpd = mask_remove(pd_input, fiber_mask, preserve_point_data=True, preserve_cell_data=False, verbose=False)

        wma.io.write_polydata(outpd, os.path.join(args.outDir, vtk_name))

