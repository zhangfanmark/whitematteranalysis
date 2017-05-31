import vtk
import slicer
import numpy
import glob
import os
import argparse

def rasterizeFibers(fiberNode, labelNode, labelValue=10, samplingDistance=0.1):
    """Trace through the given fiber bundles and
    set the corresponding pixels in the given labelNode volume"""
    #print('rasterizing...')
    rasToIJK = vtk.vtkMatrix4x4()
    labelNode.GetRASToIJKMatrix(rasToIJK)
    labelArray = slicer.util.array(labelNode.GetID())
    polyData = fiberNode.GetPolyData()
    cellCount = polyData.GetNumberOfCells()

    selection = vtk.vtkSelection()
    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
    selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)

    extractor = vtk.vtkExtractSelectedPolyDataIds()
    extractor.SetInputData(0, polyData)

    resampler = vtk.vtkPolyDataPointSampler()
    resampler.GenerateEdgePointsOn()
    resampler.GenerateVertexPointsOff()
    resampler.GenerateInteriorPointsOff()
    resampler.GenerateVerticesOff()
    resampler.SetDistance(samplingDistance)

    cellIds = vtk.vtkIdTypeArray()

    # as a compromise between memory blowup (for large fiberbundles) and not
    # eating time in the python loop, resample CELLS_TO_PROCESS at once.
    CELLS_TO_PROCESS = 100

    for startCellId in xrange(0, cellCount, CELLS_TO_PROCESS):
        # generate list of cell Ids to sample
        cellIdsAr = numpy.arange(startCellId, min(cellCount, startCellId + CELLS_TO_PROCESS))
        cellIds.SetVoidArray(cellIdsAr, cellIdsAr.size, 1)

        # update the selection node to extract those cells
        selectionNode.SetSelectionList(cellIds)
        selection.AddNode(selectionNode)
        selection.Modified()
        extractor.SetInputData(1, selection)
        extractor.Update()
        resampler.SetInputConnection(extractor.GetOutputPort())
        resampler.Update()

        # get the resampler output
        # note: to test without resampling, just use
        # `pdSubset = extractor.GetOutput()` instead
        pdSubset = resampler.GetOutput()
        pointCount = pdSubset.GetNumberOfPoints()
        points = pdSubset.GetPoints()

        for pointIndex in xrange(pointCount):
            point = points.GetPoint(pointIndex)
            ijkFloat = rasToIJK.MultiplyPoint(point + (1,))[:3]

            # skip any negative indices to avoid wrap-around
            # to the other side of the image.
            if (any(coord < 0 for coord in ijkFloat)):
                continue

            ijk = [int(round(element)) for element in ijkFloat]
            ijk.reverse()
            try:
                # paint the voxels
                labelArray[tuple(ijk)] = labelValue
            except IndexError:
                pass
        # reset the selection node
        selection.RemoveAllNodes()

    labelNode.GetImageData().Modified()
    labelNode.Modified()
    #print('finished')


def labelStatistics(labelNode):
    grayscaleNode = labelNode

    cubicMMPerVoxel = reduce(lambda x, y: x * y, labelNode.GetSpacing())
    ccPerCubicMM = 0.001

    labelStats = {}
    labelStats['Labels'] = []

    stataccum = vtk.vtkImageAccumulate()
    stataccum.SetInputConnection(labelNode.GetImageDataConnection())
    stataccum.Update()
    lo = int(stataccum.GetMin()[0])
    hi = int(stataccum.GetMax()[0])

    for i in xrange(lo, hi + 1):

        # this->SetProgress((float)i/hi);
        # std::string event_message = "Label "; std::stringstream s; s << i; event_message.append(s.str());
        # this->InvokeEvent(vtkLabelStatisticsLogic::LabelStatsOuterLoop, (void*)event_message.c_str());

        # logic copied from slicer3 LabelStatistics
        # to create the binary volume of the label
        # //logic copied from slicer2 LabelStatistics MaskStat
        # // create the binary volume of the label
        thresholder = vtk.vtkImageThreshold()
        thresholder.SetInputConnection(labelNode.GetImageDataConnection())
        thresholder.SetInValue(1)
        thresholder.SetOutValue(0)
        thresholder.ReplaceOutOn()
        thresholder.ThresholdBetween(i, i)
        thresholder.SetOutputScalarType(grayscaleNode.GetImageData().GetScalarType())
        thresholder.Update()

        # this.InvokeEvent(vtkLabelStatisticsLogic::LabelStatsInnerLoop, (void*)"0.25");

        #  use vtk's statistics class with the binary labelmap as a stencil
        stencil = vtk.vtkImageToImageStencil()
        stencil.SetInputConnection(thresholder.GetOutputPort())
        stencil.ThresholdBetween(1, 1)

        # this.InvokeEvent(vtkLabelStatisticsLogic::LabelStatsInnerLoop, (void*)"0.5")

        stat1 = vtk.vtkImageAccumulate()
        stat1.SetInputConnection(grayscaleNode.GetImageDataConnection())
        stencil.Update()
        stat1.SetStencilData(stencil.GetOutput())

        stat1.Update()

        # this.InvokeEvent(vtkLabelStatisticsLogic::LabelStatsInnerLoop, (void*)"0.75")

        if stat1.GetVoxelCount() > 0:
            # add an entry to the LabelStats list
            labelStats["Labels"].append(i)
            labelStats[i, "Index"] = i
            labelStats[i, "Count"] = stat1.GetVoxelCount()
            labelStats[i, "Volume mm^3"] = labelStats[i, "Count"] * cubicMMPerVoxel
            labelStats[i, "Volume cc"] = labelStats[i, "Volume mm^3"] * ccPerCubicMM
            labelStats[i, "Min"] = stat1.GetMin()[0]
            labelStats[i, "Max"] = stat1.GetMax()[0]
            labelStats[i, "Mean"] = stat1.GetMean()[0]
            labelStats[i, "StdDev"] = stat1.GetStandardDeviation()[0]

    return labelStats

def list_vtk_files(input_dir):
    # Find input files
    input_mask = "{0}/*.vtk".format(input_dir)
    input_mask2 = "{0}/*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_pd_fnames = sorted(input_pd_fnames)
    return(input_pd_fnames)

cwd = os.getcwd()
cluster_dir = cwd

p = cluster_dir.find('-ukftrack')
sub_id = cluster_dir[p-6:p]

label_map_path = '/data/lmi/projects/HCP/BrainMask/fs_label_mask/'+sub_id+'_wmparc_1.25mm.nrrd'
#cluster_dir = '/Users/fan/RemoteSSHFS/HCP/Tractography/UKF_2T_fs_label_mask/clusters_l40_k0800_f10k/cluster_atlas_01_00002_remove_outliers/100307-ukftrack_b3000_fsmask_421a7ad_minGA0_trans_trans.06_minFA0_trans_trans.08_seedFALimit0_trans_trans.1_with_region_trans_trans_pp_outlier_removed/'

#checkSuccessM, volumeNode = slicer.util.loadVolume(label_map_path, returnNode=True)

pd_cluster_paths = list_vtk_files(cluster_dir)

outstr = 'cluster index' + '\t' + 'count' + '\t' + 'Volume mm^3' + '\n'

for pd_cluster_path in pd_cluster_paths:
    print os.path.split(pd_cluster_path)[1],
    checkSuccessM, fiberNode = slicer.util.loadModel(pd_cluster_path, returnNode=True)
    checkSuccessM, labelNode = slicer.util.loadLabelVolume(label_map_path, returnNode=True)

    try:
        print 'rasterize...',
        rasterizeFibers(fiberNode, labelNode)
        print 'statistics...'
        labelStats = labelStatistics(labelNode)
        count = labelStats[10, "Count"]
        volume = labelStats[10, "Volume mm^3"]
    except:
        count = 0
        volume = 0

    outstr = outstr + os.path.split(pd_cluster_path)[1] + '\t' + str(count) + '\t' + str(volume) + '\n'

    slicer.mrmlScene.RemoveNode(fiberNode)
    slicer.mrmlScene.RemoveNode(labelNode)

output_file = open(os.path.join(cluster_dir, 'volume_cluster.txt'), 'w')
output_file.write(outstr)
output_file.close()