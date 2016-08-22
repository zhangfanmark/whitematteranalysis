import argparse
import os
import glob
import slicer
import vtk

parser = argparse.ArgumentParser(
    description="Compute volume statistics with Slicer.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'inputLabelmapDirectory',
    help='')
parser.add_argument(
    'inputVolumeDirectory',
    help='')
parser.add_argument(
    'outputFile',
    help='')
parser.add_argument(
    'region', type=int,
    help='')

args = parser.parse_args()

input_labelmap_dir = os.path.abspath(args.inputLabelmapDirectory)
input_volume_dir = os.path.abspath(args.inputVolumeDirectory)
output_csv = os.path.abspath(args.outputFile)
region_index = args.region

def list_files(input_dir):

    mask = "{0}/*".format(input_dir)
    fnames = glob.glob(mask)
    fnames = sorted(fnames)

    return fnames

labelmap_paths = list_files(input_labelmap_dir)
volume_paths = list_files(input_volume_dir)

keys = ("Subject", "Index", "Count", "Volume mm^3", "Volume cc", "Min", "Max", "Mean", "StdDev")
labelStats = {}
labelStats['Labels'] = []
sub_count = 0
for labelmap_path, volume_path in zip(labelmap_paths, volume_paths):
    print "======", labelmap_path, ' <TO> ', volume_path

    check_load, labelmap_node = slicer.util.loadLabelVolume(labelmap_path, {}, True)
    check_load, volume_node = slicer.util.loadVolume(volume_path, {}, True)

    ## logic as in Slicer LabelStatistics.py
    cubicMMPerVoxel = reduce(lambda x,y: x*y, labelmap_node.GetSpacing())
    ccPerCubicMM = 0.001

    stataccum = vtk.vtkImageAccumulate()
    stataccum.SetInputConnection(labelmap_node.GetImageDataConnection())
    stataccum.Update()
    lo = int(stataccum.GetMin()[0])
    hi = int(stataccum.GetMax()[0])

    if region_index > hi or region_index < lo:
        print "Error: Input region index ", region_index, "does not exist."

    i = region_index

    thresholder = vtk.vtkImageThreshold()
    thresholder.SetInputConnection(labelmap_node.GetImageDataConnection())
    thresholder.SetInValue(1)
    thresholder.SetOutValue(0)
    thresholder.ReplaceOutOn()
    thresholder.ThresholdBetween(i,i)
    thresholder.SetOutputScalarType(volume_node.GetImageData().GetScalarType())
    thresholder.Update()

    #  use vtk's statistics class with the binary labelmap as a stencil
    stencil = vtk.vtkImageToImageStencil()
    stencil.SetInputConnection(thresholder.GetOutputPort())
    stencil.ThresholdBetween(1, 1)

    stat1 = vtk.vtkImageAccumulate()
    stat1.SetInputConnection(volume_node.GetImageDataConnection())
    stencil.Update()
    stat1.SetStencilData(stencil.GetOutput())

    stat1.Update()

    if stat1.GetVoxelCount() > 0:
        # add an entry to the LabelStats list
        labelStats["Labels"].append(sub_count)
        labelStats[sub_count,"Subject"] = os.path.split(volume_path)[1]
        labelStats[sub_count,"Index"] = i
        labelStats[sub_count,"Count"] = stat1.GetVoxelCount()
        labelStats[sub_count,"Volume mm^3"] = labelStats[sub_count,"Count"] * cubicMMPerVoxel
        labelStats[sub_count,"Volume cc"] = labelStats[sub_count,"Volume mm^3"] * ccPerCubicMM
        labelStats[sub_count,"Min"] = stat1.GetMin()[0]
        labelStats[sub_count,"Max"] = stat1.GetMax()[0]
        labelStats[sub_count,"Mean"] = stat1.GetMean()[0]
        labelStats[sub_count,"StdDev"] = stat1.GetStandardDeviation()[0]
        sub_count = sub_count + 1

    slicer.mrmlScene.RemoveNode(labelmap_node)
    slicer.mrmlScene.RemoveNode(volume_node)

# output csv
csv = ""
header = ""
for k in keys[:-1]:
    header += "\"%s\"" % k + ","
header += "\"%s\"" % keys[-1] + "\n"
csv = header
for sub_count in labelStats["Labels"]:
    line = ""
    for k in keys[:-1]:
        line += str(labelStats[sub_count,k]) + ","
    line += str(labelStats[sub_count,keys[-1]]) + "\n"
    csv += line

fp = open(output_csv, "w")
fp.write(csv)
fp.close()