#!/usr/bin/env python
import argparse
import os

#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Add base path to MRML file",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'inputMRML',
    help='Input mrml file')
parser.add_argument(
    'outputDir',
    help='Output mrml file')
parser.add_argument(
    '-outputstr', type=str, default='output',
    help='Output mrml file')
parser.add_argument(
    '-basepath', type=str, default='',
    help='Base path that will be added before the tag fileName in <FiberBundleStorage>')
parser.add_argument(
    '-hierarchy', type=str, default='hierarchy',
    help='The clusters in the MRML will be assigned to one hierarchy so that multiple MRML files could be loaded together')

args = parser.parse_args()

if not os.path.isfile(args.inputMRML):
    print "<wm_fan_add_basepath_to_mrml> Error: Input MRML", args.inputMRML, "does not exist."
    exit()

outdir = args.outputDir
if not os.path.exists(outdir):
    print "<wm_fan_add_basepath_to_mrml> Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

input_file_name = os.path.split(args.inputMRML)[1]

output_file_name = input_file_name[:-5] + '_' + args.outputstr + '.mrml'

hierarchy = args.hierarchy

outstr =         '<MRML  version="Slicer4" userTags="">'
outstr += '\n' + ' <ModelDisplay'
outstr += '\n' + '  id="vtkMRMLModelDisplayNode' + hierarchy + '"  name="ModelDisplay"  hideFromEditors="true"  selectable="true"  selected="false"  color="0.5 0.5 0.5"  edgeColor="0 0 0"  selectedColor="1 0 0"  selectedAmbient="0.4"  ambient="0"  diffuse="1"  selectedSpecular="0.5"  specular="0"  power="1"  opacity="1"  pointSize="1"  lineWidth="1"  representation="2"  lighting="true"  interpolation="1"  shading="true"  visibility="true"  edgeVisibility="false"  clipping="false"  sliceIntersectionVisibility="false"  sliceIntersectionThickness="1"  frontfaceCulling="false"  backfaceCulling="true"  scalarVisibility="false"  vectorVisibility="false"  tensorVisibility="false"  interpolateTexture="false"  scalarRangeFlag="2"  autoScalarRange="true"  scalarRange="0 100"  ></ModelDisplay>'
outstr += '\n' + ' <ModelHierarchy'
outstr += '\n' + '  id="vtkMRMLModelHierarchyNode' + hierarchy + '"  name="' + hierarchy + '"  hideFromEditors="false"  selectable="true"  selected="false"  sortingValue="0"  allowMultipleChildren="true"  displayNodeID="vtkMRMLModelDisplayNode' + hierarchy + '"  expanded="true"  ></ModelHierarchy>'
outstr += '\n'

with open (args.inputMRML) as input_mrml:
    for line in input_mrml:

        if line.find('<MRML  version="Slicer4" userTags="">') >=0 or line.find('</MRML>') >=0: # ignore the fist and last lines
            continue

        k_f = line.find('fileName')
        if k_f >= 0:
            line = line[:k_f+10] + args.basepath + line[k_f+10:]

        line = line.replace("vtkMRMLFiberBundleStorageNode", 'vtkMRMLFiberBundleStorageNode' + hierarchy)
        line = line.replace("vtkMRMLFiberBundleLineDisplayNode", 'vtkMRMLFiberBundleLineDisplayNode' + hierarchy)
        line = line.replace("vtkMRMLFiberBundleTubeDisplayNode", 'vtkMRMLFiberBundleTubeDisplayNode' + hierarchy)
        line = line.replace("vtkMRMLFiberBundleGlyphDisplayNode", 'vtkMRMLFiberBundleGlyphDisplayNode' + hierarchy)
        line = line.replace("vtkMRMLFiberBundleNode", 'vtkMRMLFiberBundleNode' + hierarchy)
        line = line.replace("vtkMRMLDiffusionTensorDisplayPropertiesNode", 'vtkMRMLDiffusionTensorDisplayPropertiesNode' + hierarchy)
        line = line.replace("vtkMRMLModelHierarchyNode", 'vtkMRMLModelHierarchyNode' + hierarchy)

        outstr += line

        k_b = line.find('vtkMRMLFiberBundleNode')
        if k_b >= 0:
            k_n = line.find('name')
            vtkMRMLFiberBundleNode_id = line[k_b + len('vtkMRMLFiberBundleNode'):k_n - 3]

            outstr += ' <ModelHierarchy'
            outstr += '\n' + '  id="vtkMRMLModelHierarchyNode' + vtkMRMLFiberBundleNode_id + '"  name="ModelHierarchy"  hideFromEditors="true"  selectable="true"  selected="false"  parentNodeRef="vtkMRMLModelHierarchyNode' + hierarchy + '"  associatedNodeRef="vtkMRMLFiberBundleNode' + vtkMRMLFiberBundleNode_id + '"  sortingValue="1"  allowMultipleChildren="true"  expanded="true"  ></ModelHierarchy>'
            outstr += '\n'

outstr += '</MRML>'

output_file = open(os.path.join(args.outputDir, output_file_name), 'w')
output_file.write(outstr)
output_file.close()


print '<wm_fan_add_basepath_to_mrml> A mrml file to load atlas clusters is generated as', os.path.join(args.outputDir, output_file_name), '\n'









