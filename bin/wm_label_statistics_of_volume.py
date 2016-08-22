import argparse
import os

parser = argparse.ArgumentParser(
	description="Compute label statistics of the input volume given a labelmap with Slicer.",
	epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
parser.add_argument("-v", "--version",
	action="version", default=argparse.SUPPRESS,
	version='1.0',
	help="Show program's version number and exit")

parser.add_argument(
	'inputLabelmapDirectory',
	help='Directory of input labelmaps (e.g., freesurfer parcellation). The script will calculate label statistics of the input volume (as in LabelStatistics module in 3D Slicer) given one region of the input label map.')
parser.add_argument(
	'inputVolumeDirectory',
	help='Directory of input image volumes (e.g., T1, T2 or B0 image). Note: Image volume and the corresponding labelmap must have same dimensions.')
parser.add_argument(
	'outputFile',
	help='A CVS file that contains the output statistics.')
parser.add_argument(
	'Slicer',
	help='Path of 3D Slicer.')
parser.add_argument(
	'-r',  dest='region', type=int,
	help='Region index.')

args = parser.parse_args()

input_labelmap_dir = os.path.abspath(args.inputLabelmapDirectory)
input_volume_dir = os.path.abspath(args.inputVolumeDirectory)
if not os.path.isdir(input_labelmap_dir) or not os.path.isdir(input_volume_dir):
	print "<wm_label_statistics_of_volume> Error: Input directory", args.inputLabelmapDirectory, "or",  args.imputVolumeDirectory, "does not exist."
	exit()

slicer_path = os.path.abspath(args.Slicer)
if not os.path.exists(args.Slicer):
	print "Error: 3D Slicer", args.Slicer, "does not exist."
	exit()

output_csv = args.outputFile
region_index = args.region

print "<wm_label_statistics_of_volume> Starting label statistics."
print ""
print "=====input labelmap directory======\n", input_labelmap_dir
print "=====input volume directory======\n", input_volume_dir
print "=====output CSV file=====\n", output_csv
print "=====3D Slicer====\n", slicer_path

cmd = slicer_path + " --no-main-window --python-script $(which label_statistics_with_slicer.py) " + \
          input_labelmap_dir + " " + input_volume_dir + " " + output_csv + " " + str(region_index) + " --python-code 'slicer.app.quit()' "

os.system(cmd)










