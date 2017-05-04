#!/usr/bin/env python
import argparse
import os
import glob
import vtk
import numpy
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
    description="EndPointAnalysis (EPA) using the labels already stored in the fiber cluster vtk/vtp files.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")

parser.add_argument(
    'inputDirectory',
    help='Directory of the subject-specific fiber clustering results.')
parser.add_argument(
    'outputDirectory',
    help='The output directory will be created if it does not exist.')
parser.add_argument(
    'regionList', type=int, nargs='+',
    help='One or two regions to be analyzed, such as 8 1024')
parser.add_argument(
    '-hemi', type=str, dest="hemi", default="whole",
    help='Hemispheres or whole brain fiber clusters to be analyzed.')
parser.add_argument(
    '-fn', type=str, dest="folder_name", default="",
    help='String in the folder name of the subject-specific fiber clustering')
parser.add_argument(
    '-cp_mrml', action='store_true', dest="flag_cp_mrml",
    help='If given, a mrml will be copied to the fiber clustering folder of each subject. Otherwise, only one mrml will be generated under the input directory. ')
parser.add_argument(
    '-j', action="store", dest="numberOfJobs", type=int,
    help='Number of processors to use.')

args = parser.parse_args()

if not os.path.isdir(args.inputDirectory):
    print "Error: Input directory", args.inputDirectory, "does not exist."
    exit()

outdir = args.outputDirectory
if not os.path.exists(outdir):
    print "Output directory", outdir, "does not exist, creating it."
    os.makedirs(outdir)

if len(args.regionList) > 2:
    print " At most two regions can analyzed."
    exit()

hemi = ''
if args.hemi == "comm":
    hemi = 'tracts_commissural'
elif args.hemi == "left":
    hemi = 'tracts_left_hemisphere'
elif args.hemi == 'right':
    hemi = 'tracts_right_hemisphere'
elif args.hemi == 'whole':
    hemi = ''
else:
    print "Error: Hemisphere or whole brain: ", args.hemi, "does not exist."
    exit()

if args.numberOfJobs is not None:
    number_of_jobs = args.numberOfJobs
else:
    number_of_jobs = 1

flag_cp_mrml = args.flag_cp_mrml

output_file_path = os.path.join(args.inputDirectory, 'EPA_regions_per_cluster_in_'+hemi+'.txt')

load_previous = False
if os.path.isfile(output_file_path):
    print "\n<EndPointAnalysis> Result from a previous run is detected in", output_file_path
    print "  Analysis will be conducted with this result. To re-run, please delete this file."
    load_previous = True

if not load_previous:

    def list_fc_folders(input_dir, hemi, folder_name):
        # Find input files
        input_mask = ("{0}/*"+folder_name+"*/"+hemi).format(input_dir)
        input_fc_folders = glob.glob(input_mask)
        input_fc_folders = sorted(input_fc_folders)
        return (input_fc_folders)

    fc_folders = list_fc_folders(args.inputDirectory, hemi, args.folder_name)

    num_subjects = len(fc_folders)
    num_clusters = 0
    flag_first = True
    for fc_folder in fc_folders:
        pd_fc_paths = wma.io.list_vtk_files(fc_folder)
        if flag_first:
            num_clusters = len(pd_fc_paths)
            flag_first = False
        else:
            if num_clusters != len(pd_fc_paths):
                print "Error: All subjects should have the same number of clusters"
                print len(pd_fc_paths), 'clusters are found in', fc_folder, 'while', num_clusters, 'clusters were deteced before.'
                exit()

    def get_endpoint_labels_per_cluster(pd_fc):

        pd_fc.GetLines().InitTraversal()
        line_ptids = vtk.vtkIdList()
        inpointdata = pd_fc.GetPointData()

        point_data_array_indices = range(inpointdata.GetNumberOfArrays())
        label_array = None
        for idx in point_data_array_indices:
            label_array = inpointdata.GetArray(idx)
            if label_array.GetName().lower() == 'region_label':
                break

        eps_all = []
        num_fibers = pd_fc.GetNumberOfLines()

        for lidx in range(0, num_fibers):
            pd_fc.GetLines().GetNextCell(line_ptids)
            line_length = line_ptids.GetNumberOfIds()

            ep_1 = [0]
            ep_2 = [0]
            label_array.GetTupleValue(line_ptids.GetId(0), ep_1)
            label_array.GetTupleValue(line_ptids.GetId(line_length-1), ep_2)

            ep_1 = ep_1[0]
            ep_2 = ep_2[0]

            if ep_1 > 3000:
                ep_1 = ep_1 - 2000
            if ep_2 > 3000:
                ep_2 = ep_2 - 2000

            eps_all.append(ep_1)
            eps_all.append(ep_2)

        #eps_all = eps_all(numpy.where(eps_all >=0)) # remove the points that outside brain

        tmp = numpy.argsort(numpy.bincount(eps_all))[-2:] # The two points that mostly appear

        try:
            ep_1_label = min(tmp)
            ep_1_label_percent = numpy.sum(numpy.sum(eps_all == ep_1_label)) / float(len(eps_all))
        except:
            ep_1_label = 0
            ep_1_label_percent = 0

        try:
            ep_2_label = max(tmp)
            ep_2_label_percent = numpy.sum(numpy.sum(eps_all == ep_2_label)) / float(len(eps_all))
        except:
            ep_2_label = 0
            ep_2_label_percent = 0

        return (ep_1_label, ep_2_label, ep_1_label_percent, ep_2_label_percent)

    def get_endpoint_labels_per_subject(fc_folder, verbose=False):

        print '  Subject-specific clustering result: ', fc_folder
        pd_fc_paths = wma.io.list_vtk_files(fc_folder)

        ep_1_label_all_clusters = []
        ep_2_label_all_clusters = []
        ep_1_label_percent_all_clusters = []
        ep_2_label_percent_all_clusters = []

        # count = 0
        for pd_fc_path in pd_fc_paths:

            if verbose:
                print '  - Cluster', os.path.split(pd_fc_path)[1]

            # count = count + 1
            # if count > 20:
            #     break

            pd_fc = wma.io.read_polydata(pd_fc_path)
            ep_1_label, ep_2_label, ep_1_label_percent, ep_2_label_percent = get_endpoint_labels_per_cluster(pd_fc)

            ep_1_label_all_clusters.append(ep_1_label)
            ep_2_label_all_clusters.append(ep_2_label)
            ep_1_label_percent_all_clusters.append(ep_1_label_percent)
            ep_2_label_percent_all_clusters.append(ep_2_label_percent)

        # concatenate these four list for Parallel computing
        return numpy.concatenate((ep_1_label_all_clusters, ep_2_label_all_clusters, ep_1_label_percent_all_clusters, ep_2_label_percent_all_clusters))

    print '\n<EndPointAnalysis> A total of', num_subjects, 'subjects, in which each has', num_clusters, 'clusters'

    ep_label_all_clusters_all_subjects_concatenate = \
        Parallel(n_jobs=number_of_jobs, verbose=1)(
            delayed(get_endpoint_labels_per_subject)(fc_folder)
            for fc_folder in fc_folders)

    ep_label_all_clusters_all_subjects_concatenate = numpy.array(ep_label_all_clusters_all_subjects_concatenate)

    num_clusters_tmp = ep_label_all_clusters_all_subjects_concatenate.shape[1] / 4
    num_clusters = num_clusters_tmp

    ep_1_label_all_clusters_all_subjects = ep_label_all_clusters_all_subjects_concatenate[:, 0:num_clusters]
    ep_2_label_all_clusters_all_subjects = ep_label_all_clusters_all_subjects_concatenate[:, num_clusters:2*num_clusters]
    ep_1_label_all_clusters_all_subjects_percent = ep_label_all_clusters_all_subjects_concatenate[:, 2*num_clusters:3*num_clusters]
    ep_2_label_all_clusters_all_subjects_percent = ep_label_all_clusters_all_subjects_concatenate[:, 3*num_clusters:4*num_clusters]

    output_file = open(output_file_path, 'w')

    outstr = 'cluster index' + '\t'+ 'region-1' + '\t'+ 'region-1-percentage' + '\t'+ 'region-2' + '\t'+ 'region-2-percentage' + '\n'
    for c_idx in range(num_clusters):

        ep_1_label_per_cluster_all_subjects_tmp = ep_1_label_all_clusters_all_subjects[:, c_idx]
        ep_2_label_per_cluster_all_subjects_tmp = ep_2_label_all_clusters_all_subjects[:, c_idx]

        ep_1_label_per_cluster_all_subjects_tmp = ep_1_label_per_cluster_all_subjects_tmp.astype(numpy.int64)
        ep_2_label_per_cluster_all_subjects_tmp = ep_2_label_per_cluster_all_subjects_tmp.astype(numpy.int64)

        ep_both_label_per_cluster_all_subjects_tmp = numpy.concatenate((ep_1_label_per_cluster_all_subjects_tmp, ep_2_label_per_cluster_all_subjects_tmp))

        ep_label_occurrence = numpy.bincount(ep_both_label_per_cluster_all_subjects_tmp)

        top_two_labels = numpy.argsort(ep_label_occurrence)[-2:]

        if len(top_two_labels) == 2:
            ep_1_label_per_cluster_final = top_two_labels[1]
            ep_2_label_per_cluster_final = top_two_labels[0]

            num_subjects_with_ep_1_final = 0
            for s_idx in range(num_subjects):
                if ep_1_label_per_cluster_all_subjects_tmp[s_idx] == ep_1_label_per_cluster_final or \
                   ep_2_label_per_cluster_all_subjects_tmp[s_idx] == ep_1_label_per_cluster_final:
                    num_subjects_with_ep_1_final = num_subjects_with_ep_1_final + 1

            num_subjects_with_ep_2_final = 0
            for s_idx in range(num_subjects):
                if ep_1_label_per_cluster_all_subjects_tmp[s_idx] == ep_2_label_per_cluster_final or \
                   ep_2_label_per_cluster_all_subjects_tmp[s_idx] == ep_2_label_per_cluster_final:
                    num_subjects_with_ep_2_final = num_subjects_with_ep_2_final + 1

            ep_1_label_occurrence = ep_label_occurrence[ep_1_label_per_cluster_final]

            if ep_1_label_occurrence - num_subjects_with_ep_1_final > num_subjects_with_ep_2_final:
                ep_2_label_per_cluster_final = ep_1_label_per_cluster_final
                num_subjects_with_ep_2_final = ep_1_label_occurrence - num_subjects_with_ep_1_final

            ep_1_label_per_cluster_final_percent = num_subjects_with_ep_1_final / float(num_subjects)
            ep_2_label_per_cluster_final_percent = num_subjects_with_ep_2_final / float(num_subjects)
        else:
            ep_2_label_per_cluster_final = 0
            ep_2_label_per_cluster_final = 0
            ep_1_label_per_cluster_final_percent = 0
            ep_2_label_per_cluster_final_percent = 0


        # try:
        #     ep_1_label_per_cluster_final = numpy.argmax(numpy.bincount(ep_1_label_per_cluster_all_subjects))
        # except:
        #     ep_1_label_per_cluster_final = 0
        #
        # try:
        #     ep_2_label_per_cluster_final = numpy.argmax(numpy.bincount(ep_2_label_per_cluster_all_subjects))
        # except:
        #     ep_2_label_per_cluster_final = 0
        #
        # ep_1_label_per_cluster_final_percent = numpy.sum(numpy.sum(ep_1_label_per_cluster_final == ep_1_label_per_cluster_all_subjects)) / float(len(ep_2_label_per_cluster_all_subjects))
        # ep_2_label_per_cluster_final_percent = numpy.sum(numpy.sum(ep_2_label_per_cluster_final == ep_2_label_per_cluster_all_subjects)) / float(len(ep_2_label_per_cluster_all_subjects))

        outstr = outstr + 'cluster_{0:05d}.vtp'.format(c_idx+1) + '\t' + str(ep_1_label_per_cluster_final) + '\t' + str(ep_1_label_per_cluster_final_percent) + '\t' + str(
            ep_2_label_per_cluster_final) + '\t' + str(ep_2_label_per_cluster_final_percent) + '\n'

    output_file.write(outstr)
    output_file.close()

result_csv = numpy.genfromtxt(output_file_path, delimiter='\t')

num_clusters = result_csv.shape[0] - 1
print '\n<EndPointAnalysis> Number of clusters:', num_clusters

region_list_all = []
if len(args.regionList) == 1 and args.regionList[0] == -1:
    all_region_pairs = []
    for c_idx in range(num_clusters):
        ep_1_label_per_cluster_final = result_csv[c_idx + 1, 1]
        ep_2_label_per_cluster_final = result_csv[c_idx + 1, 3]

        if ep_1_label_per_cluster_final > ep_2_label_per_cluster_final:
            all_region_pairs.append([ep_2_label_per_cluster_final, ep_1_label_per_cluster_final])
        else:
            all_region_pairs.append([ep_1_label_per_cluster_final, ep_2_label_per_cluster_final])

    region_list_all =  numpy.array([numpy.array(x) for x in set(tuple(x) for x in all_region_pairs)])
else:
    region_list_all.append(args.regionList)


for region_list in region_list_all:
    region_list = region_list.tolist()
    result_cluster_list = []
    print '\n<EndPointAnalysis> The endpoint regions connected per cluster, calculated across all subjects.'
    print '  The percentage shows how many subjects have their endpoints connecting the detected region.'
    print "Cluster index:  region-1 (percentage)  region-2 (percentage)"
    print "------------------------------------------------------------"
    for c_idx in range(num_clusters):
        ep_1_label_per_cluster_final = result_csv[c_idx + 1, 1]
        ep_1_label_per_cluster_final_percent = result_csv[c_idx + 1, 2]
        ep_2_label_per_cluster_final = result_csv[c_idx + 1, 3]
        ep_2_label_per_cluster_final_percent = result_csv[c_idx + 1, 4]

        # print "cluster_%05d:  %8d (%10f) %8d (%10f)" % (c_idx + 1, ep_1_label_per_cluster_final, ep_1_label_per_cluster_final_percent, ep_2_label_per_cluster_final,
        #       ep_2_label_per_cluster_final_percent)

        if len(region_list) == 1:
            if ep_1_label_per_cluster_final == region_list[0] or ep_2_label_per_cluster_final == region_list[0]:
                result_cluster_list.append(c_idx)
        elif len(region_list)==2:
            if (ep_1_label_per_cluster_final == region_list[0] and ep_2_label_per_cluster_final == region_list[1]) or \
               (ep_1_label_per_cluster_final == region_list[1] and ep_2_label_per_cluster_final == region_list[0]):
                result_cluster_list.append(c_idx)

    print '\n<EndPointAnalysis> Find', len(result_cluster_list), 'cluster(s) connecting region(s)', region_list

    if len(result_cluster_list) > 0:

        outstr = ' - cluster_'
        for result_cluster in result_cluster_list:
            outstr = outstr + '{0:05d}'.format(result_cluster+1) + ', '
        print outstr[:-2]

        mrml_filename = "EPA_clusters_in_"+hemi+"_connecting_region_" + str(region_list).replace(', ', '-') + ".mrml"
        print "\n<EndPointAnalysis> A mrml file to display the result cluster(s) is generated as", mrml_filename
        print '  Copy this file to a folder containing the result and load it to 3D Slicer to display the fiber clusters.'

        cluster_polydatas = []
        for c_idx in result_cluster_list:
            cluster_polydatas.append("cluster_"+str(c_idx+1).zfill(5)+".vtp")

        number_of_files = len(cluster_polydatas)
        step = int(100 * 255.0 / (number_of_files))
        R = numpy.array(range(0, 100 * 255 + 1, step)) / 100.0
        G = numpy.abs(range(100 * -127, 100 * 128 + 1, step)) * 2.0 / 100.0
        B = numpy.array(range(100 * 255 + 1, 0, -step)) / 100.0

        colors = list()
        idx = 0
        for pd in cluster_polydatas:
            colors.append([R[idx], G[idx], B[idx]])
            idx += 1
        colors = numpy.array(colors)

        wma.mrml.write(cluster_polydatas, colors, os.path.join(outdir, mrml_filename), ratio=1.0)

        if flag_cp_mrml:
            print 'Not yet implemented.'

print '\n<EndPointAnalysis> Done! \n'