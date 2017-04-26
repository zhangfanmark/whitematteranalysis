import numpy
import whitematteranalysis as wma
import os
import vtk
import scipy

def get_label_array(input_data):

    input_data.GetLines().InitTraversal()
    num_fibers = input_data.GetNumberOfLines()
    inpoints = input_data.GetPoints()
    inpointsdata = input_data.GetPointData()
    num_points = inpoints.GetNumberOfPoints()

    print 'Total fiber number in the atlas:', num_fibers, ', point number in the atlas:', num_points

    if inpointsdata.GetNumberOfArrays() > 0:
        point_data_array_indices = range(inpointsdata.GetNumberOfArrays())
        for idx in point_data_array_indices:
            array = inpointsdata.GetArray(idx)
            if array.GetName().lower() == 'region_label':
                label_array = array

    return label_array

def _by_two_endpoints(input_data, outdir, output_fibers_per_group=True):
    label_array = get_label_array(input_data)

    num_fibers = input_data.GetNumberOfLines()
    #connectivity = numpy.zeros([num_fibers, num_fibers])

    dtype = [('ep_1', int), ('ep_2', int)]
    endpoint_regions = []

    # get the endpoint regions for each fiber
    input_data.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()
    for lidx in range(0, num_fibers):
        input_data.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        ptidx_1 = line_ptids.GetId(0)
        ptidx_2 = line_ptids.GetId(line_length - 1)
        label_1 = int(label_array.GetTuple(ptidx_1)[0])
        label_2 = int(label_array.GetTuple(ptidx_2)[0])

        label_1 = region_label(label_1)
        label_2 = region_label(label_2)

        if label_1 <= label_2:
            endpoint_regions.append((label_1, label_2))
        else:
            endpoint_regions.append((label_2, label_1))

    # sort the fibers according to the endpoint label of the first point
    endpoint_regions = numpy.array(endpoint_regions, dtype=dtype)
    sort_indices = numpy.argsort(endpoint_regions, order=['ep_1', 'ep_2'])

    # assign each fiber a group label, where clusters have the same endpoint regions are with a same label
    ep_group_label = numpy.zeros([num_fibers, 1])
    eps_previous = (-1, -1)
    label_count = 0
    for si in sort_indices:
        eps_ii = endpoint_regions[si]
        if eps_ii[0] != eps_previous[0] or eps_ii[1] != eps_previous[1]:
            label_count = label_count + 1
            ep_group_label[si] = label_count
        else:
            ep_group_label[si] = label_count
        eps_previous = eps_ii

    output_file = open(os.path.join(outdir, 'num_fibers_per_region_pair.txt'), 'w')
    outstr = 'Group_Label' + '\t' + 'Region_1' + '\t' + 'Region_2' + '\t' + 'Num_Fibers\n'


    threshold_num_fiber_per_group = 10
    print '\n<EP region connectivity> Region pair (group) that has less than', threshold_num_fiber_per_group, 'fibers will be set Group_-1 that are unclassified fibers.'
    # Assign group label to -1, if there are too few fibers in this group.
    # All the group_-1 fibers are the unclassified ones.
    for ui in numpy.unique(ep_group_label):
        group_label_subjects = numpy.where(ep_group_label == ui)[0]

        ui_mask = numpy.zeros([num_fibers, 1])
        ui_mask[group_label_subjects] = 1

        num_fb = numpy.sum(ui_mask)
        if num_fb < threshold_num_fiber_per_group:
            ep_group_label[group_label_subjects] = -1

        outstr = outstr + str(ui) + '\t' + str(endpoint_regions[group_label_subjects][0][0]) + '\t' + \
                 str(endpoint_regions[group_label_subjects][0][1]) + '\t' + str(num_fb) + '\n'

    output_file.write(outstr)
    output_file.close()

    # Output how many fibers are in each group
    output_file = open(os.path.join(outdir, 'num_fibers_per_group.txt'), 'w')
    outstr = 'index' + '\t' + 'Group_Label' + '\t' + 'Region_1' + '\t' + 'Region_2' + '\t' + 'Num_Fibers\n'

    label_idx = 1
    for ui in numpy.sort(numpy.unique(ep_group_label)):
        group_label_subjects = numpy.where(ep_group_label == ui)[0]
        # for rr in group_label_subjects:
        #     connectivity[rr, group_label_subjects] = 1

        ui_mask = numpy.zeros([num_fibers, 1])
        ui_mask[group_label_subjects] = 1

        num_fb = numpy.sum(ui_mask)
        if ui == -1:
            region_1 = -1
            region_2 = -1
        else:
            region_1 = endpoint_regions[group_label_subjects][0][0]
            region_2 = endpoint_regions[group_label_subjects][0][1]
        outstr = outstr + str(label_idx) + '\t' + str(ui) + '\t' + str(region_1) + '\t' + \
                 str(region_2) + '\t' + str(num_fb) + '\n'

        if output_fibers_per_group:
            pd_ui = wma.filter.mask(input_data, ui_mask, color=None, preserve_point_data=True, preserve_cell_data=True, verbose=False)
            wma.io.write_polydata(pd_ui, os.path.join(outdir, 'Conn_I{0:04d}'.format(int(label_idx)) + '_G{0:04d}'.format(int(ui)) + '_R{0:04d}'.format(int(region_1)) + '-{0:04d}'.format(int(region_2)) + '.vtp'))
        label_idx = label_idx + 1

    '{0:05d}'

    output_file.write(outstr)
    output_file.close()

    print '<EP region connectivity> Number of valid region pairs (groups) is', len(numpy.unique(ep_group_label))
    print '  Number of clusters (-k) should be set greater than', len(numpy.unique(ep_group_label))

    #connectivity = scipy.sparse.csr_matrix(connectivity)
    connectivity = ep_group_label

    return connectivity

def _by_one_endpoint(input_data, outdir):

    label_array = get_label_array(input_data)

    num_fibers = input_data.GetNumberOfLines()
    #connectivity = numpy.zeros([num_fibers, num_fibers])
    endpoint_regions = numpy.zeros([num_fibers, 2])

    input_data.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()
    for lidx in range(0, num_fibers):
        input_data.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        ptidx_1 = line_ptids.GetId(0)
        ptidx_2 = line_ptids.GetId(line_length-1)
        label_1 = int(label_array.GetTuple(ptidx_1)[0])
        label_2 = int(label_array.GetTuple(ptidx_2)[0])

        if label_1 >= 3000 and label_1 <= 4035:
            label_1 = label_1 - 2000

        if label_2 >= 3000 and label_2 <= 4035:
            label_2 = label_2 - 2000

        if label_1 <= label_2:
            endpoint_regions[lidx, 0] = label_1
            endpoint_regions[lidx, 1] = label_2
        else:
            endpoint_regions[lidx, 0] = label_2
            endpoint_regions[lidx, 1] = label_1

    region_label_list = numpy.sort(numpy.unique(endpoint_regions))

    output_file = open(os.path.join(outdir, 'region_number.txt'), 'w')
    outstr = 'Regions' + '\t' + 'Number of Fibers\n'

    tot_num_fb = 0
    for r_idx in region_label_list:
        region_mask = numpy.zeros([num_fibers, 1])
        region_mask[numpy.where(endpoint_regions[:, 0] == r_idx)] = 1
        region_mask[numpy.where(endpoint_regions[:, 1] == r_idx)] = 1

        # for rr in numpy.where(region_mask==1)[0]:
        #     connectivity[rr, numpy.where(region_mask==1)[0]] = 1

        num_fb = len(numpy.where(region_mask==1)[0])
        print 'Region', r_idx, ',number of fibers:', num_fb
        outstr = outstr + str(r_idx) + '\t' + str(num_fb) + '\n'
        tot_num_fb = tot_num_fb + num_fb
        pd_r_idx = wma.filter.mask(input_data, region_mask, color=None, preserve_point_data=True, preserve_cell_data=True, verbose=False)
        wma.io.write_polydata(pd_r_idx, os.path.join(outdir, 'fiber_'+str(r_idx)+'.vtp'))

    output_file.write(outstr)
    output_file.close()

    print 'Total number of fibers:', tot_num_fb, '(should be 2 *', num_fibers, ')'

    #connectivity = scipy.sparse.csr_matrix(connectivity)
    connectivity = None
    return connectivity


def region_label(label):

    left_sub_cortical_regions = [2, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28]
    right_sub_cortical_regions = [41, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]

    left_GM_cortical_regions = range(1001, 1036)
    right_GM_cortical_regions = range(2001, 2036)

    left_WM_cortical_regions = range(3001, 3036)
    right_WM_cortical_regions = range(4001, 4036)

    CC_regions = range(251, 256)
    comm_sub_cortical_regions = [16, 77, 85]

    WM_Unsegmented = [5001, 5002]

    if label in right_WM_cortical_regions:
        label = label - 3000
    elif label in left_WM_cortical_regions:
        label = label - 2000
    elif label in right_GM_cortical_regions:
        label = label - 1000
    elif label in left_GM_cortical_regions:
        label = label
    elif label in CC_regions or label in comm_sub_cortical_regions:
        label = label
    elif label in right_sub_cortical_regions:
        label = left_sub_cortical_regions[right_sub_cortical_regions.index(label)]
    elif label in left_sub_cortical_regions:
        label = label
    elif label in WM_Unsegmented:
        label = 5001
    else:
        label = 0

    return label