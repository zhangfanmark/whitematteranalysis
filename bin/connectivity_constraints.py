import numpy
import whitematteranalysis as wma
import os
import vtk
import scipy

def get_endpiont_region(input_data):

    label_array = get_label_array(input_data, verbose=False)

    num_fibers = input_data.GetNumberOfLines()
    # connectivity = numpy.zeros([num_fibers, num_fibers])

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

    return endpoint_regions


def get_along_tract_region(input_data, verbose=True):

    label_array = get_label_array(input_data, verbose=False)

    num_fibers = input_data.GetNumberOfLines()

    # loop over lines
    input_data.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()

    along_tract_regions = []
    for lidx in range(0, num_fibers):

        input_data.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        # loop over the indices that we want and get those points
        region_one_fiber = []
        for line_index in range(line_length):
            # do nearest neighbor interpolation: round index
            ptidx = line_ptids.GetId(int(round(line_index)))
            array_label = label_array.GetTuple(ptidx)[0]
            array_label = region_label(array_label)

            region_one_fiber.append(array_label)

        region_one_fiber = numpy.unique(region_one_fiber)
        region_one_fiber = region_one_fiber[numpy.where(region_one_fiber>0)]

        if verbose:
            print ' --Fiber', lidx, 'regions:', region_one_fiber

        along_tract_regions.append(region_one_fiber)

    return along_tract_regions

def get_along_tract_region_OLD(input_data, points_per_fiber=12):

    label_array = get_label_array(input_data, verbose=False)

    num_fibers = input_data.GetNumberOfLines()

    along_tract_regions = numpy.zeros((num_fibers, points_per_fiber))
    # loop over lines
    input_data.GetLines().InitTraversal()
    line_ptids = vtk.vtkIdList()

    for lidx in range(0, num_fibers):

        input_data.GetLines().GetNextCell(line_ptids)
        line_length = line_ptids.GetNumberOfIds()

        # loop over the indices that we want and get those points
        pidx = 0
        for line_index in _calculate_line_indices(line_length, points_per_fiber):
            # do nearest neighbor interpolation: round index
            ptidx = line_ptids.GetId(int(round(line_index)))
            array_label = label_array.GetTuple(ptidx)[0]
            array_label = region_label(array_label)

            along_tract_regions[lidx, pidx] = array_label

            pidx = pidx + 1

    return along_tract_regions


def _calculate_line_indices(input_line_length, output_line_length):
    """ Figure out indices for downsampling of polyline data.

    The indices include the first and last points on the line,
    plus evenly spaced points along the line.  This code figures
    out which indices we actually want from a line based on its
    length (in number of points) and the desired length.

    """

    # this is the increment between output points
    step = (input_line_length - 1.0) / (output_line_length - 1.0)

    # these are the output point indices (0-based)
    ptlist = []
    for ptidx in range(0, output_line_length):
        # print(ptidx*step)
        ptlist.append(ptidx * step)

    # test
    if __debug__:
        # this tests we output the last point on the line
        # test = ((output_line_length - 1) * step == input_line_length - 1)
        test = (round(ptidx * step) == input_line_length - 1)
        if not test:
            print "<fibers.py> ERROR: fiber numbers don't add up."
            print step
            print input_line_length
            print output_line_length
            print test
            raise AssertionError

    return ptlist

def get_label_array(input_data, verbose=False):

    input_data.GetLines().InitTraversal()
    num_fibers = input_data.GetNumberOfLines()
    inpoints = input_data.GetPoints()
    inpointsdata = input_data.GetPointData()
    num_points = inpoints.GetNumberOfPoints()

    if verbose:
        print 'Total fiber number in the atlas:', num_fibers, ', point number in the atlas:', num_points

    label_array = None
    if inpointsdata.GetNumberOfArrays() > 0:
        point_data_array_indices = range(inpointsdata.GetNumberOfArrays())
        for idx in point_data_array_indices:
            array = inpointsdata.GetArray(idx)
            if array.GetName().lower() == 'region_label':
                label_array = array

    return label_array

def get_subject_ID(input_data, verbose=False):

    input_data.GetLines().InitTraversal()
    num_fibers = input_data.GetNumberOfLines()
    inpoints = input_data.GetPoints()
    inpointsdata = input_data.GetPointData()
    num_points = inpoints.GetNumberOfPoints()

    if verbose:
        print 'Total fiber number in the atlas:', num_fibers, ', point number in the atlas:', num_points

    label_array = None
    if inpointsdata.GetNumberOfArrays() > 0:
        point_data_array_indices = range(inpointsdata.GetNumberOfArrays())
        for idx in point_data_array_indices:
            array = inpointsdata.GetArray(idx)
            if array.GetName().lower() == 'subject_id':
                label_array = array

    return label_array


def _by_two_endpoints(input_data, outdir, subject_fiber_list, output_fibers_per_group=False):
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
    group_label_per_fiber = numpy.zeros(num_fibers)
    eps_previous = (-1, -1)
    label_count = 0
    for si in sort_indices:
        eps_ii = endpoint_regions[si]
        if eps_ii[0] != eps_previous[0] or eps_ii[1] != eps_previous[1]:
            label_count = label_count + 1
            group_label_per_fiber[si] = label_count
        else:
            group_label_per_fiber[si] = label_count
        eps_previous = eps_ii

    output_file = open(os.path.join(outdir, 'num_fibers_per_region_pair.txt'), 'w')
    outstr = 'Group_Label\t' + 'Region_1\t' + 'Region_2\t' + 'Num_Fibers\t' + 'Num_Subjects\n'

    print '\n<EP region connectivity> Region pair (group) that are less commonly present will be set Group_-1 as are unclassified fibers.'
    # Assign group label to -1, if there are too few fibers in this group.
    # All the group_-1 fibers are the unclassified ones.
    num_subjects = len(numpy.unique(subject_fiber_list))
    for gl in numpy.unique(group_label_per_fiber):
        fiber_idx_per_group = numpy.where(group_label_per_fiber == gl)[0]

        num_fb = len(fiber_idx_per_group)

        subjects_with_connection_per_group = subject_fiber_list[fiber_idx_per_group]
        num_subjects_per_group = len(numpy.unique(subjects_with_connection_per_group))

        # threshold_num_fiber_per_group = 10
        # if num_fb < threshold_num_fiber_per_group:
        if num_subjects_per_group < num_subjects / float(2):
            group_label_per_fiber[fiber_idx_per_group] = -1

        outstr = outstr + str(gl) + '\t' + str(endpoint_regions[fiber_idx_per_group][0][0]) + '\t' + \
                 str(endpoint_regions[fiber_idx_per_group][0][1]) + '\t' + str(num_fb) + '\t' + str(num_subjects_per_group) + '\n'

    output_file.write(outstr)
    output_file.close()

    # Output how many fibers are in each group
    output_file = open(os.path.join(outdir, 'num_fibers_per_group_after_combining_less_common_connection.txt'), 'w')
    outstr = 'index\t' + 'Group_Label\t' + 'Region_1\t' + 'Region_2\t' + 'Num_Fibers\t' + 'Num_Subjects\n'

    label_idx = 1
    for gl in numpy.sort(numpy.unique(group_label_per_fiber)):
        fiber_idx_per_group = numpy.where(group_label_per_fiber == gl)[0]
        # for rr in fiber_idx_per_group:
        #     connectivity[rr, fiber_idx_per_group] = 1

        gl_mask = numpy.zeros([num_fibers, 1])
        gl_mask[fiber_idx_per_group] = 1

        subjects_with_connection_per_group = subject_fiber_list[fiber_idx_per_group]
        num_subjects_per_group = len(numpy.unique(subjects_with_connection_per_group))
        num_fb = numpy.sum(gl_mask)
        if gl == -1:
            region_1 = -1
            region_2 = -1
        else:
            region_1 = endpoint_regions[fiber_idx_per_group][0][0]
            region_2 = endpoint_regions[fiber_idx_per_group][0][1]
        outstr = outstr + str(label_idx) + '\t' + str(gl) + '\t' + str(region_1) + '\t' + \
                 str(region_2) + '\t' + str(num_fb) + '\t' + str(num_subjects_per_group) + '\n'

        if output_fibers_per_group:
            pd_ui = wma.filter.mask(input_data, gl_mask, color=None, preserve_point_data=True, preserve_cell_data=True, verbose=False)
            wma.io.write_polydata(pd_ui, os.path.join(outdir, 'conn{0:04d}'.format(int(label_idx)) + '_G{0:04d}'.format(int(gl)) +
                                                      '_R{0:04d}'.format(int(region_1)) + '-{0:04d}'.format(int(region_2)) + '.vtp'))
        label_idx = label_idx + 1

    output_file.write(outstr)
    output_file.close()

    print '<EP region connectivity> Number of valid region pairs (groups) is', len(numpy.unique(group_label_per_fiber))
    print '  Number of clusters (k) should be set greater than', len(numpy.unique(group_label_per_fiber))

    #connectivity = scipy.sparse.csr_matrix(connectivity)
    connectivity = group_label_per_fiber

    return (connectivity, endpoint_regions)

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
    commissural_sub_cortical_regions = [16, 77, 85]

    WM_Unsegmented = [5001, 5002]

    if label in right_WM_cortical_regions:
        label = label - 2000 - 1000
    elif label in left_WM_cortical_regions:
        label = label - 2000
    elif label in right_GM_cortical_regions:
        label = label - 1000
    elif label in left_GM_cortical_regions:
        label = label
    elif label in CC_regions or label in commissural_sub_cortical_regions:
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