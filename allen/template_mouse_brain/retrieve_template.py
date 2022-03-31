#
# Retrieve a set of template mouse brain slices from the Allen Developing Mouse
# Brain atlas.
#

import requests
import json
import collections # For ordered dictionary
import os # For directory creation

# Returns offsets keyed by image id. Better to key by section num
def getImageToImageOffsets_byid (images, x, y):
    print ('images: {0}'.format(images))
    # Convert images into an ordered dictionary
    oimages = collections.OrderedDict(sorted(images.items()))
    ## Really just want oimages.first here:
    first = True
    first2 = True
    firstimage = ''
    firstimageKey = '' # not used
    otherimages = ''
    otherimagesKey = '' # not used
    for k, v in oimages.items():
        if first:
            firstimage = '{0}'.format(v)
            firstimageKey = '{0}'.format(k)
            first = False
        else:
            if first2:
                otherimages = '{0}'.format(v)
                otherimagesKey = '{0}'.format(k)
                first2 = False
            else:
                otherimages = '{0},{1}'.format(otherimages,v)
                otherimagesKey = '{0},{1}'.format(otherimagesKey,k)

    rurl = 'http://api.brain-map.org/api/v2/image_to_image_2d/{0}.json?x={1}&y={2}&section_image_ids={3}'.format(firstimage,x,y,otherimages)
    print ('rurl: {0}'.format(rurl))
    r = requests.get (rurl)
    rjson = r.json()
    if not rjson['success']:
        print ('Failed response was: {0}'.format(rjson))
        return {}

    offsets = {}
    offsets[firstimage] = (x,y)
    for m in rjson['msg']:
        print ('m: {0}'.format(m))
        imid = m['image_sync']['section_image_id']
        x_ = m['image_sync']['x']
        y_ = m['image_sync']['y']
        offsets[imid] = (x_,y_)
        print ('set offsets[{0}] to {1}'.format(imid, offsets[imid]))

    return offsets

def getImageSliceX_byid (images):
    # Convert images into an ordered dictionary:
    oimages = collections.OrderedDict(sorted(images.items()))
    exes = {}
    for k, v in oimages.items():
        rurl = 'http://api.brain-map.org/api/v2/image_to_reference/{0}.json?x=0&y=0'.format(v)
        print ('image_to_reference rurl: {0}'.format(rurl))
        r = requests.get (rurl)
        rjson = r.json()
        print ('image_to_reference rjson: {0}'.format(rjson))
        if not rjson['success']:
            print ('Failed response was: {0}'.format(rjson))
            return {}
        # Get 'x' from rjson
        x = rjson['msg']['image_to_reference']['x']
        print ('The x value for slice {0} is {1} with image id {2}'.format (v, x, k))
        exes[v] = (k, x) # pack up the 'section_image_id' and the 'section_data_set_id'
    return exes

def getImageSliceX (images):
    slicexes = getImageSliceX_byid (images)
    print ('slicexes: {0}'.format (slicexes))
    slicexes_bynum = {}
    for num in images:
        print ('num = {0} images[{0}] = {1} slicexes[images[{0}]] = {2}'.format(num, images[num], slicexes[images[num]]))
        in_key = images[num]
        print ('in_key: {0}'.format(in_key))
        se_key = slicexes[in_key][1]
        print ('se_key: {0}'.format(se_key))
        slicexes_bynum[num] = se_key
    # Sort them into order:
    slicexes_sorted = collections.OrderedDict(sorted(slicexes_bynum.items()))
    return slicexes_sorted

# Return offsets keyed by section num
def getImageToImageOffsets (images, x, y):

    offsets = getImageToImageOffsets_byid (images, x, y)
    offsets_bynum = {}
    for num in images:
        offsets_bynum[num] = offsets['{0}'.format(images[num])]
    # Sort them into order:
    offsets_sorted = collections.OrderedDict(sorted(offsets_bynum.items()))
    return offsets_sorted

# All the work happens in here. This will download images and save them, and also
# retrieve metadata.
#
# @atlasstr The Atlas number, as a string. You find this on the Allen
# website. This number identifies a sequence of images for one particular mouse
# brain atlas.
#
# Do we downsample the image? downsample 1 divides images in two; downsample 2 makes
# image a quarter size (along a length), etc.
def retrieve (atlasstr, downsample=2, do_download_annot=0):

    # Values that are going to go into the top level of the expt json file
    section_thickness = 0.1 # HARDCODE to 10 um for atlas as I don't know how to get this out of the API.
    plane = -1

    # Do we download the ISH images? Yes, unless we're debugging
    do_download = 1

    # Create an experiment directory
    atdir = './atlas/' + atlasstr + '_annot' if do_download_annot else ''
    print ('Create directory {0}'.format(atdir))
    if not os.path.exists (atdir):
        os.makedirs (atdir)

    # The json file lives in atdir
    atjson = atdir + '/' + atlasstr + '.json'

    # We need a dictionary structure to collect together information about each slice.
    sliceinfo = {}

    # http://api.brain-map.org/api/v2/data/query.xml?criteria=model::AtlasImage,rma::criteria,[annotated$eqtrue],atlas_data_set(atlases[id$eq602630314]),alternate_images[image_type$eq'Atlas+-+Adult+Mouse'],rma::options[order$eq'sub_images.section_number'][num_rows$eqall]

    # Query SectionDataSet to obtain the section_thickness and plane information
    rurl = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::AtlasImage,rma::criteria,[annotated$eqtrue],atlas_data_set(atlases[id$eq" + atlasstr + "]),alternate_images[image_type$eq'Atlas+-+Adult+Mouse'],rma::options[order$eq'sub_images.section_number'][num_rows$eqall]"

    sliceinfo['thickness'] = 0.01 # HARDCODE

    # TO fill with key-value pairs where key is slice number and value is an object
    sliceinfo['by_slice'] = {}

    sliceinfo['pixels_per_mm'] = 0 # HARDCODE?

    sliceinfo['downsample'] = downsample

    # Obtain metadata via query interface
    #rurl = 'http://api.brain-map.org/api/v2/data/query.json' \
    #       + '?criteria=model::SectionImage,rma::criteria,[data_set_id$eq' \
    #       + atlasstr +']'

    r = requests.get (rurl)
    rjson = r.json()
    print ('rjson: {0}'.format (rjson))
    print ('Image response success: {0}'.format(rjson['success']))
    success = rjson['success']
    if not rjson['success']:
        print ('Failed response was: {0}'.format(rjson))

    # images is a dictionary of the image IDs keyed by the section_number
    images = {}
    if (success):
        num_rows = rjson['num_rows']
        for i in range(0,num_rows):
            print ("rjson['msg'][{0}]: {1}".format(i, rjson['msg'][i]))
            image_id = rjson['msg'][i]['id']
            section_num = rjson['msg'][i]['section_number']
            print ('Image {0} is section {2} and its id is {1}'.format(i, image_id, section_num))
            # Although image_id is in sliceinfo[].image_id, I create this simpler
            # dictionary, so that it can be ordered when passed to
            # getImageToImageOffsets()
            images[section_num] = image_id
            # I also put the image ID into sliceinfo, which holds more info.
            sliceinfo['by_slice'][section_num] = {}
            sliceinfo['by_slice'][section_num]['image_id'] = image_id
            sliceinfo['by_slice'][section_num]['resolution'] = rjson['msg'][i]['resolution'] # THIS sets pixels_per_mm.
            # Note that this is pixels per mm in the 'best' possible resolution:
            sliceinfo['by_slice'][section_num]['pixels_per_mm'] =  1000.0 / (float(2**downsample) * float(rjson['msg'][i]['resolution']))
            print ('resolution: {0} ppm: {1}'.format(rjson['msg'][i]['resolution'], sliceinfo['by_slice'][section_num]['pixels_per_mm']))
            if sliceinfo['by_slice'][section_num]['pixels_per_mm'] != sliceinfo['pixels_per_mm']:
                sliceinfo['pixels_per_mm'] = sliceinfo['by_slice'][section_num]['pixels_per_mm']

            # We also have width, height (of the sub-image) and image_width  and image_height (of the 'whole' image
            if section_thickness != -1:
                sliceinfo['by_slice'][section_num]['axial_position'] = section_num * section_thickness
            #slices.insert (sl_idx, {'filename': atdir + '/e{0}_{1:02d}_{2}_expr.jpg'.format(atlasstr,section_num,image_id),
            #                'x': section_num*section_thickness});
            #sl_idx = sl_idx + 1

    else:
        print ('Response failed, exiting.')
        return

    # Loop through the image IDs downloading each one.
    if do_download:
        urltail = '?downsample={0}&annotation={1}'.format(downsample, 'true' if do_download_annot else 'false')
        for im in images:
            # Can I make this get the expression version, rather than ISH?
            rurl = 'http://api.brain-map.org/api/v2/atlas_image_download/'+str(images[im])+urltail
            print ('URL: {0}'.format (rurl))
            r = requests.get (rurl, stream=True)
            filename = atdir + '/e{0}_{1:02d}_{2}.jpg'.format(atlasstr,im,images[im])
            print ('Downloading image ID: {0} to {1}'.format(images[im], filename))
            with open(filename, 'wb') as fd:
                for chunk in r.iter_content(chunk_size=128):
                    fd.write(chunk)

    # Image to image registration - lining up the slices
    # http://api.brain-map.org/api/v2/image_to_image_2d/68173101.xml?x=6208&y=2368&section_image_ids=68173103,68173105,68173107

    ofs1 = getImageToImageOffsets (images, 1000, 1000)
    # This gives us, for each non-first image, the x and y which are in the same
    # position as x,y in the first image; An offset for one pixel. If we do another
    # pixel, then we get offset and rotation.
    ofs2 = getImageToImageOffsets (images, 2000, 2000)

    for o in ofs1:
        sliceinfo['by_slice'][o]['offset1'] = ofs1[o]
        print ('ofs1: {0} {1}'.format(o, ofs1[o]))
    for o in ofs2:
        print ("Trying to set sliceinfo['by_slice'][{0}]['offset2']...".format(o))
        sliceinfo['by_slice'][o]['offset2'] = ofs2[o]
        print ('ofs2: {0} {1}'.format(o, ofs2[o]))

    exes = getImageSliceX (images)
    for ex in exes:
        print ("Trying to set sliceinfo['by_slice'][{0}]['slice_x']...".format(ex))
        sliceinfo['by_slice'][ex]['slice_x'] = exes[ex] # Need to order by
        print ('key {0}: val {1}'.format (ex, exes[ex]))

    # Now take the info info sliceinfo, and create a thing called slices, which is the
    # right format for Stalefish.
    slices = []
    sl_idx = 0
    ordered_by_slice = collections.OrderedDict(sorted(sliceinfo['by_slice'].items()))
    for sl in ordered_by_slice:
        slices.insert (sl_idx, {'filename': atdir + '/e{0}_{1:02d}_{2}.jpg'.format(atlasstr,sl,sliceinfo['by_slice'][sl]['image_id']),
                                'x': sliceinfo['by_slice'][sl]['slice_x']/1000.0, # /1000 does um to mm conversion
                                'resolution': sliceinfo['by_slice'][sl]['resolution']});
        sl_idx = sl_idx + 1

    sliceinfo['slices'] = slices

    sliceinfo['colourmodel'] = 'allen_atlas' if do_download_annot else 'greyscale'
    sliceinfo['slices_prealigned'] = True

    # Remove what we don't want in the JSON (to save it becoming cluttered)
    del sliceinfo['by_slice']

    # Copy and paste in the json created with colors.m and the program compare:
    colour_transform = json.loads('{"colour_trans" : [ 0.000000, 14.374037, 36.991335 ], "colour_rot" : [ 0.597316, 0.580850, 0.553017, 0.786215, -0.287928, -0.546775, -0.158365, 0.761387, -0.628657 ], "ellip_axes" : [ 29.768308, 6.712414 ], "luminosity_factor" : -0.004, "luminosity_cutoff" : 250.000000}')

    # Merge colour_transform into slice_info
    sliceinfo = {**sliceinfo, **colour_transform}

    # Write sliceinfo to file
    with open (atjson, 'w') as f:
        f.write (json.dumps(sliceinfo, sort_keys=True, indent=4))

    return

# Adult P56 mouse brain
# retrieve ('602630314')

import sys
if len(sys.argv) < 2:
    print ('')
    print ('Usage: {0} atlas_number [downsample_num] [download_annotations]'.format (sys.argv[0]))
    print ('')
    print (' Where atlas_number is the number you find for a given template atlas on the ')
    print (' Allen website')
    print ('')
    print (' [downsample_num] defaults to 2 and is the number of times the original (huge) image')
    print (' on the Allen server is halved in side length before being sent to you')
    print ('')
    print (' [download_annotations] should be 0 or 1 (if present; defaults to 0). If 1, ')
    print (' coloured annotation maps are downloaded rather than the greyscale MRI images.')
    print ('')
    exit(0)

enum = str(sys.argv[1])
if (len(sys.argv) > 3):
    retrieve (enum, int(sys.argv[2]), int(sys.argv[3]))
elif (len(sys.argv) > 2):
    retrieve (enum, int(sys.argv[2]))
else:
    retrieve (enum)
