#
# Retrieve a set of embryo or brain slices from the Allen Developing Mouse Brain
# atlas.
#

import requests
import collections # For ordered dictionary

def getImageToImageOffsets (images, x, y):
    # Convert images into an ordered dictionary
    oimages = collections.OrderedDict(sorted(images.items()))
    ## Really just want oimages.first here:
    first = True
    first2 = True
    firstimage = ''
    otherimages = ''
    for k, v in oimages.items():
        if first:
            firstimage = v
            first = False
        else:
            if first2:
                otherimages = '{0}'.format(v)
                first2 = False
            else:
                otherimages = '{0},{1}'.format(otherimages,v)

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
        #print ('m: {0}'.format(m))
        imid = m['image_sync']['section_image_id']
        x_ = m['image_sync']['x']
        y_ = m['image_sync']['y']
        offsets[imid] = (x_,y_)
        print ('set offsets[{0}] to {1}'.format(imid, offsets[imid]))

    return offsets

# All the work happens in here
def retrieve (exptstr):

    section_thickness = -1
    storage_directory = ''
    plane = -1

    # Query SectionDataSet to obtain the section_thickness and plane information
    rurl = 'http://api.brain-map.org/api/v2/data/SectionDataSet/' + exptstr + '.json'
    r = requests.get (rurl)
    rjson = r.json()
    print ('Response success: {0}'.format(rjson['success']))
    success = rjson['success']
    if (success):
        for m in rjson['msg'][0]:
            print ('{0} = {1}'.format(m, rjson['msg'][0][m]))
            section_thickness = rjson['msg'][0]['section_thickness']
            storage_directory = rjson['msg'][0]['storage_directory']
            # 2 is saggital (I think).
            plane = rjson['msg'][0]['plane_of_section_id']
            expression = rjson['msg'][0]['expression']
            if not expression:
                print ('Not expression experiment!')
                return
    else:
        print ('Response failed, exiting.')
        return

    # Obtain metadata via query interface
    rurl = 'http://api.brain-map.org/api/v2/data/query.json' \
           + '?criteria=model::SectionImage,rma::criteria,[data_set_id$eq' \
           + exptstr +']'

    r = requests.get (rurl)
    rjson = r.json()
    print ('Image response success: {0}'.format(rjson['success']))
    if not rjson['success']:
        print ('Failed response was: {0}'.format(rjson))

    # images is a dictionary of the image IDs keyed by the section_number
    images = {}
    if (success):
        num_rows = rjson['num_rows']
        for i in range(0,num_rows):
            image_id = rjson['msg'][i]['id']
            section_num = rjson['msg'][i]['section_number']
            print ('Image {0} is section {2} and its id is {1}'.format(i, image_id, section_num))
            images[section_num] = image_id

    else:
        print ('Response failed, exiting.')
        return

    # Loop through the image IDs downloading each one.
    do_download = 0
    do_download_expr = 1
    if do_download:
        urltail = '' # or '?downsample=4'
        for im in images:
            # Can I make this get the expression version, rather than ISH?
            rurl = 'http://api.brain-map.org/api/v2/image_download/'+str(images[im])+urltail
            r = requests.get (rurl, stream=True)
            filename = 'e{0}_{1:02d}_{2}.jpg'.format(exptstr,im,images[im])
            print ('Downloading image ID: {0} to {1}'.format(images[im], filename))
            with open(filename, 'wb') as fd:
                for chunk in r.iter_content(chunk_size=128):
                    fd.write(chunk)
    if do_download_expr:
        urltail = '?view=expression' # or '?downsample=4&view=expression'
        for im in images:
            # Can I make this get the expression version, rather than ISH? Yes, just add &view=expression
            rurl = 'http://api.brain-map.org/api/v2/image_download/'+str(images[im])+urltail
            r = requests.get (rurl, stream=True)
            filename = 'e{0}_{1:02d}_{2}_expr.jpg'.format(exptstr,im,images[im])
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
        print ('ofs1: {0}'.format(ofs1[o]))

    for o in ofs2:
        print ('ofs2: {0}'.format(ofs2[o]))

    # Values in ofs1 and ofs2 should give a scaling/rotation matrix to apply to the
    # coordinates in each frame.

    return

# Emx2 sets. Use this expt number to retrieve all image numbers and then download
# them, along with relevant metadata, such as scale.
#
# E11.5
#retrieve ('100047257');
# E13.5 100041799
#retrieve ('100041799')
# E15.5 100041837
#retrieve ('100041837')
# E18.5 100041580
retrieve ('100041580')
