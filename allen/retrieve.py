#
# Retrieve a set of embryo or brain slices from the Allen Developing Mouse Brain
# atlas.
#

# Emx2 sets. Use this expt number to retrieve all image numbers and then download
# them, along with relevant metadata, such as scale.
#
# Time  Expt number
# E11.5 100047257
# E13.5 100041799
# E15.5 100041837
# E18.5 100041580
#
import requests

# All the work happens in here
def retrieve (exptstr):

    include = 'expression,section_thickness'

    rurl = 'http://api.brain-map.org/api/v2/data/SectionDataSet/' \
           + exptstr + '.json' \
           #+ '?include=' + include # fails

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

    # http://api.brain-map.org/api/v2/data/query.xml?criteria=
    # model::SectionImage,
    # rma::criteria,[data_set_id$eq70813257]
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
            print ('Image {0} id is {1}'.format(i, image_id))
            images[section_num] = image_id

    else:
        print ('Response failed, exiting.')
        return

    for im in images:
        rurl = 'http://api.brain-map.org/api/v2/image_download/'+str(images[im])+'?downsample=2'
        r = requests.get (rurl, stream=True)
        filename = 'e{0}_{1:02d}_{2}.jpg'.format(exptstr,im,images[im])
        print ('Downloading image ID: {0} to {1}'.format(images[im], filename))
        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

    return

# E11.5
retrieve ('100047257');
# E13.5 100041799
# retrieve ('100041799')
# E15.5 100041837
# retrieve ('100041837')
