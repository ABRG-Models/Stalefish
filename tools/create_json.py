import argparse
import json
import os

def get_files(folder, suffix = '.tif'):
	last_folder = os.path.split(folder)[1]
	files = os.listdir(folder)
	files.sort()
	to_return = []
	for file in files:
		if file.endswith(suffix):
			to_return.append(os.path.join(last_folder, file))
	return to_return
def run(folder, thickness, file_to_write):
	files = get_files(folder)
	curr_x = 0
	slices = []
	for file in files:
		slices.append({"filename": file, "x": curr_x})
		curr_x+=thickness
	data =  {
		"desc_thickness" : "The thickness of a brain slice, assuming they'll all be the same.",
		"thickness" : thickness,
		"desc_pixels_per_mm" : "How many pixels in the image corresponds to 1 mm? Assumed to apply equally to horz(y) and vert(z) directions.",
	    "pixels_per_mm" : 233,
	    "desc_bg_blur" : "bg_blur_screen_proportion is multiplied by the width of the image in pixels to get a sigma for the Gaussian which is used to blur the image to subtract the background.",
	    "bg_blur_screen_proportion" : 0.1667,
	    "desc_bg_blur_subtraction_offset" : "signal_img = 255 - (original_img + (bg_blur_subtraction_offset - blurred_bg)). Thus, possible range for bg_blur_subtraction_offset is 0 to 255.",
	    "bg_blur_subtraction_offset" : 180,
	    "desc_save_per_pixel_data" : "If true, then save out the signal values and coordinates of every individual pixel in each box and freehand loop. Can lead to very large data files.",
	    "save_per_pixel_data" : True,
	    "desc_save_auto_align_data" : "Default is true. If false, then don't write coordinates in the autoaligned frame of reference into the HDF5 data file.",
	    "save_auto_align_data" : True,
	    "desc_save_landmark_align_data" : "Default is true. If false, then don't write coordinates in the landmark aligned frame of reference into the HDF5 data file.",
	    "save_landmark_align_data" : True,
	    "desc_slices" : "The brain slices in the y-z plane, where y points to the anterior, z points dorsally (and so x is the lateral dimension).",
	    "slices" : slices
	}
	with open(os.path.join(folder, file_to_write), "wt") as fp:
		json.dump(data, fp, indent = 4)
	print(file_to_write, " added to folder", folder)
def main():
    global args
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument("-f", "--folder_path",
                    help="Full path to folder which contains file to read",
                    type=str, default="../data")
    parser.add_argument("-t", "--thickness",
                    help="thickness to use",
                    type=int, default=0.1)
    parser.add_argument("-w", "--file_to_write",
                    help="JSON to write to",
                    type=str, default="test.json")
    args = parser.parse_args()
    run(args.folder_path, args.thickness, args.file_to_write)
if __name__ == "__main__":
    main()
