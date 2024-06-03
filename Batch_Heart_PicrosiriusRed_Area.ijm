//GUI creation
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Integer (label="Pyramid Level to analyze", style="slider", min=1, max=6, value= 2, stepSize=1, persist=true) level
#@ String (label="Specify Staining 1", style="text field", value= "Picrosirius Red", persist=true) CH1
#@ String (label="Specify Staining 2", style="text field", value= "Residue", persist=true) CH2
#@ Integer (label="Number of Areas to analyze", style="slider", min=1, max=6, value= 3, stepSize=1, persist=true) nROI
#@ Integer (label="ROI diameter [µm]", style="slider", min=100, max=5000, value= 1000, stepSize=1, persist=true) size

//Preparing Stage
run("Bio-Formats Macro Extensions");

setOption("BlackBackground", true);
print("\\Clear");
close("*");

if (roiManager("count")>0) {
roiManager("Deselect");
roiManager("Delete");
}

File.makeDirectory(output + "\\Manual_ROIS");
Table.create("Summary");

//Processing the folder to create previes on lower scale
processFolder_regionSelect(input);

waitForUser("Preselection Phase Finished", "Preselection loop has been completed, processing will continue in background now");

//Processing the folder again to do the complete Processing
setBatchMode(true);
processFolder_regionProcess(input);

//save Results and close all open windows
selectWindow("Summary");
saveAs("Results", output + "\\Results_Level_" + level + ".csv");
close("Results_Level_" + level + ".csv");
close("Results");
close("ROI Manager");
close("*");

setBatchMode(false);

print("\\Clear");
print("Batch processing completed");

print("            { ");
print("            { ");
print("         {   } ");
print("          }_{ { ");
print("       .-{   }   }-. ");
print("      (   }     {   ) ");
print("   |`-.._____..-'| ");
print("   |                    ;--. ");
print("   |                   (__     | ");
print("   |                    |    )  ) ");
print("   |                    |  /  / ");
print("   |                     /  /    ");
print("   |                   (  / ");
print("   |                    y ");
print("    `-.._____..-'; ");

//End of Macro


// function to scan folders/subfolders/files to find files with correct suffix and process the region creation on a lower scale
function processFolder_regionSelect(input) {
	list = getFileList(input);
	A=list.length;
	print("0 of "  + A + " files preprocessed");
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder_regionSelect(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			prepareRois(input, output, list[i]);
	}
}
// function to scan folders/subfolders/files to find files with correct suffix and do the final processing
function processFolder_regionProcess(input) {
	list = getFileList(input);
	A=list.length;
	print("0 of "  + A + " files processed");
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			rocessFolder_regionProcess(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			processFile(input, output, list[i]);
	}
}

//function to annotate ROIs in low-magnification images which load faster 
function prepareRois(input, output, file) {
	//Import&readout of pyramid level factor
	setBatchMode(true);
	run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_4");
	height2=getHeight();
	close();
	run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_5");
	height3=getHeight();
	close();
	base=round(height2/height3);
	setBatchMode(false);

	//Opening file on lower magnification
	run("Bio-Formats Importer", "open=[" + input + "\\" + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_4");
	title=getTitle();
	array1=split(title," ");
	title=array1[0];

	//check wheter the presegmentation for this file already exists
	if (File.exists(output + "\\Manual_ROIS\\" + title  +".zip")) {
	close("*");
	} else {

	//preparing file & getting statistics
	rename(title);
	Stack.getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);

	//Correct for wrong scaling in Pyramid formats
	pixelWidth2=(pixelWidth*pow(base, 3));
	pixelHeight2=(pixelHeight*pow(base, 3));
	run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" pixel_width="+pixelWidth2+" pixel_height="+pixelHeight2+" voxel_depth=3.0000000");

	//adjusting to widescreen for improved segmentation properties
	if (width<height) {
	run("Rotate 90 Degrees Right");
	}

	//Detecting total slice area
	run("Duplicate...", "title=Mask_Aorta duplicate");
	run("RGB Color");
	run("8-bit");
	getStatistics(area, mean, min, max, std, histogram);
	run("Duplicate...", "title=Mask_Aorta_backfill duplicate");

		//Compensate missing tiles
		if (min==0) {
		setThreshold(0, 1);
		run("Convert to Mask");
		run("Create Selection");
		selectWindow("Mask_Aorta (RGB)");
		run("Restore Selection");
		run("Fill", "slice");
		run("Select None");
		}

		//Compensate White Tiles
		selectWindow("Mask_Aorta (RGB)");
		run("Duplicate...", "title=Mask_Aorta_frontfill duplicate");
		setThreshold(255, 255);
		run("Convert to Mask");
		run("Create Selection");
		run("Make Inverse");

		//Tissue Detection
		selectWindow("Mask_Aorta (RGB)");
		run("Restore Selection");
		setAutoThreshold("Huang");
		run("Convert to Mask");
		run("Remove Outliers...", "radius=3 threshold=50 which=Bright");
		run("Remove Outliers...", "radius=3 threshold=50 which=Dark");
		run("Create Selection");
	
	roiManager("Add");
	roiManager("Select", 0);
	roiManager("Rename", "Section");
	close("Mask*");

	//create manual ROIs for analysis
	selectWindow(title);
	for (i = 1; i < nROI+1; i++) {
	run("Specify...", "width="+size+" height="+size+" x=5000 y=5000 slice=1 oval constrain scaled");
	setTool("oval");
	waitForUser("Create ROI #"+i+" for Analysis");
	roiManager("add");
	roiManager("Select", i);
	roiManager("rename", "ROI_"+i);
	run("Select None");
	}

	// save the RoiSet for further processing steps
	roiManager("save", output + "\\Manual_ROIS\\" + title  +".zip");
	roiManager("deselect");
	roiManager("delete");
	run("Close All");
	}
	}

// function to do the actual processing
function processFile(input, output, file) {

	// Bio-Formats Importer opens files in specified pyramid stage and gets metadata
	//getting the Pyramid stage scale
	run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_4");
	height2=getHeight();
	close();
	run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_5");
	height3=getHeight();
	close();
	base=round(height2/height3);

	//Opening the file and getting statistics
	run("Bio-Formats Importer", "open=[" + input + "\\" + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+ level);
	title=getTitle();
	array1=split(title," ");
	title=array1[0];
	rename(title);
	Stack.getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);


	//correction for widescreen 
	if (width<height) {
		run("Rotate 90 Degrees Right");
	}

	//Correct for wrong scaling in Pyramid formats
	pixelWidth2=(pixelWidth*pow(base, level-1));
	pixelHeight2=(pixelHeight*pow(base, level-1));
	run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" pixel_width="+pixelWidth2+" pixel_height="+pixelHeight2+" voxel_depth=3.0000000");

	//rescale ROIs
	roiManager("open", output + "\\Manual_ROIs\\" + title + ".zip");
	for (o=0; o<roiManager("count"); ++o) {
		roiManager("Select", o);
		ROI=Roi.getName;
		// Scale ROI
		run("Scale... ", "x="+ pow(base, 4-level)+" y="+ pow(base, 4-level));
	
		// Replace old ROI with scaled one
		roiManager("update");
		}

	//Combine ROIs to exclude empty spaces
	for (o=1; o<roiManager("count"); ++o) {
	roiManager("Select", newArray(0,1));
	roiManager("AND");
	setSelectionName("ROI_"+o);
	roiManager("add");
	roiManager("deselect");
	roiManager("Select", 1);
	roiManager("delete");
	roiManager("deselect");
	}

	//Color Deconvolution - Preset PicrosiriusRed
	run("RGB Color");
	rename(title);
	run("Colour Deconvolution2", "vectors=[User values] output=[32bit_Absorbance] simulated [r1]=0.20913965868271805 [g1]=0.8164369523939927 [b1]=0.5382297891530023 [r2]=0.21965491211132712 [g2]=0.34945497487968585 [b2]=0.9108418853550836 [r3]=0.5263468784045569 [g3]=0.598699182427183 [b3]=0.6037534700147082");
	print("\\Clear");
	print(i  + " of "  + A + " files processed");

	//Thresholding Channels
		//1st Channel
		selectWindow(title + "-(Colour_1)A");
		run("Gaussian Blur...", "sigma=0.5 scaled");

		roiManager("Select", "Section");
		setAutoThreshold("Otsu dark");
		run("Convert to Mask");

		//2nd Channel
		selectWindow(title + "-(Colour_2)A");
		run("Gaussian Blur...", "sigma=0.5 scaled");
		
		roiManager("Select", "Section");
		setAutoThreshold("Otsu dark");
		run("Convert to Mask");


	//Creating Measurements 
	selectWindow("Summary");
	run("Set Measurements...", "area area_fraction display redirect=None decimal=0");
	Table.set("Slice", i, title)
	Table.update;
	roiManager("Select", "Section");
	run("Measure");
		area=getResult("Area", 0);
		selectWindow("Summary");
		Table.set("Section_Total area [µm²]", i, area);
		run("Clear Results");

	for (o=1; o<roiManager("count"); ++o) {
	selectWindow(title + "-(Colour_1)A");
	roiManager("Select", o);
		ROI=Roi.getName;
		//measure ROI areas
		run("Measure");
		selectWindow("Summary");
		Table.set(ROI + "_Total area [µm²]", i, (getResult("Area", 0)));
		Table.set(ROI + "_" + CH1+ "_"+" relative [%]", i, getResult("%Area", 0));
		run("Clear Results");
		
		selectWindow(title + "-(Colour_2)A");
	    roiManager("Select", o);
		ROI=Roi.getName;
		run("Measure");
		selectWindow("Summary");
		Table.set(ROI + "_"+ CH2 + "_"+" relative [%]", i, getResult("%Area", 0));
		run("Clear Results");
		}

//clean up
close("*");
roiManager("deselect");
roiManager("delete");

print((i+1)  + " of "  + A + " files processed");
}

