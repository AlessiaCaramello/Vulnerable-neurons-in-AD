//create pop up dialog for selecting "automated image processing", "manual adjustment" 
Dialog.create("ImageJ processing options");
	Dialog.setInsets(0, 0, 0);
	Dialog.addMessage("Please organise your image in the following format");
	Dialog.setInsets(5, 20, 0);
	Dialog.addMessage("->Source Directory");
	Dialog.setInsets(0, 40, 0);
	Dialog.addMessage("->Subfolder(eg. mask_channels)");
	Dialog.setInsets(0, 60, 0);
	Dialog.addMessage("->ImageJ will create an output folder here");
	Dialog.setInsets(5, 0, 0);
	Dialog.addDirectory("Locate source directory", "/Users/caramea/Dropbox (UK Dementia Research Institute)/Alessia Lab Dropbox/Projects/UKDRI project/Image analysis/HistoCat:SIMPLI analysis/051022 TREM2 cohort/R code/analysis of plaques/ImageJ marco/");
	Dialog.addString("Specify subfolder name", "mask_channels",10);
	Dialog.addMessage("Specify channels to enlarge below:");
	Dialog.addString("Channel 1:","Ab",15) 
	Dialog.addNumber("Channel 1 Enlarged for:", 50, 0, 3, "um") 
Dialog.show();

SourceDir = Dialog.getString();
subfolder = Dialog.getString();
enlargemarker = Dialog.getString();
enlargesize = Dialog.getNumber();

print("Starting image analysis log...")
print(" ");
print("Enlarge "+enlargemarker+" by "+enlargesize+"um from signal perimeter");
print(" ");

//end of user setting


//list all the images to identified to process
print("Getting images located here:")
print(SourceDir+subfolder+File.separator)
print(" ");

subImageList = getFileList(SourceDir+subfolder);
print("Images located:");
Array.print(subImageList);
print(" ");

//create folder for storing processed file
out_enlarged = SourceDir+"enlarged_output_files"
File.makeDirectory(out_enlarged);
print("Plaques and ring masks will be saved here:");
print(out_enlarged);
print(" ");


//image processing starts 

Myoperation(subImageList);	

function Myoperation(subImageList) {
	
	//open and processing each images in the folder
	for (x = 0; x < subImageList.length; x++) {
		open(SourceDir+subfolder+File.separator+subImageList[x]);
		ROIname = getTitle();
		run("Duplicate...", "title="+ROIname+"_dup");
	    run("Despeckle"); 
					
		//enlarge pathology channel
		setAutoThreshold("Intermodes dark");     
		//setOption("BlackBackground", true);
		run("Analyze Particles...", "size=200-Infinity show=Masks");  
	    run("Gaussian Blur...", "sigma=5");      
	    setAutoThreshold("Otsu dark");           
	
		//if there is pathology marker present 
    	run("Create Selection");    
		run("Make Inverse");        
						
			if (selectionType() != -1) {
				run("Create Mask");        
				run("Create Selection");
						roiManager("Add");
						roiManager("Select", 0);
						roiManager("Rename", "plaques");
						run("Create Mask");
						saveAs("png", out_enlarged+File.separator+substring(ROIname, 0, lengthOf(ROIname) - 5)+"_"+enlargemarker+"_plaques.png");
						
						
						roiManager("Select", 0);
						run("Enlarge...", "enlarge=enlargesize");
						roiManager("Add");
						roiManager("Select", 1);
						roiManager("Rename", "50");
						
						roiManager("Select", newArray(0,1));
						roiManager("XOR");
						roiManager("Add");
						roiManager("Select", 2);
						roiManager("Rename", "ring");
						run("Create Mask");
						saveAs("png", out_enlarged+File.separator+substring(ROIname, 0, lengthOf(ROIname) - 5)+"_"+enlargemarker+"_ring.png");
						
						run("Close All");
						close("ROI manager");
					
					}
					//if no pathology marker present in the image
					else {
					close("*dup");
					close("*.tiff");
					print("No plaques detected in "+ROIname);
					}
			}

	
	