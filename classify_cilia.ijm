

#@ File (label = "Input folder", style = "directory") input_dir
#@ File (label = "Output directory for segmentation", style = "directory") output_dir
#@ File (label = "Cilia classifier", style = "file") cilia_classifier

#@ Integer (label = "Channel for nuclei", min=1, max=8, value = 1) nuc_channel
#@ Integer (label = "Channel for kif", min=1, max=8, value = 2) kif_channel
#@ Integer (label = "Channel for membrane", min=1, max=8, value = 3) membrane_channel
#@ Integer (label = "Channel for arl13b", min=1, max=8, value = 4) arl13b_channel

File.makeDirectory(output_dir);


file_list = getFileList(input_dir)
number_of_files = file_list.length;


for (i=0; i<number_of_files; i++) {

	close("*");
	
	
	open(input_dir + File.separator + file_list[i]);
	name = getTitle();
	rename("input");
	
	setSlice(nuc_channel);
	run("Duplicate...", "title=nuclei");
	selectImage("input");
	setSlice(kif_channel);
	run("Duplicate...", "title=cell_intensity");
	selectImage("input");
	setSlice(membrane_channel);
	run("Duplicate...", "title=cilia_membrane");
	selectImage("input");
	setSlice(arl13b_channel);
	run("Duplicate...", "title=arl13b");
	selectImage("input");
	run("Close");
	run("Merge Channels...", "c1=arl13b c2=cilia_membrane create keep");
	
	
	//Semantically Segment Cilia
	run("Segment Image With Labkit", "input=Composite segmenter_file=" + cilia_classifier + " use_gpu=true");
	rename("cilia_segmentation");
	selectImage("cilia_segmentation");
	save(output_dir + File.separator + "classified_" + name);
	
	
	}
	close("*");
	
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	
	f = File.open(output_dir + File.separator + "script_params.txt");
	print(f, "Script run on: " + year + " " + month+1 + " " + dayOfMonth + " at " + hour + ":" + minute + ":" + second);
	print(f, "Input Data: " + input_dir);
	print(f, "Output Data: " + output_dir);
	print(f, "Cilia Classifier Used: " + cilia_classifier);
	File.close(f);


	print("pipeline complete!")