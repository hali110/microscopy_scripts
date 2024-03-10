
// confocal folder
folder = "Z:\\Data";

// experiment specs
cell_line = "SUM44";
expe_date = "2022-06-16-HA";
expe_name = "SUM44_GFP_2_4_8k";
stack_siz = 21;

// import data
folder = folder + "\\" + expe_date + "\\";
open(folder + expe_name + ".nd2");
dir = folder + expe_name;
File.makeDirectory(dir);

for (n = 1; n < 96+1; n++) {
	
	// select DIC channel
	if (n < 10) {
		name = expe_name + ".nd2 - " + expe_name + ".nd2 (series 0" + toString(n) + ")";
		} else {
		name = expe_name + ".nd2 - " + expe_name + ".nd2 (series " + toString(n) + ")";
	}
	if (n < 33)
	{
		cell_num = "2k";
	}
	else if (n < 65)
	{
		cell_num = "4k";
	}
	else 
	{
		cell_num = "8K";
	}
	
	selectWindow(name);
	run("Split Channels");

	// FLU
	selectWindow("C1-" + name);
	sub_dir = dir + "\\FLU\\";
	File.makeDirectory(sub_dir);
	sub_dir = dir + "\\FLU\\MAX\\";
	File.makeDirectory(sub_dir);
	name_new = "FLU_MAX_" + expe_name + "_" + cell_num + "_series" + toString(n);

	selectWindow("C1-" + name);
	run("Z Project...", "projection=[Max Intensity] all");
	selectWindow("MAX_C1-" + name);
	run("Enhance Contrast", "saturated=0.35");
	run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
	saveAs("Tiff", sub_dir + name_new + ".tif");
	close();
	selectWindow("C1-" + name);
	close();













	
	// slice 11
	sub_dir = dir + "\\DIC\\";
	File.makeDirectory(sub_dir);
	sub_dir = dir + "\\DIC\\Slice 11\\";
	File.makeDirectory(sub_dir);
	name_new = "DIC_Slice11_" + expe_name + "_" + cell_num + "_series" + toString(n);
		
	selectWindow("C2-" + name);
	run("Duplicate...", "duplicate slices=11");
	selectWindow("C2-" + name + "-1");
	run("Enhance Contrast", "saturated=0.35");
	run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
	
	
		
	saveAs("Tiff", sub_dir + name_new + ".tif");
	close();

	// max intensity projection
	sub_dir = dir + "\\DIC\\MAX\\";
	File.makeDirectory(sub_dir);
	name_new = "DIC_MAX_" + expe_name + "_" + cell_num + "_series" + toString(n);

	selectWindow("C2-" + name);
	run("Z Project...", "projection=[Max Intensity] all");
	selectWindow("MAX_C2-" + name);
	run("Enhance Contrast", "saturated=0.35");
	run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
	
	
	
	saveAs("Tiff", sub_dir + name_new + ".tif");
	close();

	// min intensity projection
	sub_dir = dir + "\\DIC\\MIN\\";
	File.makeDirectory(sub_dir);
	name_new = "DIC_MIN_" + expe_name + "_" + cell_num + "_series" + toString(n);

	selectWindow("C2-" + name);
	run("Z Project...", "projection=[Min Intensity] all");
	selectWindow("MIN_C2-" + name);
	run("Enhance Contrast", "saturated=0.35");
	run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
	
	saveAs("Tiff", sub_dir + name_new + ".tif");
	close();

	// close stack
	selectWindow("C2-" + name);
	close();

}