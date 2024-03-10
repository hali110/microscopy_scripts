
// confocal folder
folder = "Z:/Data";

// experiment specs
cell_line = "MDA_MB_231";
cell_numb = "1K";
expe_date = "2022-10-27-HA";
expe_name = "MDA-MB-231_1mg_4mg_48hrs";
doxy = newArray(0, 0, 0, 1, 1, 1);
sample_id = newArray(1, 2, 3, 3, 2, 1);
doxy_label = newArray("1mgml", "4mgml");


// import data
folder = folder + "/" + expe_date + "/";
//for (i = 0; i < 7; i += 4) 
{

	run("Bio-Formats Importer", "open=" + folder + expe_name + ".nd2" + " color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
	// open(folder + expe_name + ".nd2");
	dir = folder + expe_name;
	File.makeDirectory(dir);
	stack_siz = nSlices;

	for (n = 1; n < sample_id.length+1; n++) 
	{
		
		// select DIC channel
		name = folder + expe_name + ".nd2 - " + expe_name + ".nd2 (series " + toString(n) + ")";
		selectWindow(name);
		run("Split Channels");

			
		
		//DIC channel 2
		// min intensity projection
		sub_dir = dir + "/DIC/";
		File.makeDirectory(sub_dir);
//		sub_dir = dir + "/DIC/Channel2";
//		File.makeDirectory(sub_dir);
		sub_dir = sub_dir + "/MIN/";
		File.makeDirectory(sub_dir);
		
		name_new = "DIC_MIN_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
	
		selectWindow("C2-" + name);
		run("Z Project...", "projection=[Min Intensity] all");
		selectWindow("MIN_C2-" + name);
		saveAs("Tiff", sub_dir + name_new + ".tif");
		run("Enhance Contrast", "saturated=0.35");
		
		
		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
		run("Colors...", "foreground=black background=white selection=yellow");
		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
		
		run("8-bit");
		saveAs("Gif", sub_dir + name_new + ".gif");
		
		close();
		
		// max intensity projection
		sub_dir = dir + "/DIC/MAX/";
		File.makeDirectory(sub_dir);
		name_new = "DIC_MAX_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
	
		selectWindow("C2-" + name);
		run("Z Project...", "projection=[Max Intensity] all");
		selectWindow("MAX_C2-" + name);
		saveAs("Tiff", sub_dir + name_new + ".tif");
		run("Enhance Contrast", "saturated=0.35");
		
		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
		run("Colors...", "foreground=black background=white selection=yellow");
		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
		
		run("8-bit");
		saveAs("Gif", sub_dir + name_new + ".gif");
		
		close();
		selectWindow("C2-" + name);
		close();

		
		
//		//DIC channel 4
//		// min intensity projection
//		sub_dir = dir + "/DIC/Channel4";
//		File.makeDirectory(sub_dir);
//		sub_dir = sub_dir + "/MIN/";
//		File.makeDirectory(sub_dir);
//		name_new = "C4_DIC_MIN_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//	
//		selectWindow("C4-" + name);
//		run("Z Project...", "projection=[Min Intensity] all");
//		selectWindow("MIN_C4-" + name);
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		run("Enhance Contrast", "saturated=0.35");
//		
//		
//		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
//		run("Colors...", "foreground=black background=white selection=yellow");
//		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
//		
//		run("8-bit");
//		saveAs("Gif", sub_dir + name_new + ".gif");
//		
//		close();
//		
//		// max intensity projection
//		sub_dir = dir + "/DIC/Channel4/MAX/";
//		File.makeDirectory(sub_dir);
//		name_new = "C4_DIC_MAX_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//	
//		selectWindow("C4-" + name);
//		run("Z Project...", "projection=[Max Intensity] all");
//		selectWindow("MAX_C4-" + name);
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		run("Enhance Contrast", "saturated=0.35");
//		
//		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
//		run("Colors...", "foreground=black background=white selection=yellow");
//		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
//		
//		run("8-bit");
//		saveAs("Gif", sub_dir + name_new + ".gif");
//		
//		close();
//		selectWindow("C4-" + name);
//		close();



	

//		// Flu
//		sub_dir = dir + "/CRM/";
//		File.makeDirectory(sub_dir);
//		sub_dir = dir + "/CRM/MAX/";
//		File.makeDirectory(sub_dir);
//		name_new = "CRM_" + "_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//		selectWindow("C1-" + name);
//		run("Z Project...", "projection=[Max Intensity] all");
//		selectWindow("MAX_C1-" + name);
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		run("Enhance Contrast", "saturated=0.35");
//
//		
//		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
//		run("Colors...", "foreground=white background=black selection=yellow");
//		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
//		
//		run("8-bit");
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		close();
//		selectWindow("C1-" + name);
//		close();

		
		
	
		// CRM max intensity projection
		sub_dir = dir + "/CRM/";
		File.makeDirectory(sub_dir);
//		sub_dir = dir + "/CRM/640/";
//		File.makeDirectory(sub_dir);
		sub_dir = dir + "/CRM/MAX/";
		File.makeDirectory(sub_dir);

		name_new = "CRM_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
	
		selectWindow("C1-" + name);
		run("Z Project...", "projection=[Max Intensity] all");
		selectWindow("MAX_C1-" + name);
		saveAs("Tiff", sub_dir + name_new + ".tif");
		run("Enhance Contrast", "saturated=0.35");

		
		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
		run("Colors...", "foreground=white background=black selection=yellow");
		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
		
		run("8-bit");
		saveAs("Gif", sub_dir + name_new + ".gif");
		close();


		// CRM Slice 11
		sub_dir = dir + "/CRM/slice11/";
		File.makeDirectory(sub_dir);
		name_new = "CRM_slice11_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);

			
		selectWindow("C1-" + name);
		run("Duplicate...", "duplicate slices=11");
		selectWindow("C1-" + name + "-1");
		saveAs("Tiff", sub_dir + name_new + ".tif");
		run("Enhance Contrast", "saturated=0.35");
		
		
		
		run("Scale Bar...", "width=500 height=10 font=50 color=Black background=None location=[Lower Right] bold overlay label");
		run("Colors...", "foreground=white background=black selection=yellow");
		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-49");
		
		run("8-bit");
		saveAs("Gif", sub_dir + name_new + ".gif");
		close();

		// Stack
		sub_dir = dir + "/CRM/STACK/";
		name_new = "CRM_stack_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
		File.makeDirectory(sub_dir);
		selectWindow("C1-" + name);
		saveAs("Tiff", sub_dir + name_new + ".tif");
		close();





//		sub_dir = dir + "/CRM/";
//		File.makeDirectory(sub_dir);
//		sub_dir = dir + "/CRM/488/";
//		File.makeDirectory(sub_dir);
//		sub_dir = dir + "/CRM/488/MAX/";
//		File.makeDirectory(sub_dir);
//
//		name_new = "CRM_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//	
//		selectWindow("C3-" + name);
//		run("Z Project...", "projection=[Max Intensity] all");
//		selectWindow("MAX_C3-" + name);
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		run("Enhance Contrast", "saturated=0.35");
//
//		
//		run("Scale Bar...", "width=500 height=10 font=50 color=White background=None location=[Lower Right] bold overlay label");
//		run("Colors...", "foreground=white background=black selection=yellow");
//		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-288");
//		
//		run("8-bit");
//		saveAs("Gif", sub_dir + name_new + ".gif");
//		close();
//
//
//		// CRM Slice 11
//		sub_dir = dir + "/CRM/488/slice11/";
//		File.makeDirectory(sub_dir);
//		name_new = "CRM_slice11_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//
//			
//		selectWindow("C3-" + name);
//		run("Duplicate...", "duplicate slices=11");
//		selectWindow("C3-" + name + "-1");
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		run("Enhance Contrast", "saturated=0.35");
//		
//		
//		
//		run("Scale Bar...", "width=500 height=10 font=50 color=Black background=None location=[Lower Right] bold overlay label");
//		run("Colors...", "foreground=white background=black selection=yellow");
//		run("Label...", "format=00:00 starting=0 interval=10 x=50 y=100 font=50 text=hours range=1-49");
//		
//		run("8-bit");
//		saveAs("Gif", sub_dir + name_new + ".gif");
//		close();
//
//		// Stack
//		sub_dir = dir + "/CRM/488/STACK/";
//		name_new = "CRM_stack_" + expe_name + "_" + doxy_label[doxy[n-1]] + "_spheroid" + toString(sample_id[n-1]);
//		File.makeDirectory(sub_dir);
//		selectWindow("C3-" + name);
//		saveAs("Tiff", sub_dir + name_new + ".tif");
//		close();



		

	}
}
