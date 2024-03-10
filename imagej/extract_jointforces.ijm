// confocal folder
//folder = "/Volumes/haider/tmr/data/Nikon AX";
folder = "Y:\\Data\\Adil"


// experiment specs
cell_line = "MCF10A_YAP5SA";
expe_date = "2023-10-19-MCF10ATAZ\\Timelapses";
expe_name = "Day";
doxy = newArray(0, 0, 0, 1, 1, 1);
sample_id = newArray(1, 2, 3, 3, 2, 1);
doxy_label = newArray("dox0", "dox1");
days = newArray(0,2,4);



// import data
folder = folder + "\\" + expe_date + "\\";
save_dir = folder + "tfm\\";
File.makeDirectory(save_dir);


for(i = 0; i < days.length; i++)
{
	
	day_dir = save_dir + "day" + toString(days[i]) + "\\";
	File.makeDirectory(day_dir);
	
	for(series = 1; series < sample_id.length+1; series++)
	{
		run("Bio-Formats Importer", "open=[" + folder + expe_name + toString(days[i]) + ".nd2]" + " color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack series_" + toString(series));
		name = expe_name + toString(days[i]) + ".nd2 - " + expe_name + toString(days[i]) + ".nd2 (series " + toString(series) + ")";
		stack_siz = nSlices;
		
		selectWindow(name);
		run("Split Channels");
		series_id = doxy_label[doxy[series-1]] + "_spheroid" + toString(sample_id[series-1]);
		series_dir = day_dir + series_id + "\\";
		File.makeDirectory(series_dir);
		
		// DIC MIN
		sub_dir = series_dir + "DIC\\";
		File.makeDirectory(sub_dir);
		name_new = "DIC_MIN_" + expe_name + "_" + doxy_label[doxy[series-1]] + "_spheroid" + toString(sample_id[series-1] + "_");
		selectWindow("C2-" + name);
		run("Z Project...", "projection=[Min Intensity] all");
		selectWindow("MIN_C2-" + name);
		run("Image Sequence... ", "dir=[" + sub_dir + "] format=TIFF name=" + name_new + " use");
		
		
		// CRM MID SCLICE
		sub_dir = series_dir + "CRM\\";
		File.makeDirectory(sub_dir);
		selectWindow("C1-" + name);
		Stack.getDimensions(w, h, c, s, f);
		midz = round(s*0.5);
		run("Duplicate...", "duplicate slices=" + toString(midz));
		selectWindow("C1-" + name + "-1");
		name_new = "CRM_slice" + toString(midz) + "_" + expe_name + toString(days[i]) + "_" + doxy_label[doxy[series-1]] + "_spheroid" + toString(sample_id[series-1] + "_");
		close("\\Others");
		run("Image Sequence... ", "dir=[" + sub_dir + "] format=TIFF name=" + name_new);
		close();
		
	}
	
}

