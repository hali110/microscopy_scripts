// data folder
folder = "Y:\\tmr_data\\Olympus MPE-RS TWIN\\FLIM\\Troubleshooting\\20221008\\Olympus";

// experiment specs
cell_name = "MCF10A_YAP5SA";
coll_conc = "3mgml";
trea_type = newArray("-D","-D","-D","+D","+D","+D");
sphe_numb = newArray(1,2,3,1,2,3);
area_numb = newArray(1,1,1,1,2,2);

// import data
for (i = 1; i < trea_type.length+1; i++) {
	
	for (j = 1; j < area_numb[i-1]+1; j++) {
		
		subfolder = folder + "\\" + coll_conc + "\\" + trea_type[i-1] + "\\spheroid" + toString(sphe_numb[i-1]) + "\\area" + toString(j);
		name_nadh = cell_name + "_" + coll_conc + "_" + trea_type[i-1] + "_spheroid" + toString(sphe_numb[i-1]) + "_740nm_stack";
		name_shg = cell_name + "_" + coll_conc + "_" + trea_type[i-1] + "_spheroid" + toString(sphe_numb[i-1]) + "_1060nm_stack";
		
		// open DIC + NADH
		open(subfolder + "\\" + name_nadh + ".oir");
		title_nadh = getTitle();
		run("Split Channels");
		selectWindow("C1-" + title_nadh);
		run("Z Project...", "projection=[Min Intensity]");
		
		// save MIN DIC image
		name_save_dic = cell_name + "_" + coll_conc + "_" + trea_type[i-1] + "_spheroid" + toString(sphe_numb[i-1]) + "_area" + toString(j) + "_DIC";
		run("Scale Bar...", "width=100 height=10 font=30 color=Black background=None location=[Upper Right] bold overlay label");
		saveAs("Tiff", folder + "\\" + name_save_dic + ".tif");
		close();
		selectWindow("C1-" + title_nadh);
		close();
		
		// open SHG
		open(subfolder + "\\" + name_shg + ".oir");
		title_shg = getTitle();
		if (i == 6 && j == 2) {
			run("Duplicate...", "duplicate range=1-9 use");
			title_shg_new = getTitle();
			selectWindow(title_shg);
			close();
			selectWindow(title_shg_new);
			rename(title_shg);
			}
		
		// merge NADH and SHG and save
		run("Merge Channels...", "c1=C2-" + title_nadh + " c2=" + title_shg + " create");
		run("Channels Tool...");
		Stack.setChannel(1);
		run("Red");
		setMinAndMax(500, 3000);
		Stack.setChannel(2);
		run("Green");
		setMinAndMax(0, 4095);
		name_save_mpm = cell_name + "_" + coll_conc + "_" + trea_type[i-1] + "_spheroid" + toString(sphe_numb[i-1]) + "_area" + toString(j) + "_MPM";
		run("Scale Bar...", "width=100 height=10 font=30 color=White background=None location=[Upper Right] bold overlay label");
		saveAs("Tiff", folder + "\\" + name_save_mpm + ".tif");
		close();

	}
	
}