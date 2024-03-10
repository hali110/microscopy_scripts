

// define parameters
data_folder = "/Users/haider/Documents/tmr/projects/flim/data/20221120/tifs";
cell_line = "MCF10A_CNTRL";
exp_type = "monolayer";

stiffs = newArray("soft", "medium", "stiff");
sampleid_per_stiff = newArray(1, 1, 1);
// sample_zsize = newArray(7, 10, 10);
num_tiles = 6;
waves = newArray("740nm", "880nm");




// function to find z slice
function num_slices(dir) 
{ 
	
	list = getFileList(dir);
	count = 0;
	
    for (i=0; i<list.length; i++) 
    {
    	
    	if (endsWith(list[i], "Intensity.tif"))
    	{
       		count++;
    	}

     }
     
     return count;
	
}

// loop through, read, concatenate, save
for (i = 0; i < stiffs.length; i++) 
{
	
	for (j = 0; j < sampleid_per_stiff[i]; j++) 
	{
		
		for (k = 0; k < waves.length; k++)
		{
			
			// define read and save location paths,names
			sample_folder = data_folder + "/" + stiffs[i] + "/" + exp_type + (j+1) + "/" + waves[k] + "/";
			sample_name = cell_line + "_" + stiffs[i] + "_" + exp_type + (j+1) +  "_" + waves[k] + "_";
			
			extraction_path = sample_folder + "/extracted_images";
			File.makeDirectory(extraction_path);
			
			// find how mnay total z slices there are in the directory
			z_slices = num_slices(sample_folder);
			
			// find how many z slices there are in each tile 
			z_per_tile = z_slices / num_tiles;
			 
			
			// concatonate into stacks
			for (tile = 1; tile <= num_tiles; tile++)
			{
				
				concat_string = "";
				for (z = 1; z <= z_per_tile; z++)
				{
	
					z_eff = ((tile - 1) * z_per_tile + z);
					open(sample_folder + sample_name + "z" + z_eff + "-Intensity.tif");
					
					concat_string = concat_string + "image" + z + "=" + sample_name + "z" + z_eff + "-Intensity.tif ";
					
				}
				
				// print("  title=" + sample_name + "stack open " + concat_string);
				run("Concatenate...", "  title=" + sample_name + "stack open " + concat_string);
				saveAs("Tiff", extraction_path + "/" + sample_name + "tile" + tile + "_stack.tif");
				// run("Z Project...", "projection=[Max Intensity]");
				// saveAs("Tiff", sample_folder + sample_name + "MAX.tif");
				
				run("Close All");
			
			}
			
			
			// run("Grid/Collection stitching", "type=[Grid: row-by-row] order=[Right & Down                ] grid_size_x=6 grid_size_y=1 tile_overlap=10 first_file_index_ii=1 directory=" + extraction_path + " file_names=" + sample_name + "tile{i}_stack.tif" + " output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");

			
		}
			
	}

}