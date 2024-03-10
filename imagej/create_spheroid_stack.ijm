

// define parameters
data_folder = "/Users/haider/Documents/tmr/data/Olympus MPE-RS TWIN/FLIM/Data/20230203/tifs";
stiff = newArray("-oxphos", "-glyco", "control");
sample_per_stiff = newArray(2,2,2);
area_per_sample = newArray(2,3, 1,3, 3,3);
cell_line = "MCF10A_CNTRL";
exp_type = "monolayer";


// function to find z slice
function num_zs(dir) 
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

for (i = 0; i < stiff.length; i++) 
{
	
	for (j = 0; j < sample_per_stiff[i]; j++) 
	{
		
		for (k = 0; k < area_per_sample[(i*2)+j]; k++) 
		{
			
			sample_folder = data_folder + "/" + stiff[i] + "/" + exp_type + (j+1) + "/area" + (k+1) + "/740nm/";
			sample_name = cell_line + "_-D_" + stiff[i] + "_" + exp_type + (j+1) + "_area" + (k+1) + "_740nm_";
			
			
			z_slices = num_zs(sample_folder);
			print(sample_folder);
			print(z_slices);
			
			concat_string = "";
			for (z = 1; z < z_slices + 1; z++)
			{

				open(sample_folder + sample_name + "z" + z + "-Intensity.tif");
				
				concat_string = concat_string + "image" + z + "=" + sample_name + "z" + z + "-Intensity.tif ";
				
			}
			
			run("Concatenate...", "  title=" + sample_name + "stack open " + concat_string);
			saveAs("Tiff", sample_folder + sample_name + "stack.tif");
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("Tiff", sample_folder + sample_name + "MAX.tif");
			
			run("Close All");
			
			
			
		}

	}

}