# Run OpenMx optimizer on a list of given models
# to be found in a specified dir

type file{}

(file min) mxmodel_processor (file dot_r, file model){
    app{
	RInvoke @filename(dot_r) @filename(model);
    }
}

file models[]<filesys_mapper;location="models/", suffix=".rdata">;  

file opt_script<single_file_mapper; file="scripts/optim.R">;

foreach m in models{
		file modmin<single_file_mapper; file=@strcat("results/", @filename(m))>;
		(modmin) = mxmodel_processor(opt_script, m);		
		}






