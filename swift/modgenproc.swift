###  generate openmx models in parallel
###  by number of connections

type file{}

(file min) mxmodel_processor (file cov_matrix, file dot_r, int numcol, int modnum){
    app{
	RInvoke @filename(dot_r) @filename(cov_matrix) numcol modnum;
    }
}

file cov_script<single_file_mapper; file="scripts/singlemodels.R">;
file covmx<single_file_mapper;file="matrices/gestspeech.cov">;

float initweight = .75;
int numcol = 4;
int totalperms[] = [0:65534];

foreach perm in totalperms{
		file modmin<single_file_mapper; file=@strcat("results/",perm,".min")>;
		(modmin) = mxmodel_processor(covmx, cov_script, numcol, perm);		
		}






