###  generate openmx models in parallel
###  by number of connections

type file{}

# produce .RData file for each model object to be tested

(file mxModels[]) mxmodel_generator (int conn, int ncol, float wt, file dot_r){
    app{
	RInvoke @filename(dot_r) conn ncol wt;
    }
}

# map input and output for generator

file modelgen<single_file_mapper;file="scripts/mxModelGen.R">;

int mxsize = 9;
int nconnections[] = [2:(mxsize-3)];
int initweight = .75;
int numcol = 3;

foreach c in nconnections {
	file potentialmodels[]<ext;exec="./modmap.py", size=mxsize, conn=c>;
	(potentialmodels) = mxmodel_generator(c, numcol, initweight, modelgen);
	}



