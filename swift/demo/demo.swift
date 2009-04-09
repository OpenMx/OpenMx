type file;
type mxMin;
type Rscript;
type dbConnect;
type mxModel{
int modnum;
int dof;
string best;
}
# -------------- atomic procedures -------------- #

(file matrix) runQuery (dbConnect dbconn, string query, Rscript calcCov){
app{
mysqlPythonDBI query @calcCov @dbconn;
}
}

(external inserted) insertMxResult (dbConnect dbconn, string query, file datafile){
app{
mysqlPythonDBI query @dbconn stdout=@filename(inserted) @datafile;
}
}

(file min) mxModelProcessor ( file cov, Rscript mxModProc, int modnum, float weight, string cond, int net){
app{
RInvoke @mxModProc @filename(cov) modnum weight cond net;
}
}

# ------------ SEM user-defined procedures --------------- #

multiNetworkSEM(string condition,dbConnect emblemdb, dbConnect semdb, int n, string net, int totalperms[]){
float initweight = .75;
file covariance<single_file_mapper;file=@strcat("matrices/net",n,"/",condition,".cov")>;
covariance = getCovariance(condition, n, net, emblemdb);
Rscript mxModProc<single_file_mapper;file="scripts/singlemodels.R">;
foreach perm in totalperms{
file modmin<single_file_mapper;file=@strcat("results/net",n,"/",condition,"_",perm,".stat")>;
modmin =  mxModelProcessor(covariance,mxModProc,perm,initweight,condition,n);
external doneflag = insertOptMod(n, semdb, condition, modmin);

}

(external ins) insertOptMod(int net, dbConnect dbconn, string cond, file modfile){
string mysqlstr = @strcat("INSERT INTO optimized_models (network,deg_of_freedom,mx_minimum, modnum, cond) VALUES (",net,",DOF,BEST,MODNUM,",cond,");");
string argList = @strcat(
" --query ", mysqlstr,
" --data ", @filename(modfile),
" --conf ", @filename(dbconn));
ins = insertMxResult(dbconn, argList, modfile);
}

(file covariance) getCovariance (string cond, int net, string rois, dbConnect dbconn){
string mysqlstr = @strcat("SELECT avg(",cond,"0B), avg(",cond,"1B), avg(",cond,"2B),",
"avg(",cond,"3B), avg(",cond,"4B), avg(",cond,"5B),",
"avg(",cond,"6B), avg(",cond,"7B), avg(",cond,"8B) ",
"FROM emblemfemlh where roi in (",rois,") group by roi ");
string argList = @strcat(
" --conf ", "user.config",
" --query ", mysqlstr,
" --r_script ", "scripts/cov.R",
" --r_swift_args ", "matrices/net",net, "/", cond);
Rscript calcCov<single_file_mapper;file="scripts/cov.R">;
covariance = runQuery(dbconn, argList, calcCov);
}

(file plotfile) plotLogLik(int net, string cond, dbConnect dbconn){
Rscript rplot<single_file_mapper;file="scripts/plotloglik.R">;
string mysqlstr = @strcat("SELECT deg_of_freedom,mx_minimum FROM optimized_models",
" where network = ",net,";");
string argList = @strcat(
" --conf ", @filename(dbconn),
" --query ", mysqlstr,
" --r_script ", "scripts/plotloglik.R",
" --r_swift_args ", @filename(plotfile));
plotfile = runQuery(dbconn, argList, rplot);
}# ----------------------- Main ------------------------ #


string condition = "gestspeech";
string networks[] = ["42, 34, 33, 60", "42, 15, 60, 80", "33, 34, 23, 80"];
dbConnect emblemdb <single_file_mapper; file="./user.config">;
dbConnect semdb <single_file_mapper; file="./user2.config">;
int totalperms[] = [1:65536];
foreach net,n in networks{
multiNetworkSEM(condition,emblemdb,semdb,n,net,totalperms);
}
