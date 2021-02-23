##### This script generates point clouds and compute distance in parallel
rm(list = ls())
#Install the packages if needed
#Set the directory to RscriptSfpca
require(foreach)
require(doParallel)
require(svMisc)
require(matlabr)

choose_Scenario = 2 ### 1, 2 ,3 , the choice must be coherent with run_from_r.m TO BE SET
#### To handle quietly CloudCompare in command line mode
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}

#### Call matlab and generate Point Clouds
setwd("../MatlabScriptEgg2D/")
run_matlab_script("run_from_r.m",verbose = FALSE,invisible = T)
setwd("../RscriptSfpca/")

#### Choose the defect
defect_names = c("OOC_WS", "OOC_MS", "OOC_IRR") ### Wide, missing, irregular (Scenario I, II, III)

root_name = defect_names[choose_Scenario]


### This is support for CloudCompare
namevec1 = NULL;
namevec2 = NULL;
namevec3 = NULL;
nomname = "../MatlabScriptEgg2D/PointClouds/Nominal_model.csv"
Ncontrol = 95
Ndef1 = 5
Ntot = Ncontrol + Ndef1
for (i in 1:Ncontrol)
{
  namevec1 = c(namevec1,paste0("../MatlabScriptEgg2D/PointClouds/In_control_",i,".csv"));
  namevec2 = c(namevec2,paste0("../CCsupport/forw",i,".csv"));
  namevec3 = c(namevec3,paste0("../CCsupport/back",i,".csv"));
}
for ( i in 1:Ndef1)
{
  namevec1 = c(namevec1,paste0("../MatlabScriptEgg2D/PointClouds/",root_name,i,".csv"));
  namevec2 = c(namevec2,paste0("../CCsupport/forw",i+Ncontrol,".csv"));
  namevec3 = c(namevec3,paste0("../CCsupport/back",i+Ncontrol,".csv"));
}
### Prepare parallel environment
cl = makeCluster(detectCores(logical = F))
registerDoParallel(cl)

#### Call CloudCompare in parallel, compute distance maps and save them
a = foreach( i = (1:Ntot)) %dopar%
  {
    quiet(system(paste0("../CloudCompareBinaries/CloudCompare -SILENT -AUTO_SAVE OFF -C_EXPORT_FMT ASC  -O ",namevec1[i]," -O ",nomname," -C2C_DIST -POP_CLOUDS -SAVE_CLOUDS FILE ",namevec2[i]," -POP_CLOUDS -O ",nomname," -O ",namevec1[i]," -C2C_DIST -POP_CLOUDS -SAVE_CLOUDS FILE ",namevec3[i]," -CLEAR_CLOUDS")))
  }

#### Put distance maps in a list
alldistances = foreach ( i = 1:Ntot)%dopar%
{
  return(list( forward = read.csv(file = namevec2[i],header = F,sep = "",colClasses = c("NULL","NULL","NULL","numeric"))[,1],backward = read.csv(file = namevec3[i],header = F,sep = "",colClasses = c("NULL","NULL","NULL","numeric"))[,1]))
}
stopCluster(cl)
#### Save the list, (almost) ready for PCA
save(alldistances,file = "alldistances.rdata")
sapply(namevec1,unlink)
sapply(namevec2, unlink)
sapply(namevec3,unlink)
