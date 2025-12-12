

ML.Tables <- list()

for(class in classes){
  depVars <- All.class.data[[variable]][[level]][[class]]
  if(variable == "ampl_or_del"){
    depVars <- depVars[,c("variable","amp.freq","del.freq","type")]
    colnames(depVars) <- c("bin","ampl_score","del_score","Type")
    ml.table <- merge(depVars, features, by = c("bin","Type"))
    ML.Tables[[variable]][[class]] <- ml.table
  }else{
    depVars <- depVars[,c("variable","Total.CNA.Freq","type")]
    colnames(depVars) <- c("bin","Total_score","Type")
    ml.table <- merge(depVars, features, by = c("bin","Type"))
    ML.Tables[[variable]][[class]] <- ml.table
    }

}

