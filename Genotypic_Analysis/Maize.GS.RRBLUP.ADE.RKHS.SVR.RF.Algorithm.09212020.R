

diploid.polyploid.5cfv <- function(phenames=NULL, diploid.Geno=NULL, m_train.pheno=NULL, species=NULL, k=NULL, cycles = NULL,
                                               number.of.folds=NULL, date=NULL){
  
  Prediction.Accuracy.By.Trait.NULL=NULL
  Intercept.By.Trait.NULL=NULL
  Slope.By.Trait.NULL=NULL
  CI.By.Trait.NULL=NULL
  
  for(l in 1:length(phenames)){ #length(phenames)#1:length(phenames)
    
    
    print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], "  !!!!!!!----------", sep = ""))
    
    Pheno1=as.matrix(m_train.pheno[,phenames[l]])
    rownames(Pheno1) <- m_train.pheno$Taxa
    colnames(Pheno1) <- phenames[l]
    
    r.gy.diplo <- matrix(NA, cycles, k)
    r.gy.RK.diplo <- matrix(NA, cycles, k)
    r.gy.SV.diplo <- matrix(NA, cycles, k)
    r.gy.RF.diplo <- matrix(NA, cycles, k)
    
    RB.Intercept <- matrix(NA, cycles, k)
    RK.Intercept <- matrix(NA, cycles, k)
    SV.Intercept <- matrix(NA, cycles, k)
    RF.Intercept <- matrix(NA, cycles, k)
    
    RB.Slope <- matrix(NA, cycles, k)
    RK.Slope <- matrix(NA, cycles, k)
    SV.Slope <- matrix(NA, cycles, k)
    RF.Slope <- matrix(NA, cycles, k)
    
    RB.CI <- matrix(NA, cycles, k)
    RK.CI <- matrix(NA, cycles, k)
    SV.CI <- matrix(NA, cycles, k)
    RF.CI <- matrix(NA, cycles, k)
    

        for(r in 1:cycles){
      
      print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles, " cycles !!!!!!!----------", sep = "")) 
      
      Pheno1.TRS1 <- data.frame(Pheno1)
      Pheno1.TRS1$id <- sample(1:k, nrow(Pheno1.TRS1), replace = TRUE)
      
    list <- 1:k
    
    
    for (f in 1:k){
      print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles," --- CrossFolds ", f," !!!!!!!----------", sep = ""))
      pheno_trainingset=NULL
      pheno_testset=NULL
      Geno_trainingset=NULL
      Geno_testset=NULL
      # remove rows with id i from dataframe to create training set
      # select rows with id i to create test set
      #set.seed()
      pheno_trainingset <- subset(Pheno1.TRS1, id %in% list[-f])
      
      #colnames(pheno_trainingset)[1] <- phenames[l]
    
      
      pheno_testset <- subset(Pheno1.TRS1, id %in% list[f])
      #colnames(pheno_testset) <- phenames[l]
      
      ##################################################
      Geno_trainingset.2 <- diploid.Geno[match(rownames(pheno_trainingset), rownames(diploid.Geno)),]
      Geno_trainingset.2 <- as.matrix(Geno_trainingset.2)
      
      
      #pheno_testset <- subset(Pheno1.TRS1, id %in% list[f])
      #colnames(pheno_testset) <- phenames[l]
      Geno_testset.2 <- diploid.Geno[match(rownames(pheno_testset), rownames(diploid.Geno)),]
      Geno_testset.2 <- as.matrix(Geno_testset.2)
      
       # Perform GS for each training size
        ww <- pheno_testset # delete data for 1/5 of the population
        
        library(rrBLUP)
        print(paste("--------- Initializing Analysis.....Trait to analyse: ", phenames[l], " --- cycle ", r, " out of ", cycles," --- CrossFolds ", f," RRBLUP GS ", sep = ""))
        
        Geno_trainingset.2.rrblup <- Geno_trainingset.2-1
        Yt=(pheno_trainingset[,1])
        rrMod.Ro.2<-mixed.solve(Yt, X=NULL, Z=Geno_trainingset.2.rrblup, K=NULL, SE=F, return.Hinv=F)
        mEff.Ro.2<-rrMod.Ro.2$u
        e.Ro.2= as.matrix(mEff.Ro.2)
        
        Geno_testset.2.rrblup <- Geno_testset.2-1
        predYv.Ro.2 = Geno_testset.2.rrblup%*%e.Ro.2
        
        predYr.Ro.2 = predYv.Ro.2[,1]+ rrMod.Ro.2$beta
        
        Y_valid.2=pheno_testset[,1]
        
        r.gy.diplo[r,f] <-  cor(predYr.Ro.2,Y_valid.2,use="complete")
        the.fitted.model.RB <- lm(ww[,1] ~ predYr.Ro.2)
        the.coefficients.RB <- c(the.fitted.model.RB$coefficients[1], the.fitted.model.RB$coefficients[2])
        RB.Intercept[r,f] <- the.coefficients.RB[1]
        RB.Slope[r,f] <- the.coefficients.RB[2]
        
        these.observed.and.predicted.phenotypic.values.RRBLUP.diplo <- data.frame(rownames(pheno_testset), Y_valid.2, predYr.Ro.2)
        rownames(these.observed.and.predicted.phenotypic.values.RRBLUP.diplo) <- NULL
        colnames(these.observed.and.predicted.phenotypic.values.RRBLUP.diplo) <- c("Taxa", "Observed.Value", "Predicted.Value")
        x.p.RRBLUP.diplo=these.observed.and.predicted.phenotypic.values.RRBLUP.diplo[,c(1,3)]
        y.o.RRBLUP.diplo=these.observed.and.predicted.phenotypic.values.RRBLUP.diplo[,c(1,2)]
        write.table(these.observed.and.predicted.phenotypic.values.RRBLUP.diplo, paste("these.observed.and.predicted.phenotypic.values.RRBLUP.diplo.Trait", phenames[l], "Validation.Set", f,"Cycle", r, "txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
        RB.CI[r,f] <- round(CI(x.p.RRBLUP.diplo,y.o.RRBLUP.diplo,s=0.2,top=T),2)
        
        # ADE
        
        Pheno.comb.1 <- rbind(pheno_trainingset, pheno_testset)
        train.comb.1 <- 1:nrow(pheno_trainingset)
        valid.comb.1 <-setdiff(1:nrow(Pheno.comb.1),train.comb.1)
        
        y.trn.1 <- Pheno.comb.1[,1] # for prediction accuracy
        
        Pheno.comb.2 <- Pheno.comb.1
        
        valid.comb.2 <- valid.comb.1
        
        y.trn.2 <- Pheno.comb.2[,1] # for prediction accuracy
        
        ww <- pheno_testset # delete data for 1/5 of the population

        y.trn.2[valid.comb.2] <- NA #valid.comb
        
        Geno.comb.2 <- rbind(Geno_trainingset.2, Geno_testset.2)
        Geno.comb.2 <- as.matrix(Geno.comb.2)
        #Geno.comb.2 <- Geno.comb.2-1
        
        
        Y_valid=pheno_testset
        
        #RKHS
        library(BGLR)
        X.2 <- Geno.comb.2#diploid.Geno[match(rownames(Geno.comb.2), rownames(diploid.Geno)),]
        #X.2 <- Geno.comb.2
        #X.2 <-scale(X.2 ,center=TRUE,scale=TRUE)
        M.2 <-tcrossprod(X.2)/ncol(X.2)
        y.trn.2[valid.comb.2] <- NA
        
        ETA.2 <-list(list(K=M.2,model='RKHS')) 
        fm.RK.2<-BGLR(y=y.trn.2,ETA=ETA.2,response_type="gaussian", nIter=12000, burnIn=2000)
        r.gy.RK.diplo[r,f] <- cor(fm.RK.2$yHat[valid.comb.2], ww[,1], use="complete")
        the.fitted.model.RK <- lm(ww[,1] ~ fm.RK.2$yHat[valid.comb.2])
        the.coefficients.RK <- c(the.fitted.model.RK$coefficients[1], the.fitted.model.RK$coefficients[2])
        RK.Intercept[r,f] <- the.coefficients.RK[1]
        RK.Slope[r,f] <- the.coefficients.RK[2]
        
        these.observed.and.predicted.phenotypic.values.RK.diplo <- data.frame(rownames(pheno_testset), ww[,1], fm.RK.2$yHat[valid.comb.2])
        rownames(these.observed.and.predicted.phenotypic.values.RK.diplo) <- NULL
        colnames(these.observed.and.predicted.phenotypic.values.RK.diplo) <- c("Taxa", "Observed.Value", "Predicted.Value")
        x.p.RK.diplo=these.observed.and.predicted.phenotypic.values.RK.diplo[,c(1,3)]
        y.o.RK.diplo=these.observed.and.predicted.phenotypic.values.RK.diplo[,c(1,2)]
        write.table(these.observed.and.predicted.phenotypic.values.RK.diplo, paste("these.observed.and.predicted.phenotypic.values.RK.diplo.Trait", phenames[l], "Validation.Set", f,"Cycle", r, "txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
        RK.CI[r,f] <- round(CI(x.p.RK.diplo,y.o.RK.diplo,s=0.2,top=T),2)
       
        #SVR
        print(paste("--------- Trait being analysed: ", phenames[l], " --cycle-- ", r, " ; GS Model: SVR"," !!!!!!!----------", sep = ""))
        library(e1071)
        
        svm.model <- svm(Geno_trainingset.2, Yt, type='eps-regression', kernel='radial', scale=FALSE)
        pred <- predict(svm.model, Geno_testset.2, decision.values = TRUE)
        r.gy.SV.diplo[r,f] <- cor(pred, ww[,1],use="complete")
        the.fitted.model.SV <- lm(ww[,1] ~ pred)
        the.coefficients.SV <- c(the.fitted.model.SV$coefficients[1], the.fitted.model.SV$coefficients[2])
        SV.Intercept[r,f] <- the.coefficients.SV[1]
        SV.Slope[r,f] <- the.coefficients.SV[2]
        
        
        these.observed.and.predicted.phenotypic.values.SV.diplo <- data.frame(rownames(ww), ww[,1], pred)
        rownames(these.observed.and.predicted.phenotypic.values.SV.diplo) <- NULL
        colnames(these.observed.and.predicted.phenotypic.values.SV.diplo) <- c("Taxa", "Observed.Value", "Predicted.Value")
        x.p.SV.diplo=these.observed.and.predicted.phenotypic.values.SV.diplo[,c(1,3)]
        y.o.SV.diplo=these.observed.and.predicted.phenotypic.values.SV.diplo[,c(1,2)]
        write.table(these.observed.and.predicted.phenotypic.values.SV.diplo, paste("these.observed.and.predicted.phenotypic.values.SV.diplo.Trait", phenames[l], "Validation.Set", f,"Cycle", r, "txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
        SV.CI[r,f] <- round(CI(x.p.SV.diplo,y.o.SV.diplo,s=0.2,top=T),2)
        

        
        print(paste("--------- Trait being analysed: ", phenames[l], " --cycle-- ", r, " ; GS Model: rForest"," !!!!!!!----------", sep = ""))
        library(randomForest)
        Yt <- as.numeric(Yt)
        rf <- randomForest(Geno_trainingset.2, Yt, tree=5000, mtry=floor(sqrt(ncol(Geno_trainingset.2))), 
          replace=TRUE, classwt=NULL, nodesize = 1, importance=FALSE, localImp=FALSE, nPerm=1,norm.votes=TRUE, do.trace=FALSE)
        
        pred.rf = predict(rf, newdata=Geno_testset.2)
        r.gy.RF.diplo[r,f] <- cor(pred.rf, ww[,1],use="complete")
        the.fitted.model.RF <- lm(ww[,1] ~ pred.rf)
        the.coefficients.RF <- c(the.fitted.model.RF$coefficients[1], the.fitted.model.RF$coefficients[2])
        RF.Intercept[r,f] <- the.coefficients.RF[1]
        RF.Slope[r,f] <- the.coefficients.RF[2]
        
        
        these.observed.and.predicted.phenotypic.values.RF.diplo <- data.frame(rownames(ww), ww[,1], pred.rf)
        rownames(these.observed.and.predicted.phenotypic.values.RF.diplo) <- NULL
        colnames(these.observed.and.predicted.phenotypic.values.RF.diplo) <- c("Taxa", "Observed.Value", "Predicted.Value")
        x.p.RF.diplo=these.observed.and.predicted.phenotypic.values.RF.diplo[,c(1,3)]
        y.o.RF.diplo=these.observed.and.predicted.phenotypic.values.RF.diplo[,c(1,2)]
        write.table(these.observed.and.predicted.phenotypic.values.RF.diplo, paste("these.observed.and.predicted.phenotypic.values.RF.diplo.Trait", phenames[l], "Validation.Set", f,"Cycle", r, "txt",sep="."), sep="\t", quote=F, row.names=F, col.names=T)
        RF.CI[r,f] <- round(CI(x.p.RF.diplo,y.o.RF.diplo,s=0.2,top=T),2)
        

    } # End of i in k cross valid
    

    
    } # End of cycle
    
          r.gy.diplo.df <- as.data.frame(r.gy.diplo)
          r.gy.diplo.df$Model <- rep("RRBLUP", nrow(r.gy.diplo.df))
          
          r.gy.RK.diplo.df <- as.data.frame(r.gy.RK.diplo)
          r.gy.RK.diplo.df$Model <- rep("RKHS", nrow(r.gy.RK.diplo.df))
          
          r.gy.SV.diplo.df <- as.data.frame(r.gy.SV.diplo)
          r.gy.SV.diplo.df$Model <- rep("SVR", nrow(r.gy.SV.diplo.df))
          
          r.gy.RF.diplo.df <- as.data.frame(r.gy.RF.diplo)
          r.gy.RF.diplo.df$Model <- rep("rForest", nrow(r.gy.RF.diplo.df))
          
          Prediction.Accuracy.By.Trait <- rbind(r.gy.diplo.df, r.gy.RK.diplo.df, r.gy.SV.diplo.df, r.gy.RF.diplo.df)
          Prediction.Accuracy.By.Trait <- as.data.frame(Prediction.Accuracy.By.Trait)
          Prediction.Accuracy.By.Trait$Trait <- rep(phenames[l], nrow(Prediction.Accuracy.By.Trait))
          write.table(Prediction.Accuracy.By.Trait, paste(species, "Prediction.RK.AD.RRBLUP.SVR.RF",phenames[l], date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
          
          Prediction.Accuracy.By.Trait.NULL <- rbind(Prediction.Accuracy.By.Trait.NULL, Prediction.Accuracy.By.Trait)
          
          RB.Intercept.df <- as.data.frame(RB.Intercept)
          RB.Intercept.df$Model <- rep("RRBLUP", nrow(RB.Intercept.df))
          
          RK.Intercept.df <- as.data.frame(RK.Intercept)
          RK.Intercept.df$Model <- rep("RKHS", nrow(RK.Intercept.df))
          
          SV.Intercept.df <- as.data.frame(SV.Intercept)
          SV.Intercept.df$Model <- rep("SVR", nrow(SV.Intercept.df))
          
          RF.Intercept.df <- as.data.frame(RF.Intercept)
          RF.Intercept.df$Model <- rep("rForest", nrow(RF.Intercept.df))
          
          Intercept.By.Trait <- rbind(RB.Intercept.df, RK.Intercept.df, SV.Intercept.df, RF.Intercept.df)
          Intercept.By.Trait <- as.data.frame(Intercept.By.Trait)
          Intercept.By.Trait$Trait <- rep(phenames[l], nrow(Intercept.By.Trait))
          Intercept.By.Trait.NULL <- rbind(Intercept.By.Trait.NULL, Intercept.By.Trait)
          write.table(Intercept.By.Trait, paste(species, "Intercept.RK.AD.RRBLUP.SVR.RF",phenames[l], date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
          
          
          
          RB.Slope.df <- as.data.frame(RB.Slope)
          RB.Slope.df$Model <- rep("RRBLUP", nrow(RB.Slope.df))
          
          RK.Slope.df <- as.data.frame(RK.Slope)
          RK.Slope.df$Model <- rep("RKHS", nrow(RK.Slope.df))
          
          SV.Slope.df <- as.data.frame(SV.Slope)
          SV.Slope.df$Model <- rep("SVR", nrow(SV.Slope.df))
          
          RF.Slope.df <- as.data.frame(RF.Slope)
          RF.Slope.df$Model <- rep("rForest", nrow(RF.Slope.df))
          
          Slope.By.Trait <- rbind(RB.Slope.df, RK.Slope.df, SV.Slope.df, RF.Slope.df)
          Slope.By.Trait <- as.data.frame(Slope.By.Trait)
          Slope.By.Trait$Trait <- rep(phenames[l], nrow(Slope.By.Trait))
          Slope.By.Trait.NULL <- rbind(Slope.By.Trait.NULL, Slope.By.Trait)
          write.table(Slope.By.Trait, paste(species, "Slope.RK.AD.RRBLUP.SVR.RF",phenames[l], date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
          
 
          RB.CI.df <- as.data.frame(RB.CI)
          RB.CI.df$Model <- rep("RRBLUP", nrow(RB.CI.df))
          
          RK.CI.df <- as.data.frame(RK.CI)
          RK.CI.df$Model <- rep("RKHS", nrow(RK.CI.df))
          
          SV.CI.df <- as.data.frame(SV.CI)
          SV.CI.df$Model <- rep("SVR", nrow(SV.CI.df))
          
          RF.CI.df <- as.data.frame(RF.CI)
          RF.CI.df$Model <- rep("rForest", nrow(RF.CI.df))
          
          CI.By.Trait <- rbind(RB.CI.df, RK.CI.df, SV.CI.df, RF.CI.df)
          CI.By.Trait <- as.data.frame(CI.By.Trait)
          CI.By.Trait$Trait <- rep(phenames[l], nrow(CI.By.Trait))
          CI.By.Trait.NULL <- rbind(CI.By.Trait.NULL, CI.By.Trait)
          
          write.table(CI.By.Trait, paste(species, "CI.RK.AD.RRBLUP.SVR.RF",phenames[l], date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
          
    
  } # End of trait
  
  print("--------- Writing Final REsults In Working Directory ------------")
  
  
  
  write.table(Prediction.Accuracy.By.Trait.NULL, paste(species, "Prediction.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(Intercept.By.Trait.NULL, paste(species, "Intercept.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(Slope.By.Trait.NULL, paste(species, "Slope.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  write.table(CI.By.Trait.NULL, paste(species, "CI.Traits", date,"csv", sep='.'), sep=",", quote=F, row.names=F, col.names=T)
  
  
} # End of library
  # 09212020 11:47 AM


