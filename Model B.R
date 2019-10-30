#Prediction Kebara 2 from N=64 H. sapiens: MODEL B
require(Morpho)
require(geomorph)
require(Rvcg)
require(ape)
require(rgl)

#Load landmark data----
Trunk_data <- read.morphologika("E:/R progress/Trunk_analyses_R/Kebara2_prediction/N64_trunks.txt")
Pelves_Kebara2 <- read.morphologika("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_pelves.txt")
Kebara2_thorax_reconstructions <- read.morphologika("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_thorax_reconstructions.txt")

#Resliding trunk data ----
#Define landmarks and semilandmarks----
fix <- c(1:134)
Curve1 <- c(1,135:147,46)
Curve2 <- c(23,148:160,68)
Curve3 <- c(3,161:173,47)
Curve4 <- c(25,174:186,69)
Curve5 <- c(5,187:199,49)
Curve6 <- c(27,200:212,71)
Curve7 <- c(7,213:225,51)
Curve8 <- c(29,226:238,73)
Curve9 <- c(9,239:251,53)
Curve10 <- c(31,252:264,75)
Curve11 <- c(11,265:277,55)
Curve12 <- c(33,278:290,77)
Curve13 <- c(13,291:303,57)
Curve14 <- c(35,304:316,79)
Curve15 <- c(15,317:329,59)
Curve16 <- c(37,330:342,81)
Curve17 <- c(17,343:355,61)
Curve18 <- c(39,356:368,83)
Curve19 <- c(19,369:381,63)
Curve20 <- c(41,382:394,85)
Curve21 <- c(21,395:407,65)
Curve22 <- c(43,408:420,87)
Curve23 <- c(2,421:433,45)
Curve24 <- c(24,434:446,67)
Curve25 <- c(4,447:459,48)
Curve26 <- c(26,460:472,70)
Curve27 <- c(6,473:485,50)
Curve28 <- c(28,486:498,72)
Curve29 <- c(8,499:511,52)
Curve30 <- c(30,512:524,74)
Curve31 <- c(10,525:537,54)
Curve32 <- c(32,538:550,76)
Curve33 <- c(12,551:563,56)
Curve34 <- c(34,564:576,78)
Curve35 <- c(14,577:589,58)
Curve36 <- c(36,590:602,80)
Curve37 <- c(16,603:615,60)
Curve38 <- c(38,616:628,82)
Curve39 <- c(18,629:641,62)
Curve40 <- c(40,642:654,84)
Curve41 <- c(20,655:667,64)
Curve42 <- c(42,668:680,86)
Curve43 <- c(22,681:693,66)
Curve44 <- c(44,694:706,88)
Curve45 <- c(89,707:722,128)
Curve46 <- c(106,723:737,108)
Curve47 <- c(107,738:752,109)
Curve48 <- c(110,753:760,112)
Curve49 <- c(111,761:768,113)
Curve50 <- c(114,769:776,118)
Curve51 <- c(115,777:784,119)
Curve52 <- c(108,785:790,124)
Curve53 <- c(109,791:796,125)
Curve54 <- c(116,797:802,122)
Curve55 <- c(117,803:808,123)
Curve56 <- c(114,809:826,120)
Curve57 <- c(115,827:844,121)
Curve58 <- c(126,845:849,118)
Curve59 <- c(127,850:854,119)
Curve60 <- c(118,855:862,133)
Curve61 <- c(119,863:870,134)

Curves <- list(Curve1,Curve2,Curve3,Curve4,Curve5,Curve6,Curve7,Curve8,Curve9,Curve10,Curve11,Curve12,Curve13,Curve14,Curve15,Curve16,Curve17,Curve18,Curve19,Curve20,Curve21,Curve22,Curve23,Curve24,Curve25,Curve26,Curve27,Curve28,Curve29,Curve30,Curve31,Curve32,Curve33,Curve34,Curve35,Curve36,Curve37,Curve38,Curve39,Curve40,Curve41,Curve42,Curve43,Curve44,Curve45,Curve46,Curve47,Curve48,Curve49,Curve50,Curve51,Curve52,Curve53,Curve54,Curve55,Curve56,Curve57,Curve58,Curve59,Curve60,Curve61)

Surfaces <- c(1:1030)[-c(fix,Curve1,Curve2,Curve3,Curve4,Curve5,Curve6,Curve7,Curve8,Curve9,Curve10,Curve11,Curve12,Curve13,Curve14,Curve15,Curve16,Curve17,Curve18,Curve19,Curve20,Curve21,Curve22,Curve23,Curve24,Curve25,Curve26,Curve27,Curve28,Curve29,Curve30,Curve31,Curve32,Curve33,Curve34,Curve35,Curve36,Curve37,Curve38,Curve39,Curve40,Curve41,Curve42,Curve43,Curve44,Curve45,Curve46,Curve47,Curve48,Curve49,Curve50,Curve51,Curve52,Curve53,Curve54,Curve55,Curve56,Curve57,Curve58,Curve59,Curve60,Curve61)]

proc_trunk <- procSym(Trunk_data)

relaxarray <- Trunk_data

for (i in 1:64) {
  print(i)
  
  relaxarray[,,i] <- relaxLM(Trunk_data[,,i],proc_trunk$mshape,SMvector=fix, deselect=TRUE, surp=Surfaces,outlines=Curves)
  
}

Todo_trunk <- relaxarray[,,1:64]

# Thorax and pelvis n=64 H. sapiens
Thorax_temp <- read.morphologika("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_temp_mlk_stern.txt")
Thorax_indices <- vcgKDtree(Trunk_data[,,36],Thorax_temp[,,1],k=1)$index
Thorax <- Todo_trunk[Thorax_indices,,]
Pelvis_temp <- read.morphologika("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Pelvis_temp_mlk.txt")
Pelvis_indices <- vcgKDtree(Todo_trunk[,,36],Pelvis_temp[,,1],k=1)$index
Pelvis <- Todo_trunk[Pelvis_indices,,]

Thorax_mesh <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Pelvis_mesh <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Pelvis_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)


# GPA DE 64 H. SAPIENS 

# Procrustes Superimposition (GPA) of the two blocks separately

Thorax_GPA_Homo <- procSym(Thorax)
Pelvis_GPA_Homo <- procSym(Pelvis)

#MODEL B: PLS Homo sapiens sample corrected for the effect of population and sexual dimorphism ----

# Thorax and pelvis pooled by pop

Thorax_pop <- name2factor(Thorax_GPA_Homo$rotated,which=4)
Pelvis_pop <- name2factor(Pelvis_GPA_Homo$rotated,which=4)

vecdata_Thorax <- vecx(Thorax_GPA_Homo$orpdata) #Transforma 3D array en 2D matrix
pooledcov_Thorax <- covW(vecdata_Thorax,Thorax_pop) ## Calcula pooled within population covariance matrix

vecdata_pelvis <- vecx(Pelvis_GPA_Homo$orpdata) #Transforma 3D array en 2D matrix
pooledcov_Pelvis <- covW(vecdata_pelvis,Pelvis_pop)

Thorax_corr_pop <- lm(vecdata_Thorax ~ Thorax_pop)$residuals ## This is the data corrected for group means
Pelvis_corr_pop <- lm(vecdata_pelvis ~ Pelvis_pop)$residuals

Thorax_corr_pop_array <- vecx(Thorax_corr_pop,revert=TRUE,lmdim=3)
Pelvis_corr_pop_array <- vecx(Pelvis_corr_pop,revert=TRUE,lmdim =3)

Thorax_ad_corrected <- sweep(Thorax_corr_pop_array,1:2,-Thorax_GPA_Homo$mshape)
Pelvis_ad_corrected <- sweep(Pelvis_corr_pop_array,1:2,-Pelvis_GPA_Homo$mshape)

Thorax_GPA_corr_pop <- procSym(Thorax_ad_corrected)
Pelvis_GPA_corr_pop <- procSym(Pelvis_ad_corrected)

# Thorax and pelvis pooled by pop and sex

Thorax_sex <- name2factor(Thorax_GPA_Homo$rotated,which=3)
Pelvis_sex <- name2factor(Pelvis_GPA_Homo$rotated,which=3)

vecdata_Thorax <- vecx(Thorax_GPA_corr_pop$orpdata) #Transforma 3D array en 2D matrix
pooledcov_Thorax <- covW(vecdata_Thorax,Thorax_sex) ## Calcula pooled within population covariance matrix

vecdata_pelvis <- vecx(Pelvis_GPA_corr_pop$orpdata) #Transforma 3D array en 2D matrix
pooledcov_Pelvis <- covW(vecdata_pelvis,Pelvis_sex)

Thorax_corr_sex <- lm(vecdata_Thorax ~ Thorax_sex)$residuals ## This is the data corrected for group means
Pelvis_corr_sex <- lm(vecdata_pelvis ~ Pelvis_sex)$residuals

Thorax_corr_sex_array <- vecx(Thorax_corr_sex,revert=TRUE,lmdim=3)
Pelvis_corr_sex_array <- vecx(Pelvis_corr_sex,revert=TRUE,lmdim =3)

Thorax_ad_corrected <- sweep(Thorax_corr_sex_array,1:2,-Thorax_GPA_corr_pop$mshape)
Pelvis_ad_corrected <- sweep(Pelvis_corr_sex_array,1:2,-Pelvis_GPA_corr_pop$mshape)

Thorax_GPA_corr_sex <- procSym(Thorax_ad_corrected)
Pelvis_GPA_corr_sex <- procSym(Pelvis_ad_corrected)

#MODEL B: PLS Homo sapiens sex-corrected ----

PLS_ModelB <- pls2B(Thorax_GPA_corr_sex$rotated,Pelvis_GPA_corr_sex$rotated, rounds=999,cv=T)
PLS_ModelB 

col_plot_sex <- c((c(rep("indianred1",4),rep("mediumseagreen",6),rep("indianred1",1),rep("mediumseagreen",2),rep("indianred1",2),rep("mediumseagreen",2),rep("indianred1",7),rep("mediumseagreen",7),rep("indianred1",3),rep("mediumseagreen",12),rep("indianred1",18))))
col_leg_sex <- c("indianred1","mediumseagreen")
plot(PLS_ModelB$Xscores[,1],PLS_ModelB$Yscores[,1], xlab="Block 1 (Thorax)",ylab="Block 2 (Pelvis)",main="TWO BLOCKS PLS ANALYSIS (1st latent variable)",pch=18,col=col_plot_sex,cex=2,ylim=c(-0.08,0.08),xlim=c(-0.12,0.12))
abline(h=0.00,v=0.00,col="black")
legend("bottomright", inset=.005, legend = c("Females","Males"), horiz=FALSE, cex=1.2,bty="n",fill=col_leg_sex, border=col_leg_sex)

Covar_ModelB <- plsCoVar(PLS_ModelB,i=1,sdx=1,sdy=1)

PLS1_ModelB_neg_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelB$x[,,1])
shade3d(PLS1_ModelB_neg_x, col="grey",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelB_pos_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelB$x[,,2])
shade3d(PLS1_ModelB_pos_x, col="grey",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelB_neg_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelB$y[,,1])
shade3d(PLS1_ModelB_neg_y,col="grey",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelB_pos_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelB$y[,,2])
shade3d(PLS1_ModelB_pos_y,col="grey",back="lines",box=FALSE,axes=FALSE,specular="white")

#Error de predicción

truth_thorax <- vecx(Thorax_GPA_corr_sex$rotated) #Convert the 3D array into a matrix (Thorax)

error <- NULL
for (i in 1:dim(PLS_ModelB$predicted.x)[3]) # loop for cross-validation
  error <- c(error,mean((truth_thorax-PLS_ModelB$predicted.x[,,i])^2))

NoCVerror <- NULL
for (i in 1:dim(PLS_ModelB$predicted.x)[3]) # loop for MSE non-validated
  NoCVerror <- c(NoCVerror,mean((truth_thorax-vecx(predictPLSfromData(PLS_ModelB,y=Pelvis_GPA_corr_sex$rotated,ncomp=i)))^2))

# Cross-validation plot----

plot(1:62,error,main="Cross-validated prediction error",pch=19,xlab="Latent Variables used for Prediction",ylab="MSE",ylim=c(0,max(c(error,NoCVerror))))
lines(1:62,error)
points(1:62,NoCVerror,col="red",pch=20)                  
lines(1:62,NoCVerror,col="red")
legend(x="topright",y=NULL,legend=c("Cross-validated", "Non-validated"), col = c(1,"red"),pch=c(19,20),lty=1)

#PREDICTIONS based on MODEL B----

#Prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model B
Kebara2RA_toGPA_ModelB <- rotonto(Pelvis_GPA_corr_sex$mshape,Pelves_Kebara2[,,1],scale = T)$yrot
Pred_Kebara2RA_ModelB <- predictPLSfromData(PLS_ModelB,y=Kebara2RA_toGPA_ModelB,ncomp=5) #number of LVs to perform the prediction
deformGrid3d(Thorax_GPA_corr_sex$mshape,Pred_Kebara2RA_ModelB)
deformGrid3d(Pred_Kebara2RA_ModelB,Pred_Kebara2RA_ModelB)
Thorax_Kebara2RA_ModelB <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2RA_ModelB)

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model B with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2RA_ModelB,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelB <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelB <- rotonto(Thorax_GPA_corr_sex$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelB <- tps3d(Thorax_Kebara2GO_mesh_ModelB,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelB)
shade3d(Thorax_Kebara2GO_ModelB,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model B with thorax reconstruction Sawyer and Maley (2005)
shade3d(Thorax_Kebara2RA_ModelB,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_Sawyer_Maley.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelB <- rotonto(Thorax_GPA_corr_sex$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelB <- tps3d(Thorax_Kebara2S_M_mesh,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelB)
shade3d(Thorax_Kebara2S_M_ModelB,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#Prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model B
Kebara2SM_toGPA_ModelB <- rotonto(Pelvis_GPA_corr_sex$mshape,Pelves_Kebara2[,,2],scale = T)$yrot
Pred_Kebara2SM_ModelB <- predictPLSfromData(PLS_ModelB,y=Kebara2SM_toGPA_ModelB,ncomp=5) #number of LVs to perform the prediction
deformGrid3d(Thorax_GPA_corr_sex$mshape,Pred_Kebara2SM_ModelB)
deformGrid3d(Pred_Kebara2SM_ModelB,Pred_Kebara2SM_ModelB)
Thorax_Kebara2SM_ModelB <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2SM_ModelB)

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model B with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2SM_ModelB,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelB <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelB <- rotonto(Thorax_GPA_corr_sex$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelB <- tps3d(Thorax_Kebara2GO_mesh_ModelB,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelB)
shade3d(Thorax_Kebara2GO_ModelB,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model B with thorax reconstruction of Sawyer and Maley (2005)
shade3d(Thorax_Kebara2SM_ModelB,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh <- vcgImport("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_Sawyer_Maley.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelB <- rotonto(Thorax_GPA_corr_sex$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelB <- tps3d(Thorax_Kebara2S_M_mesh,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelB)
shade3d(Thorax_Kebara2S_M_ModelB,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#PLS scores of predictions MODEL B
groups64thorax <- read.table("E:/R progress/Trunk_analyses_R/Kebara2_prediction/Groups_N64sapiens.txt",header=TRUE)
col_plot_sex <- c((c(rep("indianred1",4),rep("mediumseagreen",6),rep("indianred1",1),rep("mediumseagreen",2),rep("indianred1",2),rep("mediumseagreen",2),rep("indianred1",7),rep("mediumseagreen",7),rep("indianred1",3),rep("mediumseagreen",12),rep("indianred1",18))))
col_leg_sex <- c("indianred1","mediumseagreen")
plot(PLS_ModelB$Xscores[,1],PLS_ModelB$Yscores[,1], xlab="Block 1 (Thorax)",ylab="Block 2 (Pelvis)",main="TWO BLOCKS PLS ANALYSIS (1st latent variable)",pch=16,col=col_plot_sex,cex=2,ylim=c(-0.08,0.08),xlim=c(-0.12,0.12))
abline(h=0.00,v=0.00,col="black")
legend("bottomright", inset=.005, legend = c("Females","Males"), horiz=FALSE, cex=1.2,bty="n",fill=col_leg_sex, border=col_leg_sex)
col_polyg <- c("black")

if (require(car)) {
  for(ii in 1:length(levels(groups64thorax$SPECIES))){
    dataEllipse(PLS_ModelB$Xscores[groups64thorax$SPECIES==levels(groups64thorax$SPECIES)[ii],1],PLS_ModelB$Yscores[groups64thorax$SPECIES==levels(groups64thorax$SPECIES)[ii],1],add=TRUE,center.pch = FALSE,levels=c(.95), col=col_polyg[ii],fill=F, fill.alpha=0.10,robust=TRUE)}
}

PLSXscore <- getPLSscores(PLS_ModelB,y=Kebara2RA_toGPA_ModelB)
PLSYscore <- getPLSscores(PLS_ModelB,x=Pred_Kebara2RA_ModelB)
points(PLSXscore[,1],PLSYscore[,1],pch=15,col="orange",cex=0.5)

PLSXscore <- getPLSscores(PLS_ModelB,y=Kebara2SM_toGPA_ModelB)
PLSYscore <- getPLSscores(PLS_ModelB,x=Pred_Kebara2SM_ModelB)
points(PLSXscore[,1],PLSYscore[,1],pch=15,col="purple",cex=0.5)
