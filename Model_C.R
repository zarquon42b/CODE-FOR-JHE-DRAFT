#Prediction Kebara 2 from N=29 male H. sapiens: MODEL C
require(Morpho)
require(geomorph)
require(Rvcg)
require(ape)
require(rgl)

#Load landmark data----
Trunk_data <- read.morphologika("N64_trunks.txt")
Pelves_Kebara2 <- read.morphologika("Kebara2_pelves.txt")
Kebara2_thorax_reconstructions <- read.morphologika("Kebara2_thorax_reconstructions.txt")
Trunk_males <- bindArr(Trunk_data[,,5:10],Trunk_data[,,12:13],Trunk_data[,,16:17],Trunk_data[,,25:31],Trunk_data[,,35:46],along=3)
  
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

proc_trunk <- procSym(Trunk_males)

relaxarray <- Trunk_males

for (i in 1:29) {
  print(i)
  
  relaxarray[,,i] <- relaxLM(Trunk_males[,,i],proc_trunk$mshape,SMvector=fix, deselect=TRUE, surp=Surfaces,outlines=Curves)
  
}

Todo_trunk <- relaxarray[,,1:29]

# Thorax and pelvis n=29 male H. sapiens
Thorax_temp <- read.morphologika("Thorax_temp_mlk_stern.txt")
Thorax_indices <- vcgKDtree(Trunk_data[,,36],Thorax_temp[,,1],k=1)$index
Thorax_males <- Trunk_males[Thorax_indices,,]
Pelvis_temp <- read.morphologika("Pelvis_temp_mlk.txt")
Pelvis_indices <- vcgKDtree(Trunk_data[,,36],Pelvis_temp[,,1],k=1)$index
Pelvis_males <- Trunk_males[Pelvis_indices,,]

Thorax_mesh <- vcgImport("Thorax_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Pelvis_mesh <- vcgImport("Pelvis_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)

#GPA N=29 males H. sapiens

Thorax_males_GPA <- procSym(Thorax_males)
Pelvis_males_GPA <- procSym(Pelvis_males)

#MODEL C: PLS n=29 male Homo sapiens ----
PLS_ModelC <- pls2B(Thorax_males_GPA$rotated,Pelvis_males_GPA$rotated, rounds=999,cv=T)
PLS_ModelC

col_plot <- c(rep("mediumseagreen",29))
col_leg <- c("mediumseagreen")

plot(PLS_ModelC$Xscores[,1],PLS_ModelC$Yscores[,1], xlab="Block 1 (Thorax)",ylab="Block 2 (Pelvis)",main="TWO BLOCKS PLS ANALYSIS (1st latent variable)",pch=18,col=col_plot,cex=2,ylim=c(-0.08,0.08),xlim=c(-0.10,0.10))
abline(h=0.00,v=0.00,col="black")
legend("topleft", inset=.005, legend = c("Males"), horiz=FALSE, cex=1.2,bty="n",fill=col_leg, border=col_leg)

Covar_ModelC <- plsCoVar(PLS_ModelC,i=1,sdx=1,sdy=1)

PLS1_ModelC_neg_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelC$x[,,1])
shade3d(PLS1_ModelC_neg_x, col="mediumseagreen",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_pos_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelC$x[,,2])
shade3d(PLS1_ModelC_pos_x, col="mediumseagreen",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_neg_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelC$y[,,1])
shade3d(PLS1_ModelC_neg_y,col="mediumseagreen",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_pos_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelC$y[,,2])
shade3d(PLS1_ModelC_pos_y,col="mediumseagreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#Error de predicción

truth_thorax <- vecx(Thorax_males_GPA$rotated) #Convert the 3D array into a matrix (Thorax)

error <- NULL
for (i in 1:dim(PLS_ModelC$predicted.x)[3]) # loop for cross-validation
  error <- c(error,mean((truth_thorax-PLS_ModelC$predicted.x[,,i])^2))

NoCVerror <- NULL
for (i in 1:dim(PLS_ModelC$predicted.x)[3]) # loop for MSE non-validated
  NoCVerror <- c(NoCVerror,mean((truth_thorax-vecx(predictPLSfromData(PLS_ModelC,y=Pelvis_males_GPA$rotated,ncomp=i)))^2))

# Cross-validation plot----

plot(1:27,error,main="Cross-validated prediction error",pch=19,xlab="Latent Variables used for Prediction",ylab="MSE",ylim=c(0,max(c(error,NoCVerror))))
lines(1:27,error)
points(1:27,NoCVerror,col="red",pch=20)                  
lines(1:27,NoCVerror,col="red")
legend(x="topright",y=NULL,legend=c("Cross-validated", "Non-validated"), col = c(1,"red"),pch=c(19,20),lty=1)

#PREDICTIONS based on MODEL C----
#Prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model C
Kebara2RA_toGPA_ModelC <- rotonto(Pelvis_males_GPA$mshape,Pelves_Kebara2[,,1],scale = T)$yrot
Pred_Kebara2RA_ModelC <- predictPLSfromData(PLS_ModelC,y=Kebara2RA_toGPA_ModelC,ncomp=3) #number of LVs to perform the prediction
deformGrid3d(Thorax_males_GPA$mshape,Pred_Kebara2RA_ModelC)
deformGrid3d(Pred_Kebara2RA_ModelC,Pred_Kebara2RA_ModelC)
Thorax_Kebara2RA_ModelC <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2RA_ModelC)
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model C with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelC <- vcgImport("Kebara2_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelC <- rotonto(Thorax_males_GPA$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelC <- tps3d(Thorax_Kebara2GO_mesh_ModelC,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelC)
shade3d(Thorax_Kebara2GO_ModelC,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model C with thorax reconstruction Sawyer and Maley (2005)
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh_ModelC <- vcgImport("Thorax_Sawyer_Maley.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelC <- rotonto(Thorax_males_GPA$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelC <- tps3d(Thorax_Kebara2S_M_mesh_ModelC,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelC)
shade3d(Thorax_Kebara2S_M_ModelC,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#Prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model C
Kebara2SM_toGPA_ModelC <- rotonto(Pelvis_males_GPA$mshape,Pelves_Kebara2[,,2],scale = T)$yrot
Pred_Kebara2SM_ModelC <- predictPLSfromData(PLS_ModelC,y=Kebara2SM_toGPA_ModelC,ncomp=3) #number of LVs to perform the prediction
deformGrid3d(Thorax_males_GPA$mshape,Pred_Kebara2SM_ModelC)
deformGrid3d(Pred_Kebara2SM_ModelC,Pred_Kebara2SM_ModelC)
Thorax_Kebara2SM_ModelC <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2SM_ModelC)
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model C with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelC <- vcgImport("Kebara2_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelC <- rotonto(Thorax_males_GPA$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelC <- tps3d(Thorax_Kebara2GO_mesh_ModelC,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelC)
shade3d(Thorax_Kebara2GO_ModelC,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model B with thorax reconstruction of Sawyer and Maley (2005)
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh_ModelC <- vcgImport("Thorax_Sawyer_Maley.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelC <- rotonto(Thorax_males_GPA$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelC <- tps3d(Thorax_Kebara2S_M_mesh_ModelC,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelC)
shade3d(Thorax_Kebara2S_M_ModelC,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#PLS scores of predictions Model C
col_polyg <- c("black")
col_plot <- c(rep("mediumseagreen",29))
col_leg <- c("mediumseagreen")
plot(PLS_ModelC$Xscores[,1],PLS_ModelC$Yscores[,1], xlab="Block 1 (Thorax)",ylab="Block 2 (Pelvis)",main="TWO BLOCKS PLS ANALYSIS (1st latent variable)",pch=16,col=col_plot,cex=2,ylim=c(-0.08,0.08),xlim=c(-0.10,0.10))
abline(h=0.00,v=0.00,col="black")
legend("bottomright", inset=.005, legend = c("Males"), horiz=FALSE, cex=1.2,bty="n",fill=col_leg, border=col_leg)

groups29thorax <- read.table("Groups_males29.txt",header=TRUE)
if (require(car)) {
  for(ii in 1:length(levels(groups29thorax$SPECIES))){
    dataEllipse(PLS_ModelC$Xscores[groups29thorax$SPECIES==levels(groups29thorax$SPECIES)[ii],1],PLS_ModelC$Yscores[groups29thorax$SPECIES==levels(groups29thorax$SPECIES)[ii],1],add=TRUE,center.pch = FALSE,levels=c(.95), col=col_polyg[ii],fill=F, fill.alpha=0.10,robust=TRUE)}
}

#Thorax prediction Model C from Rak & Arensburg
PLSXscore <- getPLSscores(PLS_ModelC,y=Kebara2RA_toGPA_ModelC)
PLSYscore <- getPLSscores(PLS_ModelC,x=Pred_Kebara2RA_ModelC)
points(PLSXscore[,1],PLSYscore[,1],pch=15,col="orange",cex=0.5)

#Thorax prediction Model C from Sawyer & Maley
PLSXscore <- getPLSscores(PLS_ModelC,y=Kebara2SM_toGPA_ModelC)
PLSYscore <- getPLSscores(PLS_ModelC,x=Pred_Kebara2SM_ModelC)
points(PLSXscore[,1],PLSYscore[,1],pch=15,col="purple",cex=0.5)


#All the thoraces

N72thorax <- bindArr(Thorax,Pred_Kebara2RA_ModelA,Pred_Kebara2RA_ModelB,Pred_Kebara2RA_ModelC,Pred_Kebara2SM_ModelA,Pred_Kebara2SM_ModelB,Pred_Kebara2SM_ModelC,Kebara2_thorax_reconstructions[,,1],Kebara2_thorax_reconstructions[,,2],along = 3)

#Define landmarks and semilandmarks----

Fix <- c(1:105)
Curve1 <- c(1,106:118,46)
Curve2 <- c(23,119:131,68)
Curve3 <- c(3,132:144,47)
Curve4 <- c(25,145:157,69)
Curve5 <- c(5,158:170,49)
Curve6 <- c(27,171:183,71)
Curve7 <- c(7,184:196,51)
Curve8 <- c(29,197:209,73)
Curve9 <- c(9,210:222,53)
Curve10 <- c(31,223:235,75)
Curve11 <- c(11,236:248,55)
Curve12 <- c(33,249:261,77)
Curve13 <- c(13,262:274,57)
Curve14 <- c(35,275:287,79)
Curve15 <- c(15,288:300,59)
Curve16 <- c(37,301:313,81)
Curve17 <- c(17,314:326,61)
Curve18 <- c(39,327:339,83)
Curve19 <- c(19,340:352,63)
Curve20 <- c(41,353:365,85)
Curve21 <- c(21,366:378,65)
Curve22 <- c(43,379:391,87)
Curve23 <- c(2,392:404,45)
Curve24 <- c(24,405:417,67)
Curve25 <- c(4,418:430,48)
Curve26 <- c(26,431:443,70)
Curve27 <- c(6,444:456,50)
Curve28 <- c(28,457:469,72)
Curve29 <- c(8,470:482,52)
Curve30 <- c(30,483:495,74)
Curve31 <- c(10,496:508,54)
Curve32 <- c(32,509:521,76)
Curve33 <- c(12,522:534,56)
Curve34 <- c(34,535:547,78)
Curve35 <- c(14,548:560,58)
Curve36 <- c(36,561:573,80)
Curve37 <- c(16,574:586,60)
Curve38 <- c(38,587:599,82)
Curve39 <- c(18,600:612,62)
Curve40 <- c(40,613:625,84)
Curve41 <- c(20,626:638,64)
Curve42 <- c(42,639:651,86)
Curve43 <- c(22,652:664,66)
Curve44 <- c(44,665:677,88)


Curves <- list(Curve1,Curve2,Curve3,Curve4,Curve5,Curve6,Curve7,Curve8,Curve9,Curve10,Curve11,Curve12,Curve13,Curve14,Curve15,Curve16,Curve17,Curve18,Curve19,Curve20,Curve21,Curve22,Curve23,Curve24,Curve25,Curve26,Curve27,Curve28,Curve29,Curve30,Curve31,Curve32,Curve33,Curve34,Curve35,Curve36,Curve37,Curve38,Curve39,Curve40,Curve41,Curve42,Curve43,Curve44)
proc_thorax <- procSym(N72thorax)


N72thorax_reslid <- N72thorax

for (i in 1:72) {
  print(i)
  
  N72thorax_reslid[,,i] <- relaxLM(N72thorax[,,i],proc_thorax$mshape,SMvector=Fix,deselect=TRUE,outlines=Curves)
  
}

proc_N72thorax_reslid <- procSym(N72thorax_reslid[,,1:72])
groups_72 <- read.table("Groups_Homo72.txt",header=T)
col_plot <- c(rep("indianred1",64),rep("gray",6),rep("deepskyblue",1),rep("green",1))
col_polyg <- c("red","transparent")

require(car)
plot(proc_N72thorax_reslid$PCscores[,1],proc_N72thorax_reslid$PCscores[,2],pch=16,col=col_plot,cex=2)
scatter3d(proc_N72thorax_reslid$PCscores[,1],proc_N72thorax_reslid$PCscores[,2],proc_N72thorax_reslid$PCscores[,3],surface=F,groups = groups_72$SPECIES, ellipsoid = T,level=0.95,surface.col = col_polyg,grid = FALSE,point.col = col_plot,xlab = "PC1", ylab = "PC2",
          zlab = "PC3")

r2morphologika(N72thorax_reslid,"N72thorax_reslid.txt")
