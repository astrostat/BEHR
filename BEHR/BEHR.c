#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"functions.h"

#define ABS(x)   ((x)>=0   ? (x) : -(x))
#define MAX(x,y) ((x)>=(y) ? (x) :  (y))
#define MIN(x,y) ((x)>=(y) ? (y) :  (x))
#define DISPLAY 0

/*************************************************/
/* BEHR (Bayesian Estimation of Hardness Ratios) */
/*                                               */
/*       written by Taeyoung Park                */
/*                  Dept of Statistics           */
/*                  Harvard University           */
/*                  tpark@stat.harvard.edu       */
/*                                               */
/*************************************************/

int numBins,MODE,CLASS,LIMIT,tblS=0,tblH=0,numS=0,numH=0;
double X1,X2,B1,B2,C1,C2,E1,E2,K1,K2,K3,norm1,norm2,*prior1,*prior2;
double **TblPRIOR_soft,**TblPRIOR_hard;

int main(int argc, char *argv[]){
  clock_t tm_std;
  double dur_std, dur_new;        /***********************************************/
  int Gibbs;                      /* 1: Gibbs sampler; 0: Gaussian quadrature    */
  int nooutput=1;                 /* 1: do not create an output file             */
  int posterior=0;                /* 1: display parameters of gamma posterior    */
  int details=0;                  /* 1: outputs of R, HR, & C are only generated */
  int HPD=1;                      /* 1: HPD interval for the Gibbs sampler       */
  int SWITCH=0;                   /* 1: Switch the sign of the color             */
  int fix=0;                      /* 1: Fix a problem of a quad function         */
  int i,j,k;                      /*                   Indexes                   */
  int dummy=0, dum=0, helps=0;    /***********************************************/
  int dum1,dum2,dum3,dum4,dum5;   /*                                             */
  int dum6,dum7,dum8,dum9,dum10;  /*                                             */
  int dum11,dum12,dum13,dum14;    /*               Dummy variables               */
  int dum15,dum16,dum17,dum18;    /*                                             */
  int dum19,dum20,dum21,dum22;    /*                                             */
  int dum23,dum24,dum25,dum26;    /*                                             */
  int dum27;                      /***********************************************/
  int flag1, flag2, flag3, flag4; /***********************************************/ 
  int flag5, flag6, flag7, flag8; /*              flag of warnings               */ 
  int flag9, flag10=0, flag11=0;  /***********************************************/ 
  int iter;                       /*              Gibbs iterations               */
  int draws;                      /*                 Gibbs draws                 */
  int burn_in;                    /*                Burn-in draws                */
  int ercode;                     /***********************************************/
  int minrul, maxrul, ruleno;     /*       Gaussian quadrature arguments         */
  int nevals, nevals_pr=0;        /***********************************************/
  int numBins_temp;               /*     The number of bins for a parameter      */
  int rejS, rejH;                 /*    Rejection rate of the M-H algorithm      */
  float r;                        /* Acceptance probability of the M-H algorithm */
  float temp1,temp2,temp3,temp4;  /***********************************************/
  float temp5,temp6,temp7,temp8;  /*            Temporary constants              */
  float tempmeanS=0,tempmeanH=0;  /***********************************************/
  float LEVEL;                    /*            Level of error bars              */
  void  pquad_();                 /*        Gaussian quadrature function         */
  double post_();                 /*         Posterior density functions         */
  double absacc, relacc;          /***********************************************/
  double esterr;                  /*       Gaussian quadrature arguments         */
  double quads[9];                /***********************************************/
  double dens_sum;                /*              Density summation              */
  double degree;                  /*         Degree of lowering a density        */
  double soft_src, soft_bkg;      /*    Soft source & Soft background counts     */
  double hard_src, hard_bkg;      /*    Hard source & Hard background counts     */
  double ratio_soft, ratio_hard;  /*      Soft area ratio & Hard area ratio      */
  double soft_eff, hard_eff;      /*  Soft effective area & Hard effective area  */
  double soft_idx, hard_idx;      /*           Soft index & Hard index           */
  double soft_scale, hard_scale;  /*           Soft scale & Hard scale           */
  double lamS_left, lamS_right;   /*      The lower and upper bounds of lamS     */
  double lamH_left, lamH_right;   /*      The lower and upper bounds of lamH     */
  double zeta_left, zeta_right;   /*      The lower and upper bounds of zeta     */
  double R_left, R_right;         /*       The lower and upper bounds of R       */
  double C_left, C_right;         /*       The lower and upper bounds of C       */
  double omega_left, omega_right; /*     The lower and upper bounds of omega     */
  double phi_left, phi_right;     /*      The lower and upper bounds of phi      */
  double HR_left, HR_right;       /*      The lower and upper bounds of HR       */
  double beta_S, beta_H;          /*  Background counts around the source area   */
  double lam_S, lam_H;            /*         Source Poisson intensities          */
  double xi_S, xi_H;              /*       Background Poisson intensities        */
  double *PRIOR_soft;             /* Prior distributions for Poisson intensities */
  double *PRIOR_hard;             /* Prior distributions for Poisson intensities */
  double prob, *param;            /*    Success probability; Gamma parameters    */
  double *tempR, *tempHplusS;     /***********************************************/
  double *tempHminusS;            /*                                             */
  double *tempHR, *tempLogS;      /*         Temporary Monte Carlo draws         */
  double *tempLogH, *tempC;       /*                                             */
  double *tempS, *tempH;          /***********************************************/
  double X1_t,X2_t,B1_t,B2_t;     /*                                             */         
  double C1_t,C2_t,E1_t,E2_t;     /*          Temporary storage of data          */
  double norm1_t,norm2_t;         /*                                             */
  double *prior1_t, *prior2_t;    /*                                             */
  double *Q1_HR;                  /***********************************************/
  double *Q2_HR;                  /*                                             */
  double *Q3_HR;                  /*                                             */
  double *Q4_HR;                  /*                                             */
  double *Q5_HR;                  /*     Summary quantiles of Hardness ratios    */
  double *Q6_HR;                  /*                                             */
  double *Q7_HR;                  /*                                             */
  double *Q8_HR;                  /*                                             */
  double *Q9_HR;                  /***********************************************/
  double *bound_S;                /***********************************************/
  double *bound_H;                /*           Bounds for parameters             */
  double *bound_R;                /***********************************************/
  double *R,*HR,*C,*lamS,*lamH;   /*              Hardness ratios                */  
  double *LogS,*LogH,*omega,*phi; /*          Parameters of interest             */  
  double *dens_R, *dens_HR;       /***********************************************/
  double *dens_C, *dens_lamS;     /*                                             */
  double *dens_lamH, *dens_omega; /*    Densities of parameters of interest      */
  double *dens_phi;               /*                                             */
  double *dens_LogS, *dens_LogH;  /***********************************************/
   
  char key[]="=";
  char *pch;
  char softsrc[9]="softsrc=";
  char softbkg[9]="softbkg=";
  char softeff[9]="softeff=";
  char softarea[10]="softarea="; 
  char hardsrc[9]="hardsrc=";
  char hardbkg[9]="hardbkg=";
  char hardeff[9]="hardeff=";
  char hardarea[10]="hardarea="; 
  char algo[6]="algo=";
  char GibbsSampler[6]="gibbs"; 
  char Quadrature[5]="quad";
  char HPDint[5]="HPD=";
  char level[7]="level=";
  char softidx[9]="softidx=";
  char hardidx[9]="hardidx=";
  char softscl[9]="softscl=";
  char hardscl[9]="hardscl=";
  char softtbl[9]="softtbl=";
  char hardtbl[9]="hardtbl=";
  char nsim[6]="nsim="; 
  char post[6]="post="; 
  char nburnin[9]="nburnin=";
  char nbins[7]="nbins=";
  char out[8]="output=";
  char outR[9]="outputR=";
  char outHR[10]="outputHR=";
  char outC[9]="outputC=";
  char outPr[10]="outputPr=";
  char outMC[10]="outputMC=";
  char detail[9]="details=";
  char none[5]="none";
  char true[5]="true";
  char false[6]="false";
  char help[5]="help";
  char tempfile[200],tempfile1[200],tempfile2[200];  
  char infile1[200],infile2[200];
  char outfile1[200],outfile2[200],outfile3[200];
  char outfile4[200],outfile5[200],outfile6[200];

  FILE *input1,*input2;
  FILE *output1,*output2,*output3,*output4,*output5,*output6;
  
  PRIOR_soft = (double *) calloc(4, sizeof(double));
  PRIOR_hard = (double *) calloc(4, sizeof(double));

  dum1 =dum2 =dum3 =dum4 =dum5 =dum6 =dum7 =dum8 =dum9 =dum10=1;
  dum11=dum12=dum13=dum14=dum15=dum16=dum17=dum18=dum19=dum20=1;
  dum26=dum27=1;
  for(i=1;i<argc;i++){
    if(sameString(help,argv[i],4)) helps=1;
    if(sameString(softsrc,argv[i],8)){ 
      pch = strstr(argv[i], key);
      soft_src = atof(pch+1);
      if(soft_src<0){
	printf("\n#****************************************************************#");
	printf("\n# ERROR: The soft source counts should be a nonnegative integer. #");
	printf("\n#****************************************************************#\n");
	exit(1);}
      if(soft_src!=(int)soft_src){
	printf("\n#*****************************************************#");
	printf("\n# ERROR: The soft source counts should be an integer. #");
	printf("\n#*****************************************************#\n");
	exit(1);}
      dum1=0;}
    else if(sameString(hardsrc,argv[i],8)){ 
      pch = strstr(argv[i], key);
      hard_src = atof(pch+1);
      if(hard_src<0){
	printf("\n#****************************************************************#");
	printf("\n# ERROR: The hard source counts should be a nonnegative integer. #");
	printf("\n#****************************************************************#\n");
	exit(1);}
      if(hard_src!=(int)hard_src){
	printf("\n#*****************************************************#");
	printf("\n# ERROR: The hard source counts should be an integer. #");
	printf("\n#*****************************************************#\n");
	exit(1);}
      dum2=0;}
    else if(sameString(softbkg,argv[i],8)){ 
      pch = strstr(argv[i], key);
      soft_bkg = atof(pch+1);
      if(soft_bkg<0){
	printf("\n#********************************************************************#");
	printf("\n# ERROR: The soft background counts should be a nonnegative integer. #");
	printf("\n#********************************************************************#\n");
	exit(1);}
      if(soft_bkg!=(int)soft_bkg){
	printf("\n#*********************************************************#");
	printf("\n# ERROR: The soft background counts should be an integer. #");
	printf("\n#*********************************************************#\n");
	exit(1);}
      dum3=0;} 
    else if(sameString(hardbkg,argv[i],8)){ 
      pch = strstr(argv[i], key);
      hard_bkg = atof(pch+1);
      if(hard_bkg<0){
	printf("\n#********************************************************************#");
	printf("\n# ERROR: The hard background counts should be a nonnegative integer. #");
	printf("\n#********************************************************************#\n");
	exit(1);}
      if(hard_bkg!=(int)hard_bkg){
	printf("\n#*********************************************************#");
	printf("\n# ERROR: The hard background counts should be an integer. #");
	printf("\n#*********************************************************#\n");
	exit(1);}
      dum4=0;} 
    else if(sameString(softarea,argv[i],9)){ 
      pch = strstr(argv[i], key);
      ratio_soft = atof(pch+1);
      if(ratio_soft<=0){
	printf("\n#************************************************#");
	printf("\n# ERROR: The soft area ratio should be positive. #");
	printf("\n#************************************************#\n\n");
	exit(1);}
      dum5=0;} 
    else if(sameString(hardarea,argv[i],9)){ 
      pch = strstr(argv[i], key);
      ratio_hard = atof(pch+1);
      if(ratio_hard<=0){
	printf("\n#************************************************#");
	printf("\n# ERROR: The hard area ratio should be positive. #");
	printf("\n#************************************************#\n\n");
	exit(1);}
      dum6=0;}
    else if(sameString(softeff,argv[i],8)){ 
      pch = strstr(argv[i], key);
      soft_eff = atof(pch+1);
      if(soft_eff<=0){
	printf("\n#****************************************************#");
	printf("\n# ERROR: The soft effective area should be positive. #");
	printf("\n#****************************************************#\n\n");
	exit(1);}
      dum7=0;} 
    else if(sameString(hardeff,argv[i],8)){ 
      pch = strstr(argv[i], key);
      hard_eff = atof(pch+1);
      if(hard_eff<=0){
	printf("\n#****************************************************#");
	printf("\n# ERROR: The hard effective area should be positive. #");
	printf("\n#****************************************************#\n\n");
	exit(1);}
      dum8=0;} 
    else if(sameString(softidx,argv[i],8)){ 
      pch = strstr(argv[i], key);
      soft_idx = atof(pch+1);
      if(soft_idx<=0){
	 printf("\n#*******************************************#");
	 printf("\n# ERROR: The soft index should be positive. #");
	 printf("\n#*******************************************#\n");
	 exit(1);}
      dum9=0;} 
    else if(sameString(hardidx,argv[i],8)){ 
      pch = strstr(argv[i], key);
      hard_idx = atof(pch+1);
      if(hard_idx<=0){
	 printf("\n#*******************************************#");
	 printf("\n# ERROR: The hard index should be positive. #");
	 printf("\n#*******************************************#\n");
	 exit(1);}
      dum10=0;} 
    else if(sameString(level,argv[i],6)){ 
      pch = strstr(argv[i], key);
      LEVEL = atof(pch+1);
      if(LEVEL>=100 | LEVEL<=0){
	printf("\n#****************************************************************#");
	printf("\n# ERROR: The confidence level must be between 0 and 100 percent. #");
	printf("\n#****************************************************************#\n");
        exit(1);}
      dum11=0;} 
    else if(sameString(algo,argv[i],5)){ 
      pch = strstr(argv[i], key);
      if(sameString(GibbsSampler,pch+1,5))    Gibbs=1;
      else if(sameString(Quadrature,pch+1,4)) Gibbs=0;
      else{
	printf("\n#**************************************************************#");
	printf("\n# ERROR: Correctly enter which algorithm to use (i.e., algo= ).#");
	printf("\n#**************************************************************#\n\n");
	exit(1);}
      dum12=0;}    
    else if(sameString(nsim,argv[i],5)){ 
      pch = strstr(argv[i], key);
      draws = atof(pch+1);
      dum13=0;} 
    else if(sameString(nburnin,argv[i],8)){ 
      pch = strstr(argv[i], key);
      burn_in = atof(pch+1);
      dum14=0;} 
    else if(sameString(nbins,argv[i],6)){ 
      pch = strstr(argv[i], key);
      numBins_temp = atof(pch+1);
      numBins = (numBins_temp%2)==0 ? (numBins_temp+1) : (numBins_temp);
      dum15=0;}  
    else if(sameString(out,argv[i],7)){ 
      pch = strstr(argv[i], key);
      if(sameString(none,pch+1,4)) nooutput=1; /* default */
      else{
	strcpy(tempfile,pch+1);
	strcpy(outfile1,pch+1);
	strcat(outfile1,".txt");
	nooutput=0;}
      dum16=0;} 
    else if(sameString(softscl,argv[i],8)){ 
      pch = strstr(argv[i], key);
      soft_scale = atof(pch+1);
      if(soft_scale<0){
	printf("\n#***********************************************#");
	printf("\n# ERROR: The soft scale should be non-negative. #");
	printf("\n#***********************************************#\n");
	exit(1);}
      dum17=0;} 
    else if(sameString(hardscl,argv[i],8)){ 
      pch = strstr(argv[i], key);
      hard_scale = atof(pch+1);
      if(hard_scale<0){
	printf("\n#***********************************************#");
	printf("\n# ERROR: The hard scale should be non-negative. #");
	printf("\n#***********************************************#\n");
	exit(1);}
      dum18=0;} 
    else if(sameString(post,argv[i],5)){ 
      pch = strstr(argv[i], key); 
      if(sameString(true,pch+1,4)) posterior=1;
      dum19=0;}
    else if(sameString(HPDint,argv[i],4)){ 
      pch = strstr(argv[i], key); 
      if(sameString(true,pch+1,4)) HPD=1;
      else HPD=0;
      dum20=0;}
    else if(sameString(softtbl,argv[i],8)){ 
      pch = strstr(argv[i], key); 
      strcpy(tempfile1,pch+1);
      strcpy(infile1,pch+1);
      tblS=1;
      dum26=0;}
    else if(sameString(hardtbl,argv[i],8)){ 
      pch = strstr(argv[i], key); 
      strcpy(tempfile2,pch+1);
      strcpy(infile2,pch+1);
      tblH=1;
      dum27=0;}
  }

  if((tblS | tblH) && !Gibbs){
    printf("\n#****************************************************#");
    printf("\n# WARNING: The tabulated prior is not yet supported  #");
    printf("\n#          in the Gaussian quadrature mode. Instead, #");
    printf("\n#          use the Gibbs mode or the gamma prior.    #");
    printf("\n#          Now changing to the default value...      #");
    printf("\n#****************************************************#\n");
    tblS=tblH=0;
    dum26=dum27=1;
  }

  if(!dum26){
    if(!(input1=fopen(infile1,"r"))){
      printf("\nCould not open the %s file.\n\n", infile1);
      fflush(stdout);
      exit(1);
    }
    fgets(infile1,200,input1);
    sscanf(infile1,"%d",&numS);

    TblPRIOR_soft = matrix(0,numS-1,0,1);
    temp8 = 0.0;

    fgets(infile1,200,input1); i=0;
    while(fgets(infile1,200,input1)){
      if(i>=numS){
	printf("\n#******************************************************#");
	printf("\n# ERROR: The number of rows in the tabulated prior     #\n");
	printf("\n#        for soft band is larger than its input value. #\n");
	printf("\n#******************************************************#\n");
	exit(1);
      }
      sscanf(infile1,"%f %f",&temp1,&temp2);
      tempmeanS += temp1 * temp2;
      temp8     += temp2;
      TblPRIOR_soft[i][0] = temp1;  
      TblPRIOR_soft[i][1] = temp2;  
      i++;
    }
    tempmeanS /= temp8;
    fclose(input1);
  }

  if(!dum27){
    if(!(input2=fopen(infile2,"r"))){
      printf("\nCould not open the %s file.\n\n", infile2);
      fflush(stdout);
      exit(1);
    }
    fgets(infile2,200,input2);
    sscanf(infile2,"%d",&numH);

    TblPRIOR_hard = matrix(0,numH-1,0,1);
    temp8 = 0.0;

    fgets(infile2,200,input2); i=0;
    while(fgets(infile2,200,input2)){
      if(i>=numH){
	printf("\n#******************************************************#");
	printf("\n# ERROR: The number of rows in the tabulated prior     #\n");
	printf("\n#        for hard band is larger than its input value. #\n");
	printf("\n#******************************************************#\n");
	exit(1);
      }
      sscanf(infile2,"%f %f",&temp1,&temp2);
      tempmeanH += temp1 * temp2;
      temp8     += temp2;
      TblPRIOR_hard[i][0] = temp1;  
      TblPRIOR_hard[i][1] = temp2;  
      i++;
    }
    tempmeanH /= temp8;
    fclose(input2);
  }

  dummy  = dum1+dum2+dum3+dum4+dum5+dum6+dum7+dum8+dum9+dum10;
  dummy += dum11+dum12+dum13+dum14+dum15+dum16+dum17+dum18+dum19+dum20;
  dummy += dum26+dum27;
  dum    = dummy;
  if(argc==1){
    printf("\n#=============================================================================#");
    printf("\n#        BEHR (Bayesian Estimation of Hardness Ratios) ver. 12-12-2013       #");
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# Park, T., Kashyap, V., Siemiginowska, A., van Dyk, D., Zezas, A.,           #");
    printf("\n#     Heinke, C., and Wargelin, B (2006). Bayesian Estimation of Hardness     #");
    printf("\n#     Ratios: Modeling and Computations. The Astrophysical Journal, 652, 610. #");
    printf("\n#=============================================================================#");
    printf("\n# Usage: BEHR softsrc=S hardsrc=H softbkg=S_background hardbkg=H_background   #"); 
    printf("\n#        softarea=S_area hardarea=H_area softidx=S_prior hardidx=H_prior      #"); 
    printf("\n#        softscl=S_prior2 hardscl=H_prior2 softtbl=filename hardtbl=filename  #");
    printf("\n#        softeff=S_eff hardeff=H_eff level=confidence_level                   #");
    printf("\n#        algo=calculation_method nbins=num_of_bins nsim=num_of_draws          #");
    printf("\n#        nburnin=num_of_burnin_draws post=(true/false)?                       #");
    printf("\n#        post=(true/false)? HPD =(true/false)? details=(true/false)?          #");
    printf("\n#        output=filename_root outputR=(true/false)? outputHR=(true/false)?    #");
    printf("\n#        outputC=(true/false)? outputMC=(true/fale)? outputPr=(true/false)?   #");
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# BEHR computes R=S/H, C=log(S)-log(H), and HR=(H-S)/(H+S) in Poisson limit.  #");
    printf("\n#=============================================================================#\n\n");
    exit(1);
  }
  if(helps | (dum1+dum2+dum3+dum4+dum5*dum6)>0){
    printf("\n#=============================================================================#");
    printf("\n#        BEHR (Bayesian Estimation of Hardness Ratios) ver. 12-12-2013        #");
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# Park, T., Kashyap, V., Siemiginowska, A., van Dyk, D., Zezas, A.,           #");
    printf("\n#     Heinke, C., and Wargelin, B (2006). Bayesian Estimation of Hardness     #");
    printf("\n#     Ratios: Modeling and Computations. The Astrophysical Journal, 652, 610. #");
    printf("\n#=============================================================================#");
    printf("\n# computes simple ratio, R = S/H                                              #");
    printf("\n#          color, C = log10(S)-log10(H)                                       #");
    printf("\n#          fractional difference, HR = (H-S)/(H+S)                            #"); 
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# Required inputs:                                                            #");
    if(dum1)
    printf("\n#    softsrc  = source counts in soft (S) band                                #");
    if(dum2) 
    printf("\n#    hardsrc  = source counts in hard (H) band                                #");
    if(dum3) 
    printf("\n#    softbkg  = background counts in soft (S) band                            #");
    if(dum4) 
    printf("\n#    hardbkg  = background counts in hard (H) band                            #");
    if(dum5|dum6){ 
    printf("\n#    softarea = ratio of background area to source area in soft (S) band      #");
    printf("\n#               That is, softarea = (background softarea)/(source softarea)   #");
    printf("\n#    hardarea = ratio of background area to source area in hard (H) band      #");
    printf("\n#               That is, hardarea = (background hardarea)/(source hardarea)   #");
    printf("\n#             + If only one of softarea and hardarea is specified,            #");
    printf("\n#               then the same value is assumed for the other.                 #");}
    if(dum7|dum8+dum9|dum10+dum11+dum12+dum13+dum14+dum15+dum16>0){
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# Optional inputs:                                                            #");
    if(dum7|dum8){
    printf("\n#    softeff  = effective area in soft (S) band   (hardeff or 1.0 by default) #");
    printf("\n#    hardeff  = effective area in hard (H) band   (softeff or 1.0 by default) #");}
    if(dum9|dum10){
    printf("\n#    softidx  = index of prior on S     (positive; hardidx or 0.5 by default) #");
    printf("\n#    hardidx  = index of prior on H     (positive; softidx or 0.5 by default) #");}
    if(dum17|dum18){
    printf("\n#    softscl  = scale of prior on S (non-negative; hardscl or 0.0 by default) #");
    printf("\n#    hardscl  = scale of prior on H (non-negative; softscl or 0.0 by default) #");}
    if(dum26|dum27){
    printf("\n#    softtbl  = tabulated prior for soft band             (unused by default) #");
    printf("\n#    hardtbl  = tabulated prior for hard band             (unused by default) #");}
    if(dum11)
    printf("\n#    level    = confidence level at which to report error   (68.0 by default) #");
    printf("\n#    details  = compute various ratios (true/false)?       (false by default) #");
    if(dum12){ 
    printf("\n#    algo     = calculation method, either the Gibbs sampler (gibbs)          #");
    printf("\n#               or efficient numerical integration (quad)  (gibbs by default) #");}
    if(dum13 && (Gibbs | dum12)){ 
    printf("\n#    nsim     = number of draws if algo=gibbs              (10000 by default) #");
    printf("\n#    nburnin  = number of burn-in draws if algo=gibbs       (5000 by default) #");}
    if(Gibbs && dum20){
    printf("\n#    HPD      = computes a HPD interval when algo=gibbs (true/false)?         #");
    printf("\n#               if false, computes an equal-tail interval   (true by default) #");}
    if(dum19)
    printf("\n#    post     = if true, displays a gamma posterior distributions             #");
    if(dum15 && !Gibbs)
    printf("\n#    nbins    = number of bins in integration if algo=quad   (500 by default) #");
    if(dum16){ 
    printf("\n#-----------------------------------------------------------------------------#");
    printf("\n# Output files:                                                               #");
    printf("\n#    output   = root of filename in which to place output   (BEHR by default) #");
    printf("\n#               output will be placed in the file, output.txt                 #");
    printf("\n#    outputR  = if true, writes output for R in output_R.txt                  #");
    printf("\n#    outputC  = if true, writes output for C in output_C.txt                  #");
    printf("\n#    outputHR = if true, writes output for HR in output_HR.txt                #");
    printf("\n#    outputMC = if true, writes Monte Carlo draws for lamS and lamH           #");
    printf("\n#               to output_draws.txt when algo=gibbs                           #");
    printf("\n#    outputPr = if true, writes the probability distributions for             #");
    printf("\n#               R, HR, C, lamS, and lamH to output_prob.txt when algo=quad    #");}}
    printf("\n#=============================================================================#\n\n");
    exit(1);
  }
  if( dum5  && !dum6)  ratio_soft = ratio_hard;
  if(!dum5  &&  dum6)  ratio_hard = ratio_soft;
  if( dum7  && !dum8)  soft_eff   = hard_eff;
  if(!dum7  &&  dum8)  hard_eff   = soft_eff;
  if( dum7  &&  dum8)  soft_eff   = hard_eff = 1.0;
  if( dum9  && !dum10) soft_idx   = hard_idx;
  if(!dum9  &&  dum10) hard_idx   = soft_idx;
  if( dum9  &&  dum10) soft_idx   = hard_idx = 0.5;
  if( dum17 && !dum18) soft_scale = hard_scale;
  if(!dum17 &&  dum18) hard_scale = soft_scale;
  if( dum17 &&  dum18) soft_scale = hard_scale = 0.0;
  if( dum11) LEVEL=68.0;
  if(!dum11 && LEVEL<1.0){
    printf("\n#*********************************************#");
    printf("\n# BEHR converts level = %4g into %2g percent. *",LEVEL,LEVEL*100);  
    printf("\n#*********************************************#\n");
    LEVEL *= 100;
  }
  if(dum12) Gibbs=1;

  PRIOR_soft[0] = (double)soft_idx;      /*     Psi_S1      */
  PRIOR_hard[0] = (double)hard_idx;      /*     Psi_H1      */
  PRIOR_soft[1] = (double)soft_scale;    /*     Psi_S2      */
  PRIOR_hard[1] = (double)hard_scale;    /*     Psi_H2      */
  PRIOR_soft[2] = PRIOR_hard[2]=0.5;     /* Psi_S3 & Psi_H3 */
  PRIOR_soft[3] = PRIOR_hard[3]=0.0;     /* Psi_S4 & Psi_H4 */

  if(tblS | tblH){
    if(!dum9 | !dum10 | !dum17 | !dum18){ 
      printf("\n#***********************************************************#");
      printf("\n# NOTE: Tabulated priors for the soft band and/or hard band #");
      printf("\n#       will be used to compute posteriors.                 #");
      printf("\n#***********************************************************#\n");
    }
  }
  
  if(Gibbs && !dum13 && !dum14 && draws<burn_in){
    printf("\n#***************************************************#");
    printf("\n# ERROR: number of draws > number of burn-in draws. #");
    printf("\n#***************************************************#\n\n");
    exit(1);
  }
  if( Gibbs && dum13) draws = (dum14) ? 10000 : MAX((int)(burn_in*2),10000);
  if( Gibbs && dum14) burn_in = MIN((int)(draws/2),5000);
  if(!Gibbs && dum15) numBins = 501;
  if(!Gibbs && numBins<100){
    printf("\n#***************************************************************#");
    printf("\n# WARNING: Increase the number of bins to get better estimates. #");
    printf("\n#***************************************************************#\n");
  }
  if(dum16) strcpy(outfile1,"BEHR.txt");
  
  dum+=6; dum21=dum22=dum23=dum24=dum25=1;
  for(i=1;i<argc;i++){
    if(sameString(outR,argv[i],8)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)){
	if(nooutput) strcpy(outfile2,"BEHR");
	else         strcpy(outfile2,tempfile);
	strcat(outfile2,"_R.txt");
	dum21=0;}}
    else if(sameString(outHR,argv[i],9)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)){
	if(nooutput) strcpy(outfile3,"BEHR");
	else         strcpy(outfile3,tempfile);
	strcat(outfile3,"_HR.txt");
	dum22=0;}}
    else if(sameString(outC,argv[i],8)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)){
	if(nooutput) strcpy(outfile4,"BEHR");
	else         strcpy(outfile4,tempfile);
	strcat(outfile4,"_C.txt");
	dum23=0;}}
    else if(sameString(outMC,argv[i],9)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)){
	if(nooutput) strcpy(outfile5,"BEHR");
	else         strcpy(outfile5,tempfile);
	strcat(outfile5,"_draws.txt");
	dum24=0;}}
    else if(sameString(outPr,argv[i],9)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)){
	if(nooutput) strcpy(outfile6,"BEHR");
	else         strcpy(outfile6,tempfile);
	strcat(outfile6,"_prob.txt");
	dum25=0;}}
    else if(sameString(detail,argv[i],8)){ 
      pch = strstr(argv[i], key); dum--;
      if(sameString(true,pch+1,4)) details=1;
    }
  }

  /* argc-1        : net number of arguments    */
  /* dum-3+2*Gibbs : number of unused arguments */
  /* 24+2*Gibbs+2  : total number of arguments  */

  if((argc-1) > (28-dum)){
    printf("\n#********************************************#");
    if(argc==(28-dum))
    printf("\n# ERROR: There exists an illegal argument.   #");
    if(argc>(28-dum))
    printf("\n# ERROR: There exist illegal arguments.      #"); 
    printf("\n#        Please check out your inputs.       #");
    printf("\n#********************************************#\n");
    exit(1);
  } 
  if(!nooutput)
    if(!(output1=fopen(outfile1,"w"))){
      printf("\n Could not open the %s file.\n", outfile1);
      exit(1);}
  if(!dum21)
    if(!(output2=fopen(outfile2,"w"))){
      printf("\n Could not open the %s file.\n", outfile2);
      exit(1);}
  if(!dum22)
    if(!(output3=fopen(outfile3,"w"))){
      printf("\n Could not open the %s file.\n", outfile3);
      exit(1);}
  if(!dum23)
    if(!(output4=fopen(outfile4,"w"))){
      printf("\n Could not open the %s file.\n", outfile4);
      exit(1);}
  if(!dum24 && Gibbs)
    if(!(output5=fopen(outfile5,"w"))){
      printf("\n Could not open the %s file.\n", outfile5);
      exit(1);}
  if(!dum25 && !Gibbs)
    if(!(output6=fopen(outfile6,"w"))){
      printf("\n Could not open the %s file.\n", outfile6);
      exit(1);}
  if(!dum24 && Gibbs && DISPLAY) fprintf(output5,"lamS\t lamH\n"); 

            printf("\n#==========================================#");
            printf("\n#           BEHR ver. 12-12-2013           #");
            printf("\n#==========================================#");
  if(Gibbs) printf("\n#      The Gibbs sampler is running.       #");  
  else      printf("\n#    The Gaussian quadrature is running.   #");  
            printf("\n#------------------------------------------#");
            printf("\n# softsrc  = %8g  hardsrc  = %8g #",soft_src,hard_src);
            printf("\n# softbkg  = %8g  hardbkg  = %8g #",soft_bkg,hard_bkg);
	    printf("\n# softarea = %8g  hardarea = %8g #",ratio_soft,ratio_hard);
            printf("\n#------------------------------------------#");
	    printf("\n# softeff  = %8g  hardeff  = %8g #",soft_eff,hard_eff);
  if(!tblS) printf("\n# softidx  = %8g  softscl  = %8g #",soft_idx,soft_scale);
  else      printf("\n# softtbl  = %20s          #",tempfile1);
  if(!tblH) printf("\n# hardidx  = %8g  hardscl  = %8g #",hard_idx,hard_scale);
  else      printf("\n# hardtbl  = %20s          #",tempfile2);
            printf("\n# level    = %8g  post     = %8s #",LEVEL,(!dum19?true:false));
  if(Gibbs) printf("\n# nsim     = %8d  nburnin  = %8d #",draws,burn_in);   
  else      printf("\n# nbins    = %8d                      #",numBins-1); 
            printf("\n#------------------------------------------#");
  if(!dum24 && Gibbs)
            printf("\n#     Monte Carlo draws will be saved.     #"); 
  if(!dum25 && !Gibbs) 
            printf("\n#   Probability densities will be saved.   #"); 
  if(HPD)   printf("\n#     A HPD interval will be computed.     #");
  else      printf("\n# An equal-tail interval will be computed. #");
  if(!nooutput)
            printf("\n# The outputs are placed in %13s. #",outfile1);
  if(!dum21)printf("\n# The output file for R will be generated. #"); 
  if(!dum22)printf("\n# The output file for C will be generated. #"); 
  if(!dum23)printf("\n# The output file for HR will be generated.#"); 
            printf("\n#==========================================#\n"); 
  
  temp1 = soft_src-soft_bkg/ratio_soft;
  temp2 = hard_src-hard_bkg/ratio_hard;
  flag1 = flag2 = flag3 = flag4 = flag5 = flag6 = flag7 = flag8 = 0;
  
  if(!Gibbs){ /* Gaussian quadrature */
    if((temp1>100 && PRIOR_soft[0]>=1.0) | (temp2>100 && PRIOR_hard[0]>=1.0)){
      printf("\n#********************************************************#");
      printf("\n# WARNING: The Gaussian quadrature is inefficient.       #");
      printf("\n#          Recommend using the Gibbs option, algo=gibbs. #");
      printf("\n#********************************************************#\n");
    }
  }
  
  X1 = soft_src;
  X2 = hard_src;
  B1 = soft_bkg;
  B2 = hard_bkg;
  C1 = ratio_soft;
  C2 = ratio_hard; 
  E1 = soft_eff;
  E2 = hard_eff; 
  prior1 = PRIOR_soft;
  prior2 = PRIOR_hard;
  if(posterior) param = (double *)calloc(4, sizeof(double));

  tm_std = clock(); 
  if(Gibbs){  /* begin Gibbs */
    /************************************/
    /*********   Gibbs Sampler   ********/
    /************************************/
    if(!tblS) lam_S = 1.0; /* Starting Value of lam_S */
    else      lam_S = tempmeanS;
    if(!tblH) lam_H = 1.0; /* Starting Value of lam_H */
    else      lam_H = tempmeanH;

    xi_S  = 1.0; /* Starting Value of xi_S  */
    xi_H  = 1.0; /* Starting Value of xi_H  */

    tempR         = (double *)calloc((int)draws-burn_in, sizeof(double)); 
    tempHR        = (double *)calloc((int)draws-burn_in, sizeof(double)); 
    tempC         = (double *)calloc((int)draws-burn_in, sizeof(double)); 
    if(!dum21 | details | posterior){
      tempS       = (double *)calloc((int)draws-burn_in, sizeof(double)); 
      tempH       = (double *)calloc((int)draws-burn_in, sizeof(double));}
    if(!dum22 | details){
      tempHplusS  = (double *)calloc((int)draws-burn_in, sizeof(double)); 
      tempHminusS = (double *)calloc((int)draws-burn_in, sizeof(double));}
    if(!dum23 | details){
      tempLogS    = (double *)calloc((int)draws-burn_in, sizeof(double)); 
      tempLogH    = (double *)calloc((int)draws-burn_in, sizeof(double));}

    for(flag9=0,k=0,iter=1;iter<=draws;iter++){	
      if(tblS) temp1 = lam_S;
      if(tblH) temp2 = lam_H;
  
      /********************************************/
      /*********   DRAW THE MISSING DATA   ********/
      /********************************************/
      beta_S = (double)ignbin(X1,xi_S/(lam_S+xi_S));
      beta_H = (double)ignbin(X2,xi_H/(lam_H+xi_H));
      
      /*********************************************/
      /************   DRAW PARAMETERS   ************/
      /*********************************************/
      lam_S  = gengam(E1+prior1[1],X1-beta_S+prior1[0]);
      lam_H  = gengam(E2+prior2[1],X2-beta_H+prior2[0]); 
      xi_S   = gengam(E1+C1*E1+prior1[3],B1+beta_S+prior1[2]);
      xi_H   = gengam(E2+C2*E2+prior2[3],B2+beta_H+prior2[2]); 
    
      if(tblS){
	lam_S = gengam(E1,X1-beta_S+1.0);
	rejS  = 0;
	/* M-H algorithm for lamS */ 
	while((TblPRIOR_soft[0][0] > lam_S | lam_S > TblPRIOR_soft[numS-1][0]) && rejS < 100){
	  lam_S  = gengam(E1,X1-beta_S+1.0);
	  rejS++;
	}
	if(rejS==100 | flag10>=iter){
	  tblS = 0;
	  printf("\n#******************************************************************#");
	  printf("\n# WARNING: THE TABULATED PRIOR FOR THE 'SOFT' BAND CONTRADICTS THE #");
	  printf("\n#          LIKELIHOOD OF THE DATA, I.E., THE MLE IS OUTSIDE THE    #");
	  printf("\n#          RANGE OF THE TABULATED PRIOR. NOW BEHR IS FINISHING THE #");
	  printf("\n#          CALCULATION WITH ITS DEFAULT GAMMA PRIOR...             #");
	  printf("\n#******************************************************************#\n");
	}
        else{
	  r = findprob(lam_S,numS,TblPRIOR_soft) - findprob(temp1,numS,TblPRIOR_soft);
	  if(ranf()>exp(r)){ lam_S = temp1; flag10++; }/* reject the draw */ 
	}
      }
      if(tblH){
	lam_H = gengam(E2,X2-beta_H+1.0);
        rejH  = 0;
	/* M-H algorithm for lamH */ 
	while((TblPRIOR_hard[0][0] > lam_H | lam_H > TblPRIOR_hard[numH-1][0]) && rejH < 100){
	  lam_H  = gengam(E2,X2-beta_H+1.0);
	  rejH++;
	}
	if(rejH==100 | flag11>=iter){
	  tblH = 0;
	  printf("\n#******************************************************************#");
	  printf("\n# WARNING: THE TABULATED PRIOR FOR THE 'HARD' BAND CONTRADICTS THE #");
	  printf("\n#          LIKELIHOOD OF THE DATA, I.E., THE MLE IS OUTSIDE THE    #");
	  printf("\n#          RANGE OF THE TABULATED PRIOR. NOW BEHR IS FINISHING THE #");
	  printf("\n#          CALCULATION WITH ITS DEFAULT GAMMA PRIOR...             #");
	  printf("\n#******************************************************************#\n");
	}
	else{
	  r = findprob(lam_H,numH,TblPRIOR_hard) - findprob(temp2,numH,TblPRIOR_hard);
	  if(ranf()>exp(r)){ lam_H = temp2; flag11++; }/* reject the draw */ 
	}
      }

      if(iter>burn_in){
	tempR[k]         = lam_S/lam_H;
	tempHR[k]        = (1.0-lam_S/lam_H)/(1.0+lam_S/lam_H);
	tempC[k]         = log10(lam_S)-log10(lam_H);
	if(!dum21 | details | posterior){
	  tempS[k]       = lam_S;
	  tempH[k]       = lam_H;}
	if(!dum22 | details){
	  tempHminusS[k] = lam_H-lam_S;
	  tempHplusS[k]  = lam_H+lam_S;}
	if(!dum23 | details){
	  tempLogS[k]    = log10(lam_S);
	  tempLogH[k]    = log10(lam_H);}
	if(!dum24){
	  fprintf(output5,"%g\t%g\n",lam_S,lam_H);
	  fflush(output5);}
       
	if(k==draws-burn_in-1){ /* begin writing the output of intervals */
	    Q1_HR = (double *) calloc(6, sizeof(double));
	    Q2_HR = (double *) calloc(6, sizeof(double));
	    Q3_HR = (double *) calloc(6, sizeof(double));
  	    summary_fn(tempR,  draws-burn_in, LEVEL, HPD, Q1_HR); free(tempR);
	    summary_fn(tempHR, draws-burn_in, LEVEL, HPD, Q2_HR); free(tempHR);
	    summary_fn(tempC,  draws-burn_in, LEVEL, HPD, Q3_HR); free(tempC);

	  if(!dum21 | details | posterior){
	    Q4_HR = (double *) calloc(6, sizeof(double));
	    Q5_HR = (double *) calloc(6, sizeof(double));
	    summary_fn(tempS, draws-burn_in, LEVEL, HPD, Q4_HR);  free(tempS);
	    summary_fn(tempH, draws-burn_in, LEVEL, HPD, Q5_HR);  free(tempH);}

	  if(!dum22 | details){
	    Q6_HR = (double *) calloc(6, sizeof(double));
	    Q7_HR = (double *) calloc(6, sizeof(double));
	    summary_fn(tempHminusS, draws-burn_in, LEVEL, HPD, Q6_HR); free(tempHminusS);
	    summary_fn(tempHplusS,  draws-burn_in, LEVEL, HPD, Q7_HR); free(tempHplusS);}

	  if(!dum23 | details){
	    Q8_HR = (double *) calloc(6, sizeof(double));
	    Q9_HR = (double *) calloc(6, sizeof(double));
	    summary_fn(tempLogS, draws-burn_in, LEVEL, HPD, Q8_HR); free(tempLogS);
	    summary_fn(tempLogH, draws-burn_in, LEVEL, HPD, Q9_HR); free(tempLogH);}
	 
	  flag1 = flag2 = flag3 = flag4 = flag5 = flag6 = flag7 = flag8 = flag9 = 0;

	  if(Q1_HR[1]>1e+5)                          flag1 = 1;
	  if(Q1_HR[3]<    1e-5 && Q1_HR[4]<    1e-5) flag1 = 1;
	  if(Q2_HR[3]> 0.99999 && Q2_HR[4]> 0.99999) flag2 = 1;
	  if(Q2_HR[3]<-0.99999 && Q2_HR[4]<-0.99999) flag2 = 1;  	  
	  if(Q3_HR[1]>1e+5|Q3_HR[1]<-1e+5)           flag3 = 1;
	  if(details){
	    if(Q8_HR[1]<-1e+5|Q8_HR[1]>1e+5)         flag7 = 1;
	    if(Q9_HR[1]<-1e+5|Q9_HR[1]>1e+5)         flag8 = 1;
	  }

	  if(flag1+flag2+flag3+flag4+flag5+flag6+flag7+flag8+flag9>0 ){
	    printf("\n#************************************************************#");
	    printf("\n# WARNING: The Monte Carlo simulations may be unstable.      #");
	    printf("\n#          Watch out for the hardness ratios marked by '#.'  #");
	    printf("\n#          Use the quadrature option for stable results.     #");
	    printf("\n#************************************************************#\n");
	  }
	
	  printf("\t\tMode\t\tMean\t\tMedian\t\tLower Bound\tUpper Bound\n"); 
	  if(!nooutput){ 
	    fprintf(output1,"\t\tMode\t\tMean\t\tMedian\t\tLower Bound\tUpper Bound\n");
	    fflush(output1);}
	  if(flag1) printf("#");
	  printf("(S/H)\t\t");     for(j=0;j<5;j++) printf("%f\t",Q1_HR[j]);
	  flag2 ? printf("\n#") : printf("\n"); 
	  printf("(H-S)/(H+S)\t"); for(j=0;j<5;j++) printf("%f\t",Q2_HR[j]);
	  flag3 ? printf("\n#") : printf("\n"); 
	  printf("log10(S/H)\t");  for(j=0;j<5;j++) printf("%f\t",Q3_HR[j]);
	  
	  if(!dum21 | details){
	    flag4 ? printf("\n#") : printf("\n"); 
	    printf("S\t\t");       for(j=0;j<5;j++) printf("%f\t",Q4_HR[j]);
	    flag5 ? printf("\n#") : printf("\n"); 
	    printf("H\t\t");       for(j=0;j<5;j++) printf("%f\t",Q5_HR[j]);}
	  
	  if(!dum22 | details){
	    flag6 ? printf("\n#") : printf("\n"); 
	    printf("(H-S)\t\t");   for(j=0;j<5;j++) printf("%f\t",Q6_HR[j]); 
	    flag7 ? printf("\n#") : printf("\n"); 
	    printf("(H+S)\t\t");   for(j=0;j<5;j++) printf("%f\t",Q7_HR[j]);}
	  
	  if(!dum23 | details){ 
	    flag8 ? printf("\n#") : printf("\n"); 
	    printf("log10(S)\t");  for(j=0;j<5;j++) printf("%f\t",Q8_HR[j]);
	    flag9 ? printf("\n#") : printf("\n"); 
	    printf("log10(H)\t");  for(j=0;j<5;j++) printf("%f\t",Q9_HR[j]);}
    	  
	  if(!nooutput){
	    if(flag1) fprintf(output1,"#"); 
	    fprintf(output1,"(S/H)\t\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q1_HR[j]);
	    fprintf(output1,"\n(H-S)/(H+S)\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q2_HR[j]);
	    flag3 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"log10(S/H)\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q3_HR[j]);
	    
	    if(!dum21 | details){
	    flag4 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"S\t\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q4_HR[j]);
	    flag5 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"H\t\t");
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q5_HR[j]);}
	    
	    if(!dum22 | details){
	    flag6 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"(H-S)\t\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q6_HR[j]); 
	    flag7 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"(H+S)\t\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q7_HR[j]);}
	    
	    if(!dum23 | details){ 
	    flag8 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"log10(S)\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q8_HR[j]);
	    flag9 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
	    fprintf(output1,"log10(H)\t"); 
	    for(j=0;j<5;j++) fprintf(output1,"%f\t",Q9_HR[j]);}
          fprintf(output1,"\n"); 
	    fflush(output1);
	  }
	  if(!dum21){
	    for(j=0;j<5;j++) fprintf(output2,"%f\t",Q1_HR[j]);
	    fprintf(output2,"%f\t%f\n",Q4_HR[0],Q5_HR[0]);
	    fflush(output2);
	  }
	  if(!dum22){
	    for(j=0;j<5;j++) fprintf(output3,"%f\t",Q2_HR[j]);
	    fprintf(output3,"%f\t%f\n",Q6_HR[0],Q7_HR[0]);
	    fflush(output3);
	  }	  
	  if(!dum23){
	    for(j=0;j<5;j++) fprintf(output4,"%f\t",Q3_HR[j]);
	    fprintf(output4,"%f\t%f\n",Q8_HR[0],Q9_HR[0]);
	    fflush(output4);
	  }
	}
	k++;
      }
    }
  } /* end Gibbs */
  else{ /* begin Gaussian quadrature */
    /****************************************/
    /*********  Gaussian quadrature  ********/
    /****************************************/     
    dens_R       = (double *)calloc((int) numBins,sizeof(double));     
    dens_HR      = (double *)calloc((int) numBins,sizeof(double));
    dens_C       = (double *)calloc((int) numBins,sizeof(double));

    if(!dum21 | !dum25 | !dum6 | details | posterior){
      dens_lamS  = (double *)calloc((int) numBins,sizeof(double));     
      dens_lamH  = (double *)calloc((int) numBins,sizeof(double));}
    if(!dum22 | details){
      dens_omega = (double *)calloc((int) numBins,sizeof(double)); 
      dens_phi   = (double *)calloc((int) numBins,sizeof(double));}
    if(!dum23 | details){
      dens_LogS  = (double *)calloc((int) numBins,sizeof(double));     
      dens_LogH  = (double *)calloc((int) numBins,sizeof(double));}

    absacc = 0.0;
    relacc = 1.0e-5;
    minrul = 3;
    maxrul = 9;

    lam_S = MAX(1.0e-8,X1-B1/(C1*E1));
    lam_H = MAX(1.0e-10,X2-B2/(C2*E2));
    temp8 = lam_S/lam_H; /* R */

    if(lam_S > lam_H) SWITCH = 1;

    temp1 = sqrt(X1+0.75)+1; /* hat(sig)_S   */
    temp2 = sqrt(X2+0.75)+1; /* hat(sig)_H   */
    temp3 = sqrt(B1+0.75)+1; /* hat(sig)_B.S */
    temp4 = sqrt(B2+0.75)+1; /* hat(sig)_B.H */

    temp5 = ABS(temp8)*sqrt((pow(temp1,2)+pow(temp3/(C1*E1),2))
			    /pow(X1-B1/(C1*E1),2)+(pow(temp2,2)
                            +pow(temp4/(C2*E2),2))/pow(X2-B2/(C2*E2),2));
    temp6 = 2*sqrt(pow(X2-B2/(C2*E2),2)*(pow(temp1,2)+pow(temp3/(C1*E1),2))
		   +pow(X1-B1/(C1*E1),2)*(pow(temp2,2)+pow(temp4/(C2*E2),2)))
                   /pow(X2-B2/(C2*E2)+X1-B1/(C1*E1),2);
    temp7 = temp5/(ABS(temp8)*log(10));

    lamS_left  = MAX(1.0e-100,lam_S-5*temp1);
    lamS_right = lam_S+20*temp1;
    lamH_left  = MAX(1.0e-100,lam_H-5*temp2);
    lamH_right = lam_H+20*temp2;
    R_left     = MAX(1.0e-10,temp8-5*temp5);
    
    if(tblS){
      lamS_left  = TblPRIOR_soft[0][0];
      lamS_right = TblPRIOR_soft[numS-1][0];
      for(i=0,lam_S=0;i<numS;i++) lam_S += TblPRIOR_soft[i][0]*TblPRIOR_soft[i][1];
    }
    if(tblH){
      lamH_left  = TblPRIOR_hard[0][0];
      lamH_right = TblPRIOR_hard[numS-1][0];
      for(i=0,lam_H=0;i<numH;i++) lam_H += TblPRIOR_hard[i][0]*TblPRIOR_hard[i][1];
    }
    if((tblS | tblH) && (lam_S > lam_H)) SWITCH = 1;
    

    if((temp8+5*temp5)>1e+3) R_right = 1e+3; 
    else if(temp5<1) R_right = (temp8+5)*MAX(1,1.8-PRIOR_hard[0]);
    else R_right = (temp8+5*temp5)*MAX(1,1.8-PRIOR_hard[0]);

    if(lam_S>=10 && lam_H>=10){ 
	  if(lam_S>=50 && lam_H>=50){
	    C_left  = log10(temp8)-2.0;
		C_right = log10(temp8)+2.0;
		if(lam_S > lam_H) SWITCH = 0;
	  }
	  else{ /* 1.5 + exp(0.5 * (1 - idx)) */
		C_left  = -1.5-exp(0.5*(1-MAX(0.1,MIN(soft_idx,hard_idx))));
        C_right =  1.5+exp(0.5*(1-MAX(0.1,MIN(soft_idx,hard_idx))));
	  }
    }
    else if(lam_S>=1 && lam_H>=1){ /* 4 + exp(2 * (1 - idx)) */
      C_left  = -4-exp(2*(1-MAX(0.1,MIN(soft_idx,hard_idx))));
      C_right =  4+exp(2*(1-MAX(0.1,MIN(soft_idx,hard_idx))));
    }
    else{ /* 20 * sqrt(1 - idx) + 5 */
      C_left  = -20*sqrt(1-MAX(0.1,MIN(soft_idx,hard_idx)))-5;
      C_right =  20*sqrt(1-MAX(0.1,MIN(soft_idx,hard_idx)))+5;
    }

    omega_left  = lamS_left  + lamH_left;  
    omega_right = lamS_right + lamH_right;      
    phi_left    = lamH_left  - lamS_right;  
    phi_right   = lamH_right - lamS_left; 
    HR_left     = -1.0;
    HR_right    =  1.0;
              
    R  = (double *)calloc((int) numBins,sizeof(double));     
    C  = (double *)calloc((int) numBins,sizeof(double));
    HR = (double *)calloc((int) numBins,sizeof(double));

    ppoints(numBins,R_left,R_right,R);
    ppoints(numBins,C_left,C_right,C);
    segment(numBins,HR_left,HR_right,HR); 

    if(!dum21 | !dum25 | !dum6 | details | posterior){
      lamS  = (double *)calloc((int) numBins,sizeof(double));
      lamH  = (double *)calloc((int) numBins,sizeof(double));
      ppoints(numBins,lamS_left,lamS_right,lamS);
      ppoints(numBins,lamH_left,lamH_right,lamH);}
    if(!dum22 | details){
      omega = (double *)calloc((int) numBins,sizeof(double));
      phi   = (double *)calloc((int) numBins,sizeof(double));
      ppoints(numBins,omega_left,omega_right,omega);
      ppoints(numBins,phi_left,phi_right,phi);}
    if(!dum23 | details){ 
      LogS  = (double *)calloc((int) numBins,sizeof(double));
      LogH  = (double *)calloc((int) numBins,sizeof(double));
      ppoints(numBins,log10(lamS_left)-2,log10(lamS_right)+2,LogS);
      ppoints(numBins,log10(lamH_left)-2,log10(lamH_right)+2,LogH);}
    norm1   = maximum(lamS_left,lamS_right,X1,B1,C1,E1,prior1);
    norm2   = maximum(lamH_left,lamH_right,X2,B2,C2,E2,prior2);

    if(SWITCH){
      X1_t = X1;         X2_t = X2;
      B1_t = B1;         B2_t = B2;
      C1_t = C1;         C2_t = C2;
      E1_t = E1;         E2_t = E2;
      norm1_t = norm1;   norm2_t = norm2;
      prior1_t = prior1; prior2_t = prior2;

      X1 = X2;           X2 = X1_t;
      B1 = B2;           B2 = B1_t;
      C1 = C2;           C2 = C1_t;
      E1 = E2;           E2 = E1_t;
      norm1 = norm2;     norm2 = norm1_t;
      prior1 = prior2;   prior2 = prior1_t;
    }

    if(lam_S>3 && lam_H>3){
      if(SWITCH) zeta_left = -2 / sqrt(0.5*soft_idx);
      else       zeta_left = -2 / sqrt(0.5*hard_idx); 
    }
    else if(lam_S<=3 && lam_H<=3) zeta_left = -100;
    else zeta_left = -30;

    if((SWITCH && lam_S>100) | (!SWITCH && lam_H>100)) zeta_right = 7;
    else if(SWITCH) zeta_right = MIN(5,MAX(2,log(lamS_right))); 
    else            zeta_right = MIN(5,MAX(2,log(lamH_right))); 

    for(k=0;k<numBins;k++){
      /******* computing log10(S/H) *******/
      K1   = pow(10,C[k]);
      K2   = 1.0;
      K3   = pow(10,C[k])*pow(log(10),2);
      MODE = 3;
      pquad_(post_, &zeta_left, &zeta_right, &absacc, &relacc, 
	     &minrul, &maxrul, &ruleno, quads, &nevals, &ercode);
      if(SWITCH) dens_C[numBins-k-1] = quads[ruleno-1];
      else       dens_C[k]           = quads[ruleno-1];
    }

    if(SWITCH){
      X1 = X1_t;         X2 = X2_t;
      B1 = B1_t;         B2 = B2_t;
      C1 = C1_t;         C2 = C2_t;
      E1 = E1_t;         E2 = E2_t;
      norm1 = norm1_t;   norm2 = norm2_t;
      prior1 = prior1_t; prior2 = prior2_t;
    }

    for(k=0;k<numBins;k++){
      /********** computing S/H ***********/
      K1 = R[k];
      K2 = 1.0;
      K3 = 1.0;
      MODE = 1;
      pquad_(post_, &lamH_left, &lamH_right, &absacc, &relacc, 
	     &minrul, &maxrul, &ruleno, quads, &nevals, &ercode);
      dens_R[k] = quads[ruleno-1];

      /******* computing (H-S)/(H+S) ******/
      K1   = (1-HR[k])/2.0;
      K2   = (1+HR[k])/2.0;
      K3   = 1/2.0;
      MODE = 1;
      pquad_(post_, &omega_left, &omega_right, &absacc, &relacc, 
	     &minrul, &maxrul, &ruleno, quads, &nevals, &ercode);
      dens_HR[k] = quads[ruleno-1];
      
      if(DISPLAY){
	esterr = fabs(quads[ruleno - 2] - dens_HR[k]) / dens_HR[k];
	printf ("\nEstimated integral = %22.15e\n", dens_HR[k]);
	printf ("Estimated relative error = %10.3e\n", esterr);
	printf ("Number of function evaluations = %4d\n", nevals);
	printf ("Error code = %2d\n", ercode);
      }

      if(fix){
	dens_HR[k-1] = (dens_HR[k-2] + dens_HR[k])/2;
	fix = 0;
      }
      if(nevals < nevals_pr) fix = 1;     
      nevals_pr = nevals;
      
      if(!dum21 | !dum25 | details | posterior){
      /*********** computing S ************/
      K1 = 1.0;
      dens_lamS[k] = dens_S(lamS[k]);
     
      /*********** computing H ************/
      K1 = 1.0;
      dens_lamH[k] = dens_H(lamH[k]);}
     
      if(!dum22 | details){
      /********* computing (H+S) **********/
      K1   = omega[k]/2.0;
      K2   = omega[k]/2.0;
      K3   = omega[k]/2.0;
      MODE = 0;
      pquad_(post_, &HR_left, &HR_right, &absacc, &relacc, 
	     &minrul, &maxrul, &ruleno, quads, &nevals, &ercode);
      dens_omega[k] = quads[ruleno-1];
      
      /********* computing (H-S) **********/
      K1   = phi[k];
      K2   = 1.0;
      K3   = 1.0;
      MODE = 2;
      lamH_left = MAX(0.0, phi[k]);
      pquad_(post_, &lamH_left, &lamH_right, &absacc, &relacc, 
	     &minrul, &maxrul, &ruleno, quads, &nevals, &ercode);
      dens_phi[k] = quads[ruleno-1];}
      
      if(!dum23 | details){      
      /********* computing LogS ***********/
      K1 = pow(10,LogS[k])*log(10);
      dens_LogS[k] = dens_S(pow(10,LogS[k]));
      
      /********* computing LogH ***********/
      K1 = pow(10,LogH[k])*log(10);
      dens_LogH[k] = dens_H(pow(10,LogH[k]));}  
    }
 
    dum7=dum8=dum9=dum10=dum11=dum12=dum13=dum14=dum15=dum16=dum17=dum18=0;

    Q1_HR = (double *) calloc(6, sizeof(double));
    Q2_HR = (double *) calloc(6, sizeof(double));
    Q3_HR = (double *) calloc(6, sizeof(double));

    CLASS=1; HPD_fn(dens_R,R,LEVEL,0,Q1_HR);         if(LIMIT) dum7 =1;
    CLASS=2; HPD_fn(dens_HR,HR,LEVEL,0,Q2_HR); 
    CLASS=0; HPD_fn(dens_C,C,LEVEL,1,Q3_HR);         if(LIMIT) dum8 =1; 

    if(Q1_HR[0]==R[0])      Q1_HR[0]= 0;           /*  R has a mode at 0      */
    if(Q1_HR[3]==R[0])     {Q1_HR[3]= 0; dum16=1;} /*  R has a lower bound  0 */
    if(Q2_HR[3]<-9.9989e-1){Q2_HR[3]=-1; dum15=1;} /* HR has a lower bound -1 */
    if(Q2_HR[4]> 9.9989e-1){Q2_HR[4]= 1; dum15=1;} /* HR has a upper bound  1 */      
    if(Q2_HR[0]<-9.9989e-1) Q2_HR[0]=-1;           /* HR has a mode at -1     */
    if(Q2_HR[0]> 9.9989e-1) Q2_HR[0]= 1;           /* HR has a mode at  1     */
    if(Q2_HR[0]<0.0 && Q2_HR[0]>-1e-10) Q2_HR[0]=0;/* HR has a mode at  1     */
    
    if(!dum21 | !dum25 | details | posterior){
      Q4_HR = (double *) calloc(6, sizeof(double));    
      Q5_HR = (double *) calloc(6, sizeof(double));
      
      CLASS=1; HPD_fn(dens_lamS,lamS,LEVEL,0,Q4_HR);   if(LIMIT) dum9 =1; 
      CLASS=1; HPD_fn(dens_lamH,lamH,LEVEL,0,Q5_HR);   if(LIMIT) dum10=1;

      if(Q4_HR[0]==lamS[0])      Q4_HR[0]= 0;           /*  S has a mode at 0      */
      if(Q4_HR[3]==lamS[0])     {Q4_HR[3]= 0; dum17=1;} /*  S has a lower bound  0 */
      if(Q5_HR[0]==lamH[0])      Q5_HR[0]= 0;           /*  H has a mode at 0      */
      if(Q5_HR[3]==lamH[0])     {Q5_HR[3]= 0; dum18=1;} /*  H has a lower bound  0 */		
    } 

    if(!dum22 | details){
      Q6_HR = (double *) calloc(6, sizeof(double));
      Q7_HR = (double *) calloc(6, sizeof(double));

      CLASS=0; HPD_fn(dens_phi,phi,LEVEL,0,Q6_HR);     if(LIMIT) dum11=1; 
      CLASS=1; HPD_fn(dens_omega,omega,LEVEL,0,Q7_HR); if(LIMIT) dum12=1;}  

    if(!dum23 | details){
      Q8_HR = (double *) calloc(6, sizeof(double));
      Q9_HR = (double *) calloc(6, sizeof(double));

      CLASS=0; HPD_fn(dens_LogS,LogS,LEVEL,1,Q8_HR);   if(LIMIT) dum13=1; 
      CLASS=0; HPD_fn(dens_LogH,LogH,LEVEL,1,Q9_HR);   if(LIMIT) dum14=1;} 
  
    if(dum7+dum8+dum9+dum10+dum11+dum12+dum13+dum14>0){
                printf("#******************************************************#");
      if(dum7+dum8+dum9+dum10+dum11+dum12+dum13+dum14>1)      
	        printf("\n# WARNING: Due to numerical constraints, the limits of #");
      else if(dum7+dum8+dum9+dum10+dum11+dum12+dum13+dum14==1)      
		printf("\n# WARNING: Due to numerical constraints, the limit of  #");
      if(dum7){ printf("\n#          (S/H)                                       #"); flag1=1;}
      if(dum8){ printf("\n#          log10(S/H)                                  #"); flag2=1;}
      if(dum9){ printf("\n#          S                                           #"); flag3=1;}
      if(dum10){printf("\n#          H                                           #"); flag4=1;}
      if(dum11){printf("\n#          H-S                                         #"); flag5=1;}
      if(dum12){printf("\n#          H+S                                         #"); flag6=1;}
      if(dum13){printf("\n#          log10(S)                                    #"); flag7=1;}
      if(dum14){printf("\n#          log10(H)                                    #"); flag8=1;}
      if(dum7+dum8+dum9+dum10+dum11+dum12+dum13+dum14>1)   
     	        printf("\n#          are confined to pre-determined range.       #");
      else if(dum7+dum8+dum9+dum10+dum11+dum12+dum13+dum14==1) 
                printf("\n#          is confined to pre-determined range.        #");
                printf("\n#******************************************************#\n");
    }

    printf("\t\tMode\t\tMean\t\tMedian\t\tLower Bound\tUpper Bound\n"); 
    if(!nooutput){ 
      fprintf(output1,"\t\tMode\t\tMean\t\tMedian\t\tLower Bound\tUpper Bound\n");
      fflush(output1);
    }

    if(flag1) printf("#"); 
    printf("(S/H)\t\t");     for(j=0;j<5;j++) printf("%f\t",Q1_HR[j]); 
    printf("\n(H-S)/(H+S)\t"); for(j=0;j<5;j++) printf("%f\t",Q2_HR[j]);
    flag2 ? printf("\n#") : printf("\n"); 
    printf("log10(S/H)\t");  for(j=0;j<5;j++) printf("%f\t",Q3_HR[j]);
    
    if(!dum21 | !dum25 | details){
    flag3 ? printf("\n#") : printf("\n"); 
    printf("S\t\t");         for(j=0;j<5;j++) printf("%f\t",Q4_HR[j]);
    flag4 ? printf("\n#") : printf("\n"); 
    printf("H\t\t");         for(j=0;j<5;j++) printf("%f\t",Q5_HR[j]);}
    
    if(!dum22 | details){
    flag5 ? printf("\n#") : printf("\n");  
    printf("(H-S)\t\t");     for(j=0;j<5;j++) printf("%f\t",Q6_HR[j]); 
    flag6 ? printf("\n#") : printf("\n");  
    printf("(H+S)\t\t");     for(j=0;j<5;j++) printf("%f\t",Q7_HR[j]);}
    
    if(!dum23 | details){
    flag7 ? printf("\n#") : printf("\n");      
    printf("log10(S)\t");    for(j=0;j<5;j++) printf("%f\t",Q8_HR[j]);
    flag8 ? printf("\n#") : printf("\n");
    printf("log10(H)\t");    for(j=0;j<5;j++) printf("%f\t",Q9_HR[j]);}

    if(!nooutput){
      if(flag1) fprintf(output1,"#"); 
      fprintf(output1,"(S/H)\t\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q1_HR[j]);
      fprintf(output1,"\n(H-S)/(H+S)\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q2_HR[j]);
      flag2 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"log10(S/H)\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q3_HR[j]);
      
      if(!dum21 | !dum25 | details){
      flag3 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"S\t\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q4_HR[j]);
      flag4 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"H\t\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q5_HR[j]);}
      
      if(!dum22 | details){	
      flag5 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"(H-S)\t\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q6_HR[j]); 
      flag6 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"(H+S)\t\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q7_HR[j]);}
      
      if(!dum23 | details){
      flag7 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"log10(S)\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q8_HR[j]);
      flag8 ? fprintf(output1,"\n#") : fprintf(output1,"\n"); 
      fprintf(output1,"log10(H)\t"); 
      for(j=0;j<5;j++) fprintf(output1,"%f\t",Q9_HR[j]);} 
      fprintf(output1,"\n");
      fflush(output1);
    }
    if(!dum21){
      for(j=0;j<5;j++) fprintf(output2,"%f\t",Q1_HR[j]);
      fprintf(output2,"%f\t%f\n",Q4_HR[0],Q5_HR[0]);
      fflush(output2);
    }
    if(!dum22){
      for(j=0;j<5;j++) fprintf(output3,"%f\t",Q2_HR[j]);
      fprintf(output3,"%f\t%f\n",Q6_HR[0],Q7_HR[0]);
      fflush(output3);
    }	  
    if(!dum23){
      for(j=0;j<5;j++) fprintf(output4,"%f\t",Q3_HR[j]);
      fprintf(output4,"%f\t%f\n",Q8_HR[0],Q9_HR[0]);
      fflush(output4);
    }
    if(!dum25){
      if(DISPLAY) 
	fprintf(output6,"R\t densR\t HR\t densHR\t C\t densC\t S\t densS\t H\t densH\n"); 
      for(k=0;k<numBins;k++){ 
	fprintf(output6,"%f %14g\t",R[k],dens_R[k]);
	fprintf(output6,"%f %14g\t",HR[k],dens_HR[k]);
	fprintf(output6,"%f %14g\t",C[k],dens_C[k]);
	fprintf(output6,"%f %14g\t",lamS[k],dens_lamS[k]);
	fprintf(output6,"%f %14g\n",lamH[k],dens_lamH[k]);
      }
      fflush(output6);
    }
  } /* end Gaussian quadrature */
  
  if(posterior){
    param[0] = pow(Q4_HR[1],2)/Q4_HR[5];
    param[1] = pow(Q5_HR[1],2)/Q5_HR[5];
    param[2] = Q4_HR[1]/Q4_HR[5];
    param[3] = Q5_HR[1]/Q5_HR[5];
  }

  dur_std=(double)(clock()-tm_std)/CLOCKS_PER_SEC;
  printf("\n#==========================================#\n");
  if(posterior){
    printf("# Posterior distributions of lamS and lamH #\n");
    printf("#------------------------------------------#\n");
    printf("# X~Gamma(a,b) where E(X)=a/b, V(X)=a/b^2  #\n");
    printf("#------------------------------------------#\n");
    printf("# p(lamS|data) = Gamma(%8g,%9g) #\n",param[0],param[2]);
    printf("# p(lamH|data) = Gamma(%8g,%9g) #\n",param[1],param[3]);
    printf("#------------------------------------------#\n");}
  if(dur_std<60.0) printf("# Total amount of running time =%6g sec #",dur_std); 
  else printf("# Total amount of running time =%6g min #",dur_std/60); 
  printf("\n#==========================================#\n");
  exit(0);
} /* main */
