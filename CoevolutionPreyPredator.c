/*  Created by Weini Huang on 10/24/17.
 *  Coevolution Predator-Prey Dynamics + Competitive Lotka-Volterra Dynamics
 *  Evolution of prey species in a system with one fixed predator type (set up the mutation rate of predator as 0)
 *  Copyright 2012 __WeiniHuang__. All rights reserved.
 *  Diversity:ShannonIndex+SimpsonIndex and number of types
 *  g k value in a certain range between [0.02,1.0) and [0.01,Kmax)
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mtwist.h"
#define pi 3.1415926



double Prey[100][2];// 1.g-values of the prey species 2.population size of all prey species
double Prey2[100][2];//renew the Prey[][] every 99 mutations of prey or predator 
double TimeBirthPrey[100];//index as prey index, prey reproduce
double TimeDeathPrey[100];//index as prey index, prey dies
double TimePredation_Prey[100][100];//index as prey and predator index, prey dies
double TimePredation_Predator[100][100];//index as prey and predator index, prey dies, predator reproduces
double TimeDeathPredator[100];//index as predator index, predator dies
double Predator[100][2];// 1.k-values of the predator species 2.population size of all predator species
double Predator2[100][2];//renew Predator[][] every 99 mutations 
int CT[100];//the index of current prey types
int Num_CT;//Number of current types
int Num_PreyHistory;//Number of all specis in the history, need to reset to be Num_CT_Prey every 99 mutations
int CT_Predator[100];//the index of current prey types
int Num_CT_Predator;//Number of current types
int Num_PredatorHistory;//Number of all specis in the history, need to reset to be Num_CT_Predator every 99 mutations
double TimeResCompPrey[100][100]; //index as prey index, resource competition, the diagonal number is competition between individuals of the same prey species; TimeResCompPrey[i][j]; indicate the competition of i and j, which results i survivies and j dies.
double ShannonIndex_Prey;
double ShannonIndex_Predator;
double SimpsonIndex_Prey;
double SimpsonIndex_Predator;
int NumberNoiseFree_Prey; //number of prey types whose abuandances are above the threshold of Noise - percentage
int NumberNoiseFree_Predator;//number of predator types whose abuandances are above the threshold of Noise - percentage
int NumberNoiseFree2_Prey; //number of prey types whose abuandances are above the threshold of Noise - absolute number
int NumberNoiseFree2_Predator;//number of predator types whose abuandances are above the threshold of Noise - absolute number
double Threshold_Noise; // The threshold of the noise, like 1%, 2%,...
double Threshold2_Noise;// The threshold of the noise, like absolute number 10, 20,...

//generate random numbers for the payoff entries
double RandomGauss(double mean, double variance)
{     double z;
    double rnd1;
    rnd1 = mt_drand();
    double rnd2;
    rnd2 = mt_drand();
    double x,y;
    x=sqrt(-2*log(rnd1));
    y=cos(2*pi*rnd2);
    z=mean+x*y*sqrt(variance);// the mean of the number of this function is zero
    return z;
}

//Extinction, remove from CT
void Extinct_Prey(int x) // removes extinct type x
{
	int y=x; // x is the index in CT[]
	//delete in CT[]
	int i;
	i=CT[x];
	//printf("extinction of prey type %d, x=%d, PreyNum=%lf,Num_CT=%d\n", CT[x],x,Prey[i][1],Num_CT);
	for(;CT[y]!=1000;y++)
	CT[y]=CT[y+1];
	Num_CT--;
	
}

void Extinct_Predator(int x) // removes extinct type x
{
	int y=x; // x is the index in CT_Predator[]
	//delete in CT_Predator[]
	int i;
	i=CT_Predator[x];
	//printf("extinction of predator type %d, x=%d, PredatorNum=%lf,Num_CT_Predator=%d\n", CT_Predator[x],x,Predator[i][1],Num_CT_Predator);
	for(;CT_Predator[y]!=1000;y++)
	CT_Predator[y]=CT_Predator[y+1];
	Num_CT_Predator--;
}

//Mutation of Prey
void MutatePrey(int x)
{
	int y=x;//x is already the index in prey[][];
	
	Num_CT++;
	Num_PreyHistory++;
	
	//add the new type in CT[]
	CT[Num_CT-1]=Num_PreyHistory;  //maybe it is better to be : CT[Num_CT-1]=Num_PreyHistory-1; but it will not influence the results
	//add the new type in prey[][]
	Prey[Num_PreyHistory][1]=1.;
	Prey[Num_PreyHistory][0]=mt_drand()*0.98+0.02;//give g-value between [0.02,1)
	//printf("muNum=%d,g=%lf\n",Num_PreyHistory,Prey[Num_PreyHistory][0]);   
}

//Mutation of Predator
void MutatePredator(int x, double k_max)
{
	int y=x; //x is the index in predator[][];
	Num_CT_Predator++;
	Num_PredatorHistory++;
	CT_Predator[Num_CT_Predator-1]=Num_PredatorHistory;
	Predator[Num_PredatorHistory][1]=1.;
	Predator[Num_PredatorHistory][0]=mt_drand()*(k_max-0.01)+0.01;//give k-value of the predator mutate type
}

//Growth Function of Prey
double GrowthFunction(double g) //how fast a prey grow
{
	double Growth;
	Growth=g;
	return Growth;
}

//Predation Function based on Prey and predator
double DefenceFunction(double g, double k, double m, double p, double k_max)//how fast a prey can be consumed by a predator
{
	double Defence;
	double powervalue;
	powervalue=(m*k)/k_max;
	Defence=pow(g,powervalue)*p;
	return Defence;
}


//reactions from reation rate to reation time, and the reaction happening next is the one with the shortest time
double ReactionRateToTime(double ReactionRate)
{
    double timeR;
    timeR=(1./ReactionRate)*log(1./(1.-mt_drand())); //mt_drand() gives random numbers between [0,1)
    return timeR;
}


            
int main (int argc, const char*argv[]) //see the definition of argv[] below 

{   
	mt_goodseed(); //Automatically seeds from /dev/random or system time (can be slow)
	
	double Prey_b,Prey_d,Predator_d,Predation_p;
	double Realization;
	double ShortestTime;
	double Time;
	int ReactionID_X;
	int ReactionID_Y;
	int BirthDeathIndex;//for the same index of reactionsIDs, 1 for birth, 0 for death
	double mu_Prey,mu_Predator;
	double rnd1,rnd2;
	int n,m,count,i,j;
	int PrintIndex;
	int Notice;//notice once that all prey speceis went extinct
	double power;// the power of g for the DefenceFunction
	Notice=0;
	double Resource_c;// resource competition rate for all prey individuals
	int ResComp_index;//will be used as an index to caculate rates in the resource_competition 
	int ResComp_index_i;
	int CompetitionIndex;//1 for resource competition, 0 for no competition events
	double k_max;// the max value of the rate of predator species can reprudce by consuming, between 0 and 1; k_max=0.5 means that by eating 2 preys the predator can reproduce one offspring
	double defaultg;
	double defaultk;
	
	
	for(n=0;n<100;n++)
	   {   
		Prey[n][0]=0.;
	    Prey[n][1]=0.;
	  }
	  
  	for(n=0;n<100;n++)
  	   {   
  		Predator[n][0]=0.;
  	    Predator[n][1]=0.;
  	  }
	
	Realization=atoi(argv[1]); //how long the program will be run
	
	//argv[]
	//prey
    Prey[0][0]=atof(argv[2]);//g value of prey species 0
	Prey[0][1]=atof(argv[3]);//initial number of prey 0
    Prey_b=atof(argv[4]);//basic birth rate of prey
    Prey_d=atof(argv[5]);//basic death rate of prey
	
	Num_CT=1;
	Num_PreyHistory=1;
	
	//predator
	Predator[0][1]=atof(argv[6]);//initial number of predator
	Predator_d=atof(argv[7]);
	Predation_p=atof(argv[8]);
	Predator[0][0]=atof(argv[9]);//how fast the first predator reproduce according to predation, k belongs to (0,Kmax)
	Num_CT_Predator=1;
	Num_PredatorHistory=1; 
	
	mu_Prey=atof(argv[10]);
	power=atof(argv[11]); // the power index of the predator function. m 
	Resource_c=atof(argv[12]);
	k_max=atof(argv[13]); // the max value of the rate of predator species can reprudce by consuming, between 0 and 1; k_max=0.5 means that by eating 2 preys the predator can reproduce one offspring
	mu_Predator=atof(argv[14]);
	int TimeToStop=atoi(argv[15]);//time to stop the programm
	Threshold_Noise=atof(argv[16]);//Threshold of the noise to identify number of types
	Threshold2_Noise=atof(argv[17]);//Threshold of the noise to identify number of types
	int NumPrintForGKValues=atoi(argv[18]);//how many time points to print the g, k values 
	int FileNum=atoi(argv[19]); // the index of the realisation
	
	double TimeGapToPrintForGKValues;// TimeToStop/NumPrintForGKValues
	TimeGapToPrintForGKValues=(double)TimeToStop/(double)NumPrintForGKValues; // print g k values every xxx time
	int NumPrinted_GK=1;
	
    //initialize CT[], current prey types
    int l;
    for(l=0;l<100;l++) {CT[l]=1000;}
	CT[0]=0;
	
	 //initialize CT[], current predator types
    for(l=0;l<100;l++) {CT_Predator[l]=1000;}
	CT_Predator[0]=0;
	
	//define output files
	char arq2[300],arq3[300],arq4[300],arq5[300];
	
	FILE *fp2,*fp3,*fp4,*fp5;
	
	if(mu_Predator==0.)
	sprintf(arq2,"/home/outputDirectoryForEvolutionOfOnlyThePrey/Time_Prey_Predator_TotalAbundance_m%g_File%d.txt",power,FileNum);
	
	if(mu_Predator>0.)
	sprintf(arq2,"/home/outputDirectoryForCoevolution/Time_Prey_Predator_TotalAbundance_m%g_File%d.txt",power,FileNum);
	
	fp2=fopen(arq2,"a");
	if(fp2== NULL)
	{
	fprintf(stderr,"Can't open input file in list!\n");
	exit(1);
	}
	
	
	if(mu_Predator==0.)
	sprintf(arq3,"/home/outputDirectoryForEvolutionOfOnlyThePrey/Gvalues/Time_Prey_Gvalues_longruns_m%g_File%d.txt",power,FileNum);
	
	if(mu_Predator>0.)
	sprintf(arq3,"/home/outputDirectoryForCoevolution/Gvalues/Time_Prey_Gvalues_longruns_m%g_File%d.txt",power,FileNum);
	fp3=fopen(arq3,"a");
	if(fp3== NULL)
	{
	fprintf(stderr,"Can't open input file in list!\n");
	exit(1);
	}
	
	if(mu_Predator==0.)
	sprintf(arq4,"/home/outputDirectoryForEvolutionOfOnlyThePrey/Kvalues/Time_Predator_Kvalues_longruns_m%g_File%d.txt",power,FileNum);
	
	if(mu_Predator>0.)
	sprintf(arq4,"/home/outputDirectoryForCoevolution/Kvalues/Time_Predator_Kvalues_longruns_m%g_File%d.txt",power,FileNum);
	fp4=fopen(arq4,"a");
	if(fp4== NULL)
	{
	fprintf(stderr,"Can't open input file in list!\n");
	exit(1);
	}
	
	
	
	if(mu_Predator==0.)
	sprintf(arq5,"/home/outputDirectoryForEvolutionOfOnlyThePrey/Diversity/Time_Prey_Predator_Diversity_m%g_File%d.txt",power,FileNum);
	
	if(mu_Predator>0.)
	sprintf(arq5,"/home/outputDirectoryForCoevolution/Diversity/Time_Prey_Predator_Diversity_m%g_File%d.txt",power,FileNum);
	
	fp5=fopen(arq5,"a");
	if(fp5== NULL)
	{
	fprintf(stderr,"Can't open input file in list!\n");
	exit(1);
	}
	

	
	
	Time=0;
	PrintIndex=0;
	double SumPrey;
	double SumPredator;
	
	double TestSumPrey;
	double TestSumPredator;
	
	//print the data when time =0.
    SumPrey=0.;
	SumPredator=0.;
	fprintf(fp2,"%lf\t",Time);
	for(m=0;m<100;m++) 
		{ 
			i=CT[m];
			if(CT[m]!=1000) SumPrey=SumPrey+Prey[i][1];
			else break;
		}	
		
	for(m=0;m<100;m++) 
		{
				i=CT_Predator[m];
				if(CT_Predator[m]!=1000) SumPredator=SumPredator+Predator[i][1];
				else break;
		 }
	fprintf(fp2,"%lf\t%lf\n",SumPrey,SumPredator);
 	
	
	ShannonIndex_Prey=0;
	ShannonIndex_Predator=0;
	SimpsonIndex_Prey=0;
	SimpsonIndex_Predator=0;
	NumberNoiseFree_Prey=0;
	NumberNoiseFree_Predator=0;
	NumberNoiseFree2_Prey=0;
	NumberNoiseFree2_Predator=0;
	TestSumPrey=0;
	TestSumPredator=0;
	fprintf(fp5,"%lf\t",Time);
	
	for(m=0;m<100;m++)
	{
	if(CT[m]!=1000)
	 {
	  i=CT[m];
	  ShannonIndex_Prey=ShannonIndex_Prey+(Prey[i][1]/SumPrey)*log(Prey[i][1]/SumPrey);
	  SimpsonIndex_Prey=SimpsonIndex_Prey+(Prey[i][1]/SumPrey)*(Prey[i][1]/SumPrey);
	  if((Prey[i][1]/SumPrey)>Threshold_Noise) NumberNoiseFree_Prey++;
	  if(Prey[i][1]>Threshold2_Noise) NumberNoiseFree2_Prey++;
	  TestSumPrey=TestSumPrey+Prey[i][1];
	 }
	else break;
	}   
    if(TestSumPrey!=SumPrey) printf("Mistake in counting prey types: SumPrey=%lf,TestSumPredy=%lf.\n",SumPrey,TestSumPrey);
    fprintf(fp5,"%d\t%d\t%d\t%lf\t%lf\t",m,NumberNoiseFree_Prey,NumberNoiseFree2_Prey,-1.0*ShannonIndex_Prey,1.-SimpsonIndex_Prey);		

	for(m=0;m<100;m++) 
	{
   if(CT_Predator[m]!=1000) 
	   {
		i=CT_Predator[m];
		ShannonIndex_Predator=ShannonIndex_Predator+(Predator[i][1]/SumPredator)*log(Predator[i][1]/SumPredator);
		SimpsonIndex_Predator=SimpsonIndex_Predator+(Predator[i][1]/SumPredator)*(Predator[i][1]/SumPredator);
		if((Predator[i][1]/SumPredator)>Threshold_Noise) NumberNoiseFree_Predator++;
		if(Predator[i][1]>Threshold2_Noise/2.) NumberNoiseFree2_Predator++;
		TestSumPredator=TestSumPredator+Predator[i][1];
	   }
	else break;
	}
	if(TestSumPredator!=SumPredator) printf("Mistake in counting predator types: SumPredator=%lf,TestSumPredator=%lf.\n",SumPredator,TestSumPredator);
	fprintf(fp5,"%d\t%d\t%d\t%lf\t%lf\n",m,NumberNoiseFree_Predator,NumberNoiseFree2_Predator,-1.0*ShannonIndex_Predator,1.-SimpsonIndex_Predator);	
	
	
	//start the population dynamics
	n=0;
	while(Time<=TimeToStop)
	{
	n++;
	while ((Num_PreyHistory<99)&&(Num_PredatorHistory<99&&Time<=TimeToStop)) //when the number of prey or predator species reaches 100, renew all the matrix
	{
		
		//print the data to output files
		PrintIndex++;
		if(PrintIndex==1000)
		{
		PrintIndex=0;
		
	    SumPrey=0.;
		SumPredator=0.;
		fprintf(fp2,"%lf\t",Time);
		for(m=0;m<100;m++) 
			{ 
				i=CT[m];
				if(CT[m]!=1000) SumPrey=SumPrey+Prey[i][1];
				else break;
			}	
			
			for(m=0;m<100;m++) 
				{
					i=CT_Predator[m];
					if(CT_Predator[m]!=1000) SumPredator=SumPredator+Predator[i][1];
					else break;
				}
		fprintf(fp2,"%lf\t%lf\n",SumPrey,SumPredator);	
		
		ShannonIndex_Prey=0;
		ShannonIndex_Predator=0;
		SimpsonIndex_Prey=0;
		SimpsonIndex_Predator=0;
		NumberNoiseFree_Prey=0;
		NumberNoiseFree_Predator=0;
		NumberNoiseFree2_Prey=0;
		NumberNoiseFree2_Predator=0;
		TestSumPrey=0;
		TestSumPredator=0;
		fprintf(fp5,"%lf\t",Time);
		for(m=0;m<100;m++)
		{
		if(CT[m]!=1000)
		 {
		  i=CT[m];
		  ShannonIndex_Prey=ShannonIndex_Prey+(Prey[i][1]/SumPrey)*log(Prey[i][1]/SumPrey);
		  SimpsonIndex_Prey=SimpsonIndex_Prey+(Prey[i][1]/SumPrey)*(Prey[i][1]/SumPrey);
		  if((Prey[i][1]/SumPrey)>Threshold_Noise) NumberNoiseFree_Prey++;
		  if(Prey[i][1]>Threshold2_Noise) NumberNoiseFree2_Prey++;
		  TestSumPrey=TestSumPrey+Prey[i][1];
		 }
		else break;
		}   
	    if(TestSumPrey!=SumPrey) printf("Mistake in counting prey types: SumPrey=%lf,TestSumPredy=%lf,m=%d,Num_CT=%d.\n",SumPrey,TestSumPrey,m,Num_CT);
	    fprintf(fp5,"%d\t%d\t%d\t%lf\t%lf\t",m,NumberNoiseFree_Prey,NumberNoiseFree2_Prey,-1.0*ShannonIndex_Prey,1.-SimpsonIndex_Prey);		
	
		for(m=0;m<100;m++) 
		{
	   if(CT_Predator[m]!=1000) 
		   {
			i=CT_Predator[m];
			ShannonIndex_Predator=ShannonIndex_Predator+(Predator[i][1]/SumPredator)*log(Predator[i][1]/SumPredator);
			SimpsonIndex_Predator=SimpsonIndex_Predator+(Predator[i][1]/SumPredator)*(Predator[i][1]/SumPredator);
			if((Predator[i][1]/SumPredator)>Threshold_Noise) NumberNoiseFree_Predator++;
			if(Predator[i][1]>Threshold2_Noise/2.) NumberNoiseFree2_Predator++;
			TestSumPredator=TestSumPredator+Predator[i][1];
		   }
		else break;
		}
		if(TestSumPredator!=SumPredator) printf("Mistake in counting predator types: SumPredator=%lf,TestSumPredator=%lf,m=%d,Num_CT_Predator=%d.\n",SumPredator,TestSumPredator,m,Num_CT_Predator);
		fprintf(fp5,"%d\t%d\t%d\t%lf\t%lf\n",m,NumberNoiseFree_Predator,NumberNoiseFree2_Predator,-1.0*ShannonIndex_Predator,1.-SimpsonIndex_Predator);	
		
		//only print in certain time points, defined by TimeGapToPrintForGKValues
		if(Time>=TimeGapToPrintForGKValues*NumPrinted_GK-1.)
		{//Begin-print gk values
		NumPrinted_GK++;//increase the print to the next time point
		fprintf(fp3,"%lf\t",Time);
		for(m=0;m<100;m++)
		{
			if(CT[m]!=1000) { i=CT[m]; fprintf(fp3,"%lf\t%lf\t",Prey[i][0],Prey[i][1]/SumPrey);}
			else break;
		}
		fprintf(fp3,"\n");
	
		fprintf(fp4,"%lf\t",Time);
		for(m=0;m<100;m++)
		{
			if(CT_Predator[m]!=1000) { i=CT_Predator[m];fprintf(fp4,"%lf\t%lf\t",Predator[i][0],Predator[i][1]/SumPredator);}
			else break;
		}
	   fprintf(fp4,"\n");
       }//End - print gk values	
	   
	    }// End of if(PrintIndex==1000)
		
		
		//do some check of the program 
		for(count=0;CT[count]!=1000;count++);//figure out the number of prey species and give it to count
		if(count!=Num_CT) 
			{
			 i=CT[count-1];		
			 printf("\n The number of the current prey types are counted wrong. count=%d,index=%d,Indi=%lf,Num_CT=%d, TimeSteps=%d (Mistake in main function?)\n", count,i,Prey[i][1],Num_CT,n);
		     break;}
			 
	 	for(count=0;CT_Predator[count]!=1000;count++);//figure out the number of predator species and give it to count
	 	if(count!=Num_CT_Predator) 
	 			{
	 			 i=CT_Predator[count-1];
	 			 printf("\n The number of the current predator types are counted wrong. count=%d,index=%d,Indi=%lf,Num_CT_Predator=%d, TimeSteps=%d (Mistake in main function?)\n", count,i,Predator[i][1],Num_CT_Predator,n);
	 		     break;}
				 
		
		
		//Define the reaction rates
		for(m=0;m<Num_CT;m++)
		{
			i=CT[m];
			TimeBirthPrey[m]=ReactionRateToTime(Prey[i][1]*GrowthFunction(Prey[i][0])*Prey_b);
			TimeDeathPrey[m]=ReactionRateToTime(Prey[i][1]*Prey_d);
			for (l=0;l<Num_CT_Predator;l++) //Prey type i, predator type j are the indexs in Prey[][] and Predator[][]
			{
			j=CT_Predator[l];
			//reaction rate of the predation without the reproduction of the predator
			TimePredation_Prey[m][l]=ReactionRateToTime(Prey[i][1]*Predator[j][1]*DefenceFunction(Prey[i][0],Predator[j][0],power,Predation_p,k_max)*(1-Predator[j][0]));
			//Reaction rate of the predation with the reproduction of the predator
			TimePredation_Predator[m][l]=ReactionRateToTime(Prey[i][1]*Predator[j][1]*DefenceFunction(Prey[i][0],Predator[j][0],power,Predation_p,k_max)*Predator[j][0]);
		    }
			//Resource Competition for preys
			for(ResComp_index=0;ResComp_index<Num_CT;ResComp_index++)
			{	
			ResComp_index_i=CT[ResComp_index];
			TimeResCompPrey[m][ResComp_index]=ReactionRateToTime(Prey[i][1]*Prey[ResComp_index_i][1]*Resource_c);
			//printf("TimeResCompPrey[m][ResComp_index]= %lf \n", TimeResCompPrey[m][ResComp_index]);
			}	
		}
		
		for(l=0;l<Num_CT_Predator;l++)
			{ 
				j=CT_Predator[l];
				TimeDeathPredator[l]=ReactionRateToTime(Predator[j][1]*Predator_d);
		    }
		
		//Pick up the shortest time for the reaction
		ShortestTime=TimeBirthPrey[0];
        ReactionID_X=0;
        ReactionID_Y=1000;
        BirthDeathIndex=1;//1 for birth, 0 for death
		CompetitionIndex=0; //1 for competition, 0 for no competition happens
		
		for(m=1;m<Num_CT;m++)//Birth of Prey
			{if(ShortestTime>TimeBirthPrey[m]) {ShortestTime=TimeBirthPrey[m];ReactionID_X=m; ReactionID_Y=1000;BirthDeathIndex=1;}}
		for(m=0;m<Num_CT;m++)//natural death of prey
			{if(ShortestTime>TimeDeathPrey[m]) {ShortestTime=TimeDeathPrey[m];ReactionID_X=m; ReactionID_Y=1000;BirthDeathIndex=0;}}
		for(m=0;m<Num_CT;m++)//predation, prey die
			for(l=0;l<Num_CT_Predator;l++)
		    {if(ShortestTime>TimePredation_Prey[m][l]) {ShortestTime=TimePredation_Prey[m][l];ReactionID_X=m;ReactionID_Y=1000;BirthDeathIndex=0;}}
		for(m=0;m<Num_CT;m++)//predation,prey die, predator reproduce
			for(l=0;l<Num_CT_Predator;l++)
		    {if(ShortestTime>TimePredation_Predator[m][l]) {ShortestTime=TimePredation_Predator[m][l];ReactionID_X=m;ReactionID_Y=l;BirthDeathIndex=0;}}
		for(l=0;l<Num_CT_Predator;l++)//natural death of predator
		   {if(ShortestTime>TimeDeathPredator[l]) {ShortestTime=TimeDeathPredator[l];ReactionID_X=1000;ReactionID_Y=l;}}
				
		//competition
		for(m=0;m<Num_CT;m++)//competition of preys, prey die
			for(ResComp_index=0;ResComp_index<Num_CT;ResComp_index++)
		{if(ShortestTime>TimeResCompPrey[m][ResComp_index]) {ShortestTime=TimeResCompPrey[m][ResComp_index];ReactionID_X=ResComp_index;ReactionID_Y=1000;CompetitionIndex=1;}}
		
		
				
		//change population size according to the reaction with the shortest Time
		Time=ShortestTime+Time;
		//  birth, death, competition of prey species, predation without reproduction
		if(ReactionID_Y==1000)
			{
				i=CT[ReactionID_X];
				if(CompetitionIndex==1)Prey[i][1]--; //death because of competition
				if(BirthDeathIndex==1&&CompetitionIndex==0)//birth
					 { 
						 
						 rnd1=mt_drand();
						 if (rnd1<mu_Prey) 
							 {
								// printf("\n Mutation number %d and now the time step is %d, and the time is %lf.\n",Num_PreyHistory,n,Time);
								 MutatePrey(i);//give the index of the parent type
							 }
						 
						 else Prey[i][1]++; //reproduction of prey
					 }
				if(BirthDeathIndex==0&&CompetitionIndex==0) Prey[i][1]--;  // natural death, or eaten by predator
			 }
			 
		// predation with reproduction	or nautral death of predator 
		 if(ReactionID_Y!=1000)
		 {
			 j=CT_Predator[ReactionID_Y];
			 if(ReactionID_X!=1000) // begin : prey dies, predator reproduces
			 {
				 i=CT[ReactionID_X];
				 Prey[i][1]--; //one prey dies
				 rnd2=mt_drand();
				 if (rnd2<mu_Predator)
					 { MutatePredator(j,k_max);} //one predator mutants
				 else Predator[j][1]++; //one predator reproduces
			 }//end :prey dies, predator reproduces
			 else Predator[j][1]--; //natural death of predator
		 }
	     
		 //every time step one action happens, so there will be only at most one type going extinct in one time step.
		 for(m=0;m<Num_CT;m++)
			 {i=CT[m];if(Prey[i][1]==0){Extinct_Prey(m); break;}}
		 
		 for(m=0;m<Num_CT_Predator;m++)
			 {
				 //test
				 i=CT_Predator[m];
				 if(i==1000) 
				 {
				 printf("the count of Num_CT_Predator is %d, which is wrong\n",Num_CT_Predator); 
				 for(l=0;l<Num_CT_Predator;l++) printf("CT_Predator[%d]=%d\n",m,CT_Predator[m]);
			     }
				 if(Predator[i][1]==0){Extinct_Predator(m); break;}
			 }
		 
			 //break the program when both prey species and predator species go extinct
		 if(Num_CT==0)
			  { 
				if(Notice==0)
				{
				printf("all prey species went extinct\n"); Notice=1;} 
		        if(Num_CT_Predator==0) {printf("all predator species went extinct\n");break;}
			  }
			  
          //count the current prey and predator types, just for the same reason. It is not necessary
         for(m=0;CT[m]!=1000;m++); Num_CT=m;
		 for(m=0;CT_Predator[m]!=1000;m++); Num_CT_Predator=m;
		 
	}//	while ((Num_PreyHistory<100)&&(Num_PredatorHistory<100))
	
	// to renew the Prey[][], Predator[][] every 99 mutations
	for(m=0;m<100;m++)
	   {   
		Prey2[m][0]=0.;
	    Prey2[m][1]=0.;
	  }
	 
	for(m=0;CT[m]!=1000;m++); Num_CT=m;
	
	for(m=0;m<Num_CT;m++)
	{
		i=CT[m];
		Prey2[m][0]=Prey[i][0];
		Prey2[m][1]=Prey[i][1];
	}
	
	for(m=0;m<Num_CT;m++) CT[m]=m;
	for(m=0;m<Num_CT;m++)
	{Prey[m][0]=Prey2[m][0];Prey[m][1]=Prey2[m][1];}
	
	  
  	for(m=0;m<100;m++)
  	   {   
  		Predator2[m][0]=0.;
  	    Predator2[m][1]=0.;
  	  }
	for(m=0;CT_Predator[m]!=1000;m++); Num_CT_Predator=m;
	for(m=0;m<Num_CT_Predator;m++)  
	{
		i=CT_Predator[m];
		Predator2[m][0]=Predator[i][0];
		Predator2[m][1]=Predator[i][1];
	}
	
	for(m=0;m<Num_CT_Predator;m++) CT_Predator[m]=m;
	for(m=0;m<Num_CT_Predator;m++)
		{Predator[m][0]=Predator2[m][0];Predator[m][1]=Predator2[m][1];}
	Num_PreyHistory=Num_CT;
	Num_PredatorHistory=Num_CT_Predator;
	//printf("n=%d\n",n);
    
   }//while (n<10000000)
   //printf("Out of while(n<**)\n");
	return 0;
}