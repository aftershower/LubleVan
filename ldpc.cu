
//LDPC decoder
//#define N 19942
//#define K 10003
//#define L 9939

#define N 8176
#define K 7156
#define L 1020

//#define N 9968
//#define K 4984
//#define L 4984

#define TOTAL_LENGTH 100
#define Eb_No 0.5
//#define FILENAME "PCM.5000.10000"

/*
#define N 2048
#define K 1024
#define L 1024
*/

#define ITERATION 200
#define MAX_Le 10.0
#define MIN_Le 0.0
#define _HH	printf("Have been Here!\n")
#define CLIP 12.0
#define min(x,y) ((x<=y) ? x:y)
#define max(x,y) ((x>=y) ? x:y)
#define sign(x) ((x <= 0) ? (-1.0):(1.0))
#define absd(x) ((x <= 0)? (-1.0*x):(x))
#define absi(x) ((x <= 0)? (-1*x):(x))

#include<stdio.h>
#include<math.h>
#include<memory.h>
#include<stdlib.h>


/****************************************************************************/
/*  RANDOM.H                                                                */
/*  Includes functions which simulate 
/*                                                                          */
/*    float ran2(long *idum) - generates random numbers uniformily         */
/*                              distributed over (0,1).  Random sequence is */
/*                              initialized by calling ran2 with *idum a    */
/*                              negative number.                            */
/*                                                                          */
/*    float gasdev(long *idum) - generates Gaussian random numbers, N(0,1) */
/*                              utilizing the function ran2.  Initialized   */
/*                              the same way as ran2.                       */
/*                                                                          */
/*    int irbit(unsigned long *iseed) - generates a pseudo-random noise     */
/*					See files "pn_gen#.h"		 */
/*    float poidev(float xm, long *idum) - generates an integer value     */
/*				that is a random deviate drawn from a 	    */
/* 				Poisson distribution of mean xm, using      */
/*				ran2(idum) as a source of the uniform dev.  */
/*									    */
/*    float gammln(float xx) - returns the value of the log of the gamma  */
/*                           	function evaluated at xx for xx>0           */
/*				(used for poidev()).			    */
/*    									    */
/*    float expdev(long *idum) - returns an exponentially distributed,     */
/* 				positive, random deviate of unit mean.      */
/*				To change mean to lambda, multiply the      */
/*				result of the call by lambda.		    */
/*                                                                          */
/*   Note:  requires that math.h be included !!!!!!!!!!!!!!!!!!!            */
/****************************************************************************/

/****************************************************************************/
/* function:  float ran2(long *idum)                                       */
/****************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IM1     2147483563
#define IM2     2147483399
#define AM      (1.0/IM1)
#define IMM1    (IM1-1)
#define IA1     40014
#define IA2     40692
#define IQ1     53668
#define IQ2     52744
#define IR1     12211
#define IR2     3791
#define NTAB    32
#define NDIV    (1+IMM1/NTAB)
#define EPS     1.2e-7
#define RNMX    (1.0-EPS)
#define PI	3.1415927

float ran2(long *idum)                 /////////////////产生均匀分布的随机数。。
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0)
	{
		if (-(*idum) < 1)
			*idum = 1;
		else
			*idum = -(*idum);
		idum2 = (*idum);
		for (j=NTAB+7 ; j>=0 ; j--)
		{
			k = (*idum)/IQ1;
			*idum = IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0)
				*idum +=IM1;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ2;                
	*idum = IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0)
		*idum += IM1;
	k = idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0)
		idum2 += IM2;
	j = iy/NDIV;
	iy = iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1)
		iy += IMM1;
	if ( (temp=AM*iy) > RNMX)
		return RNMX;
	else
		return temp;
}


/****************************************************************************/
/* function:  float gasdev(long *idum)                                     */
/****************************************************************************/

float gasdev(long *idum)               ///////////高斯噪声 
{
	float ran2(long *idum);
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;

	if (iset == 0)
	{
		do {
			v1 = 2.0*ran2(idum)-1.0;
			v2 = 2.0*ran2(idum)-1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else {
		iset = 0;
		return gset;
	}
}

/*****************************************************************************/
/* function:  float poidev(float xm, long *idum)                           */
/*****************************************************************************/

float poidev(float xm, long *idum)            ///////////////柏松分布。。。

{
	float gammln(float xx);
	float ran2(long *idum);
	static float sq, alxm, g, oldm=(-1.0);
	float em, t, y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran2(idum);
		} while (t>g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran2(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2(idum) > t);
	}
	return em;
}


/*****************************************************************************/
/* function:  float gammln(float xx)                                       */
/*****************************************************************************/

float gammln(float xx)                 			//////////计算对数值？？？？？？？？？           
{
	float x,y,tmp,ser;
	static float cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0; j<=5; j++) 
		ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

/*****************************************************************************
 * function:  float expdev(long *idum)					     *
 *****************************************************************************/
 
 float expdev(long *idum)                   //////////////指数分布？/////////
 {
 	float dum;
 	
 	do
 		dum=ran2(idum);
 	while (dum==0.0);
 	return -log(dum);
 }


typedef struct RECORD1
{
int no_elem;
int *pos;
int *pos2;
int count;
} record1;

record1 row[L],col[N];
float stdd,probzero,nvi;



int parity_check(unsigned short c[])        ////////////////////检验
//unsigned short c[];
{
int i,j,k,sum;
for(i=0;i<=L-1;i++)
	{
	sum = 0;
	for (j=0;j<=row[i].no_elem-1;j++)
		{
		sum = (sum + c[row[i].pos[j]])%2;
		}
	if (sum != 0)
		return(0);
	}
return(1);
}
		
float boxminus(float x,float y)
{
float z;

z = log((exp(x)-exp(y))/(1-exp(x+y)));                 ///////////////某个函数。。

return(z);
}

float psifunc(float x)
{
return(log(tanh(absd(x/2))));             /////////////双曲正切函数的对数。。。（对数似然用？）
}

		
float boxplus( float x,float  y)
{
float z;

z = log((exp(x)+exp(y))/(1+exp(x+y)));           ///////////某个函数  类似于boxminus（x，y）。。。。


return(z);
}

void main()                       //////////////main函数
{
int i,j,k,l,np,m,n,p,q,simpt;
unsigned short d[K],c[N],dec[N];
int count,N_ITER,iter,itemp[100];               ////////////////////n=1000
int cf,pkt_err,detner,undetner,nber;
float detber,undetber,totber,fer;
long seed;
float ccia[]={Eb_No};
float cci,cdr,const1,**Lc,**Le,Lch[N],Lf[N],temp,r[N],tmp;
float tmpmag1,tmpsign1,tmpmag2,tmpsign2,tmparr[N];
FILE *Gmat,*Hmat,*fpres,*fptemp1,*fptemp2;
char tstr[80],tc='0';
int checkH[L];
FILE *fp;
int nozero;
char FILENAME[100];

for (i=0;i<L;i++) checkH[i]=0;	

         /////////////////fopen

sprintf(FILENAME,"PCM.%d.%d", K,N);           //////////////sprintf返回值：字符串长度
Hmat = fopen("D:\\Code\\C\\asdfsadfg\\asdfsadfg\\PCM.8176.7156","r");

fgets(tstr,80,Hmat);                                //////////////fgets从流读取n-1字符除非读完行参数s来接收字符串成功则返回tstr指针否则返回NULL
for(i=0;i<=N-1;i++)
        {
        fscanf(Hmat,"<%d> ",&l);                                //////////////fscanf返回值：返回实际被转换并赋值的输入项的数目。返回成功读入的个数
        col[i].no_elem = l;
        col[i].pos = (int *) calloc(col[i].no_elem,sizeof(int)); /////calloc   分配l个sizeof单元  并返回地址
        for (j=0;j<=col[i].no_elem-1;j++)
                {
                fscanf(Hmat,"%d  ",&l); 
                col[i].pos[j] = l;
//		fprintf(fp,"[%d] ",col[i].pos[j]);
				checkH[l]=1;
                (row[l].no_elem)++;
                }
        fscanf(Hmat," \n",tc);
        }

nozero=1;
for (i=0;i<L;i++) 
{ if (checkH[i]==0)
	{ printf("i=%d \n", i);
		nozero=0;
		}	
}
if (nozero) printf("nozero!!!!\n");
printf("Successfully read H matrix \n");
fflush(stdout);           //////////////////////fflush(stdout)刷新标准输出缓冲区把输出缓冲区里东西打印标准输出设备上
/*///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////////////////*/
for(i=0;i<=L-1;i++)
        {
        row[i].pos = (int *) calloc(row[i].no_elem,sizeof(int));
        row[i].pos2 = (int *) calloc(row[i].no_elem,sizeof(int));
        row[i].count = 0;
        }

for(i=0;i<=N-1;i++)
        {
        for(j=0;j<=col[i].no_elem-1;j++)
                {
                k = col[i].pos[j];
                row[k].pos[row[k].count] = i;
                row[k].pos2[row[k].count] = j;
                /*if (row[k].count >= row[k].no_elem)
                        printf("something wrong row index exceeds \n");*/
                (row[k].count)++;
                }
	}
	/*    可能是在产生矩阵///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////*/
printf("successfully formed row and col positions \n");

Le = (float **)calloc(N,sizeof(float *));
Lc = (float **)calloc(N,sizeof(float *));
for(i=0;i<=N-1;i++)
	{
	Le[i] = (float *)calloc(col[i].no_elem,sizeof(float));
	Lc[i] = (float *)calloc(col[i].no_elem,sizeof(float));                                          ///////////////？？？
	if ((Le[i] == NULL)||(Lc[i] == NULL))
		 printf("problem in allocated memory for %d \n",i);
		 system("pause");
	}


seed = -108; N_ITER = ITERATION;

for(simpt=0;simpt<=0;simpt++)
	{
	cci = ccia[simpt]-10.0*log10((float)(N)/(float)(K));
	stdd = sqrt(0.5/(pow(10.0,(float)cci/10.0))); 
	
	printf("std = %f \n",stdd);
	nvi = 1.0/(2.0*stdd*stdd);const1 = -2.0/(stdd*stdd);               //////////////const1
	probzero = exp(-1.0*nvi);
	

	pkt_err = 0;undetner=0;detner=0;
	
	for(np=1;((np<=TOTAL_LENGTH) && (pkt_err <= 50));np++)
		{
		for(i=0;i<=K-1;i++)
			{
			/*d[i] = (unsigned short)(ran2(&seed)*2.0);*/
			d[i] = 0;
			}
		
		/*encode(d,c,G);*/

		for(i=0;i<=N-1;i++)
			c[i]=0;
	

		if (parity_check(c) != 1) 
			printf("Something wrong parity check is not satisfied \n"); 
		else
			printf("Parity check satisfied \n");
		
		for (i=0;i<=N-1;i++)
			{
			r[i] = (float)(2*c[i]-1)+gasdev(&seed)*stdd;            ///////////////加噪声
			}
	
		for (i=0;i<=N-1;i++)
			{
			Lch[i] = const1*r[i];
			for(j=0;j<=col[i].no_elem-1;j++)
				{
				Lc[i][j] = Lch[i];
				}
			}

		cf = 0;
		printf("Eb/No = %f packet %d iter = ",ccia[simpt],np);
		_HH;
		for (iter=1;((iter<=N_ITER) && (cf != 1));iter++)              //////////////迭代译码？
			{
//			printf("[%d] ",iter);
//_HH;
			for (j=0;j<=L-1;j++)															/////for1main
				{
				
				for(l=0;l<=row[j].no_elem-1;l++)									//for1of1					
					{
					q = row[j].pos2[l];
					p = row[j].pos[l];
					tmparr[l] = psifunc(Lc[p][q]);              /////////////////双曲正切函数的对数
					}
				tmpmag1 = tmparr[0]; 

				tmpsign1 = sign(Lc[row[j].pos[0]][row[j].pos2[0]]);
				
				for(l=1;l<=row[j].no_elem-1;l++)								///for2of1
                                	{
                                        p = row[j].pos[l];
					q = row[j].pos2[l];
					tmpmag1 = tmpmag1 + tmparr[l];
					tmpsign1 = tmpsign1*sign(Lc[p][q]);
					if (tmpsign1<-CLIP) tmpsign1=-CLIP;
					if (tmpsign1>CLIP) tmpsign1=CLIP;
                                        }
			//	_HH;
				for(i=0;i<=row[j].no_elem-1;i++)							/////for3of1
					{
                                        p = row[j].pos[i];
					q = row[j].pos2[i];
					tmpmag2 = tmpmag1 - tmparr[i];
					tmpsign2 = tmpsign1*sign(Lc[p][q]);

					Le[p][q] = tmpsign2*(-psifunc(tmpmag2));
					if (Le[p][q]<-CLIP) Le[p][q]=-CLIP;
					if (Le[p][q]>CLIP) Le[p][q]=CLIP;
					}
//					_HH;
				}				
				
			for(i=0;i<=N-1;i++)                          ////for cycle1
				{
				Lf[i] = Lch[i];
				for(j=0;j<=col[i].no_elem-1;j++)
                                        {
					Lf[i] += Le[i][j];										////checkNode1  lf=N   Le=N*no_elem
					}
				for(j=0;j<=col[i].no_elem-1;j++)
                                        {
					Lc[i][j] = Lf[i] - Le[i][j];											////checkNode2
					}
				}

			for (i=0;i<=N-1;i++)                   //////////////for2
				{
				if (Lf[i] > 0.0) dec[i] = 0;
					else dec[i] = 1;
				Lf[i] = Lf[i]-Lch[i];
				}
//		printf("Lf[i]=%e\n",Lf[229]);
                nber = 0;
            for (i=L;i<=N-1;i++) 						/////////////for3
                        {
                        if (dec[i] != d[i-L]) nber++;
                        }

//		printf("After %d iterarions there are %d errors \n",iter,nber);
				
			cf = parity_check(dec);
			if (cf == 1) printf("\n parity check satisfied during %d iter can stop iterations \n",iter);
			}/* for iter = */

		if (cf != 1) pkt_err++;
		nber = 0;
		for (i=L;i<=N-1;i++)
			{
			if (dec[i] != d[i-L]) nber++;
			}
		if (cf != 1) detner = detner + nber;
		else
		undetner = undetner + nber; 

		undetber = (float)undetner/(float)((np)*K);
		detber = (float)detner/(float)((np)*K);
		totber = (float)undetber + detber;
		fer = (float)pkt_err/(float)(np);

		printf("Eb/No = %f K = %d N = %d L = %d \n",ccia[simpt],K,N,L);
		printf("Detected BER = %e Undetected BER = %e Total BER = %e ",detber,undetber,totber);
		printf("FER = %e \n \n",fer);
		}/* for np = */

		undetber = (float)undetner/(float)((np-1)*K);
		detber = (float)detner/(float)((np-1)*K);
		totber = (float)undetber + detber;
		fer = (float)pkt_err/(float)(np-1);

		fpres =	fopen("RESULTS_91","a");
		fprintf(fpres,"Eb/No = %f K = %d N = %d L = %d \n",ccia[simpt],K,N,L);
		fprintf(fpres,"Detected BER = %e Undetected BER = %e Total BER = %e ",detber,undetber,totber);
		fprintf(fpres,"FER = %e \n \n",fer);
		fclose(fpres);

	} /* for simpt */ 
} /* for main */
