/***************************************************************

METROPOLIS Monte Carlo Simulation for self-avoiding bead-spring
polymer in 3D
Beads can be in one of the three states (s): +1 or 0
+1  : Nucleosome bound
0   : Nucleosome unbound


***************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define Sqr(xx)     ((xx) * (xx))
#define Cube(xx)    ((xx) * (xx) * (xx))
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))

//Random number generator
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}

double Energy(double xo,double yo,double zo,int ss, int i);

/********************************************************************
VARIABLES:
N                       Number of beads
m                       mass
k                       spring constant
x,y,z                   Position of particle
xn,yn,zn                new position of particle
xo,yo,zo                old position of particle
dr                      random displacement factor
E,Eo,En                 current, old, new potential energy
ks                      spring constant
r0                      equilibrium bond length
rg                      radius of gyration

********************************************************************/

double *x,*y,*z;
int *s,*frz;
int *tau;
double alpha0,alpha1;
double sigma,eps,ecut,rc,rc2;
int N;
double Js,Jw;
double ks,kc,mu;
double r0;
double en_lj,en_sp;


int main()
{
char file[64];
double m=1.0;
int i,j,k,id,count=0;
int n,nstep,naccept=0;
int nf;
double xr,yr,zr,rr;
double xn,yn,zn;
double xo,yo,zo;
double cx,cy,cz;
double dr,d,dd;
double E,Eo,En;
int sn,so;
int taun,tauo;
int nrepeat;
double vend;
double dE,dE_b;
double theta,phi;
double T;
double beta;
double r2i,r6i;
double rg;
double rnum;
double R,sum_R=0.0;
long idum;
FILE *ip,*op;

nstep=1000000;
Jw=0.2;
Js=2.0;
mu=0.4;

nrepeat=10;

N=200;
sigma = 1.0;
eps = 1.0;
dr = sigma/2.0;
rc = pow(2.0,1.0/6.0)*sigma;
rc2 = Sqr(rc);
r0 = sigma*1.0;
ks = 100.0;
T = 1.0;
beta  = 1.0/T;
srand(time(NULL));
idum=-955430+(rand()%1000);
r2i=Sqr(sigma)/rc2;
r6i=r2i*r2i*r2i;
ecut=4.0*eps*r6i*(r6i-1.0);

//Memory allocation
AllocMem(x,N,double);
AllocMem(y,N,double);
AllocMem(z,N,double);
AllocMem(s,N,int);
AllocMem(frz,N,int);
AllocMem(tau,N,int);

//Generate initial random configuration 
x[0]=0.0;
y[0]=0.0;
z[0]=0.0;
    for(i=1;i<N;i++)
    {
    d=0.0;
        while(d<sigma)
        {
        xr=2.0*ran2(&idum)-1.0;
        yr=2.0*ran2(&idum)-1.0;
        zr=2.0*ran2(&idum)-1.0;
        dd=sqrt(Sqr(xr)+Sqr(yr)+Sqr(zr));
        xr=xr/dd;
        yr=yr/dd;
        zr=zr/dd;
        x[i]=x[i-1]+xr*r0;
        y[i]=y[i-1]+yr*r0;
        z[i]=z[i-1]+zr*r0;
            for(j=0;j<i;j++)
            {
                if(j==0)
                {
                xr=x[i]-x[j];
                yr=y[i]-y[j];
                zr=z[i]-z[j];
                d=Sqr(xr)+Sqr(yr)+Sqr(zr);
                }
                else
                {
                xr=x[i]-x[j];
                yr=y[i]-y[j];
                zr=z[i]-z[j];
                dd=Sqr(xr)+Sqr(yr)+Sqr(zr);
                }
                if(dd<d)d=dd;
            }
        }
    }

//Generate initial random modifications s[i]

    for(i=0;i<N;i++)
    {
        if(i%nrepeat==0)
        tau[i]=1;
        else 
        tau[i]=0;
    s[i]=0;
    frz[i]=0;
    }

    //Monte Carlo steps
    for(n=0;n<nstep;n++)
    {
        //Dupmp the positions
        if(n%10000==0)
        {
        sprintf(file,"nr%d_Js_%de-1_Jw_%de-1_mu_%de-1.xyz",nrepeat,(int)(Js*10),(int)(Jw*10.0),(int)(mu*10.0));
        op=fopen(file,"a");
        fprintf(op,"%d\n%d\n",N,n);
            for(i=0;i<N;i++)
            {
                fprintf(op,"%d %lf %lf %lf\n",s[i],x[i],y[i],z[i]);
            }
        fclose(op);
        }

        for(i=0;i<N;i++)
        {
        //Move 1: Change position of the bead
        id=(int)(ran2(&idum)*N)%N;
        xo=x[id];
        yo=y[id];
        zo=z[id];
        Eo=Energy(xo,yo,zo,s[id],id);
        xr=(2.0*ran2(&idum)-1.0);
        yr=(2.0*ran2(&idum)-1.0);
        zr=(2.0*ran2(&idum)-1.0);
        xn=x[id]+xr*dr;
        yn=y[id]+yr*dr;
        zn=z[id]+zr*dr;
        En=Energy(xn,yn,zn,s[id],id);
        dE=En-Eo;
        dE_b=dE*beta;
            if(dE<=0.0)
            {
            E=E+dE;
            x[id]=xn;
            y[id]=yn;
            z[id]=zn;
            }
            else if(ran2(&idum)<exp(-dE_b))
            {
            E=E+dE;
            x[id]=xn;
            y[id]=yn;
            z[id]=zn;
            }

        //Move 2: flip the state (s) 
	
        id=(int)(ran2(&idum)*N)%N;
        so=s[id];
        Eo=Energy(x[id],y[id],z[id],so,id);
            if(s[id]==1)
            {
            sn=0;
            }
            else
            {
            sn=1;
            }
        En=Energy(x[id],y[id],z[id],sn,id);
        dE=En-Eo;
        dE_b=dE*beta;
            if(dE<=0.0)
            {
            E=E+dE;
            s[id]=sn;
            }
            else if(ran2(&idum)<exp(-dE_b))
            {
            E=E+dE;
            s[id]=sn;
            }
	
        }
        if(n%100000==0)
        printf("Step = %d \n",n);
    }
   
return 0;
}


double Energy(double xo,double yo,double zo,int ss,int i)
{
double en=0.0;
double xr,yr,zr,r2,r2i,r6i;
int j,k;
double tauij;
en_sp=0.0;
en_lj=0.0;

    for(j=0;j<N;j++)
    {
        if(j!=i)
        {
        xr=xo-x[j];
        yr=yo-y[j];
        zr=zo-z[j];
        r2=Sqr(xr)+Sqr(yr)+Sqr(zr);
            //LJ potential 
            if(r2<rc2)
            {
            r2i=Sqr(sigma)/r2;
            r6i=r2i*r2i*r2i;
            en=en+4.0*eps*r6i*(r6i-1.0)-ecut;
            en_lj=en_lj+4.0*eps*r6i*(r6i-1.0)-ecut;
            }
            //Spring potential
            if(j==i-1)
            {
            en=en+0.5*ks*Sqr((sqrt(r2)-r0));
            en_sp=en_sp+0.5*ks*Sqr((sqrt(r2)-r0));
            }
            if(j==i+1)
            {
            en=en+0.5*ks*Sqr((sqrt(r2)-r0));
            en_sp=en_sp+0.5*ks*Sqr((sqrt(r2)-r0));
            }
            tauij=tau[i]*tau[j];
            //Modification potential
            if(r2<Sqr(1.5*sigma))
            {
            en=en-Js*tauij*((double)(ss*s[j])) -Jw*(1.0-tauij)*((double)(ss*s[j]));
            }
        }
    }
en=en+mu*(double)ss;
    
return en;

}


