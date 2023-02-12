#include<stdio.h>
#include<stdlib.h>
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer
#include<math.h>

int main()
{
int i,j,k;
int ii,jj,kk,ll;
int n;
int N;
N=200;
int s[N];
int t,nn,clust_id[N];
int clust_size[N];
int size_dist[N];
int n_ens[5000],ncount;
int count;
int tmax,tmin;
double rgmax=0;
int nframes;
int n_sys;
int Jw,Js,mu;
int nrepeat;
int flag;
int sum;
int bin;
int s0;
double pp;
double x[N],y[N],z[N];
double xr,yr,zr,rr;
double cx,cy,cz;
double rg;
double sigma,rc,rc2;
FILE *ip,*op;
char file[100];
sigma=1.0;
rc=1.5*sigma;
rc2=rc*rc;
n_sys=40;
mu=4;
Jw=2;
Js=2;
mu=4;
nrepeat=10;
nframes=1000;
count=0;
sprintf(file,"nr%d_Js_%de-1_Jw_%de-1_mu_%de-1.xyz",nrepeat,Js*10,Jw,mu);
ip=fopen(file,"r");
ncount=0;
count=0;
rg=0.0;
    while(ncount<(N+2)*nframes)
    {
    fscanf(ip,"%d\n%d\n",&nn,&t);
    ncount=ncount+N+2;
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
        clust_id[i]=i;
        }
    rr=0.0;
    cx=0.0;
    cy=0.0;
    cz=0.0;
        for(i=0;i<N;i++)
        {
        cx=cx+x[i];
        cy=cy+y[i];
        cz=cz+z[i];
        }
    cx=cx/(double)N;
    cy=cy/(double)N;
    cz=cz/(double)N;
        for(i=0;i<N;i++)
        {
	    xr=x[i]-cx;
	    yr=y[i]-cy;
	    zr=z[i]-cz;
	    rr=rr+xr*xr+yr*yr+zr*zr;
        }
    rr=rr/(double)N;
    rr=sqrt(rr);
        if(ncount>(N+2)*100)
        {
        rg+=rr;
        count++;
        }
    }
fclose(ip);
sprintf(file,"%d_perc_rg.txt",100/nrepeat);
op=fopen(file,"a");
fclose(op);
return 0;
}
