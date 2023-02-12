#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main()
{
int i,j,k;
int ii,jj,kk,nn;
int N;
int t,n,count=0;
double x[1000],y[1000],z[1000];
double ss[400],s[400];
int nt;
double sl,sr;
int a1,J,mu;
int ecount=0;
double dx;
double xr,yr,zr,rr;
double alpha_0,alpha_1;
double sum_s,sum_t;
char file[90],element[20];
FILE *ip,*op;
N=200;
  J=2;
  mu=4;

sprintf(file,"rho_avg_J_%de-1_mu_%de-1.dat",J,mu);
op=fopen(file,"w");
  for(i=0;i<400;i++)
    ss[i]=0;

  count=0;
  ecount=0;
  sprintf(file,"posJs_2_Jw_%de-1_mu_%de-1.xyz",J,mu);
  ip=fopen(file,"r");
    while(count!=(N+2)*1000)
    {
    count=count+N+2;
    fscanf(ip,"%d\n%d\n",&n,&t);
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&j,&x[i],&y[i],&z[i]);
            if(count>(N+2)*50)
            {
            ss[i]=ss[i]+(double)j;
            }
        }
        if(count>(N+2)*100)
        ecount++;
    }
  fclose(ip);
  
  for(i=0;i<N;i++)
  {
  ss[i]=ss[i]/(double)ecount;
  }

  for(i=0;i<N;i++)
  fprintf(op,"%d %lf\n",i,ss[i]);
fclose(op);

return 0;
}
