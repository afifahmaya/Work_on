#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{   
    double T=1., dt = 0.01, t, h=0.1;
    int Nt = (int)(T/Nt);
    int Nh = (int)(1./h);
    double x,y,z,u,v,w;

    FILE *plot;

    plot = fopen("convergence.txt","w");

    for(int i=1;i<=Nt;i++){
        t= i*dt;
        fprintf(plot," \" t = %.3f \" \n",t);
        for(int j=1;j<=Nh;j++){
            x = j*h;
            for(int k=1;k<=Nh;k++){
                y = k*h;
                for(int l=1;l<=Nh;l++){
                    z = k*h;
                    u= -cos(x)*sin(y)*cos(z)*exp(-2*t);
                    v= sin(x)*cos(y)*cos(z)*exp(-2*t);
                    w=0;
                    fprintf(plot, "%.3f %.3f %.3f \n",u,v,w);
                }
            }
        }
        fprintf(plot,"\n");
    }
    fclose(plot);

    return(0);
}