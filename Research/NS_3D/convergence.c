#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
    double t,dt,x,y,z,h,u1,u2,u3;
    int Nt,Nh;

    dt=0.01;
    h=0.1;
    Nt=(int)1./dt;
    Nh=(int)1./h;

    printf("Nt = %d", Nt);

    FILE *plot;
    plot = fopen("trial.txt","w+");

    for(int i=1; i<=Nt; i++){
        t = i*dt;
        fprintf(plot,"\"t=%f \"\n",t);
        for(int j=1; j<=Nh; j++){
            x = j*h;
            for(int k=1; k<=Nh; k++){
                y = k*h;
                for(int l=1; l<=Nh; l++){
                    z = l*h;
                    u1 = -cos(x)*sin(y)*cos(z)*exp(-2*t);
                    u2 = sin(x)*cos(y)*cos(z)*exp(-2*t);
                    u3 = 0.;
                    fprintf(plot,"%f \t %f \t %f \n",u1,u2,u3);
                }
            }
        }
        fprintf(plot,"\n\n");
    }
}