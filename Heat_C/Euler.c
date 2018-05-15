#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
    int i;
    double y, yy, h, t, er, k1, k2, k3, k4;

    FILE *euler;
    FILE *rk;
    char filename1[30];
    char filename2[30];

    for(int dev=4;dev<=128;dev=dev*2)
    {
        h = 1./dev;
        
        sprintf(filename1, "euler_%.3d.txt",dev);
        sprintf(filename2, "rk_%.3d.txt",dev);
        
        euler = fopen(filename1,"w");
        rk = fopen(filename2, "w");

        y = 0.04; //initial condition    
        
        fprintf(euler,"%.3f %.3f %.3f %.3f \n", 0., y, y, 0.);
        fprintf(rk,"%.3f %.3f %.3f %.3f \n", 0., y, y, 0.);

        //Using Euler
        for(i=1;i<=dev;i++)
        {
            t=i*h;
            y = 5*h*(2-y)*y + y;
            yy = 2./(1+49*exp(-10*t));
            er = fabs(yy-y);
            fprintf(euler,"%.3f %.3f %.3f %.3f \n", t, y, yy, er);
        }

        //Using RK
        y = 0.04; //initial condition
        for(i=1;i<=dev;i++)
        {
            t=i*h;
            k1 = 5*(2-y)*y;
            k2 = 5*(2-(y+(k1*h/2.)))*(y+(k1*h/2.));
            k3 = 5*(2-(y+(k2*h/2.)))*(y+(k2*h/2.));
            k4 = 5*(2-(y+(k3*h)))*(y+(k3*h));
            y = h*(k1 + 2*k2 + 2*k3 + k4)/6. + y;
            yy = 2./(1+49*exp(-10*t));
            er = fabs(yy-y);
            fprintf(rk,"%.3f %.3f %.3f %.3f \n", t, y, yy, er);
        }

        //plot
        FILE *pipe = popen("gnuplot", "w");
        fprintf(pipe, "reset \n");
        fprintf(pipe, "set terminal png \n");
        fprintf(pipe, "set output 'plot%.3d.png'\n", dev);
        fprintf(pipe, "set title 'Euler and Runge-Kutta for devider %.3d' \n", dev);
        fprintf(pipe, "set xrange [0:1] \n");
        fprintf(pipe, "set xlabel 't' \n");
        fprintf(pipe, "set yrange [0:4] \n");
        fprintf(pipe, "set ylabel 'y' \n");
        fprintf(pipe, "plot '%s' using 1:2 w linespoint title 'Euler y numeric', '%s' using 1:2 w linespoint title 'RK y numeric', '%s' using 1:3 w linespoint title 'Exact y' , '%s' using 1:4 w linespoint title 'Euler differences' , '%s' using 1:4 w linespoint title 'RK differences'\n", filename1, filename2, filename1, filename1, filename2);
        fclose(pipe);
    }

    fclose(euler);
    fclose(rk);
}