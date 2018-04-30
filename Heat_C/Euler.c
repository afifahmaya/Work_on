#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
    int i, devider = 10;
    double y, yy, h =1./devider;

    FILE *euler;
    FILE *error;

    y = 0.04; //initial condition
    euler = fopen("euler.txt", "w");
    error = fopen("euler_er.txt","w");
    fprintf(euler,"%.3f %.3f \n", 0., y);
    fprintf(error,"%.3f %.3f \n", 0., y);
    for(i=1;i<=devider;i++){
        y = y + h*5*(2-y)*y;
        fprintf(euler,"%.3f %.3f \n", h*i, y);
    }
    fclose(euler);

}