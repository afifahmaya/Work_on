#include<stdio.h>
#include<stdlib.h>
#include<math.h>


void main()
{
    double a,n=10;
    double A[n+1][n+1],x[n+1],b[n+1],AA[n+1][n+2];
    int i,j;

    printf("Solving system of linear equation\n");                                                                               //calculate value of c0, c1 and c2
      for(i=1;i<=n+1;i++){
          for(j=1;j<=n+1;j++)                                                       //build the system of linear equation A
            {
                if(i==j) A[i][j]=2;
                else if(abs(i-j)==1) A[i][j]=-1;
                else A[i][j]=0;
                AA[i][j]=A[i][j];                                                   //AA is going to be processed
            }
        b[i]=1;
        AA[i][n+2]=b[i];
      }
      
      for ( i=1 ; i<=n+1 ; i++ )                                                     //Gauss Elimination process
      {
        if ( AA[i][i]==0 )                                                          //checking if the diagonal is 0
        {
          for ( j=i+1 ; j<=n+2 ; j++ )                                                //if yes, it is swapped with other row that is not zero
          {
            if ( AA[j][i]!=0 ) 
            {
              for ( d=1 ; d<=n+1 ; d++ )                                            //swapping the row
              {
                a = AA[i][d] ;
                AA[i][d] = AA[j][d] ;
                AA[i][d] = a ;
              }
            }
          }
        }
        a = AA[i][i];
        if(AA[i][i]!=0)                                                              //Devide row with diagonal
        {                                                                            //such that the diagonal is 1
          for (j=1;j<=n+2;j++)
          {
            AA[i][j] = AA[i][j]/a ;
          }
          for ( c=j+1; c<=n+2 ; c++)
          {
            a = AA[c][j];
            for ( d=1 ; d<=n+2 ; d++)
            {
              AA[c][d]=AA[c][d]-(a*AA[j][d]);
            }
          }
        }                   //baru sampe sini
      }
      for ( j=3 ; j>=1 ; j--)                                                         //Finding varphi by result of Gauss Elimination over AA
      {
        varphi[i][j]=AA[j][4];
        for(c=j+1 ; c<=3 ; c++)
        {
          varphi[i][j] = varphi[i][j] - (AA[j][c]*varphi[i][c]);
        }
      }
      printf("Using Gauss Elimination with solution : %.3f %.3f %.3f \n\n", varphi[i][1],varphi[i][2],varphi[i][3]); //showing the result
      for(j=1;j<=3;j++)                                                               //calculate local & global matrix A
      {
        locB[j]=0.;                                                                   //initialization of local B
        for(c=1;c<=3;c++)
        {
          a = (varphi[j][2]*varphi[c][2]+varphi[j][3]*varphi[c][3]);                  //Matrix A
          locA[j][c] = area*a;                                                        //local matrix A
          A[elnp[k][j]][elnp[k][c]] = A[elnp[k][j]][elnp[k][c]]+locA[j][c];           //global matrix A
          if(j==c) locB[j] = locB[j]+f*(area/6.);                                     //Matrix B
          else locB[j] = locB[j]+f*(area/12.);                                        //local matrix B
        }
        B[elnp[k][j]] = B[elnp[k][j]]+locB[j];                                        //global matrix B
        locB[j]=0.;
      }
}
