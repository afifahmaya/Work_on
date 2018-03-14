#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"./functions.c"

double norm_L2(int **elnp, double *area , double *v, int ne)
{
  double norm, a1, result=0.;
  int i,j,k;
  for(k=1;k<=ne;k++)
  {
    norm =0.;
    //printf("=====element %d(%d,%d,%d)=====\n",k,elnp[k][1],elnp[k][2],elnp[k][3]);
    for(i=1;i<=3;i++)
    {
      a1=0.;
      for(j=1;j<=3;j++)
      {
        if(i==j) 
        {
          //printf("a1 = %.5f ",a1);
          a1 = a1 + 2*v[elnp[k][j]];
          //printf("+ 2*%.5f(%d) = %.5f\n",v[elnp[k][j]],elnp[k][j],a1);
        }
        else 
        {
          //printf("a1 = %.5f ",a1);
          a1 = a1 + v[elnp[k][j]];
          //printf("+ %.5f = %.5f\n",v[elnp[k][j]],a1);
        }
      }
      //printf("norm = %.5f ",norm);
      norm = norm + a1*v[elnp[k][i]];
      //printf("+ %.5f*%.5f = %.5f\n",a1,v[elnp[k][i]],norm);
    }
    result = result + area[k]*norm/12.;
    //printf("result = %f\n", result);
  }
  result = sqrt(result);
  printf("\n norm L2 = %.5f \n", result);
  return result;
}

double norm_H(int **elnp, double *area , double **c, double *v, int ne)
{
  double norm =0., v1, v2, result;
  int i,j,k;
  for(k=1;k<=ne;k++)
  {
    v1=0.;v2=0.;
    //printf("=====element %d(%d,%d,%d)=====\n",k,elnp[k][1],elnp[k][2],elnp[k][3]);
    //printf("c = %.3f %.3f %.3f %.3f %.3f %.3f \n", c[k][1], c[k][2], c[k][3], c[k][4], c[k][5], c[k][6]);
    for(j=1;j<=3;j++)
    {
      v1 = v1 + v[elnp[k][j]]*c[k][2*j-1];
      v2 = v2 + v[elnp[k][j]]*c[k][2*j];
      //printf("v1 += v[%d]*c[%d] = %f*%f = %.5f \n", elnp[k][j],2*j-1,v[elnp[k][j]],c[k][2*j-1],v1);
      //printf("v2 += v[%d]*c[%d] = %f*%f = %.5f \n", elnp[k][j],2*j,v[elnp[k][j]],c[k][2*j],v2);
    }
    //printf("norm = %.5f + (%.5f^2+%.5f^2)*%.5f = ", norm, v1,v2,area[k]);
    norm = norm + (v1*v1 + v2*v2)*area[k];
    //printf("%.5f\n", norm);
  }
  result = sqrt(norm);
  printf("\n norm H = %.5f \n", result);
  return result;
}

int main(int argc, char *argv[])
{
  double **npxy,**varphi,**A,**locA,*locB,*B,**AA, *f,*exact, 
        **nabla, a ,*area, sumarea, a_1,b_1,c_1,d_1, xmax, max,
        **M_1, *M_2, *MB, nu, **MA, t, T, l, L2_norm_CG, H_norm_CG,
        maxL2, maxH;                       //DECLARATION
  int DIM=2,np,ne,nb,**elnp,**fanp,dummy,i,j,k,c,d,jmax,n;                            //declaration of integer variable used
  double **aa ;                                                                       //declaration used for showing result
  double *r, *p, *Ap, *x;                                                             //Declaration for CG Method
  double alpha, beta, rr, rnew, pAp;
  double *v;

  FILE *fp;
  FILE *sol_CG;
  FILE *exactt;
  FILE *norm;
  char filename[30];
  
  if((fp=fopen(argv[1],"r"))==NULL)                                                   //CHECK if file can be opened
  {
    printf("Can't open file: npxy.dat.\n"); 
    exit(0);
  }

  fscanf(fp,"%d %d %d",&np,&ne,&nb);                                                  //read first line : nodal point, element, & boundary

  npxy = dmatrix(1,np,1,DIM);                                                         //ALLOC
  elnp = imatrix(1,ne,1,DIM+1);
  fanp = imatrix(1,nb,1,DIM);
  varphi = dmatrix(1,3,1,3);
  A = dmatrix(1,np,1,np);
  AA = dmatrix(1,3,1,3+1);
  B = dvector(1,np);
  M_1 = dmatrix(1,np,1,np);
  M_2 = dvector(1,np);
  MA = dmatrix(1,np,1,np);
  MB =  dvector(1,np);
  x = dvector(1,np);
  locA = dmatrix(1,3,1,3);
  locB = dvector(1,3);
  f = dvector(1,np);
  exact = dvector(1,np);
  area = dvector(1,ne);
  
  aa = dmatrix(1,3,1,3);                                                              //allocation used for showing result

  r = dvector(1,np);                                                                  //allocation used for CG Method
  p = dvector(1,np);
  Ap = dvector(1,np);

  v = dvector(1,np);                                                                  //allocation used for calculate norm
  nabla = dmatrix(1,ne,1,6);

  printf("np = %d, ne = %d, nb = %d\n\n", np, ne, nb);                                //show number of nodal points, elements, & nodes in boundary

  for(i=1;i<=(np+ne+nb);i++)                                                          //READ from .msh file
  {
    if (i<=np) fscanf(fp,"%lf %lf %d",&npxy[i][1],&npxy[i][2],&dummy);                //read the coordinate x1,x2
    else if (i<=(np+ne)) 
      fscanf(fp,"%d %d %d %d",&elnp[i-np][1],&elnp[i-np][2],&elnp[i-np][3],&dummy);   //read element of each triangle
    else fscanf(fp,"%d %d %d",&fanp[i-(np+ne)][1],&fanp[i-(np+ne)][2],&dummy);        //read the boundary point
  }

  T=1.;                                                                               //max time
  t=(4./nb)*(4./nb)*32;                                                                             //delta t
  printf("delta t = %.3f \n",t);
  nu=1.;                                                                          //nu is positive constant
  for(k=1;k<=np;k++) x[k]= 2*sin(M_PI*npxy[k][1])*sin(M_PI*npxy[k][2]);                                                         //initial u for domain at t=0

  for(k=1;k<=ne;k++)
  {
    a_1 = npxy[elnp[k][2]][1]-npxy[elnp[k][1]][1];                                    //MEASURE of each triangle         
    b_1 = npxy[elnp[k][2]][2]-npxy[elnp[k][1]][2];
    c_1 = npxy[elnp[k][3]][1]-npxy[elnp[k][1]][1];
    d_1 = npxy[elnp[k][3]][2]-npxy[elnp[k][1]][2];
    area[k] = 0.5*(a_1*d_1-b_1*c_1);
    sumarea += area[k];
    //printf("area >> %.5f \n", area[k]);
  }

  for(k=1;k<=ne;k++)                                                                  //FIND c0,c1,c2 using Gauss Elimination
  {
    //printf("----------TRIANGLE %d---------\n",k);
    for(i=1;i<=3;i++)                                                                 //varphi_i (x) = c0 + (c1 x1) + (c2 x2)
    {
      //printf("Solving system of linear equation %d \n", i);
      for(j=1;j<=3;j++)                                                               //build the system of linear equation AA
      {
        AA[j][1]=1.;
        AA[j][2]=npxy[elnp[k][j]][1];
        AA[j][3]=npxy[elnp[k][j]][2];
        aa[j][1]=AA[j][1]; aa[j][2]=AA[j][2]; aa[j][3]=AA[j][3];                      //variable used to show the result
      }
      for(j=1;j<=3;j++)                                                               //put the right hand side as AA[][4]
      {
        if(i==j) AA[j][4]=1.;
        else AA[j][4]=0;
        aa[j][4]=AA[j][4];                                                            //variable used to show the result
        //printf(" %.3f %.3f %.3f %.3f \n", aa[j][1],aa[j][2],aa[j][3],aa[j][4]);     //show the system to be solve to find c0,c1,c2
      }
      for ( j=1 ; j<=3 ; j++ )                                                        //Gauss Elimination process
      {
        if ( AA[j][j]==0 )                                                            //checking if the diagonal is 0
        {
          for ( c=j+1 ; c<=3 ; c++ )                                                  //if yes, it is swapped with other row that is not zero
          {
            if ( AA[c][j]!=0 ) 
            {
              for ( d=1 ; d<=3+1 ; d++ )                                              //swapping the row
              {
                a = AA[j][d] ;
                AA[j][d] = AA[c][d] ;
                AA[c][d] = a ;
              }
            }
          }
        }
        a = AA[j][j];
        if(AA[j][j]!=0)                                                              //Devide row with diagonal
        {                                                                            //such that the diagonal is 1
          for ( c=1 ; c<=3+1 ; c++ )
          {
            AA[j][c] = AA[j][c]/a ;
          }
        }
        if(j!=3)                                                                      //Elimination process
        {                                                                             //such that every entries below diagonal is 0
          for ( c=j+1; c<=3 ; c++)
          {
            a = AA[c][j];
            for ( d=1 ; d<=3+1 ; d++)
            {
              AA[c][d]=AA[c][d]-(a*AA[j][d]);
            }
          }
        }
      }
      for ( j=3 ; j>=1 ; j--)                                                         //Finding varphi by result of Gauss Elimination over AA
      {
        varphi[i][j]=AA[j][4];
        for(c=j+1 ; c<=3 ; c++)
        {
          varphi[i][j] = varphi[i][j] - (AA[j][c]*varphi[i][c]);
        }
      }
      //printf("Using Gauss Elimination with solution : %.3f %.3f %.3f \n\n", varphi[i][1],varphi[i][2],varphi[i][3]); //showing the result  
    }
    d=1;
    for(i=1;i<=3;i++)                                                                 //saving the c1 & c2 as nabla, used for calculate L2 & H norm
    {
      for(j=2;j<=3;j++)
      {
        nabla[k][d] = varphi[i][j];
        d=d+1;
      }
    }
    //printf("nabla - %d = %.3f %.3f %.3f %.3f %.3f %.3f \n",k,nabla[k][1],nabla[k][2],nabla[k][3],nabla[k][4],nabla[k][5],nabla[k][6]);
  }

  for(k=1;k<=ne;k++)                                                                //calculation for each triangle k
  {
    //printf("\n------TRIANGLE %d (%d,%d,%d)-------\n", k, elnp[k][1], elnp[k][2], elnp[k][3]); 
    for(j=1;j<=3;j++)                                                               //calculate local & global matrix A
    {                                                                  
      for(c=1;c<=3;c++)
      {
        a = nabla[k][2*j-1]*nabla[k][2*c-1]+nabla[k][2*j]*nabla[k][2*c];
        locA[j][c] = area[k]*a*nu;                                                  //local matrix A
        //printf("for (j,c) = (%d,%d) loc A = %.3f * %.3f = %.3f \n", j, c, area[k], a, locA[j][c]);
        A[elnp[k][j]][elnp[k][c]] = A[elnp[k][j]][elnp[k][c]]+ locA[j][c] ;         //global matrix A
        //printf("FIX Global A[%d][%d] = %.3f \n",elnp[k][j],elnp[k][c],A[elnp[k][j]][elnp[k][c]]);
        if(j==c) 
        {
          M_1[elnp[k][j]][elnp[k][c]] = M_1[elnp[k][j]][elnp[k][c]]+area[k]/(t*6.) ;                            //Matrix M LHS
          //printf("M(%d,%d) = %.3f ",elnp[k][j],elnp[k][c],M_1[elnp[k][j]][elnp[k][c]]);
        }
        else 
        {
          M_1[elnp[k][j]][elnp[k][c]] = M_1[elnp[k][j]][elnp[k][c]]+area[k]/(t*12.) ;                           //Matrix M LHS
          //printf("M(%d,%d) = %.3f ",elnp[k][j],elnp[k][c],M_1[elnp[k][j]][elnp[k][c]]);
        }
      }
    }
  }

  for(k=1;k<=np;k++)
  {
    for(i=1;i<=np;i++)
    {
      MA[k][i] = A[k][i] + M_1[k][i];
      //printf("MA(%d,%d) = %.3f + %.3f = %.3f\n",k,i,A[k][i],M_1[k][i],MA[k][i]);
    }
  }
  /*
  printf("===============================\n");
  printf("MATRIX A WE WANT TO SOLVE\n");
  printf("===============================\n");
  for(k=1;k<=np;k++)                                                                  //show matrix A and B
  {
    //printf("line %d \n",k );
    printf("[ ");
    for(i=1;i<=np;i++) printf("%.2f ", A[k][i]);
    printf(" ] + [ ");
    for(i=1;i<=np;i++) printf("%.2f ", M_1[k][i]);
    printf(" ] \n = \n [ ");
    for(i=1;i<=np;i++) printf("%.2f ", MA[k][i]);
    printf(" ]\n");
  }
  */
  for(k=1;k<=nb;k++)                                                                   //input the BOUNDARY CONDITION x=0
  {
    for(j=1;j<=np;j++)
    {
      MA[fanp[k][1]][j]=0.;                                                            //the entry value of A is 1 in diagonal of boundary poins
      MA[j][fanp[k][1]]=0.;                                                            //and 0 for other row and column
    }
    MA[fanp[k][1]][fanp[k][1]]=1.;
  }
  /*
  printf("=====================================\n");
  printf("FINAL MATRIX A WE WANT TO SOLVE\n");
  printf("=====================================\n");
  for(k=1;k<=np;k++)                                                                  //show matrix A
  {
    for(i=1;i<=np;i++) printf(" %.2f ", MA[k][i]);
    printf("\n");
  }   
  */
  for(l=1.;l<=floor(T/t);l++)                                                  //iteration for each timestep
  {
    printf(" t = %.3f \n", t*l);
    for(k=1;k<=np;k++)
    {
      exact[k] = (2+sin(M_PI*t*l))*sin(M_PI*npxy[k][1])*sin(M_PI*npxy[k][2]);                             //exact solution used to find error : q=sin(pi*x)*sin(pi*y)
      f[k] = (M_PI*cos(M_PI*t*l)+nu*2*M_PI*M_PI*(2+sin(M_PI*t*l)))*sin(M_PI*npxy[k][1])*sin(M_PI*npxy[k][2]) ;
      //printf("f[%d] = %.3f \n",k,f[k]);
      if(exact[k]>=max) max=exact[k];
      B[k]=0.;
      M_2[k]=0.;
    }
    
    sprintf (filename, "hh_exact%.2f.txt", l);					                              //save data for every step
    exactt = fopen(filename,"w");                                               //print 3D nodes to file
    for( i=1 ; i<=ne ; i++)
    {
      for(j=1;j<=3;j++)
      {
        fprintf(exactt,"%f %f %f \n", npxy[elnp[i][j]][1], npxy[elnp[i][j]][2], exact[elnp[i][j]]);
      }
      fprintf(exactt,"%f %f %f \n", npxy[elnp[i][1]][1], npxy[elnp[i][1]][2], exact[elnp[i][1]]);
      fprintf(exactt,"\n\n");
    }
    fclose(exactt);
    
    //printf("================TIME STEP - %.2f = %.3f until %.2f=================\n",l,t*l, floor(T/t));
    for(k=1;k<=ne;k++)                                                                //calculation for each triangle k
    {
      //printf("\n------TRIANGLE %d (%d,%d,%d)-------\n", k, elnp[k][1], elnp[k][2], elnp[k][3]); 
      for(j=1;j<=3;j++)                                                               //calculate local & global matrix A
      {
        for(c=1;c<=3;c++)
        {
          if(j==c) 
          {
            M_2[elnp[k][j]] = M_2[elnp[k][j]]+(x[elnp[k][j]]*area[k])/(t*6.) ;                          //Matrix M RHS
            locB[c] = f[elnp[k][j]]*area[k]/6.;                                       //Matrix B
            //printf("ke  k=%d j=%d >> loc B[%d] = %.5f = %.5f(f)*%.5f \n", k, j, c, locB[c], f[elnp[k][j]], area[k]/6.);
          }
          else 
          {
            M_2[elnp[k][j]] = M_2[elnp[k][j]]+(x[elnp[k][j]]*area[k])/(t*12.) ;                          //Matrix M RHS
            locB[c] = f[elnp[k][j]]*area[k]/12.;                                      //local matrix B
            //printf("ke k=%d j=%d >> loc B[%d] = %.5f = %.5f(f)*%.5f \n",k, j, c, locB[c], f[elnp[k][j]], area[k]/12.);
          }
          B[elnp[k][j]] = B[elnp[k][j]]+locB[c];
          //printf("FIX Global B[%d] = %.5f \n", elnp[k][j], B[elnp[k][j]]);
        }
      }
    }

    for(k=1;k<=np;k++)
    {
      MB[k] = B[k] + M_2[k];
      //printf("MB(%d,%d) = %.3f + %.3f = %.3f\n",k,i,B[k][i],M_2[k][i],MB[k][i]);
    }
    /*
    printf("===============================\n");
    printf("MATRIX B WE WANT TO SOLVE\n");
    printf("===============================\n");
    for(k=1;k<=np;k++)                                                                  //show matrix A and B
    {
      //printf("line %d \n",k );
      printf("[%.2f] + [%.2f] = [%.2f] \n", B[k], M_2[k], MB[k]);
    }
    */
  
    for(k=1;k<=nb;k++)                                                                   //input the BOUNDARY CONDITION x=0
    {
      MB[fanp[k][1]]=0.;
      //printf("%d. Line-%d become 0 \n",k,fanp[k][1]);
    }
    /*
    printf("===============================\n");
    printf("FINAL MATRIX WE WANT TO SOLVE\n");
    printf("===============================\n");
    for(k=1;k<=np;k++)                                                                  //show matrix A
    {
      for(i=1;i<=np;i++) printf(" %.2f ", MA[k][i]);
      printf(" = %.2f \n", MB[k]);
    }
    */
                                                                                   //CONJUGATE GRADIENT METHOD
    j=0;                                                                                //number of step
    jmax = 1000;

    for(k=1;k<=np;k++)                                                                  //first iteration
    {
      a=0.;
      for(c=1;c<=np;c++)                                                                //r = b - Ax
      {
        a = a + MA[k][c]*x[c];
      }
      r[k]=MB[k]-a;
      p[k]=r[k];
      rr = rr + r[k]*r[k];
    }

    double rrnol=rr;
    double eps=0.000000000001;                                                          //termination epsilon
  
    while(j<jmax && rr>(eps*eps)*rrnol)                                                 //iteration process
    {
      for(k=1;k<=np;k++)
      {
        for(c=1;c<=np;c++)
        {
          Ap[k] = Ap[k] + MA[k][c]*p[c]; 
        }
        pAp = pAp + p[k]*Ap[k]; 
      }
      alpha = rr/pAp;                                                                   //value of alpha
      for(k=1;k<=np;k++)
      {
        x[k] = x[k] + (alpha*p[k]);
      }
      if(j%50==0)
      {
        for(k=1;k<=np;k++)
        {
          a=0.;
          for(c=1;c<=np;c++)                                                            //r = b - Ax
          {
            a = a + MA[k][c]*x[c];
          }
          r[k] = MB[k] - a ;
        }
      }
      else
      {
        for(k=1;k<=np;k++)
        {
          r[k] = r[k] - (alpha*Ap[k]);
        }
      }
      for(k=1;k<=np;k++)
      {
        rnew = rnew + r[k]*r[k]; 
      }
      beta = rnew/rr;                                                                   //value of beta
      for(k=1;k<=np;k++)
      {
        p[k] = r[k] + (beta * p[k]);
        Ap[k]=0.;
      }
      rr = rnew;
      j = j+1;                                                                          //number of step
      //printf("step - %d with norm square = %f\n", j, rr);
      pAp=0.;
      rnew = 0.;
    }

    for(k=1;k<=np;k++)
    {
      if(x[k]>=xmax) xmax = x[k];
    }

    //printf("SOLUTION OF \t CG \t| EXACT\n");
    for(k=1;k<=np;k++)
    {
      //printf("x[%d] = %.5f \t %.5f\n", k, x[k], exact[k]);
    }
    for(k=1;k<=np;k++)
    {
      v[k]=x[k]-exact[k];
      //printf("v[%d] = %.5f \n",k,v[k]);
    }
    
    L2_norm_CG = norm_L2(elnp,area,v,ne);
    if(L2_norm_CG>=maxL2) maxL2 = L2_norm_CG;
    H_norm_CG = norm_H(elnp,area,nabla,v,ne);
    if(H_norm_CG>=maxH) maxH = H_norm_CG;
  
    printf(" \n Sum of area = %f \n", sumarea);
    printf(" Maximum value of x using CG Method = %f \n", xmax);
    printf(" Maximum value of x in exact solution = %f \n", max);
    xmax =0.; max =0.;
    
    sprintf (filename, "hh_sol%.2f.txt", l);					                              //save data for every step
    sol_CG = fopen(filename,"w");                                               //print 3D nodes to file
    for( i=1 ; i<=ne ; i++)
    {
      for(j=1;j<=3;j++)
      {
        fprintf(sol_CG,"%f %f %f \n", npxy[elnp[i][j]][1], npxy[elnp[i][j]][2], x[elnp[i][j]]);
      }
      fprintf(sol_CG,"%f %f %f \n", npxy[elnp[i][1]][1], npxy[elnp[i][1]][2], x[elnp[i][1]]);
      fprintf(sol_CG,"\n\n");
    }
    fclose(sol_CG);
  }
  /*
  sprintf (filename, "norm_hh.txt");					                              //save data for every step
  norm = fopen(filename,"a");
  fprintf(norm, "%f %f %f \n", 4./nb , maxL2, maxH);
  fclose(norm);
  */
  fclose(fp);
  
  
  free_dmatrix(npxy,1,np,1,DIM);                                                      //free
  free_imatrix(elnp,1,ne,1,DIM+1);
  free_imatrix(fanp,1,nb,1,DIM);
  free_dmatrix(varphi,1,3,1,3);
  free_dmatrix(A,1,np,1,np);
  free_dmatrix(AA,1,3,1,3+1);
  free_dvector(B,1);
  free_dmatrix(M_1,1,np,1,np);
  free_dvector(M_2,1);
  free_dmatrix(MA,1,np,1,np);
  free_dvector(MB,1);
  free_dvector(x,1);
  free_dmatrix(locA,1,3,1,3);
  free_dvector(locB,1);
  free_dvector(f,1);
  free_dvector(exact,1);
  free_dvector(area,1);
  free_dmatrix(nabla,1,ne,1,9);

  free_dmatrix(aa,1,3,1,3);                                                           //free used for showing the result

  free_dvector(r,1);                                                                  //free used for CG method
  free_dvector(p,1);
  free_dvector(Ap,1);

  free_dvector(v,1);
  return 0;
}
