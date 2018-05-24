////
#define ep1 1.0
#define ep2 1.0
#define ep3 1.0
#define ep4 1.0
#define ep5 1.0
#define ep6 1.0
#define be1 1.0
#define be2 1.0
#define be3 1.0
#define be4 1.0
#define be5 1.0
#define be6 1.0
double sign(double s){ if(s<0) {return -1.0;} else {return 1.0;} }
double vphi(double a, double ep, double be){ return pow(a*a + ep, be); }

// for u_theta \noteq 0, swirl
double u1_e(int DIM,double *xy,double t){ char s[]="u1_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];
   r=sqrt(x*x+y*y); theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);
   rho=vphi(r,ep3,-be3)*vphi(z,ep4,be4);
   ur=sign(z)*rho*uz;
   ut=vphi(r,ep5,-be5)*vphi(z,ep6,-be6);
   return ur*cos(theta)-ut*sin(theta);}
 else{ exit_func(s); return -1.0;}
}
double u2_e(int DIM,double *xy,double t){ char s[]="u2_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];
   r=sqrt(x*x+y*y); theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);
   rho=vphi(r,ep3,-be3)*vphi(z,ep4,be4);
   ur=sign(z)*rho*uz;
   ut=vphi(r,ep5,-be5)*vphi(z,ep6,-be6);
   return ur*sin(theta)+ut*cos(theta);}
 else{ exit_func(s); return -1.0;}
}
/*
// for u_theta=0, no swirl
double u1_e(int DIM,double *xy,double t){ char s[]="u1_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];r=sqrt(x*x+y*y);theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);rho=vphi(r,ep3, -be3) *
vphi(z,ep4,be4);ur=sign(z)*rho*uz;ut=0.0; return
ur*cos(theta)-ut*sin(theta);}
 else{ exit_func(s); return -1.0;}
}
double u2_e(int DIM,double *xy,double t){ char s[]="u2_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];r=sqrt(x*x+y*y);theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);rho=vphi(r,ep3,-be3) *
vphi(z,ep4,be4);ur=sign(z)*rho*uz;ut=0.0; return
ur*sin(theta)+ut*cos(theta);}
 else{ exit_func(s); return -1.0;}
}
*/

/*
// for rho=0
double u1_e(int DIM,double *xy,double t){ char s[]="u1_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];r=sqrt(x*x+y*y);theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);rho=0.0;ur=sign(z)*rho*uz;ut=vphi(r,ep5,-be5)*vphi(z,ep6,-be6);
return ur*cos(theta)-ut*sin(theta);}
 else{ exit_func(s); return -1.0;}
}
double u2_e(int DIM,double *xy,double t){ char s[]="u2_e"; double
x,y,z,ur,ut,uz,rho,r,theta;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];r=sqrt(x*x+y*y);theta=atan2(y, x);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);rho=0.0;ur=sign(z)*rho*uz;ut=vphi(r,ep5,-be5)*vphi(z,ep6,-be6);
return ur*sin(theta)+ut*cos(theta);}
 else{ exit_func(s); return -1.0;}
}
*/

// u3 is common for both swirl and no-swirl
double u3_e(int DIM,double *xy,double t){ char s[]="u3_e"; double x,y,z,uz,r;
 if(DIM==2){ exit_func(s); return -1.0;}
 else if(DIM==3){x=xy[0]; y=xy[1]; z=xy[2];
   r=sqrt(x*x+y*y);
   uz=vphi(r,ep1,-be1)*vphi(z,ep2,-be2);
   return uz;}
 else{ exit_func(s); return -1.0;}}
double g1(int DIM,double *xy,double t){return 0.0;}
double g2(int DIM,double *xy,double t){return 0.0;}
double g3(int DIM,double *xy,double t){return 0.0;}
double f1(int DIM,double *xy,double t,double nu){return 0.0;}
double f2(int DIM,double *xy,double t,double nu){return 0.0;}
double f3(int DIM,double *xy,double t,double nu){return 0.0;}
double u10(int DIM,double *xy){ return u1_e(DIM,xy,0.0);}
double u20(int DIM,double *xy){ return u2_e(DIM,xy,0.0);}
double u30(int DIM,double *xy){ return u3_e(DIM,xy,0.0);}
double p_e(int DIM,double *xy,double t){ char s[]="p_e"; exit_func(s);
return -1.0;}
double tau1(int DIM,double *xy, double *nrml,double t,double nu){ char
s[]="tau1"; exit_func(s); return -1.0;}
double tau2(int DIM,double *xy, double *nrml,double t,double nu){ char
s[]="tau2"; exit_func(s); return -1.0;}
double tau3(int DIM,double *xy, double *nrml,double t,double nu){ char
s[]="tau3"; exit_func(s); return -1.0;}
////