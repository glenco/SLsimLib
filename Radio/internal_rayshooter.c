/*
 * internal_rayshooter.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */


/*void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag){*/
void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
  /* i_points need to be already linked to s_points */
  double x_rescale[2];
  static double zs_old=-1,ro;
  long i;

  if( zsource != zs_old){
    ro=4*pi*pow(sigma_lens/2.99792e5,2)*angDist(0,zlens)*angDist(zlens,zsource)/angDist(0,zsource)/(1+zlens);
    zs_old=zsource;
  }

#pragma omp parallel for private(x_rescale)
  for(i=0;i<Npoints;++i){

    if(zsource <= zlens){
      i_points[i].image->x[0]=i_points[i].x[0];
      i_points[i].image->x[1]=i_points[i].x[1];
      i_points[i].kappa=0.0;
      i_points[i].gamma[0]=0.0; i_points[i].gamma[1]=0.0;
      i_points[i].invmag=1.0;

    }else{


      x_rescale[0]=i_points[i].x[0]/ro;
      x_rescale[1]=i_points[i].x[1]/ro;

      alphaNSIE(i_points[i].image->x,x_rescale,axis_ratio,core_size/ro,orientation);
      if(!kappa_off){
    	  gammaNSIE(i_points[i].gamma,x_rescale,axis_ratio,core_size/ro,orientation);
    	  i_points[i].kappa=kappaNSIE(x_rescale,axis_ratio,core_size/ro,orientation);
      }
/*       i_points[i].image->x[0]*=ro; */
/*       i_points[i].image->x[1]*=ro; */

      i_points[i].image->x[0]=i_points[i].x[0] - ro*i_points[i].image->x[0];
      i_points[i].image->x[1]=i_points[i].x[1] - ro*i_points[i].image->x[1];

    /*   printf("x = %e %e alpha = %e %e   kappa=%e  gamma = %e %e ro=%e\n",x[0],x[1],alpha[0],alpha[1] */
/* 	 ,i_points[i].kappa,i_points[i].gamma[0],i_points[i].gamma[1],ro); */

      i_points[i].invmag=(1-i_points[i].kappa)*(1-i_points[i].kappa)
	- i_points[i].gamma[0]*i_points[i].gamma[0] - i_points[i].gamma[1]*i_points[i].gamma[1];
    }

    i_points[i].image->invmag=i_points[i].invmag;
    i_points[i].image->kappa=i_points[i].kappa;
    i_points[i].image->gamma[0]=i_points[i].gamma[0];
    i_points[i].image->gamma[1]=i_points[i].gamma[1];
  }
}

void in_source(double *y_source,ListHndl sourcelist){
  unsigned long i,Npoints;
  double dx[2],theta;
  Point *point;

  Npoints=sourcelist->Npoints;
  MoveToTopList(sourcelist);
  for(i=0;i<Npoints;++i){
    dx[0]=sourcelist->current->x[0]-y_source[0];
    dx[1]=sourcelist->current->x[1]-y_source[1];
    theta=atan2(dx[1],dx[0]);
    if( (dx[0]*dx[0]+dx[1]*dx[1]) >
	amax*amax/( pow(amax*cos(theta+angle)/amin,2) + pow(sin(theta+angle),2) ) ){
      point=TakeOutCurrent(sourcelist);
      free(point);
/*     if( sqrt(dx[0]*dx[0]+dx[1]*dx[1]) > amax ){ */
/*       point=TakeOutCurrent(sourcelist); */
/*       free(point); */
    }
    MoveDownList(sourcelist);
  }
}
