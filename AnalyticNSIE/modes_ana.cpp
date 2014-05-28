/** 
 *  Created on: March 28, 2014
 *      Author: D. Leier
 */

#include "lens_halos.h"



PosType LensHalo::modfunc(int modnumber, PosType my_slope, PosType my_fratio){
    if(modnumber==0){return 1.0;};
    PosType a=amodfunc(modnumber, my_fratio);
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType ans=a+b*log10(my_slope);
    //std::cout << "modn: " << modnumber << "slope: " << my_slope << " a: "<< a << " b: " << b << " c: " << cmodfunc(my_fratio) << " d: " << dmodfunc(my_fratio) << " " << ans << std::endl;
    if(modnumber % 4 == 0 && modnumber > 4){return -1.0*pow(10,ans);};
    if(modnumber==4){ans+=(cmodfunc(my_fratio)/(dmodfunc(my_fratio)+my_slope)); return -1.0*pow(10.0,ans);};
	return 0.0;
}

PosType LensHalo::hi_order_modfunc(PosType x, int modnumber, PosType my_slope, PosType my_fratio){
    int k=modnumber/2;
    double da=dmod(x,modnumber,my_slope,my_fratio), dda=ddmod(x,modnumber,my_slope,my_fratio);
    PosType dans=(my_slope*my_slope-k*k)/(my_slope*my_slope-k*k+x*(2.*my_slope+1.)*da+x*x*(dda+da*da));
    if(modnumber==0){return 1.0;};
    PosType a=amodfunc(modnumber, my_fratio);
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType ans=a+b*log10(my_slope);
    //std::cout << "DANS: " << dans<< " " << 2.*my_slope+1. <<" " << da << std::endl;
    //std::cout << "modn: " << modnumber << "slope: " << my_slope << " a: "<< a << " b: " << b << " c: " << cmodfunc(my_fratio) << " d: " << dmodfunc(my_fratio) << " " << ans << std::endl;
    if(modnumber % 4 == 0 && modnumber > 4){return -1.0*pow(10.0,ans)*dans;};
    if(modnumber==4){ans+=(cmodfunc(my_fratio)/(dmodfunc(my_fratio)+my_slope)); return -1.0*pow(10.0,ans)*dans;};
	return 0.0;
}

PosType LensHalo::amodfunc(int modnumber, PosType my_fratio){
    //PosType ans=mod_params[0][modnumber][0]+mod_params[0][modnumber][1]*my_fratio;
    PosType ans=analModes(0,modnumber,0)+analModes(0,modnumber,1)*my_fratio;
    if(modnumber % 4 == 0 && modnumber > 4){
        //ans+=(mod_params[0][modnumber][2]/(mod_params[0][modnumber][3]-my_fratio)+mod_params[0][modnumber][4]*pow(my_fratio,mod_params[0][modnumber][5]));
        ans+=(analModes(0,modnumber,2)/(analModes(0,modnumber,3)-my_fratio)+analModes(0,modnumber,4)*pow(my_fratio,analModes(0,modnumber,5)));
        return -1.0*ans;
    };
    if(modnumber==4){
        //ans+=mod_params[0][modnumber][2]*log10(mod_params[0][modnumber][3]*my_fratio/pow(1-my_fratio,mod_params[0][modnumber][4]));
        ans+=analModes(0,modnumber,2)*log10(analModes(0,modnumber,3)*my_fratio/pow(1-my_fratio,analModes(0,modnumber,4)));
        return -1.0*ans;
    };
	return 0.0;
}

PosType LensHalo::bmodfunc(int modnumber, PosType my_fratio){
    //PosType ans=mod_params[1][modnumber][0]+mod_params[1][modnumber][1]*my_fratio;
    PosType ans=analModes(1,modnumber,0)+analModes(1,modnumber,1)*my_fratio;
    if(modnumber % 4 == 0 && modnumber > 4){
        //ans+=(mod_params[1][modnumber][2]*exp(mod_params[1][modnumber][3]*pow(my_fratio,mod_params[1][modnumber][4]))+mod_params[1][modnumber][5]/my_fratio);
        ans+=(analModes(1,modnumber,2)*exp(analModes(1,modnumber,3)*pow(my_fratio,analModes(1,modnumber,4)))+analModes(1,modnumber,5)/my_fratio);
        return ans;
    };
    if(modnumber==4){
        //ans+=mod_params[1][modnumber][2]/(pow(my_fratio,mod_params[1][modnumber][3])-mod_params[1][modnumber][4])+mod_params[1][modnumber][5]/pow(my_fratio,mod_params[1][modnumber][6]);
        ans+=analModes(1,modnumber,2)/(pow(my_fratio,analModes(1,modnumber,3))+analModes(1,modnumber,4));
        //std::cout << "JAJAJA " << analModes(1,modnumber,1) << " " << analModes(1,modnumber,2)  << std::endl;
        return ans;
    };
	return 0.0;
}

PosType LensHalo::cmodfunc(PosType my_fratio){
    return 0.0079/(0.407+6.909*my_fratio)-0.454;
}

PosType LensHalo::dmodfunc(PosType my_fratio){
    return 0.002/(0.028+0.709*my_fratio)-2.261;
}


PosType LensHalo::dmod(PosType x, int modnumber, PosType my_slope, PosType my_fratio){
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType d=dmodfunc(my_fratio);
    PosType c=cmodfunc(my_fratio);
    PosType db=dbfunction(x); //dbnum(x)
    PosType ans=b*db/(my_slope*log(10));
    if(modnumber==4){
        ans+=-c/(d+my_slope)/(d+my_slope)*db;
    }
    return ans;
}

PosType LensHalo::dlnmod_dr(PosType x, int modnumber, PosType my_slope, PosType my_fratio){
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType d=dmodfunc(my_fratio);
    PosType c=cmodfunc(my_fratio);
    PosType dbdr=dbfunction(x); //dbnum(x)
    PosType ans=dbdr*(b/my_slope);
    if(modnumber==4){
       ans+=-1.0*((c*log(10))/(d+my_slope)/(d+my_slope))*dbdr;
    }
    return ans;
}
        
        
PosType LensHalo::ddmod(PosType x, int modnumber, PosType my_slope, PosType my_fratio){
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType d=dmodfunc(my_fratio);
    PosType c=cmodfunc(my_fratio);
    PosType db=dbfunction(x);
    PosType ddb=ddbfunction(x);
    PosType ans=b*(ddb/(my_slope*log(10))-(db*db/(1/my_slope/my_slope/log(10))));
    if(modnumber==4){
        ans+=2*c/(d+my_slope)/(d+my_slope)/(d+my_slope)*db*db-c/(d+my_slope)/(d+my_slope)*ddb;
    }
    return ans;
}

PosType LensHalo::ddlnmod_dr(PosType x, int modnumber, PosType my_slope, PosType my_fratio){
    PosType b=bmodfunc(modnumber, my_fratio);
    PosType d=dmodfunc(my_fratio);
    PosType c=cmodfunc(my_fratio);
    PosType db=dbfunction(x); //dbnum(x)
    PosType ddb=ddbfunction(x); //ddbnum(x);
    PosType ans=(ddb-db*db/my_slope)*b/my_slope;
    if(modnumber==4){
        ans+=(ddb-db*db*2./(d+my_slope))*(-1.0*c*log(10))/(d+my_slope)/(d+my_slope);
    }
    return ans;
}


// fit parameters to be used in afunc and bfunc for different i
//mfunc_param [0][][] is parameter a, [1][][] is parameter b, [][i][] is the mode number 0==4,1==8 etc., [][][p_i], p_0=p_1 in eqns. (59,60)

PosType LensHalo::analModes(int ab, int mn, int pn){
    //std::cout << mn << std::endl;
    PosType arr[2][7][7] = {{{5.45400234e-01, 2.00623559e-01, 1.74273034e-01, 1.86242095e+02, 5.70675374e+00,0,0},
         {0.98519072, 1.79227955, 0.19333932, 1.05911801, 1.80285544, 175.19794791,0},
         {1.27286496, 2.52757113, 0.35561051, 1.07921878, 0.40027514, 37.47696142,0},
         {3.90009539, 4.59744288, -8.40249184, 4.13420228, 3.31751963, 6.35692121,0},
         {4.31129234, 7.55980323, -2.70971355, 1.20304148, 10.04502808, 4.77684945,0},
         {28.27243869, 23.49608644, -43.13187974, 1.65596222, 22.72500104, 3.1378364,0},
         {57.01022684, 33.55782426, -114.00434863, 2.08922268, 22.96852996, 2.53164822,0}},
        {{7.41256977e+00,3.35544981e-02,-4.99797570e+01,-9.33421144e-02,1.01892453e+01,0,0},
         {-7.15392458e+01, -1.93808978e-02, 7.48623502e+01, -1.70286774e-02, 4.20781335e+02, 3.66600276e-04,0},
         {-2.67369155e+00, -3.43984326e-02, 6.08826051e+00, -4.03323814e-01, 7.02617948e+01, 3.95509769e-04,0},
         {1.93036483e+00, -3.94871873e-02, 1.56471642e+00, -5.54299042e+00, 3.62847996e+01, 4.40377288e-04,0},
         {2.06365018e+00, -4.79733978e-02, 1.49994453e+00, -1.21602730e+01, 2.28348079e+01, 4.49107825e-04,0},
         {2.06019572e+00, -4.22588623e-02, 1.55824013e+00, -1.62854932e+01, 1.57266141e+01, 4.93672097e-04,0},
         {2.09778383e+00, -7.99553494e-02, 1.58048481e+00, -3.27089042e+01, 1.36378731e+01, 4.02453963e-04,0}}};
    
    return arr[ab][int(mn/4)-1][pn];
}

/*PosType LensHalo::betaCorr(PosType theta,PosType bfunc){
    PosType a=1.05*bfunc-1.011;
    PosType b=0.2484-0.03362*bfunc;
    PosType c=0.00425+0.085*bfunc;
    return a+pow(cos(theta)*cos(theta)+b*sin(theta)*sin(theta),c);
}*/

PosType LensHalo::betaCorr(PosType thet, PosType bfunc){
    PosType a,b,c;
    PosType theta=thet+pi/2.;
    a=0.03910245;
    b=0.21480206;
    c=0.08913886;
    return bfunc*(a+pow((cos(theta)*cos(theta)+b*sin(theta)*sin(theta)),b));
}



