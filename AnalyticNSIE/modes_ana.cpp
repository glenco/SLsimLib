/** 
 *  Created on: March 28, 2014
 *      Author: D. Leier
 */

#include "lens_halos.h"


/* 
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
 */


// fit parameters to be used in afunc and bfunc for different i
//mfunc_param [0][][] is parameter a, [1][][] is parameter b, [][i][] is the mode number 0==4,1==8 etc., [][][p_i], p_0=p_1 in eqns. (59,60)

/*PosType LensHalo::analModes(int ab, int mn, int pn){
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
 */


void LensHalo::analModes(int modnumber, PosType my_beta, PosType q, PosType amod[3]){
    PosType x=my_beta;
    PosType xsq=x*x;
    PosType qsq=q*q;
    PosType qqq=qsq*qsq;
    if(modnumber==4){
        PosType p[4];
        p[0]=0.999247407765171+1.57671597806174*qsq+0.226801415239749*qqq-0.00113150933379816/(0.00982744492049579+q)-1.92526147723216*q-0.876057800660322*qsq*q;
        p[1]=(0.256689163217493-0.349094565459393*qsq)/(0.240187114931842+7.04774689369728*q+127.355470876729*qqq*qsq+74.7649585603172*qqq+40.6443593959714*qsq*q +7.04774689369728*qsq);
        p[2]=(0.869382229417186*q-0.807125199076924)/(1.13660761698993+45.2956858226937*q+2121.94004854273*qqq*qqq+421.206428582793*qsq*q-38.4857728532176*qsq);
        p[3]=(0.000640486669983527-0.00104518529303727*q)/(0.00710198457626786+0.800055447745725*q+120.581051723995*qqq*q+32.2823375097222*qsq*q+1.50311757016824*qsq-21.4518837566564*qqq);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*x*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq;
    }
    if(modnumber==8){
        PosType p[5];
        p[0]=0.497555478681077+0.000204902074485823/x+3.50663151993971*qsq+2.34317902193599*qqq-1.94862571121286*q-0.628046645263927*qqq*q-3.77143120591481*qsq*q;
        p[1]=2.48928468360676+6.903227095891*q+0.365066682291406*log(q)-2.54875770684561*sin(q)-1.50679862465651*qsq-3.24349713567861*sqrt(3.13261876573299*q);
        p[2]=(q-0.921090631057856)/(7382.43442477471*qqq-0.352662146083676-42.8982972220118*q-4358.26871100734*qsq*q-121695.680019394*qqq*qsq);
        p[3]=(0.0277306792378124*q-0.0217671131740949)/(0.00997700501985252+q+747.18250982068*qqq*qsq*q+43.9903595800984*qsq*q+1.61849149117731*qsq);
        p[4]=26.0020256464656/(57.5695837127945+4624.73063241391*q+44981192.5548355*qqq*qsq*q+2056746.30474371*qqq+58135.4015130137*qsq);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq+p[4]*xsq*xsq*x;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*xsq*x+p[4]*5.*xsq*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq+p[4]*20.*xsq*x;
    }
    if(modnumber==12){
        PosType p[5];
        p[0]=(1.37787428901026*q-0.630373108405854-0.754474021902366*qsq)/(-1.7724943276976-8.84252195876267*q-2.94086666623792*qsq-28.7963682928423*qsq*q);
        p[1]= (16.5298835055148*q + 55.2192594512767*qqq)/(pow(277.645516542271,q+pow(q,4.26447869683269)+sqrt(q)) - pow(46.2928262754497*q,0.0612738432254949/sqrt(q)));
        p[2]=0.0285580878225644/(0.0100292334359904+q+337.859555992645*qqq*qsq+55.3098246205842*qsq*q+49.8287458342673*qqq*q+6.05384360506144*qqq+6.05384360506144*qsq);
        p[3]=(q-qsq)/(0.00226252608467631+2244.83274325001*qqq+48.8491549267231*qsq-682.187594099543*qqq*q);
        p[4]=0.021571960010229/(0.0333135143092527+5.96741147805435*q+4467561.92719737*qqq*qqq+1169.22725067438*qsq*q);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq+p[4]*xsq*xsq*x;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*xsq*x+p[4]*5.*xsq*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq+p[4]*20.*xsq*x;
    }
    if(modnumber==16){
        PosType p[5];
        p[0]=(0.0332506690007129+0.0507849727193692*qsq-0.0818798140066402*q)/(0.117773213911455+1.00258035623967*q+23.6760013441453*qqq*qqq*qsq*q+ 7.87516159052667*qsq*q);
        //p[0]=(0.0341321710217401+0.0539589185062629*qsq-0.0854561071707757*q)/(0.120919700383776+1.01676574135655*q+7.44144054666771*qsq*q+0.0850587108046311*qsq);
        //p[1]=0.004951452353*q/(0.0004819509283+0.00596250089042225*q+133.347137165264*qqq*qqq+4.45414660698931*qqq+0.227069236590851*qsq);
        p[1]=(0.370843148866421*q-0.00313460015313055)/(0.00580285397911041+q*exp(9.3730127916388*q));
        p[2]=0.0511355873898476/(0.0176742950717808+1.54692258021449*q+15325.1426627127*qqq*qqq*qsq+565.977896268815*qqq+21.574807463678*qsq);
        p[3]=-0.00214038765774794/(0.000704011874446295+0.127676507519154*q+8717.84655997113*qqq*qsq*q+22.3344032381803*qsq*q);
        p[4]= 0.00319327126262356/(0.00480055091351625+1.00925649375844*q+101539659.312983*qqq*qqq*q+385.154374016887*qsq*q);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq+p[4]*xsq*xsq*x;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*xsq*x+p[4]*5.*xsq*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq+p[4]*20.*xsq*x;
    }
    if(modnumber==20){
        PosType p[5];
        p[0]=0.00216929246186679/(0.00874532408794626+0.150239670550345*q+661.308273994221*qqq*qqq*q+17.2984134591348*qqq+0.377245995454281*qsq);
        p[1]=(0.00690072719794233+1.55926022949841*qsq-qsq*q-0.623885933520025*q)/(-2.55838861366553*q-176.802304238127*qsq*q-3217.71548020031*qqq*qsq*q);
        p[2]=0.0250561178602656/(0.00677502952553192+q+2747.01604139518*qqq*qsq+114.481847626077*qsq*q+4.20620916094827*qsq);
        p[3]=(0.0269794302997856-0.0792906388766586*q)/(-0.0124943423636934-1.09649706499358*q-26.2275198235134*qsq-12717.1851237396*qqq*q);
        p[4]=(3.84628930462709e-6-3.31250089838063e-5*q)/(4.85021837022965e-6+0.0013811919166733*q+90.5084746032337*qqq*qsq*q+0.302838155475293*qsq*q-0.00153415114620723*qsq);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq+p[4]*xsq*xsq*x;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*xsq*x+p[4]*5.*xsq*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq+p[4]*20.*xsq*x;
    }
    if(modnumber==24){
        PosType p[5];
        p[0]=(3.82661873893811e-5+7.56042377625885e-5*qsq-0.000109078486936465*q)/(0.000175126592499607+0.00341872621572824*q+0.0815844981491867*qsq*q-0.00451293919087288*qsq);
        p[1]=(x-0.0152700264815539)/(0.0595436013358412+1256796.89494238*qqq*qqq*qsq+15685.1646694688*qqq*q+103.93281463956*qsq);
        p[2]=0.0154722687102909/(0.00340798484062739+0.739971972927153*q+6733.07078343369*qqq*qsq+133.921631245673*qsq*q-890.541218846465*qqq*q);
        p[3]=0.0137409320162211/(-0.00417866808304444-1.01379599610013*q-424.692345463305*qsq*q-11608361.2829203*qqq*qqq);
        p[4]=1.03789139769077/(0.237233723446576+594.210136818616*q+5.14323465264109e21*qqq*qqq*qqq*qqq*q+9246244.7451883*qqq);
        amod[0]=p[0]*x+p[1]*xsq+p[2]*xsq*x+p[3]*xsq*xsq+p[4]*xsq*xsq*x;
        amod[1]=p[0]+p[1]*2.*x+p[2]*3.*xsq+p[3]*4.*xsq*x+p[4]*5.*xsq*xsq;
        amod[2]=p[1]*2.+p[2]*6.*x+p[3]*12.*xsq+p[4]*20.*xsq*x;
    }
    if(modnumber==28){
        std::cout << "There are yet no analytic models for Fourier modes with array index 28 or higher." << std::endl;
        amod[0]=0.;
        amod[1]=0.;
        amod[2]=0.;
    }
    if(modnumber%4!=0){
        amod[0]=0.;
        amod[1]=0.;
        amod[2]=0.;
    }
    amod[0]=0.;
    amod[1]=0.;
    amod[2]=0.;
}

