#include "solver.h"

Solver::Solver(double Beta, int Nbath_max, int Nw, int NwED) : 
    Beta(Beta), Nbath_max(Nbath_max), Nw(Nw), NwED(NwED) { 

        n_bath_max=-1; NFock=0; 
        H=NULL; Psi=NULL; E=NULL; cuPsi=NULL; cdPsi=NULL; nud=NULL; 

        t2u=new double[Nbath_max]; t2d=new double[Nbath_max]; 
        
        eu=new double[Nbath_max]; ed=new double[Nbath_max]; 
        for1(l, Nbath_max) {
            eu[l]=0.1*l; 
            ed[l]=0.1*l;
        };

        Iw=new complex [Nw]; Iw2= new double [Nw];   
        for1(w, Nw){
            Iw[w]=I*(2.*w+1)*Pi/Beta; 
            Iw2[w]=-sqr((2.*w+1)*Pi/Beta); 
        } 

        gu=new complex [Nw]; gd=new complex [Nw]; 
        for1(w, Nw) {
            gu[w]=1./Iw[w]; gd[w]=1./Iw[w];
        } 

        sigmau=new complex [Nw]; sigmad=new complex [Nw]; 
        for1(w, Nw) {
            sigmau[w]=0; sigmad[w]=0;
        } 
}


Solver::~Solver()
{
    for1(j, NFock) {
        delete [] H[j]; delete [] Psi [j]; 
        delete [] cuPsi[j]; delete [] cdPsi[j]; 
    } 

    delete [] H; delete [] Psi; 
    delete [] cuPsi; delete [] cdPsi;  
    delete [] E; delete [] nud; 
    delete [] t2u; delete [] t2d;delete [] eu;delete [] ed; 
    delete [] gu; delete [] gd;

    delete [] Iw; delete [] Iw2;
}

double Solver::adjust_bath_levels(complex * Delta_w,   double * t2_bath, double * e_bath)
{
    double s, s_min=1e6; 
    int i_min=0;
    bool t2positive=false;
    static double * e_bath_min=new double [Nbath_max];
    static int Ntrials=10000;  
    static double de=1, *e_prev=new double [Nbath_max]; 
    
    for1(l, n_bath) e_prev[l]=e_bath[l];
 
    if (rnd()<.5) 
        for1(j, n_bath) 
            e_bath[j]=de*real(rnd_gauss2());
    else 
        for1(j, n_bath) 
            e_bath[j]=e_prev[j]+.3*de*real(rnd_gauss2());
    
    s_min=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
    for1(l, n_bath) e_bath_min[l]=e_bath[l];
     
    
    for1 (i, Ntrials)
    {
        if (rnd()<.5) 
            for1(j, n_bath) e_bath[j]=de*real(rnd_gauss2()); 
        else 
            for1(j, n_bath) e_bath[j]=e_prev[j]+.3*de*real(rnd_gauss2());
         
        s=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
        if (!(s<0) && s<.9999*s_min) {
            s_min=s; i_min=i;
            for1(j, n_bath) e_bath_min[j]=e_bath[j]; 
            
            t2positive=true; 
            for1(l, n_bath) t2positive=(t2positive && t2_bath[l]>0);
        }       

        for1 (ii, 50)
        {
            static double de=.001, * s1=new double[Nbath_max], ** s2=new_double2(Nbath_max, Nbath_max); 
            for1(j, n_bath) 
            {
                e_bath[j]+=de;  
                double sp=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 

                e_bath[j]-=2.*de; 
                double sm=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
                e_bath[j]+=de; s1[j]=(sp-sm)/(2.*de); s2[j][j]=(sp+sm-2.*s)/(de*de);
            }
            for1(j1, n_bath)    
            {
                for(int j2=j1+1; j2<n_bath; j2++)
                    {
                        e_bath[j1]+=de; e_bath[j2]+=de; 
                        double sxy=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
                        e_bath[j1]-=de; e_bath[j2]-=de;
                        
                        e_bath[j1]-=de; e_bath[j2]-=de; 
                        double sxym=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
                        e_bath[j1]+=de; e_bath[j2]+=de;
                        
                        s2[j1][j2]=(sxy+sxym-2.*s-s2[j1][j1]*de*de - s2[j2][j2]*de*de)/(2.*de*de);
                        s2[j2][j1]=s2[j1][j2];
                    }
            }
        
            if (div_left(s2,s1,n_bath))  
                for1(j, n_bath) e_bath[j]-=s1[j];  
            
            double s0=s;       
            s=adjust_bath_couplings(Delta_w, t2_bath, e_bath);  
            if (s>=s0) 
                break;    
            
            if (!(s<0) && s<.9999*s_min) {
                s_min=s; i_min=i;
                for1(j, n_bath) e_bath_min[j]=e_bath[j]; 
                 
                t2positive=true; 
                for1(l, n_bath) t2positive=(t2positive && t2_bath[l]>0);
            } 
        }
        
         
        
        if (s_min<1e-10 && t2positive) 
            break;  
        if (!(s<0) && s<1e-10 && t2positive) {
            s_min=s; i_min=i;
            for1(j, n_bath) e_bath_min[j]=e_bath[j]; 
            break;
        } 
    } 
    
    for1(j, n_bath) e_bath[j]=e_bath_min[j];  
    adjust_bath_couplings(Delta_w, t2_bath, e_bath);

    if (!t2positive) {
        std::cout<<"!!! adjust_bath_levels led to unphysical bath !!!\n"<<std::flush; 
        for1(l, n_bath) std::cout<<"e="<<e_bath[l]<<"  t2="<<t2_bath[l]<<"\n" << std::flush;
    }
    
    double de_max=0; 
    for1(j, n_bath) {
        if (de_max<abs(e_bath_min[j])) 
            de_max=abs(e_bath_min[j]); 
    }
    
    Ntrials=.9*Ntrials+i_min; 
    de=.9*de+.1*de_max; 
    if (Ntrials<300) 
        Ntrials=300; 
    
    for1(j, n_bath) e_bath[j]=e_bath_min[j];   
    
    return adjust_bath_couplings(Delta_w, t2_bath, e_bath);
}


double Solver::adjust_bath_couplings(complex * Delta_w, double  * t2_bath, double * e_bath)
{
    double a=0.0, a2=0.0, a3=0.0; 
    static double ** A=new_double2(Nbath_max, Nbath_max);

    for1(i,n_bath) {
        for(int j=i; j<n_bath; j++) {
            a=0; 
            for1(w, NwED) {
                a+=(e_bath[i]/(sqr(e_bath[i])-Iw2[w])+e_bath[j]/(sqr(e_bath[j])-Iw2[w]))/(e_bath[i]+e_bath[j]);
                if (w+1==NwED/2) a2=a; 
                if (w+1==3*NwED/4) a3=a;
            } 
            A[i][j]=8.*a-9.*a3+2.*a2;  
            A[j][i]=A[i][j];
        }
    }
    
    for1(i, n_bath) {
        a=0; 
        for1(w, NwED) {
            a+=real(conj(Delta_w[w])/(Iw[w]-e_bath[i])); 
            if (w+1==NwED/2) a2=a; 
            if (w+1==3*NwED/4) a3=a;
        } 
        t2_bath[i]=8.*a-9.*a3+2.*a2; 
    }
        
    if (!div_left(A,t2_bath,n_bath)) 
        return 100; 
    
    a=0;
    for1(w, NwED) {
        complex r=Delta_w[w]; 
        for1(i, n_bath) r-=t2_bath[i]/(Iw[w]-e_bath[i]);  
        a+=norm2(r);
        if (w+1==NwED/2) a2=a; 
        if (w+1==3*NwED/4) a3=a;
    }

    return 2.*(8.*a-9.*a3+2.*a2)/Beta;
}


double Solver::init(double U, complex * Delta_up, complex * Delta_down,  int N_bath, double h_loc, double mu_loc)
{   
    n_bath=N_bath;
    if (n_bath>n_bath_max)
    {
        for1(j, NFock) {
            delete [] H[j]; delete [] Psi [j]; 
            delete [] cuPsi[j]; delete [] cdPsi[j]; 
        } 
        delete [] H; delete [] Psi; 
        delete [] cuPsi; delete [] cdPsi; 
        delete [] E; delete [] nud;

        n_bath_max = n_bath;  
        NFock = 1 << (2*n_bath_max+2);
        
        H=new_double2(NFock, NFock); Psi=new_double2(NFock, NFock); 
        cuPsi=new_double2(NFock, NFock); cdPsi=new_double2(NFock, NFock); 
        E=new double[NFock]; nud=new int[NFock];

        for1(j, NFock) {
            nud[j]=0; 
            for1(l, n_bath+1) {
                nud[j]+=n(j,l)+(n_bath+2)*n(j,l+n_bath+1); 
            }
        } 
    } 

    tolerance=adjust_bath_levels(Delta_up, t2u, eu)+adjust_bath_levels(Delta_down, t2d, ed); 
    
    for2(i,j, NFock) {
        H[i][j]=0; cuPsi[i][j]=0; cdPsi[i][j]=0;  
    }
    
    for1(j, NFock)
    {
        H[j][j] += -mu_loc * (n(j,0) + n(j, n_bath+1)) -
                    h_loc * (n(j,0) - n(j, n_bath+1)) +
                    U * ((n(j,0)-.5) * (n(j,n_bath+1)-.5));
        
        for1(l, n_bath) {
            H[j][j] += n(j, l+1)*eu[l] + n(j, n_bath+l+2)*ed[l];
        }

        for1(l, n_bath)
        {
            int nu=0, nd=0; 
            for1(l2,l) {
                nu+=n(j,l2+1); nd+=n(j,n_bath+l2+2);
            }
            if (n(j,0)==0  && n(j, l+1)==1) {
                int j1=j+1-(1<<(l+1));                        
                
                H[j][j1] += sqrt(t2u[l]) * (1-2*(nu%2));
                H[j1][j] += H[j][j1];
            }
            if (n(j,n_bath+1)==0  && n(j, l+n_bath+2)==1) {
                int j1=j+(1<<(n_bath+1))-(1<<(l+n_bath+2));   
                
                H[j][j1] += sqrt(t2d[l]) * (1-2*(nd%2));
                H[j1][j] += H[j][j1];
            }
        }
    }
    

    EigenBlock(H, Psi, E, NFock, nud, sqr(n_bath+2));
    
    
    Z=0; 
    for1(j, NFock) Z+=expl(-Beta*E[j]);  

    long double lnZ=logl(Z); 
    for1(l, n_bath) {
        lnZ -= log(1+exp(-Beta*eu[l])) + log(1+exp(-Beta*ed[l]));
    }

 
    for1(j, NFock) {
        if (n(j,0)==1) 
            for1(j1, NFock) {
                cuPsi[j1][j-1]=Psi[j1][j];
            }
        if (n(j, n_bath+1)==1) {
            double f=1-2*((nud[j]%(n_bath+2))%2); 
            for1(j1, NFock) {
                cdPsi[j1][j-(1<<(n_bath+1))]=f*Psi[j1][j];  
            } 
        }
    }
        
    for1 (w, Nw) {
        gu[w]=0; gd[w]=0;
    }

    for2(j1, j2, NFock) 
    {
        if (nud[j1]==nud[j2]-1)     
        {
            double d=0; 
            for1(j3, NFock) {
                d += Psi[j1][j3]*cuPsi[j2][j3]; 
            }
            double x = d*d*( (expl(-Beta*E[j1]) + expl(-Beta*E[j2]))/Z ); 
        
            for1(w, Nw) {
                gu[w] += x / (Iw[w]+E[j1]-E[j2]);
            }
        }

        if (nud[j1]==nud[j2]-(n_bath+2) )     
        {
            double d=0; 
            for1(j3, NFock) {
                d += Psi[j1][j3]*cdPsi[j2][j3]; 
            }
            double x = d*d*( (expl(-Beta*E[j1]) + expl(-Beta*E[j2]))/Z ); 
        
            for1(w, Nw) {
                gd[w] += x/(Iw[w]+E[j1]-E[j2]);
            }
        }

    }

    for1(w, Nw) {
        sigmau[w] = Iw[w] - 1./gu[w] + h_loc + mu_loc; 
        sigmad[w] = Iw[w] - 1./gd[w] - h_loc + mu_loc; 

        for1(l, n_bath) {
            sigmau[w] -= t2u[l]/(Iw[w]-eu[l]); 
            sigmad[w] -= t2d[l]/(Iw[w]-ed[l]);
        } 
    }
    
    lnZ -= Beta * mu_loc;
    
    return lnZ;
}

double Solver::double_occupancy() const {
    double D = 0;
    for1(m, NFock) {
        double Dm = 0;
        for1(j, NFock)
            Dm += sqr(Psi[m][j]) * n(j, 0) * n(j, n_bath + 1);
        D += expl(-Beta * E[m]) * Dm;
    }
    D /= Z;
    return D;
}

double Solver::mean_occupancy() const {
    double N = 0;
    for1(m, NFock) {
        double Nm = 0;
        for1(j, NFock)
            Nm += sqr(Psi[m][j]) * (n(j, 0) + n(j, n_bath + 1));
        N += expl(-Beta * E[m]) * Nm;
    }
    N /= Z;
    return N;
}

double Solver::mean_n_up() const {
    double N = 0;
    for1(m, NFock) {
        double Nm = 0;
        for1(j, NFock)
            Nm += sqr(Psi[m][j]) * n(j, 0);
        N += expl(-Beta * E[m]) * Nm;
    }
    N /= Z;
    return N;
}

double Solver::mean_n_down() const {
    double N = 0;
    for1(m, NFock) {
        double Nm = 0;
        for1(j, NFock)
            Nm += sqr(Psi[m][j]) * n(j, n_bath + 1);
        N += expl(-Beta * E[m]) * Nm;
    }
    N /= Z;
    return N;
}
