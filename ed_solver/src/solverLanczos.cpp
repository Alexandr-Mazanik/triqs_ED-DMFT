#include "solverLanczos.h"

SolverLanczos::SolverLanczos(double Beta, int Nbath_max, int Nw, int NwED) : 
    Beta(Beta), Nbath_max(Nbath_max), Nw(Nw), NwED(NwED) { 

        t2u.setZero(Nbath_max); t2d.setZero(Nbath_max); 
        
        eu.resize(Nbath_max); ed.resize(Nbath_max); 
        for1(l, Nbath_max) {
            eu[l]=0.1*l; 
            ed[l]=0.1*l;
        };
  
        Iw.resize(Nw); Iw2.resize(Nw);
        for1(w, Nw){
            Iw(w) = I*(2.*w+1)*Pi/Beta; 
            Iw2(w) = -sqr((2*w+1)*Pi/Beta); 
        } 

        gu.resize(Nw); gd.resize(Nw);
        sigmau.resize(Nw); sigmad.resize(Nw);
}

double SolverLanczos::adjust_bath_levels(const VectorXcd& Delta_w, VectorXd& t2_bath, VectorXd& e_bath)
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
        throw std::runtime_error("Unphysical bath: invalid hybridization");
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


double SolverLanczos::adjust_bath_couplings(const VectorXcd& Delta_w, VectorXd& t2_bath, VectorXd& e_bath)
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
        
    if (!div_left(A, t2_bath.data(), n_bath)) 
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


double SolverLanczos::init(double U, 
                           const VectorXcd& Delta_up, const VectorXcd& Delta_down,  
                           int N_bath, int Nm, int Nkr, double h_loc, double mu_loc) {
    n_bath = N_bath;
    if (n_bath > n_bath_max) {

        n_bath_max=n_bath;  
        NFock = 1 << (2*n_bath_max+2);

        nud.setZero(NFock);
        for1(j, NFock) {
            for1(l, n_bath + 1) 
                nud(j) += n(j,l) + (n_bath+2) * n(j, l+n_bath+1);
        }

        H.resize(NFock, NFock);
        H.reserve(Eigen::VectorXi::Constant(NFock, 1 + 2*(n_bath+1)));
    } 

    tolerance = adjust_bath_levels(Delta_up, t2u, eu) + adjust_bath_levels(Delta_down, t2d, ed); 

    std::vector<Eigen::Triplet<double>> H_elements;

    for1(j, NFock) {
        double diag = 0.0;
        diag += -mu_loc * (n(j,0) + n(j, n_bath+1)) -
                        h_loc * (n(j,0) - n(j, n_bath+1)) +
                        U * ((n(j,0)-.5) * (n(j, n_bath+1)-.5));
        for1(l, n_bath) {
            diag += n(j, l+1) * eu[l] + n(j, n_bath+l+2) * ed[l];
        }
        H_elements.push_back(Eigen::Triplet<double>(j, j, diag));

        for1(l, n_bath) {
            int nu=0, nd=0; 
            for1(l2,l) {
                nu+=n(j,l2+1); nd+=n(j,n_bath+l2+2);
            }
            if (n(j,0)==0  && n(j, l+1)==1) {
                int j1=j+1-(1<<(l+1));                        

                double value = sqrt(t2u[l]) * (1-2*(nu%2));
                H_elements.push_back(Eigen::Triplet<double>(j, j1, value));
                H_elements.push_back(Eigen::Triplet<double>(j1, j, value));
            }
            if (n(j,n_bath+1)==0  && n(j, l+n_bath+2)==1) {
                int j1=j+(1<<(n_bath+1))-(1<<(l+n_bath+2));   
                
                double value = sqrt(t2d[l]) * (1-2*(nd%2));
                H_elements.push_back(Eigen::Triplet<double>(j, j1, value));
                H_elements.push_back(Eigen::Triplet<double>(j1, j, value));
            }
        }
    }
    H.setFromTriplets(H_elements.begin(), H_elements.end());
    H.makeCompressed(); 

    Spectra::SparseSymMatProd<double> op(H);
    int ncv;
    if (Nm < 10) 
        ncv = std::min(30, NFock);
    else
        ncv = std::min(3 * Nm, NFock);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, Nm, ncv);
    
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::SmallestAlge, 10000, 1e-10);
    if (eigs.info() != Spectra::CompInfo::Successful) {
        std::cerr << "Spectra convergence info: " << static_cast<int>(eigs.info()) << std::endl;
        std::cerr << "Converged eigenvalues: " << nconv << " out of " << Nm << std::endl;
        throw std::runtime_error("Spectra has not converged.");
    } 

    E = eigs.eigenvalues();
    Psi = eigs.eigenvectors();

    Z = 0; 
    for1(j, Nm) 
        Z += expl(-Beta * E(j));  

    long double lnZ = logl(Z); 
    for1(l, n_bath) {
        lnZ -= log(1+exp(-Beta*eu[l])) + log(1+exp(-Beta*ed[l]));
    }

    gu.setZero(Nw); gd.setZero(Nw);
    sigmau.setZero(Nw); sigmad.setZero(Nw);

    for1(m, Nm) {
        const VectorXd psi_m = Psi.col(m);
        const double E_m = E(m);
        
        int block_idx = get_block_idx(psi_m);
        if (block_idx == -1) {
            std::cout << "Skipping invalid eigenvector " << m << "\n" << std::flush;
            continue;  
        }

        std::vector<int> ib;
        std::unordered_map<int, int> ibu_p, ibu_h, ibd_p, ibd_h;
        int ku_p = 0, ku_h = 0, kd_p = 0, kd_h = 0;
        for1(j, NFock) {
            if (nud(j) == block_idx) {
                ib.push_back(j);
            }
            if (nud(j) == block_idx + 1) {
                ibu_p[j] = ku_p;
                ku_p++;
            }
            if (nud(j) == block_idx - 1) {
                ibu_h[j] = ku_h;
                ku_h++;                
            }
            if (nud(j) == block_idx + (n_bath + 2)) {
                ibd_p[j] = kd_p;
                kd_p++;
            }
            if (nud(j) == block_idx - (n_bath + 2)) {
                ibd_h[j] = kd_h;
                kd_h++;
            }
        }

        VectorXd v0u_p, v0u_h, v0d_p, v0d_h;
        v0u_p.setZero(ku_p); v0u_h.setZero(ku_h);
        v0d_p.setZero(kd_p); v0d_h.setZero(kd_h);

        for (size_t i = 0; i < ib.size(); ++i) {
            if (n(ib[i], 0) == 0)
                v0u_p.coeffRef(ibu_p.at(ib[i] ^ (1 << 0))) += psi_m[ib[i]];
            else 
                v0u_h.coeffRef(ibu_h.at(ib[i] ^ (1 << 0))) += psi_m[ib[i]];
            if (n(ib[i], n_bath + 1) == 0)
                v0d_p.coeffRef(ibd_p.at(ib[i] ^ (1 << (n_bath+1)))) += fermion_sign(ib[i], n_bath + 1) * psi_m[ib[i]];
            else 
                v0d_h.coeffRef(ibd_h.at(ib[i] ^ (1 << (n_bath+1)))) += fermion_sign(ib[i], n_bath + 1) * psi_m[ib[i]];
        }

        double nu = v0u_h.squaredNorm();
        double nd = v0d_h.squaredNorm();


        v0u_h *= (sqrt(nu) > 1e-12) ? 1.0 / sqrt(nu) : 0.0;
        v0u_p *= (sqrt(1 - nu) > 1e-12) ? 1.0 / sqrt(1 - nu) : 0.0;
        v0d_h *= (sqrt(nd) > 1e-12) ? 1.0 / sqrt(nd) : 0.0;
        v0d_p *= (sqrt(1 - nd) > 1e-12) ? 1.0 / sqrt(1 - nd) : 0.0;

        std::vector<Eigen::Triplet<double>> Hu_p_elem, Hu_h_elem, Hd_p_elem, Hd_h_elem;
        for (Eigen::Index k = 0; k < H.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
                int row = it.row();
                int col = it.col();
                double val = it.value();

                auto blu_p_row = ibu_p.find(row); auto blu_p_col = ibu_p.find(col);
                auto blu_h_row = ibu_h.find(row); auto blu_h_col = ibu_h.find(col);
                auto bld_p_row = ibd_p.find(row); auto bld_p_col = ibd_p.find(col);
                auto bld_h_row = ibd_h.find(row); auto bld_h_col = ibd_h.find(col);

                if (blu_p_row != ibu_p.end() && blu_p_col != ibu_p.end()) 
                    Hu_p_elem.emplace_back(blu_p_row->second, blu_p_col->second, val);
                if (blu_h_row != ibu_h.end() && blu_h_col != ibu_h.end()) 
                    Hu_h_elem.emplace_back(blu_h_row->second, blu_h_col->second, val);
                if (bld_p_row != ibd_p.end() && bld_p_col != ibd_p.end()) 
                    Hd_p_elem.emplace_back(bld_p_row->second, bld_p_col->second, val);
                if (bld_h_row != ibd_h.end() && bld_h_col != ibd_h.end()) 
                    Hd_h_elem.emplace_back(bld_h_row->second, bld_h_col->second, val);
            }
        }

        SparseMatrixXd Hu_p(ku_p, ku_p), Hu_h(ku_h, ku_h);
        SparseMatrixXd Hd_p(kd_p, kd_p), Hd_h(kd_h, kd_h);
        Hu_p.setFromTriplets(Hu_p_elem.begin(), Hu_p_elem.end());
        Hu_h.setFromTriplets(Hu_h_elem.begin(), Hu_h_elem.end());
        Hd_p.setFromTriplets(Hd_p_elem.begin(), Hd_p_elem.end());
        Hd_h.setFromTriplets(Hd_h_elem.begin(), Hd_h_elem.end());
        Hu_p.makeCompressed(); Hu_h.makeCompressed();
        Hd_p.makeCompressed(); Hd_h.makeCompressed();

        auto coeff_u_p = LanczosTridiag(Hu_p, v0u_p, Nkr);
        auto coeff_u_h = LanczosTridiag(Hu_h, v0u_h, Nkr);
        auto coeff_d_p = LanczosTridiag(Hd_p, v0d_p, Nkr);
        auto coeff_d_h = LanczosTridiag(Hd_h, v0d_h, Nkr);

        VectorXcd gu_p = VectorXcd::Zero(Nw), gu_h = VectorXcd::Zero(Nw);
        VectorXcd gd_p = VectorXcd::Zero(Nw), gd_h = VectorXcd::Zero(Nw);
        for (int i = coeff_u_p.krylov_dim - 1; i >= 0; --i) {
            double a_i = coeff_u_p.a_vec(i);
            double b_i = coeff_u_p.b_vec(i);
            VectorXcd d = Iw.array() + E_m - a_i - sqr(b_i) * gu_p.array();
            gu_p = d.array().inverse();
        }
        for (int i = coeff_u_h.krylov_dim - 1; i >= 0; --i) {
            double a_i = coeff_u_h.a_vec(i);
            double b_i = coeff_u_h.b_vec(i);
            VectorXcd d = Iw.array() - E_m + a_i - sqr(b_i) * gu_h.array();
            gu_h = d.array().inverse();
        }
        for (int i = coeff_d_p.krylov_dim - 1; i >= 0; --i) {
            double a_i = coeff_d_p.a_vec(i);
            double b_i = coeff_d_p.b_vec(i);
            VectorXcd d = Iw.array() + E_m - a_i - sqr(b_i) * gd_p.array();
            gd_p = d.array().inverse();
        }
        for (int i = coeff_d_h.krylov_dim - 1; i >= 0; --i) {
            double a_i = coeff_d_h.a_vec(i);
            double b_i = coeff_d_h.b_vec(i);
            VectorXcd d = Iw.array() - E_m + a_i - sqr(b_i) * gd_h.array();
            gd_h = d.array().inverse();
        }

        gu_h *= nu; gd_h *= nd;
        gu_p *= (1 - nu); gd_p *= (1 - nd);
        
        gu += expl(-Beta * E_m) * (gu_p + gu_h);
        gd += expl(-Beta * E_m) * (gd_p + gd_h);
    }

    gu /= Z; gd /= Z;

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

int SolverLanczos::get_block_idx(const VectorXd& psi) const {
    const double tol = 1e-6; 

    if (psi.norm() < tol)
        return -1;  

    for1(j, NFock) {
        if (std::abs(psi(j)) > tol) 
            return nud(j);
    }

    return -1;
}
