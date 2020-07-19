#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <limits>
#include <vector>
#include <valarray>
#include <cmath>

// Change these values according to your QE calculation
const int n_k = 27;
const int n_band = 144;
const int n_valence = 40;
const double alat = 8.4059; // a.u.
const double tpiba=M_PI*2/alat;
const std::valarray<double> b1 = {1,0,0};
const std::valarray<double> b2 = {0,1,0};
const std::valarray<double> b3 = {0,0,1};

// binning settings
const double dq = 0.02;
const double dE = 0.1; //eV
const int Ebins=500;
const int qbins=900;

// G cutoff
const int gcutoff = 2;
// -----------------------------------------------------


std::vector< std::valarray<double> > read_klist(std::string fname);
std::vector< std::vector<double> >  read_Ek(std::string fname);
std::vector< std::vector< std::valarray<double> > >  read_glist(std::string fname);
std::vector<std::vector< std::vector< std::vector<int> > > > build_index(const std::vector< std::vector< std::valarray<double> > >& g_list);

int main(int argc, char *argv[]){
    int ki=0;
    int kpi=0;
    if(argc>1){
        ki= atoi(argv[1]);
        //kpi= atoi(argv[2]);
    }


    auto k_list = read_klist("klist.dat");
    auto E_ki   = read_Ek("Ek.dat");
    auto g_list = read_glist("glist.dat");

    
    std::vector<std::vector<double>> f_crystal(qbins,std::vector<double>(Ebins,0));

    std::valarray<double> q_vec;
    double q,E;
    int qi,Ei;

    // for(int ki=0; ki<n_k; ++ki){
    for(int kpi=0; kpi<n_k; ++kpi){
    std::cout<< ki << kpi << std::endl;
    for(int i=0; i<n_valence; ++i){
    for(int ip=n_valence; ip<n_band; ++ip){
    for( auto gp : g_list[kpi] ){
    	if (std::abs(gp[0]) > gcutoff || std::abs(gp[1]) > gcutoff || std::abs(gp[2])> gcutoff){
            continue;
        }
        q_vec = k_list[kpi]-k_list[ki]+(gp[0]*b1+gp[1]*b2+gp[2]*b3);
        q = tpiba*sqrt( (q_vec*q_vec).sum() );
        
        E = (E_ki[kpi][ip]-Escissor-E_ki[ki][i]);

        qi = int(q/dq);
        Ei = int(E/dE);

        if(qi >= qbins || Ei >=Ebins){
            continue;
        }
      
        f_crystal[qi][Ei] += 1;
    }
    }
    }
    }
    // }

    // write to file
    std::string output_fname = "jdos_" + std::to_string(ki)+".dat";
    std::ofstream f_crystal_fs(output_fname);
    for(qi=0; qi<qbins;++qi){
        for(Ei=0;Ei<Ebins;++Ei){
            f_crystal_fs << f_crystal[qi][Ei] << ' ';
        }
        f_crystal_fs << '\n';
    }

    return 0;
}


std::vector< std::valarray<double> > read_klist(std::string fname){
    std::vector< std::valarray<double> > k_list;    
    std::ifstream k_list_fs(fname);

    if(!k_list_fs.is_open()){
        std::cout << "Failed to open klist file!" << std::endl;
        exit(-1);
    }

    std::valarray<double> k_vec (3);
    for(int ki = 0; ki<n_k; ++ki){        
        k_list_fs >> k_vec[0] >> k_vec[1] >> k_vec[2];        
        k_list.push_back(k_vec);
    }
    return k_list;
}

std::vector< std::vector<double> >  read_Ek(std::string fname){
    std::vector< std::vector<double> > E_ki;    
    std::ifstream E_ki_fs(fname);

    if(!E_ki_fs.is_open()){
        std::cout << "Failed to open Ek file!" << std::endl;
        exit(-1);
    }

    std::vector<double> E_i;
    double E;
    for(int ki = 0; ki<n_k; ++ki){ 
        E_i.clear();
        for(int i=0; i<n_band; ++i){
            E_ki_fs >>  E;   
            E_i.push_back(E);
        }
        E_ki.push_back(E_i);
    }
    return E_ki;
}

std::vector< std::vector< std::valarray<double> > >  read_glist(std::string fname){
    std::vector< std::vector< std::valarray<double> > > g_list;
    std::vector< std::valarray<double> > g_list_k;
    std::ifstream g_list_fs(fname);

    if(!g_list_fs.is_open()){
        std::cout << "Failed to open glist file!" << std::endl;
        exit(-1);
    }

    std::valarray<double> g_vec (3);
    
    double n_g;
    for(int ki = 0; ki<n_k; ++ki){        
        g_list_fs >> n_g;
        g_list_k.clear();
        for(int gi=0; gi<n_g; ++gi){
            g_list_fs >> g_vec[0] >> g_vec[1] >> g_vec[2];
            g_list_k.push_back(g_vec);
        }
        g_list.push_back(g_list_k);
    }
    return g_list;
}

