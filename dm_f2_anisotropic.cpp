#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <complex>
#include <limits>
#include <vector>
#include <valarray>
#include <cmath>

const int g_max = 20;
int n_k, n_band, n_valence, Ebins, qbins, qbin0;
double alat, tpiba, dq, dE;

std::vector< std::valarray<double> > read_klist(std::string fname);
std::vector< std::vector<double> >  read_Ek(std::string fname);
std::vector< std::vector< std::valarray<double> > >  read_glist(std::string fname);
std::vector< std::vector< std::vector<std::complex<double> > > > read_ui(std::string fname, const std::vector< std::vector< std::valarray<double> > >& g_list);
std::vector<std::vector< std::vector< std::vector<int> > > > build_index(const std::vector< std::vector< std::valarray<double> > >& g_list);
std::map<std::string,double> read_config(std::string fname,std::map<std::string,std::string>& datfile, std::map<std::string,std::valarray<double>>& basis,std::string& calc_name);


int main(int argc, char *argv[]){
    int ki=0;
    //int kpi=0;
    if(argc<3){
        std::cout << "Usage: dmf2 [Path to config file] [Index of k-point]" << std::endl;
        exit(0);
    }

    std::string calc_name = "f2";
    std::map<std::string,std::string> datfile = {
        {"Ek","Ek.dat"},
        {"klist","klist.dat"},
        {"glist","glist.dat"},
        {"ui","ui.dat"}
    };

    std::map<std::string,std::valarray<double>> basis = {
        {"b1",{1,0,0}},
        {"b2",{0,1,0}},
        {"b3",{0,0,1}}
    };
    
    auto config = read_config(argv[1],datfile,basis,calc_name);

    ki= atoi(argv[2]);

    n_k = (int) round(config["n_k"]);
    n_band = (int) round(config["n_band"]);
    n_valence = (int) round(config["n_valence"]);
    alat = config["alat"];
    tpiba = M_PI*2/alat;
    Ebins=(int) round(config["Ebins"]);
    qbins=(int) round(config["qbins"]);
    qbin0 = (qbins-1)/2;
    dq=config["dq"];
    dE=config["dE"];    

    auto k_list = read_klist(datfile["klist"]);
    auto E_ki   = read_Ek(datfile["Ek"]);
    auto g_list = read_glist(datfile["glist"]);
    auto u_kgi = read_ui(datfile["ui"],g_list);
    auto g_index = build_index(g_list); // for fast lookup 
    
    std::vector<std::vector<std::vector<std::vector<double>>>> f_crystal(qbins,std::vector<std::vector<std::vector<double>>>(qbins,std::vector<std::vector<double>>(qbins,std::vector<double>(Ebins,0))));
    
    std::valarray<double> q_vec;
    double E;
    int qix,qiy,qiz,Ei;

    int gpi;
    std::complex<double> f,x;
    
    // for(int ki=0; ki<n_k; ++ki){
    for(int kpi=0; kpi<n_k; ++kpi){
        
        std::cout << std::to_string(kpi+1) << '/' << std::to_string(n_k) << std::endl;
    for(auto gp : g_list[kpi] ){
        q_vec = tpiba*(gp[0]*basis["b1"]+gp[1]*basis["b2"]+gp[2]*basis["b3"])+k_list[kpi]-k_list[ki];
        qix = int(q_vec[0]/dq);
        qiy = int(q_vec[1]/dq);
        qiz = int(q_vec[2]/dq);
        if( std::abs(qix) >= qbins/2 || std::abs(qiy) >= qbins/2 || std::abs(qiz) >=qbins/2 ){
            continue;
        }
    for(int i=0; i<n_valence; ++i){
    for(int ip=n_valence; ip<n_band; ++ip){
        E = (E_ki[kpi][ip]-E_ki[ki][i]);
        Ei = int(E/dE);
        if (Ei >=Ebins){
            continue;
        }


        f=0;
        for(int gi=0; gi<g_list[ki].size(); ++gi){
            auto g_sum = g_list[ki][gi] + gp;
            if((std::abs(g_sum)).max()>=g_max){
                continue;
            }
            gpi = g_index[kpi][int(g_sum[0])+g_max][int(g_sum[1])+g_max][int(g_sum[2])+g_max];
            if(gpi<0){
                continue;
            }
            f += u_kgi[ki][gi][i]*std::conj(u_kgi[kpi][gpi][ip]);            
        }      
        f_crystal[qbin0+qix][qbin0+qiy][qbin0+qiz][Ei] += std::pow( std::real(std::abs(f)) ,2);        
    }
    }
    }
    }
    // }
    
    // write to file
    std::string output_fname = calc_name+"_aniso_" + std::to_string(ki)+".dat";
    std::ofstream f_crystal_fs(output_fname);
    for(qix=0; qix<qbins;++qix){
        for(qiy=0; qiy<qbins;++qiy){
            for(qiz=0; qiz<qbins;++qiz){
                for(Ei=0;Ei<Ebins;++Ei){
                    f_crystal_fs << f_crystal[qix][qiy][qiz][Ei] << ' ';
                }
                f_crystal_fs << '\n';
            }
        }
    }
    f_crystal_fs.close();
    return 0;
}


std::vector<std::string> split(std::string text, char delim) {
    std::string line;
    std::vector<std::string> vec;
    std::stringstream ss(text);
    while(std::getline(ss, line, delim)) {
        vec.push_back(line);
    }
    return vec;
}

std::map<std::string,double> read_config(std::string fname,std::map<std::string,std::string>& datfile, std::map<std::string,std::valarray<double>>& basis,std::string& calc_name){
    std::map<std::string,double> output;
    std::ifstream config_fs(fname);

    // Config file format
    // name=silicon_k333
    // Ek=Ek.dat
    // klist=klist.dat
    // glist=glist.dat
    // ui=ui.dat
    // n_k=27
    // g_max=20
    // n_band=144
    // n_valence=40
    // alat=8.4059 
    // b1=1 0 0
    // b2=0 1 0
    // b3=0 0 1
    // dq=0.2
    // dE=0.1 
    // Ebins=50
    // qbins=51

    if(!config_fs.is_open()){
        std::cout << "Failed to open config file "<< fname << std::endl;
        exit(-1);
    }

    std::string line;
    while(std::getline(config_fs,line))
    {
        std::istringstream line_is(line);
        std::string key;
        if(std::getline(line_is,key,'=')){
            std::string value;
            if(std::getline(line_is,value)){
                if(key=="Ek" || key=="glist" || key=="klist"||key=="ui"){
                    datfile[key] = value;
                }else if(key=="b1" || key=="b2" || key=="b3"){
                    std::vector<std::string> vec = split(value,' ');
                    if(vec.size() < 3){
                        std::cout << "Basis vector not enough components" << std::endl;
                        exit(-1);
                    }
                    basis[key] = {std::stod(vec[0]),std::stod(vec[1]),std::stod(vec[2])};
                }else if(key=="name"){
                    calc_name = value;
                }else{
                    output[key]=std::stod(value);
                }                
            }            
        }
    }

    return output;
}

std::vector< std::valarray<double> > read_klist(std::string fname){
    std::vector< std::valarray<double> > k_list;    
    std::ifstream k_list_fs(fname);

    if(!k_list_fs.is_open()){
        std::cout << "Failed to open klist file" << fname << std::endl;
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
        std::cout << "Failed to open Ek file"<< fname << std::endl;
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
        std::cout << "Failed to open glist file"<<fname << std::endl;
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

std::vector< std::vector< std::vector<std::complex<double> > > > read_ui(std::string fname, const std::vector< std::vector< std::valarray<double> > >& g_list){
    std::vector< std::vector< std::vector<std::complex<double> > > > u_kgi;
    std::vector< std::vector<std::complex<double> > > u_gi;
    std::vector<std::complex<double> > u_i;
    std::complex<double> u;

    std::ifstream u_kgi_fs(fname);
    if(!u_kgi_fs.is_open()){
        std::cout << "Failed to open ui file" << fname << std::endl;
        exit(-1);
    }

    for(int ki=0; ki<n_k; ++ki){
        u_gi.clear();
        for(int gi=0; gi<g_list[ki].size(); ++gi){
            u_i.clear();
            for(int i=0; i<n_band; ++i){
                u_kgi_fs >> u;
                u_i.push_back(u);
            }
            u_gi.push_back(u_i);
        }
        u_kgi.push_back(u_gi);
    }
    return u_kgi;
}

std::vector<std::vector< std::vector< std::vector<int> > > > build_index(const std::vector< std::vector< std::valarray<double> > >& g_list){
    std::vector<std::vector< std::vector< std::vector<int> > > > g_index(n_k, std::vector< std::vector< std::vector<int> > >(2*g_max, std::vector< std::vector<int> >(2*g_max,std::vector<int>(2*g_max,-1))));
    
    int g_ind,gx,gy,gz;
    for(int ki=0; ki<n_k; ++ki){
        for(int gi=0; gi<g_list[ki].size();++gi){     
            gx = int(g_list[ki][gi][0]);
            gy = int(g_list[ki][gi][1]);
            gz = int(g_list[ki][gi][2]);
            if(abs(gx)>= g_max||abs(gy)>=g_max || abs(gz)>=g_max){
                continue;
            }    
            g_index[ki][int(g_list[ki][gi][0])+g_max][int(g_list[ki][gi][1])+g_max][int(g_list[ki][gi][2])+g_max] = gi;
        }
    }
    return g_index;
}
