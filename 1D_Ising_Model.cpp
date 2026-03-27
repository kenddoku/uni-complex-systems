#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <typeinfo>
#include <ctime>
#include <sstream>
#include <iomanip>

double magnetization(const int L, const std::vector<std::string> &chain, const std::map<std::string, double> &spin_to_val)
{
    double result {0};

    for(int i=0; i<L; i++)
    {
        result += spin_to_val.at(chain[i]);
    }

    return result;
}

double local_energy(const int L, const std::string config, const std::map<std::string, double> &spin_to_val)
{
    double result {0};

    for(size_t i=0; i<config.size()-1; ++i)
    {
        result -= spin_to_val.at(std::string(1, config[i])) * spin_to_val.at(std::string(1, config[(i+1)%L]));
    }

    return result;
}

double energy(const int L, const std::vector<std::string> &chain, const std::map<std::string, double> &spin_to_val)
{
    double result {0};

    for(int i=0; i<L; ++i)
    {
        result -= spin_to_val.at(chain[i]) * spin_to_val.at(chain[(i+1)%L]); // %L jest potrzebne, żeby ostatni wezel przemnożyc przez pierwszy wezel (periodyczne warunki)
    }

    return result;
}

int main()
{
    srand(time(NULL));

    std::map<std::string, double> spin_to_val {};
    spin_to_val["d"] = -1.;
    spin_to_val["u"] = +1.;

    std::map<double, std::string> val_to_spin {};
    val_to_spin[-1.] = "d";
    val_to_spin[+1.] = "u";

    std::vector<std::string> spin_L {"d", "u"};
    std::vector<std::string> spin_R {"d", "u"};
    std::vector<std::string> spin_C {"d", "u"};

    const int L {1000};
    const int MCS {(int)1E6};
    const int tau {(int)1E5};

    std::vector<double> over_T_val {0.5, 1.5, 2.5, 4.0};

    for(double over_T : over_T_val)
    {
        // *** Tworzenie łańcucha losowo ustawionych spinów
        std::vector<std::string> chain {};
        for(int i=0; i<L; ++i)
        {
            if(rand()/(double)RAND_MAX <= 0.5) chain.push_back("d");
            else chain.push_back("u");
        }

        // *** Tablicowanie prawdopodobienst zmiany spinu w zaleznosci od stanu i temperatury
        std::map<std::string, double> p{};

        for(std::string s_L : spin_L)
        {
            for(std::string s_R : spin_R)
            {
                for(std::string s_C : spin_C)
                {
                    std::string config {s_L + s_C + s_R};
                    double sigma_L {spin_to_val[s_L]};
                    double sigma_C {spin_to_val[s_C]};
                    double sigma_R {spin_to_val[s_R]};

                    double delta_E {(sigma_C - (-1)*sigma_C) * (sigma_L + sigma_R)};

                    p[config] = std::min(1., std::exp((-1)*delta_E*over_T));
                }
            }
        }

        // *** Wypisywanie prawdopodobienstw dla danych przejsc w zaleznosi od temperatury
        std::cout << " 1/T = " << over_T << " *****************************************************\n\n";
        for(const auto& pair : p)
        {
            std::string config = pair.first;
            std::string c {config[1]};
            std::string c_prime {};

            if(c == std::string("d")) c_prime = "u";
            else c_prime = "d";

            //std::cout << "Initial state: " << pair.first << ", p = " << pair.second << "\n";
            std::cout << pair.first << " --> " << config[0] + c_prime + config[2] << ", p = " << pair.second << "\n";
        }
        std::cout << "\n";

        std::ostringstream oss;
        oss << "magnet_T" << std::fixed << std::setprecision(1) << over_T << ".csv";
        std::string filename {oss.str()};
        oss.str("");
        oss.clear();

        std::ofstream file;
        file.open(filename);
        if(!file.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

        oss << "energy_T" << std::fixed << std::setprecision(1) << over_T << ".csv";
        filename = oss.str();

        std::ofstream file2;
        file2.open(filename);
        if(!file2.is_open()) std::cerr << "ERROR opening" << filename << std::endl;

        for(int t=0; t<tau; ++t)
        {
            int nod {rand() % L};
            std::string config {chain[(L + nod - 1) % L] + chain[(L + nod) % L] + chain[(L + nod + 1) % L]};

            if(rand()/(double)RAND_MAX <= p.at(config))
            {
                chain[nod] = val_to_spin.at((-1)*spin_to_val.at(chain[nod]));
            }

            if(t % 100 == 0)
            {
                file << magnetization(L, chain, spin_to_val) / L << "\n";
                file2 << energy(L, chain, spin_to_val) / L << "\n";
            }
        }

        file.close();
        file2.close();
    }

    // *** ZADANIE 3 ------------------------------------------------------------
    over_T_val = {};
    
    for(double T = 0.1; T<4.1; T+=0.1)
    {
        over_T_val.push_back(T);
    }

    std::string filename {"e.csv"};
    std::ofstream file;

    file.open(filename);
    if(!file.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

    filename = "c.csv";
    std::ofstream file2;

    file2.open(filename);
    if(!file2.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

    for(double over_T : over_T_val)
    {
        // *** Tworzenie łańcucha losowo ustawionych spinów
        std::vector<std::string> chain {};
        for(int i=0; i<L; ++i)
        {
            if(rand()/(double)RAND_MAX <= 0.5) chain.push_back("d");
            else chain.push_back("u");
        }

        // *** Tablicowanie prawdopodobienst zmiany spinu w zaleznosci od stanu i temperatury
        std::map<std::string, double> p{};

        for(std::string s_L : spin_L)
        {
            for(std::string s_R : spin_R)
            {
                for(std::string s_C : spin_C)
                {
                    std::string config {s_L + s_C + s_R};
                    double sigma_L {spin_to_val[s_L]};
                    double sigma_C {spin_to_val[s_C]};
                    double sigma_R {spin_to_val[s_R]};

                    double delta_E {(sigma_C - (-1)*sigma_C) * (sigma_L + sigma_R)};

                    p[config] = std::min(1., std::exp((-1)*delta_E*over_T));
                }
            }
        }

        double M {0};
        double E {0};
        double E2 {0};

        //double temp_magnet {};
        //double temp_energy {};
        double chain_E {};
        double loc_E {};
        double loc_E_prime {};

        for(int t=0; t<MCS; ++t)
        {
            int nod {rand() % L};
            std::string config {chain[(L + nod - 1) % L] + chain[(L + nod) % L] + chain[(L + nod + 1) % L]};
            
            if(t >= MCS - tau)
            {
                if(t == MCS - tau)
                {
                    chain_E = energy(L, chain, spin_to_val);
                }
                loc_E = local_energy(L, config, spin_to_val);
            }

            if(rand()/(double)RAND_MAX <= p.at(config))
            {
                chain[nod] = val_to_spin.at((-1)*spin_to_val.at(chain[nod]));
            }

            if(t >= MCS - tau)
            {
                config = chain[(L + nod - 1) % L] + chain[(L + nod) % L] + chain[(L + nod + 1) % L];
                loc_E_prime = local_energy(L, config, spin_to_val);
                chain_E = chain_E - loc_E + loc_E_prime;
                //temp_magnet = magnetization(L, chain, spin_to_val);
                //temp_energy = energy(L, chain, spin_to_val);
                //M += temp_magnet;
                E += chain_E;
                E2 += chain_E * chain_E;
            }
        }

        //double M_mean {M / tau};
        double E_mean {E / tau};
        double E2_mean {E2 / tau};

        //double m {M_mean / L};
        double e {E_mean / L};
        double c {(E2_mean - E_mean * E_mean) * over_T * over_T / L};

        file << over_T << ", " << e << "\n";
        file2 << over_T << ", " << c << "\n";

        std::cout << (over_T - over_T_val[0]) / (over_T_val[over_T_val.size() - 1] - over_T_val[0]) * 100 << "%" << std::endl;
    }

    file.close();
    file2.close();
    
    return 0;
}