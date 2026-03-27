#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <cmath>
#include <random>
#include <boost/numeric/odeint.hpp>

double PI = std::numbers::pi;

void omega_initialization(std::vector<double>& omega_arr, std::mt19937& mt)
{
    std::normal_distribution<double> normDist(0, 0.5);

    for(size_t i=0; i<omega_arr.size(); ++i)
    {
        omega_arr[i] = normDist(mt);
    }
}

void phase_initialization(std::vector<double>& phi_arr, std::mt19937& mt)
{
    std::uniform_real_distribution<double> uniDist(0, 2*PI);

    for(size_t i=0; i<phi_arr.size(); ++i)
    {
        phi_arr[i] = uniDist(mt);
    }
}

// Function that calculates r and psi (phase average)
// std::pair<double, double> order_parameter(std::vector<double>& phi_arr)
// {
//     std::complex<double> sum(0.0, 0.0);

//     for(double phi : phi_arr)
//     {
//         sum += std::polar(1.0, phi);
//     }

//     sum /= (double)phi_arr.size();
//     double r = std::abs(sum);
//     double psi = std::arg(sum);

//     return {r, psi};
// }

class Kuramoto
{
    int N;
    double K;
    std::vector<double> omega_arr;

public:
    Kuramoto(int N_, double K_, std::vector<double>& omega_arr_) : N(N_), K(K_), omega_arr(omega_arr_) {}

    std::pair<double, double> order_parameter(const std::vector<double>& phi_arr)
    {
        std::complex<double> sum(0.0, 0.0);

        for(double phi : phi_arr)
        {
            sum += std::polar(1.0, phi);
        }

        sum /= (double)phi_arr.size();
        double r = std::abs(sum);
        double psi = std::arg(sum);

        return {r, psi};
    }

    void operator() (const std::vector<double>& phi_arr, std::vector<double>& dphi_dt, double t)
    {
        auto [r, psi] = order_parameter(phi_arr);

        for(int i=0; i<N; ++i)
        {
            dphi_dt[i] = omega_arr[i] + K*r*std::sin(psi - phi_arr[i]);
        }
    }
};

class Observer 
{
    std::ofstream file;

public:
    Observer(const std::string& filename)
    {
        file.open(filename);
    }

    void operator() (const std::vector<double>& phi_arr, double t)
    {
        file << t;
        for(double phi : phi_arr)
        {
            phi = fmod(phi, 2*PI);
            if(phi < 0) phi += 2*PI;
            file << ", " << phi;
        }
        file << '\n';
    }

    ~Observer()
    {
        if(file.is_open()) file.close();
    }
};

int main()
{
    std::random_device rd;
    std::mt19937 mt(rd());

    double t0 = 0.0;        // Initial time of the simulation
    double t1 = 1000.0;     // Final time of the simulation
    double dt = 0.01;       // Time step

    std::vector<int> N_arr {10, 20, 50};
    std::vector<double> K_arr {};

    for(double k=0.25; k<5; k+=0.5)
    {
        K_arr.push_back(k);
    }

    for(int N : N_arr)
    {
        for(double K : K_arr)
        {
            std::vector<double> omega_arr(N, 0);
            std::vector<double> phi_arr(N, 0);

            omega_initialization(omega_arr, mt);
            phase_initialization(phi_arr, mt);

            std::ostringstream oss;
            oss << "output_N" << N << "_K" << std::fixed << std::setprecision(2) << K << ".csv";
            std::string filename = oss.str();
            // oss.str("");
            // oss.clear();
            //std::ofstream file(filename);
            Observer obs(filename);

            Kuramoto kmoto(N, K, omega_arr);
            boost::numeric::odeint::integrate(kmoto, phi_arr, t0, t1, dt, std::ref(obs));

            //file.close();
        }
    }

    return 0;
}