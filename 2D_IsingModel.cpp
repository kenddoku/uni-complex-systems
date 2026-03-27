#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <random>

#define N 100 // Rozmiar siatki w jendym wymiarze
#define MCS 100000 // Liczba kroków procedury MonteCarlo
#define tau 10000 // Liczba kroków MC po których przeprowadzone będzie uśrednianie

int vicinity_indx(const int xi, const int yi, const std::vector<std::vector<int>>& lattice)
{
    int result {0};
    int temp_x {xi + N};
    int temp_y {yi + N};

    result += lattice[(temp_x-1)%N][yi] + lattice[(temp_x+1)%N][yi] + lattice[xi][(temp_y-1)%N] + lattice[xi][(temp_y+1)%N];

    switch(result)
    {
        case +4:
            return 0;
        case +2:
            return 1;
        case +0:
            return 2;
        case -2:
            return 3;
        case -4:
            return 4;
        default:
            std::cerr << "Error in vicinity_indx() function" << std::endl;
            return 0;
    }
}

double magnetization(const std::vector<std::vector<int>>& lattice)
{
    double result {0};
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            result += lattice[i][j];
        }
    }

    return result;
}

void MC_step(std::vector<std::vector<int>>& lattice, const std::vector<std::vector<double>>& p, std::mt19937& gen)
{
    std::uniform_int_distribution<int> coord_dist(0, N-1);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    for(int t=0; t<N*N; ++t)
    {
        int xi {coord_dist(gen)};
        int yi {coord_dist(gen)};
        int spin_indx {lattice[xi][yi]==-1 ? 0 : 1};

        if(prob_dist(gen) <= p[spin_indx][vicinity_indx(xi, yi, lattice)])
        {
            lattice[xi][yi] *= (-1);
        }
    }
}

// * NADAJMY INDEKSY 'a' i 'b' KONKRETNYM STANOM JAKIE MOŻEMY NAPOTKAĆ
// * Jeżeli spin jest do dołu (góry) to a = 0 (a = 1)
// * Jeżeli suma sąsiadów +4, to b = 0
// * Jeżeli suma sąsiadów +2, to b = 1
// * Jeżeli suma sąsiadów +0, to b = 2
// * Jeżeli suma sąsiadów -2, to b = 3
// * Jeżeli suma sąsiadów -4, to b = 4
// * Przy tych oznaczeniach p-stwo zmiany stanu to p[a][b]
void prob_calc(const double T, std::vector<std::vector<double>>& p)
{
    std::vector<int> spin_c {-1, 1};
    std::vector<int> vicinity {+4, +2, +0, -2, -4};
    
    for(int a=0; a<2; a++)
    {
        for(int b=0; b<5; b++)
        {
            double delta_E {};
            delta_E = (-1)*((-1)*spin_c[a] * vicinity[b]) - (-1)*(spin_c[a] * vicinity[b]);

            p[a][b] = std::min(1., std::exp((-1)*delta_E / T));
        }
    }
}

void prob_print(const double T, const std::vector<std::vector<double>>& p)
{
    std::vector<int> spin_c {-1, 1};
    std::vector<int> vicinity {+4, +2, +0, -2, -4};

    std::cout << "T = " << T << "\n";
    for(int a=0; a<2; a++)
    {
        for(int b=0; b<5; b++)
        {
            std::cout << "SPIN = " << spin_c[a] << ",\tVICINITY = " << vicinity[b] << ",\tPROB = " << p[a][b] << "\n";
        }
    }
    std::cout << "\n";
}

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> T_val {0.5, 2.5, 4.0};

    // *** EXERCISE 1 & 2 -----------------------------------------------------------------------------

    for(double T : T_val)
    {
        std::vector<std::vector<double>> p(2, std::vector<double>(5, 0)); // Tablica prawdopodobienstw przejscia

        prob_calc(T, p);
        prob_print(T, p);

        std::vector<std::vector<int>> lattice(N, std::vector<int>(N, 1)); // Tworzymy siatke spinów do góry (1)

        std::string filename;
        std::ostringstream oss;
        oss << "mag_dist_T" << std::fixed << std::setprecision(1) << T << ".txt";
        filename = oss.str();
        oss.str("");
        oss.clear();
        std::ofstream file;
        file.open(filename);
        if(!file.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

        oss << "mag_evolution_T" << std::fixed << std::setprecision(1) << T << ".txt";
        filename = oss.str();
        oss.str("");
        oss.clear();
        std::ofstream file_mag;
        file_mag.open(filename);
        if(!file_mag.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

        double mag {};

        for(int t=0; t<MCS; ++t)
        {
            MC_step(lattice, p, gen);
            mag = magnetization(lattice) / N / N;
            file_mag << t << "\t" << mag << "\n";
        }

        for(int xi=0; xi<N; ++xi)
        {
            for(int yi=0; yi<N; ++yi)
            {
                file << xi << "\t" << yi << "\t" << lattice[xi][yi] << "\n";
            }
        }

        file.close();
        file_mag.close();
    }

    // *** EXERCISE 4 -----------------------------------------------------------------------------

    std::string filename;
    std::ostringstream oss;
    oss << "T_mag_chi.txt";
    filename = oss.str();
    oss.str("");
    oss.clear();

    std::ofstream file;
    file.open(filename);
    if(!file.is_open()) std::cerr << "ERROR opening " << filename << std::endl;

    for(double T = 1.4; T<=3.0; T+=0.05)
    {
        std::vector<std::vector<double>> p(2, std::vector<double>(5, 0)); // Tablica prawdopodobienstw przejscia

        prob_calc(T, p);
        prob_print(T, p);

        std::vector<std::vector<int>> lattice(N, std::vector<int>(N, 1)); // Tworzymy siatke spinów do góry (1)

        double mag_average {0}; // Zmienna przechowująca wartośc średniej magnetyzacji z ostatnich tau kroków
        double mag2_average {0}; // Zmienna przechowująca wartośc średniej kwadratu magnetyzacji z ostatnich tau kroków
        double mag_temp {0};

        for(int t=0; t<MCS; ++t)
        {
            MC_step(lattice, p, gen);

            if(t >= MCS - tau)
            {
                mag_temp = magnetization(lattice);
                mag_average += mag_temp;
                mag2_average += mag_temp*mag_temp;
            }
        }

        mag_average = mag_average / tau;
        mag2_average = mag2_average / tau;

        double m_T  {mag_average / N / N};
        double chi {(1/T)*(1/N/N) * (mag2_average - mag_average*mag_average)}; // Zmienna przechowująca wartośc podatnosci magnetycznej

        file << T << "\t" << m_T << "\t" << chi << "\n";
    }

    file.close();
    return 0;
}