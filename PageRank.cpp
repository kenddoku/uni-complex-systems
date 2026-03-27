#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>
#include <sstream>
#include <vector>

void display_2D(std::vector<std::vector<double>>& arr)
{
    int rows = arr.size();
    int cols = arr[0].size();

    for(int i=0; i<rows; ++i)
    {
        for(int j=0; j<cols-1; ++j)
        {
            std::cout << std::setw(10) << arr[i][j] << ", ";
        }
        std::cout << std::setw(10) << arr[i][cols-1] << '\n';
    }
}

void display_1D(std::vector<double>& vec)
{
    double sum = 0.0;

    for(double elem : vec)
    {
        sum += elem;
    }

    for(double elem : vec)
    {
        std::cout << std::setw(5) << elem / sum << ", ";
    }
    std::cout << '\n';
}

std::vector<double> vector_times_matrix(const std::vector<double>& vec, const std::vector<std::vector<double>>& M)
{
    int n = vec.size();            
    int m = M[0].size();           
    std::vector<double> result(m, 0.0);

    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            result[j] += vec[i] * M[i][j];
        }
    }

    return result;
}

class Network
{
    std::vector<std::vector<int>> adj_matrix0;
    std::vector<std::vector<int>> adj_matrix;
    std::vector<double> prob_visit; // vector holding info about probability of visiting n-th node (n=0, 1, 2, ..., N-1)
    std::vector<double> prob_visit0;
    int N;  // number of nodes
    double eps {1E-6};  // error toleration

public:
    Network(std::vector<std::vector<int>>& adj_matrix_)
    {
        adj_matrix0 = adj_matrix_;
        adj_matrix = adj_matrix_;

        prob_visit0 = std::vector<double>((int)adj_matrix.size(), 0);
        prob_visit = std::vector<double>((int)adj_matrix.size(), 0);

        N = (int)prob_visit.size();
    }

    void random_walk(std::mt19937& mt)
    {
        std::uniform_int_distribution<int> pos_gen(0, N-1);
        int n0 = pos_gen(mt);   // initial node i.e. node that we start random walk from
        double err {1.2*eps};   // value greater than eps initially to ensure at least one while loop iteration

        std::vector<double> prob_visit_old(N, 0);
        int nstep = 1;
        prob_visit_old[n0] = 1;
        prob_visit[n0] = 1;

        while(err > eps)
        {
            nstep++;

            std::vector<int> n1_vec {};
            for(int k=0; k<N; ++k)
            {
                if(adj_matrix[n0][k] == 1) n1_vec.push_back(k);
            }

            if(n1_vec.size() == 0) n1_vec.push_back(n0);

            std::uniform_int_distribution<int> node_index_gen(0, n1_vec.size()-1);
            int n1 = n1_vec[node_index_gen(mt)];
            prob_visit[n1]++;

            err = 0;
            for(int k=0; k<N; ++k)
            {
                err += std::abs(prob_visit[k]/nstep - prob_visit_old[k]/(nstep-1));
            }

            n0 = n1;
            prob_visit_old = prob_visit;
        }

        std::cout << "RANDOM WALK -----------------------------------------------------------------\n";
        display_1D(prob_visit);
        std::cout << "nsteps = " << nstep << "\nERROR = " << err << "\n\n";

        // Going back to initial conditions of the network
        adj_matrix = adj_matrix0;
        prob_visit = prob_visit0;
    }

    void random_walk_telep(std::mt19937& mt)
    {
        std::uniform_int_distribution<int> pos_gen(0, N-1);
        std::uniform_real_distribution<double> telep_gen(0.0, 1.0);
        double p_telep = 0.15; // probability of teleportation

        int n0 = pos_gen(mt);   // initial node i.e. node that we start random walk from
        double err {1.2*eps};   // value greater than eps initially to ensure at least one while loop iteration

        std::vector<double> prob_visit_old(N, 0);
        int nstep = 1;
        prob_visit_old[n0] = 1;
        prob_visit[n0] = 1;

        while(err > eps)
        {
            nstep++;

            int is_connected = 0;
            for(int i=0; i<N; ++i)
            {
                is_connected += adj_matrix[n0][i];
            }

            if(telep_gen(mt) < p_telep ||  is_connected==0)
            {
                int n1;
            
                n1 = pos_gen(mt);
                
                prob_visit[n1]++;

                err = 0;
                for(int k=0; k<N; ++k)
                {
                    err += std::abs(prob_visit[k]/nstep - prob_visit_old[k]/(nstep-1));
                }

                n0 = n1;
                prob_visit_old = prob_visit;
            }
            else 
            {
                std::vector<int> n1_vec {};
                for(int k=0; k<N; ++k)
                {
                    if(adj_matrix[n0][k] == 1) n1_vec.push_back(k);
                }

                std::uniform_int_distribution<int> node_index_gen(0, n1_vec.size()-1);
                int n1 = n1_vec[node_index_gen(mt)];
                prob_visit[n1]++;

                err = 0;
                for(int k=0; k<N; ++k)
                {
                    err += std::abs(prob_visit[k]/nstep - prob_visit_old[k]/(nstep-1));
                }

                n0 = n1;
                prob_visit_old = prob_visit;
            }
        }

        std::cout << "RANDOM WALK WITH TELEPORTATION ----------------------------------------------\n";
        display_1D(prob_visit);
        std::cout << "nsteps = " << nstep << "\nERROR = " << err << "\n\n";

        adj_matrix = adj_matrix0;
        prob_visit = prob_visit0;
    }

    void page_rank(std::mt19937& mt)
    {
        std::vector<std::vector<double>> A(N, std::vector<double>(N, 0)); // Matrix of probabilities of transition from i-th node to j-th node

        for(int i=0; i<N; ++i)
        {
            int nconnections {0};
            for(int j=0; j<N; ++j)
            {
                nconnections += adj_matrix[i][j];
            }

            if(nconnections == 0) nconnections = 1;

            for(int j=0; j<N; ++j)
            {
                A[i][j] = (double)adj_matrix[i][j] / (double)nconnections;
            }
        }

        for(int i=0; i<N; ++i)
        {
            prob_visit[i] = 1/(double)N;
        }

        std::vector<double> prob_visit_old = prob_visit;
        double err {1.2*N*eps};
        int nstep = 0;

        while(err > N * eps)
        {
            nstep++;

            prob_visit = vector_times_matrix(prob_visit, A);

            err = 0;
            for(int i=0; i<N; ++i)
            {
                err += std::abs(prob_visit[i] - prob_visit_old[i]);
            }

            prob_visit_old = prob_visit;
        }

        std::cout << "PAGE RANK -------------------------------------------------------------------\n";
        display_1D(prob_visit);
        std::cout << "nsteps = " << nstep << "\nERROR = " << err << "\n\n";

        adj_matrix = adj_matrix0;
        prob_visit = prob_visit0;
    }

    void page_rank_telep(std::mt19937& mt)
    {
        std::vector<std::vector<double>> A(N, std::vector<double>(N, 0)); // Matrix of probabilities of transition from i-th node to j-th node
        std::vector<std::vector<double>> B(N, std::vector<double>(N, 1/(double)N ));
        std::vector<std::vector<double>> M(N, std::vector<double>(N, 0));
        double p_telep {0.15};

        for(int i=0; i<N; ++i)
        {
            int nconnections {0};
            for(int j=0; j<N; ++j)
            {
                nconnections += adj_matrix[i][j];
            }

            if(nconnections == 0) nconnections = 1;

            for(int j=0; j<N; ++j)
            {
                A[i][j] = (double)adj_matrix[i][j] / (double)nconnections;
            }
        }

        for(int i=0; i<N; ++i)
        {
            for(int j=0; j<N; ++j)
            {
                M[i][j] = (1-p_telep)*A[i][j] + p_telep*B[i][j];
            }
        }

        for(int i=0; i<N; ++i)
        {
            prob_visit[i] = 1/(double)N;
        }

        std::vector<double> prob_visit_old = prob_visit;
        double err {1.2*N*eps};
        int nstep = 0;

        while(err > N * eps)
        {
            nstep++;

            prob_visit = vector_times_matrix(prob_visit, M);

            err = 0;
            for(int i=0; i<N; ++i)
            {
                err += std::abs(prob_visit[i] - prob_visit_old[i]);
            }

            prob_visit_old = prob_visit;
        }

        std::cout << "PAGE RANK TELEPORTATION -----------------------------------------------------\n";
        display_1D(prob_visit);
        std::cout << "nsteps = " << nstep << "\nERROR = " << err << "\n\n";

        adj_matrix = adj_matrix0;
        prob_visit = prob_visit0;
    }
};

// *** MAIN FUNCTION ------------------------------------------------------------------------------
int main()
{
    std::random_device rd;
    std::mt19937 mt(rd());

    // adjacency matrix for the first example network
    std::vector<std::vector<int>> adj_matrix1 {
        {0, 1, 1, 1},
        {1, 0, 0, 0},
        {1, 1, 0, 0},
        {0, 1, 1, 0}
    };

    // adjacency matrix for the second example network
    std::vector<std::vector<int>> adj_matrix2 {
        {0, 1, 0, 0},
        {1, 0, 0, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0}
    };

    // adjacency matrix for the third example network
    std::vector<std::vector<int>> adj_matrix3 {
        {0, 1, 0, 0},
        {0, 0, 0, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0}
    };

    // Creating three networks --------------------------------------------------------------------
    Network N1(adj_matrix1);
    Network N2(adj_matrix2);
    Network N3(adj_matrix3);

    // *** EXERCISE 1 -----------------------------------------------------------------------------
    N1.random_walk(mt);
    N2.random_walk(mt);
    N3.random_walk(mt);

    // *** EXERCISE 2 -----------------------------------------------------------------------------
    N1.random_walk_telep(mt);
    N2.random_walk_telep(mt);
    N3.random_walk_telep(mt);

    // *** EXERCISE 3 -----------------------------------------------------------------------------
    N1.page_rank(mt);
    N2.page_rank(mt);
    N3.page_rank(mt);

    // *** EXERCISE 4 -----------------------------------------------------------------------------
    N1.page_rank_telep(mt);
    N2.page_rank_telep(mt);
    N3.page_rank_telep(mt);

    return 0;
}