#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

int sign(int a)
{
    return (a >= 0) - (a < 0);
}

void show_pattern(std::vector<int> patt)
{
    int n = std::sqrt((int)patt.size());
    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            if(patt[i*n+j] == -1)
                std::cout << "\033[34m#\033[0m ";
            else
                std::cout << "\033[31m#\033[0m ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

void print_weights(const std::vector<std::vector<int>>& W)
{
    int rows = (int)W.size();
    int cols = (int)W[0].size();

    for(int i=0; i<rows; ++i)
    {
        for(int j=0; j<cols; ++j)
        {
            std::cout << std::setw(2) << W[i][j] << " ";
        }
        std::cout << '\n';
    }
}

std::vector<int> remember(std::vector<int> input, const std::vector<std::vector<int>>& W)
{
    int n = (int)input.size();
    std::vector<int> output = input;

    bool to_continue = true;
    while(to_continue)
    {
        for(int i=0; i<n; ++i)
        {
            int val {0};
            for(int j=0; j<n; ++j)
            {
                val += W[i][j] * output[j];
            }
            output[i] = sign(val);
        }

        int is_equal {0};
        for(int i=0; i<n; ++i)
        {
            if(output[i] == input[i]) is_equal++;
        }
        if(is_equal == n) to_continue = false;

        for(int i=0; i<n; ++i)
        {
            input[i] = output[i];
        }
    }   

    return output;
}

int main()
{
    int N {25};
    std::vector<std::vector<int>> Patterns(3);

    std::vector<int> pattern1 { // T
        +1, +1, +1, +1, +1,
        -1, -1, +1, -1, -1,
        -1, -1, +1, -1, -1,
        -1, -1, +1, -1, -1,
        -1, -1, +1, -1, -1
    };
    Patterns[0] = pattern1;

    std::vector<int> pattern2 { // H
        +1, -1, -1, -1, +1,
        +1, -1, -1, -1, +1,
        +1, +1, +1, +1, +1,
        +1, -1, -1, -1, +1,
        +1, -1, -1, -1, +1
    };
    Patterns[1] = pattern2;

    std::vector<int> pattern3 { // A
        -1, -1, +1, -1, -1,
        -1, +1, -1, +1, -1,
        -1, +1, -1, +1, -1,
        +1, +1, +1, +1, +1,
        +1, -1, -1, -1, +1
    };
    Patterns[2] = pattern3;

    std::cout << "Initial training set:\n";
    for(auto pattern : Patterns)
    {
        show_pattern(pattern);
    }

    std::cout << "---------------------\n\n";

    // Training
    std::vector<std::vector<int>> W(N, std::vector<int>(N, 0));

    for(auto pattern : Patterns)
    {
        for(int i=0; i<N; ++i)
        {
            for(int j=0; j<N; ++j)
            {
                if(i==j) continue;
                W[i][j] += pattern[i] * pattern[j];
            }
        }
    }

    //print_weights(W);

    // Creating distorted patterns to test

    std::vector<int> distorted_pattern1 { // distorted T
        +1, +1, +1, -1, +1,
        -1, -1, +1, -1, -1,
        -1, -1, -1, -1, -1,
        +1, -1, -1, -1, -1,
        -1, -1, +1, -1, -1
    };

    std::vector<int> distorted_pattern2 { // distorted H
        +1, -1, -1, -1, +1,
        +1, -1, -1, -1, +1,
        +1, -1, +1, +1, -1,
        +1, -1, -1, -1, -1,
        +1, -1, -1, -1, +1
    };

    std::vector<int> distorted_pattern3 { // distorted A
        -1, -1, -1, -1, -1,
        -1, +1, -1, +1, -1,
        -1, +1, -1, +1, -1,
        -1, -1, -1, -1, -1,
        +1, -1, -1, -1, +1
    };

    std::cout << "Testing -------------\n\n";

    std::cout << "Distorted T:\n";
    show_pattern(distorted_pattern1);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern1, W));

    std::cout << "Distorted H:\n";
    show_pattern(distorted_pattern2);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern2, W));

    std::cout << "Distorted A:\n";
    show_pattern(distorted_pattern3);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern3, W));

    // Adding 4th pattern
    std::vector<int> pattern4 {
        +1, +1, +1, +1, +1,
        +1, -1, -1, -1, -1,
        +1, +1, +1, +1, -1,
        +1, -1, -1, -1, -1,
        +1, +1, +1, +1, +1,
    };
    Patterns.push_back(pattern4);

    // Adding 5th pattern
    std::vector<int> pattern5 {
        +1, -1, -1, -1, +1,
        -1, +1, -1, +1, -1,
        -1, -1, +1, -1, -1,
        -1, +1, -1, +1, -1,
        +1, -1, -1, -1, +1
    };
    Patterns.push_back(pattern5);

    // Reseting weights matrix
    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            W[i][j] = 0;
        }
    }

    for(auto pattern : Patterns)
    {
        for(int i=0; i<N; ++i)
        {
            for(int j=0; j<N; ++j)
            {
                if(i==j) continue;
                W[i][j] += pattern[i] * pattern[j];
            }
        }
    }

    //print_weights(W);

    std::cout << "Testing after adding the 4th and 5th pattern -----\n\n";
    std::cout << "Pattern 4:\n";
    show_pattern(pattern4);

    std::cout << "Pattern 5:\n";
    show_pattern(pattern5);

    std::cout << "Distorted T:\n";
    show_pattern(distorted_pattern1);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern1, W));

    std::cout << "Distorted H:\n";
    show_pattern(distorted_pattern2);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern2, W));

    std::cout << "Distorted A:\n";
    show_pattern(distorted_pattern3);
    std::cout << "Retrieved pattern:\n";
    show_pattern(remember(distorted_pattern3, W));


}