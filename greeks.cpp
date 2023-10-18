#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>

// Creates the tree array in which the option price is calculated
std::vector<std::vector<double>> BOX(int rows, int cols)
{
    std::vector<std::vector<double>> res;
    std::vector<double> temp;
    for(int i = 0; i < rows; ++i){
        temp.clear();
        for(int j = 0; j < cols; ++j){
            temp.push_back(0.0);
        }
        res.push_back(temp);
    }
    return res;
}

// Prints the tree
void PRINTM(std::vector<std::vector<double>> x)
{
    std::cout << "\nOPTIONS TRINOMIAL TREE\n" << std::endl;

    for(auto & i : x){
        for(auto & j : i){
            std::cout << std::setprecision(3) << j << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Computes the factor in which the stationary price changes
auto m = []()
{
    return 1.0;
};

// Computes the factor in which the price goes down
auto d = [](double u)
{
    return 1.0 / u;
};

// Computes the factor in which the price goes up
auto u = [](double v, double dt)
{
    return exp(v*sqrt(2.0*dt));
};

// Probability of price down
auto pd = [](double r, double q, double v, double dt)
{
    double a = exp((r-q)*dt/2.0);
    double e0 = exp(-v*sqrt(dt/2.0));
    double e1 = exp(v*sqrt(dt/2.0));
    return pow((e1 - a)/(e1 - e0), 2);
};

// Probability of price up
auto pu = [](double r, double q, double v, double dt)
{
    double a = exp((r-q)*dt/2.0);
    double e0 = exp(-v*sqrt(dt/2.0));
    double e1 = exp(v*sqrt(dt/2.0));
    return pow((a - e0)/(e1 - e0), 2);
};

// Probability of price staying the same
auto pm = [](double PU, double PD)
{
    return 1.0 - (PU + PD);
};

// Option Tree
double PRICE(double S, double K, double r, double q, double v, double t, int nodes, std::string optype)
{
    double dt = t / (double) nodes;
    // Calculate all functions
    double U = u(v, dt);
    double D = d(U);
    double M = m();

    double PU = pu(r, q, v, dt);
    double PD = pd(r, q, v, dt);
    double PM = pm(PU, PD);
    
    // Set parameters for grid
    int rows = 4*nodes + 2;
    int cols = nodes + 1;
    int split = rows / 2 - 1;
    int srows = rows / 2;

    double option_price = 0;
    std::vector<std::vector<double>> Tree = BOX(rows, cols);

    // Set the initial stock price into the tree vector
    Tree[split][0] = S;

    // Forward Propigation
    for(int j = 0; j < cols; ++j){
        for(int i = 1; i < cols - j; ++i){
            // Upper price simulation
            Tree[split - i*2][j + i] = Tree[split - (i-1)*2][j + i-1]*U;
            Tree[split][j + i] = Tree[split][j + i-1]*M;
            Tree[split + i*2][j + i] = Tree[split + (i-1)*2][j + i-1]*D;
        }
    }

    // Payoff Computation
    int c = 1; // Row index
    for(int i = 0; i < srows; ++i){
        if(optype == "call"){
            // Call option payoff
            Tree[c][cols - 1] = fmax(Tree[c - 1][cols - 1] - K, 0);
        } else {
            // Put option payoff
            Tree[c][cols - 1] = fmax(K - Tree[c - 1][cols - 1], 0);
        }
        c += 2;
    }

    // Backpropigation to find option price
    int f = 3;
    double A, B, C;
    double E = exp(-r*dt);
    for(int i = cols - 2; i >= 0; --i){
        for(int j = f; j < rows - f + 1; ++j){
            if(j % 2 != 0){
                A = Tree[j - 2][i + 1];
                B = Tree[j][i + 1];
                C = Tree[j + 2][i + 1];
                if(optype == "call"){
                    Tree[j][i] = fmax(E*(A*PU + B*PM + C*PD), Tree[j - 1][i] - K);
                } else {
                    Tree[j][i] = fmax(E*(A*PU + B*PM + C*PD), K - Tree[j - 1][i]);
                }
            }
        }
        f += 2;
    }

    
    return Tree[split + 1][0];
}

int main()
{
    // Define Stock Inputs
    double S = 100;
    double K = 105; 
    double r = 0.05;
    double q = 0.025;
    double v = 0.20;
    double t = 15.0/365.0;
    std::string optype = "call";
    int nodes = 30;



    std::cout << "\nStock Price: " << S << std::endl;
    std::cout << "Strike Price: " << K << std::endl;
    std::cout << "RiskFree Rate: " << r << std::endl;
    std::cout << "Dividend Yield: " << q << std::endl;
    std::cout << "Volatility: " << v << std::endl;
    std::cout << "Expiration: " << t << std::endl;
    std::cout << "Option Type: " << optype << std::endl;
   
    // Declare greek variables
    double delta, gamma, theta, vega, rho, volga, vanna;

    // Calculate delta
    double dS = 0.01*S;
    delta = (PRICE(S+dS, K, r, q, v, t, nodes, optype) - PRICE(S-dS, K, r, q, v, t, nodes, optype))/(2.0*dS);
    
    // Calculate gamma
    gamma = (PRICE(S+dS, K, r, q, v, t, nodes, optype) - 2.0*PRICE(S, K, r, q, v, t, nodes, optype) + PRICE(S-dS, K, r, q, v, t, nodes, optype))/pow(dS, 2);
    
    // Calculate theta
    double dT = 1.0/365.0;
    theta = -(PRICE(S, K, r, q, v, t+dT, nodes, optype) - PRICE(S, K, r, q, v, t, nodes, optype))/dT;

    // Calculate vega
    double dV = 0.01;
    vega = (PRICE(S, K, r, q, v+dV, t, nodes, optype) - PRICE(S, K, r, q, v-dV, t, nodes, optype))/(2.0*dV);
    vega /= 100;

    // Calculate rho
    double dR = 0.01;
    rho = (PRICE(S, K, r+dR, q, v, t, nodes, optype) - PRICE(S, K, r-dR, q, v, t, nodes, optype))/(2.0*dR);
    rho /= 100;

    // Calculate volga
    volga = (PRICE(S, K, r, q, v+dV, t, nodes, optype) - 2*PRICE(S, K, r, q, v, t, nodes, optype) + PRICE(S, K, r, q, v-dV, t, nodes, optype))/pow(dV, 2);
    volga /= 100;

    // Calculate vanna
    vanna = (PRICE(S+dS, K, r, q, v+dV, t, nodes, optype) - PRICE(S+dS, K, r, q, v, t, nodes, optype) - PRICE(S, K, r, q, v+dV, t, nodes, optype) + PRICE(S, K, r, q, v, t, nodes, optype))/(dS*dV);
    vanna /= 100;

    std::cout << std::endl;
    std::cout << "Delta: " << delta << std::endl;
    std::cout << "Gamma: " << gamma << std::endl;
    std::cout << "Theta: " << theta << std::endl;
    std::cout << "Vega: " << vega << std::endl;
    std::cout << "Volga: " << volga << std::endl;
    std::cout << "Vanna: " << vanna << std::endl;
    std::cout << "Rho: " << rho << std::endl;
    std::cout << std::endl;

    return 0;
}

