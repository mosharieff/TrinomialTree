# OptionsCalculator
In this repository there are three different C++ scripts regarding Option Pricing 

## Sources
I got the idea to create this program from the paper "Pricing Options Using Trinomial Trees" by Paul Clifford, Yan Wang, and Oleg Zaboronski.

The paper is located at https://warwick.ac.uk/fac/sci/maths/people/staff/oleg_zaboronski/fm/trinomial_tree_2009.pdf

## Compiling and Running

### Compiling

```sh
g++ -o tree tree.cpp -std=c++11
g++ -o iv iv.cpp -std=c++11
g++ -o greeks greeks.cpp -std=c++11
```
### Running
```sh
./tree
./iv
./greeks
```

## Programs Included

### Trinomial Tree
In the file labeled as ```tree.cpp```, generates a Trinomial Tree based on various inputs into the model and calculates the call or put option price. The tree simulates a stock price rising, staying the same, or falling in a set of steps. The payoff is then calculated at the last column. Finally the payoffs are multiplied by probabilites and discounted backwards until the option price is computed.

### Implied Volatility Calculator
In the file labeled as ```iv.cpp```, the Implied Volatility is computed from a given set of inputs including a sample option price. The Implied Volatility is then checked to see if it is correct by being inputted into the standard trinomial tree function.

### Options Greeks
In the file labeled as ```greeks.cpp```, the greeks for the option are calculated using the trinomial tree. The greeks include Delta, Gamma, Theta, Vega, & Rho. Each value is incremented by a percent or day or dollar and is calculated accordingly.

