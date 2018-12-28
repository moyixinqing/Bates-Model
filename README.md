# Jump-Process #

## Build status ##

[![Build Status](http://img.shields.io/travis/badges/badgerbadgerbadger.svg?style=flat-square)](https://github.com/moyixinqing/Greeks)

## Features ##
- [Numpy](http://www.numpy.org) 
- [SciPy](https://www.scipy.org) in Python

## Topics ##
### Bates_1996_Table_1_reproduction ###
Souece Code: [Bates_1996_Table_1_reproduction](https://github.com/moyixinqing/Jump-Process/tree/master/Bates_1996_Table_1_reproduction)

Table 1 of Bates (1996), first entries in each row only <br/>
European prices <br/>

| Set | K=38   | K=39   | K=40   | K=41   | K=42   |
|-----|--------|--------|--------|--------|--------|
| 0   | 0.3744 | 0.6617 | 1.074  | 1.6175 | 2.2832 |
| 1   | 0.575  | 0.902  | 1.3344 | 1.874  | 2.5147 |
| 2   | 0.3693 | 0.6485 | 1.0564 | 1.6015 | 2.2736 |
| 3   | 0.3688 | 0.658  | 1.0737 | 1.6207 | 2.289  |
| 4   | 0.3565 | 0.6194 | 1.0181 | 1.5665 | 2.2518 |

### Bates_Call_Price_Using_FFT ###
Souece Code: [Bates_Call_Price_Using_FFT](https://github.com/moyixinqing/Jump-Process/tree/master/Bates_Call_Price_Using_FFT)

|   | Strike  | Exact  | FFT Trapz | FFT Simp | #Error Tr | #Error Sim |
|---|---------|--------|-----------|----------|-----------|------------|
| 0 | 41.4102 | 9.5305 | 9.5306    | 9.5306   | 0.0007    | 0.0007     |
| 1 | 44.0956 | 7.3771 | 7.3769    | 7.3769   | -0.0024   | -0.0024    |
| 2 | 46.9551 | 5.3330 | 5.3331    | 5.3331   | 0.0019    | 0.0019     |
| 3 | 50.0000 | 3.5120 | 3.5119    | 3.5119   | -0.0015   | -0.0015    |
| 4 | 53.2424 | 2.0412 | 2.0411    | 2.0411   | -0.0015   | -0.0015    |
| 5 | 56.6950 | 1.0236 | 1.0237    | 1.0237   | 0.0099    | 0.0099     |
| 6 | 60.3716 | 0.4539 | 0.4538    | 0.4538   | -0.0183   | -0.0183    |

| Item                                | Value               |
|-------------------------------------|---------------------|
| Integration increment               | 0.02441400000000000 |
| Log strike increment                | 0.06283200000000000 |
| Trapezoidal FFT mean absolute error | 0.00519165869547100 |
| Simpsons FFT mean absolute error    | 0.00519165829025704 |

### Bates_model_DJIA_parameter_estimation ###
Souece Code: [Bates_model_DJIA_parameter_estimation](https://github.com/moyixinqing/Jump-Process/tree/master/Bates_model_DJIA_parameter_estimation)

![alt text](https://github.com/moyixinqing/Jump-Process/blob/master/Bates_model_DJIA_parameter_estimation/Market%20and%20Bates%20Implied%20Volatilities%20from%20Puts%20on%20DIA.jpg)

### Bates_Effect_of_jump_parameters_on_RND ###
Souece Code: [Bates_Effect_of_jump_parameters_on_RND](https://github.com/moyixinqing/Jump-Process/tree/master/Bates_Effect_of_jump_parameters_on_RND)

![alt text](https://github.com/moyixinqing/Jump-Process/blob/master/Bates_Effect_of_jump_parameters_on_RND/Implied%20and%20Local%20Volatility.jpg)

![alt text](https://github.com/moyixinqing/Jump-Process/blob/master/Bates_Effect_of_jump_parameters_on_RND/Effect%20of%20Jump%20Parameters%20on%20Risk%20Netural%20Densities.jpg)

### Bates_Call_Price_Using_Simulation ###
Souece Code: [Bates_Call_Price_Using_Simulation](https://github.com/moyixinqing/Jump-Process/tree/master/Bates_Call_Price_Using_Simulation)

[Euler and Milstein Discretization](https://github.com/moyixinqing/Jump-Process/blob/master/Bates_Call_Price_Using_Simulation/Euler%20and%20Milstein%20Discretization.pdf)

Bates price by Euler or Milstein Simulation simulation <br/>

| Method Price | Price  | DollarError |
|--------------|--------|-------------|
| Closed Form  | 4.7049 | n/a         |
| Simulation   | 4.6957 | 0.0091      |


