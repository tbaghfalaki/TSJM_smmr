### TSJM package
On this page, you can access the R code accompanying our paper titled "A Two-stage Joint Modeling Approach for Multiple Longitudinal Markers and Time-to-event Data" authored by Baghfalaki, Hashemi, Helmer and Jacqmin-Gadda (2025).

### TSJM packege 
To install the latest development version of the **TSJM** package and explore usage examples, visit:  
👉 [https://github.com/tbaghfalaki/TSJM](https://github.com/tbaghfalaki/TSJM)

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/TSJM")
```
This will seamlessly fetch and install the most up-to-date version of TSJM for your use.


Note: Before installing this package, please install the *parallelsugar* package using the following command:

```
  # install.packages("devtools")
  devtools::install_github('nathanvan/parallelsugar')
```

### Example Usage
We consider two methods here. The first method is the full Bayesian approach, where both stages of the analysis are performed using Bayesian techniques. As a result, the second stage requires more computational time. Both stages can be easily executed using the TSJM package. Examples from the full Bayesian approach in TSJM can be found in the following:
- Parameter estimation: This analysis is presented [here](/Exam1.md)
- Dynamic prediction: This analysis is presented [here](/Exam2.md)

The second method uses the Rubin formula to estimate the standard deviation. We apply this approach to solve two examples: 
- Example 1: Without considring explanatroty variables  [here](/Exam3.md)
- Example 2: By considring explanatroty variables [here](/Exam4.md)


### Reference 
Baghfalaki, T., Hashemi, R., Helmer, C. & Jacqmin-Gadda, H. (2024). A Two-stage Joint Modeling Approach for Multiple Longitudinal Markers and Time-to-event Data. Submitted.
