# doConcrete

doConcrete

## Usage

``` r
doConcrete(ConcreteArgs)

# S3 method for class 'ConcreteEst'
print(x, ...)

# S3 method for class 'ConcreteEst'
plot(x, convergence = FALSE, gweights = TRUE, ask = FALSE, ...)

# S3 method for class 'ConcreteOut'
print(x, ...)
```

## Arguments

- ConcreteArgs:

  "ConcreteArgs" object : output of formatArguments()

- x:

  a ConcreteOut object

- ...:

  additional arguments to be passed into print methods

- convergence:

  logical: plot the PnEIC norms for each TMLE small update step

- gweights:

  logical: plot the densities of the intervention-related nuisance
  weights for each intervention

- ask:

  logical: whether or not to prompt for user input before displaying
  plots

## Value

Object with S3 class `"ConcreteEst"`. The fitted object stores TMLE
convergence metadata, including the stopping rule, update trace, and
empirical EIC norm trajectory. Use
[`getTmleDiagnostics()`](https://blind-contours.github.io/concrete/reference/getTmleDiagnostics.md)
to inspect these diagnostics.

## Functions

- `print(ConcreteEst)`: print.ConcreteEst print method for "ConcreteEst"
  class

- `plot(ConcreteEst)`: plot.ConcreteEst plot method for "ConcreteEst"
  class

- `print(ConcreteOut)`: print.ConcreteOut print method for "ConcreteOut"
  class

## Examples

``` r
library(data.table)
library(concrete)

data <- as.data.table(survival::pbc)
data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
data[, trt := sample(0:1, nrow(data), TRUE)]
#>         id  time status   trt      age    sex
#>      <int> <int>  <int> <int>    <num> <fctr>
#>   1:     1   400      2     0 58.76523      f
#>   2:     2  4500      0     0 56.44627      f
#>   3:     3  1012      2     0 70.07255      m
#>   4:     4  1925      2     0 54.74059      f
#>   5:     5  1504      1     0 38.10541      f
#>  ---                                         
#> 196:   196  2363      0     0 57.04038      f
#> 197:   197  2365      0     1 44.62697      f
#> 198:   198  2357      0     1 35.79740      f
#> 199:   199  1592      0     1 40.71732      f
#> 200:   200  2318      0     1 32.23272      f

# formatArguments() returns correctly formatted arguments for doConcrete()

concrete.args <- formatArguments(DataTable = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 TargetTime = 2500,
                                 TargetEvent = c(1, 2),
                                 Intervention = makeITT(),
                                 CVArg = list(V = 2))

# doConcrete() returns tmle (and g-formula plug-in) estimates of targeted risks
concrete.est <- doConcrete(concrete.args)
#> 
#> Estimating Treatment Propensity:
#> 
#> Estimating Hazards:
#> Warning: Loglik converged before variable  1,3 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  3 ; coefficient may be infinite. 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=0  2500     2  0.065164882      0.009123345 0.065164882
#> 2:    A=1  2500     2  0.050489833      0.008559467 0.050489833
#> 3:    A=1  2500     1 -0.008810427      0.001669846 0.008810427
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.009123345      7.142652           Inf  7.14 relative
#> 2:                0  0.008559467      5.898712           Inf  5.90 relative
#> 3:                0  0.001669846      5.276192           Inf  5.28 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.08326815
#> 
#> Starting TMLE Update:
#> Problem dimension: k = 4 
#> Using standard universal LFM approach (iterative small steps)
#> Starting step 1 with update epsilon = 0.1 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.007857597      0.001668246 0.007857597
#> 2:    A=0  2500     2  0.031032672      0.009080694 0.031032672
#> 3:    A=1  2500     2  0.027514698      0.008536386 0.027514698
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001668246      4.710094           Inf  4.71 relative
#> 2:                0  0.009080694      3.417434           Inf  3.42 relative
#> 3:                0  0.008536386      3.223225           Inf  3.22 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.04331261
#> Starting step 2 with update epsilon = 0.1 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.006659947      0.001666426 0.006659947
#> 2:    A=0  2500     1  0.010278816      0.004770394 0.010278816
#> 3:    A=1  2500     2  0.002478026      0.008528344 0.002478026
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001666426     3.9965460           Inf  4.00 relative
#> 2:                0  0.004770394     2.1547101           Inf  2.15 relative
#> 3:                0  0.008528344     0.2905636           Inf  0.29 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.01251422
#> Starting step 3 with update epsilon = 0.1 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.004977405      0.001664270 0.004977405
#> 2:    A=1  2500     2 -0.005910836      0.008529625 0.005910836
#> 3:    A=0  2500     2  0.005110688      0.009070745 0.005110688
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664270     2.9907441           Inf  2.99 relative
#> 2:                0  0.008529625     0.6929772           Inf  0.69 relative
#> 3:                0  0.009070745     0.5634254           Inf  0.56 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.009330346
#> Starting step 4 with update epsilon = 0.1 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 4 with update epsilon = 0.05 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 4 with update epsilon = 0.025 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0047684797      0.001664046 0.0047684797
#> 2:    A=0  2500     1  0.0014107092      0.004755201 0.0014107092
#> 3:    A=0  2500     2 -0.0009671491      0.009071217 0.0009671491
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664046     2.8655944           Inf  2.87 relative
#> 2:                0  0.004755201     0.2966666           Inf  0.30 relative
#> 3:                0  0.009071217     0.1066173           Inf  0.11 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.005070701
#> Starting step 5 with update epsilon = 0.025 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0041599532      0.001663394 0.0041599532
#> 2:    A=0  2500     2  0.0015251116      0.009070901 0.0015251116
#> 3:    A=0  2500     1  0.0003560477      0.004753881 0.0003560477
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001663394    2.50088268           Inf  2.50 relative
#> 2:                0  0.009070901    0.16813231           Inf  0.17 relative
#> 3:                0  0.004753881    0.07489621           Inf  0.07 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004464216
#> Starting step 6 with update epsilon = 0.025 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0036111580      0.001662857 0.0036111580
#> 2:    A=0  2500     2 -0.0022592566      0.009071464 0.0022592566
#> 3:    A=0  2500     1  0.0005336252      0.004753824 0.0005336252
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662857     2.1716591           Inf  2.17 relative
#> 2:                0  0.009071464     0.2490509           Inf  0.25 relative
#> 3:                0  0.004753824     0.1122518           Inf  0.11 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004303761
#> Starting step 7 with update epsilon = 0.025 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 7 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0033565271      0.001662621 0.0033565271
#> 2:    A=0  2500     2  0.0007891780      0.009070979 0.0007891780
#> 3:    A=1  2500     2 -0.0001286102      0.008528439 0.0001286102
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662621    2.01881634           Inf  2.02 relative
#> 2:                0  0.009070979    0.08700032           Inf  0.09 relative
#> 3:                0  0.008528439    0.01508016           Inf  0.02 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003450563
#> Starting step 8 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0030793463      0.001662377 0.0030793463
#> 2:    A=0  2500     2 -0.0005031201      0.009071152 0.0005031201
#> 3:    A=0  2500     1  0.0001554357      0.004753432 0.0001554357
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662377    1.85237506           Inf  1.85 relative
#> 2:                0  0.009071152    0.05546375           Inf  0.06 relative
#> 3:                0  0.004753432    0.03269967           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003124308
#> Starting step 9 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.802624e-03      0.001662145 2.802624e-03
#> 2:    A=0  2500     2  4.397499e-04      0.009071021 4.397499e-04
#> 3:    A=1  2500     2 -7.122961e-05      0.008528416 7.122961e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662145   1.686149042           Inf  1.69 relative
#> 2:                0  0.009071021   0.048478549           Inf  0.05 relative
#> 3:                0  0.008528416   0.008352033           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002837864
#> Starting step 10 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.532937e-03      0.001661930 2.532937e-03
#> 2:    A=0  2500     2 -4.428983e-04      0.009071143 4.428983e-04
#> 3:    A=0  2500     1  8.646205e-05      0.004753340 8.646205e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661930    1.52409338           Inf  1.52 relative
#> 2:                0  0.009071143    0.04882497           Inf  0.05 relative
#> 3:                0  0.004753340    0.01818975           Inf  0.02 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002573131
#> Starting step 11 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.269526e-03      0.001661731 2.269526e-03
#> 2:    A=0  2500     2  5.528891e-04      0.009071007 5.528891e-04
#> 3:    A=1  2500     2 -5.339599e-05      0.008528397 5.339599e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661731   1.365760479           Inf  1.37 relative
#> 2:                0  0.009071007   0.060951238           Inf  0.06 relative
#> 3:                0  0.008528397   0.006260964           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002337542
#> Starting step 12 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0020163864      0.001661549 0.0020163864
#> 2:    A=0  2500     2 -0.0008071269      0.009071199 0.0008071269
#> 3:    A=0  2500     1  0.0001200647      0.004753358 0.0001200647
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661549    1.21355837           Inf  1.21 relative
#> 2:                0  0.009071199    0.08897687           Inf  0.09 relative
#> 3:                0  0.004753358    0.02525892           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002175417
#> Starting step 13 with update epsilon = 0.0125 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 13 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.897209e-03      0.001661467 1.897209e-03
#> 2:    A=0  2500     2  2.610421e-04      0.009071044 2.610421e-04
#> 3:    A=0  2500     1 -3.495704e-05      0.004753228 3.495704e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661467   1.141888446           Inf  1.14 relative
#> 2:                0  0.009071044   0.028777512           Inf  0.03 relative
#> 3:                0  0.004753228   0.007354378           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001915757
#> Starting step 14 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.771575e-03      0.001661382 1.771575e-03
#> 2:    A=0  2500     2 -1.308206e-04      0.009071097 1.308206e-04
#> 3:    A=1  2500     2 -3.219348e-05      0.008528379 3.219348e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661382   1.066326106           Inf  1.07 relative
#> 2:                0  0.009071097   0.014421700           Inf  0.01 relative
#> 3:                0  0.008528379   0.003774865           Inf  0.00 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001776809
#> Starting step 15 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.646529e-03      0.001661301 1.646529e-03
#> 2:    A=0  2500     2  8.137837e-05      0.009071068 8.137837e-05
#> 3:    A=1  2500     2 -3.005180e-05      0.008528375 3.005180e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661301   0.991108248           Inf  0.99 relative
#> 2:                0  0.009071068   0.008971200           Inf  0.01 relative
#> 3:                0  0.008528375   0.003523743           Inf  0.00 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001648847
```
