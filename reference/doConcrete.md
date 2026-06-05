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
#> 1:    A=1  2500     1 -0.007503576      0.001668183 0.007503576
#> 2:    A=0  2500     1  0.012855822      0.004800477 0.012855822
#> 3:    A=1  2500     2 -0.003757434      0.008589528 0.003757434
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001668183     4.4980545           Inf  4.50 relative
#> 2:                0  0.004800477     2.6780301           Inf  2.68 relative
#> 3:                0  0.008589528     0.4374436           Inf  0.44 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.01557269
#> 
#> Starting TMLE Update:
#> Problem dimension: k = 4 
#> Using standard universal LFM approach (iterative small steps)
#> Starting step 1 with update epsilon = 0.1 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.006298007      0.001666380 0.006298007
#> 2:    A=0  2500     2  0.008705936      0.009176792 0.008705936
#> 3:    A=0  2500     1  0.003574282      0.004787829 0.003574282
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001666380     3.7794531           Inf  3.78 relative
#> 2:                0  0.009176792     0.9486906           Inf  0.95 relative
#> 3:                0  0.004787829     0.7465350           Inf  0.75 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.01267847
#> Starting step 2 with update epsilon = 0.1 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 2 with update epsilon = 0.05 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.005369281      0.001665150 0.005369281
#> 2:    A=0  2500     2 -0.007152786      0.009186461 0.007152786
#> 3:    A=0  2500     1  0.003644437      0.004788008 0.003644437
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001665150     3.2245029           Inf  3.22 relative
#> 2:                0  0.009186461     0.7786226           Inf  0.78 relative
#> 3:                0  0.004788008     0.7611593           Inf  0.76 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.01043122
#> Starting step 3 with update epsilon = 0.05 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 3 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.005110785      0.001664829 0.005110785
#> 2:    A=0  2500     1  0.001794214      0.004785959 0.001794214
#> 3:    A=0  2500     2  0.001459226      0.009180357 0.001459226
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664829     3.0698555           Inf  3.07 relative
#> 2:                0  0.004785959     0.3748912           Inf  0.37 relative
#> 3:                0  0.009180357     0.1589509           Inf  0.16 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.005611189
#> Starting step 4 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.004516020      0.001664131 0.004516020
#> 2:    A=0  2500     1  0.001183582      0.004785364 0.001183582
#> 3:    A=0  2500     2 -0.001255921      0.009182059 0.001255921
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664131     2.7137400           Inf  2.71 relative
#> 2:                0  0.004785364     0.2473337           Inf  0.25 relative
#> 3:                0  0.009182059     0.1367799           Inf  0.14 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004835683
#> Starting step 5 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.003932219      0.001663499 0.003932219
#> 2:    A=0  2500     2  0.002154573      0.009179947 0.002154573
#> 3:    A=0  2500     1  0.000125501      0.004784312 0.000125501
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001663499    2.36382404           Inf  2.36 relative
#> 2:                0  0.009179947    0.23470435           Inf  0.23 relative
#> 3:                0  0.004784312    0.02623178           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004486477
#> Starting step 6 with update epsilon = 0.025 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 6 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0036671302      0.001663229 0.0036671302
#> 2:    A=0  2500     1  0.0003914089      0.004784588 0.0003914089
#> 3:    A=0  2500     2 -0.0007021947      0.009181692 0.0007021947
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001663229    2.20482604           Inf  2.20 relative
#> 2:                0  0.004784588    0.08180618           Inf  0.08 relative
#> 3:                0  0.009181692    0.07647771           Inf  0.08 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003755086
#> Starting step 7 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0033781711      0.001662947 0.0033781711
#> 2:    A=0  2500     2  0.0004825112      0.009180940 0.0004825112
#> 3:    A=0  2500     1  0.0001154230      0.004784317 0.0001154230
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662947    2.03143684           Inf  2.03 relative
#> 2:                0  0.009180940    0.05255575           Inf  0.05 relative
#> 3:                0  0.004784317    0.02412527           Inf  0.02 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003415184
#> Starting step 8 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0030922760      0.001662680 0.0030922760
#> 2:    A=0  2500     2 -0.0003410745      0.009181458 0.0003410745
#> 3:    A=0  2500     1  0.0001555877      0.004784362 0.0001555877
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662680    1.85981418           Inf  1.86 relative
#> 2:                0  0.009181458    0.03714819           Inf  0.04 relative
#> 3:                0  0.004784362    0.03252006           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003115591
#> Starting step 9 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.812110e-03      0.001662431 2.812110e-03
#> 2:    A=0  2500     2  3.455930e-04      0.009181024 3.455930e-04
#> 3:    A=1  2500     2 -5.741848e-05      0.008589298 5.741848e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662431   1.691565067           Inf  1.69 relative
#> 2:                0  0.009181024   0.037642100           Inf  0.04 relative
#> 3:                0  0.008589298   0.006684886           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002833866
#> Starting step 10 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.538560e-03      0.001662199 2.538560e-03
#> 2:    A=0  2500     2 -3.817379e-04      0.009181484 3.817379e-04
#> 3:    A=0  2500     1  8.334083e-05      0.004784295 8.334083e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662199    1.52722948           Inf  1.53 relative
#> 2:                0  0.009181484    0.04157693           Inf  0.04 relative
#> 3:                0  0.004784295    0.01741967           Inf  0.02 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002568949
#> Starting step 11 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -2.272520e-03      0.001661985 2.272520e-03
#> 2:    A=0  2500     2  5.274996e-04      0.009180912 5.274996e-04
#> 3:    A=1  2500     2 -4.399059e-05      0.008589296 4.399059e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661985   1.367352909           Inf  1.37 relative
#> 2:                0  0.009180912   0.057456127           Inf  0.06 relative
#> 3:                0  0.008589296   0.005121559           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002334154
#> Starting step 12 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0020164794      0.001661789 0.0020164794
#> 2:    A=0  2500     2 -0.0008420416      0.009181782 0.0008420416
#> 3:    A=0  2500     1  0.0001238755      0.004784337 0.0001238755
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661789    1.21343872           Inf  1.21 relative
#> 2:                0  0.009181782    0.09170786           Inf  0.09 relative
#> 3:                0  0.004784337    0.02589189           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002189066
#> Starting step 13 with update epsilon = 0.0125 
#> Update increased the active convergence objective, halving OneStepEps
#> Starting step 13 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.897407e-03      0.001661702 1.897407e-03
#> 2:    A=0  2500     2  3.260843e-04      0.009181036 3.260843e-04
#> 3:    A=0  2500     1 -4.252794e-05      0.004784171 4.252794e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661702   1.141845895           Inf  1.14 relative
#> 2:                0  0.009181036   0.035517158           Inf  0.04 relative
#> 3:                0  0.004784171   0.008889302           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001926018
#> Starting step 14 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.771510e-03      0.001661611 1.771510e-03
#> 2:    A=0  2500     2 -1.874402e-04      0.009181359 1.874402e-04
#> 3:    A=0  2500     1  2.862341e-05      0.004784242 2.862341e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661611   1.066139970           Inf  1.07 relative
#> 2:                0  0.009181359   0.020415300           Inf  0.02 relative
#> 3:                0  0.004784242   0.005982852           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001781927
#> Starting step 15 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -1.645956e-03      0.001661524 1.645956e-03
#> 2:    A=0  2500     2  1.322667e-04      0.009181157 1.322667e-04
#> 3:    A=1  2500     2 -2.989555e-05      0.008589294 2.989555e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661524    0.99063073           Inf  0.99 relative
#> 2:                0  0.009181157    0.01440632           Inf  0.01 relative
#> 3:                0  0.008589294    0.00348056           Inf  0.00 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001651624
```
