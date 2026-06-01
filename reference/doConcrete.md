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
#>       Trt  Time Event       PnEIC RelativeCriteria   AbsPnEIC AbsoluteCriteria
#>    <char> <num> <num>       <num>            <num>      <num>            <num>
#> 1:    A=0  2500    -1 -0.07292871      0.009466211 0.07292871                0
#> 2:    A=0  2500     2  0.06516488      0.009123345 0.06516488                0
#> 3:    A=1  2500     2  0.05048983      0.008559467 0.05048983                0
#>    StopCriteria RelativeRatio AbsoluteRatio ratio StopRule StopAbsTol
#>           <num>         <num>         <num> <num>   <char>      <num>
#> 1:  0.009466211      7.704108           Inf  7.70 relative          0
#> 2:  0.009123345      7.142652           Inf  7.14 relative          0
#> 3:  0.008559467      5.898712           Inf  5.90 relative          0
#> Norm PnEIC = 0.08326815
#> 
#> Starting TMLE Update:
#> Problem dimension: k = 4 
#> Using standard universal LFM approach (iterative small steps)
#> Starting step 1 with update epsilon = 0.1 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.007857597      0.001668246 0.007857597
#> 2:    A=0  2500    -1 -0.040736023      0.009413281 0.040736023
#> 3:    A=0  2500     2  0.031032672      0.009080694 0.031032672
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001668246      4.710094           Inf  4.71 relative
#> 2:                0  0.009413281      4.327505           Inf  4.33 relative
#> 3:                0  0.009080694      3.417434           Inf  3.42 relative
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
#> 3:    A=0  2500    -1 -0.009603478      0.009387906 0.009603478
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001666426      3.996546           Inf  4.00 relative
#> 2:                0  0.004770394      2.154710           Inf  2.15 relative
#> 3:                0  0.009387906      1.022963           Inf  1.02 relative
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
#> 2:    A=1  2500    -1  0.010888241      0.008591127 0.010888241
#> 3:    A=1  2500     2 -0.005910836      0.008529625 0.005910836
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664270     2.9907441           Inf  2.99 relative
#> 2:                0  0.008591127     1.2673822           Inf  1.27 relative
#> 3:                0  0.008529625     0.6929772           Inf  0.69 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.009330346
#> Starting step 4 with update epsilon = 0.1 
#> Update increased ||PnEIC||, halving OneStepEps
#> Starting step 4 with update epsilon = 0.05 
#> Update increased ||PnEIC||, halving OneStepEps
#> Starting step 4 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.004768480      0.001664046 0.004768480
#> 2:    A=1  2500    -1  0.004549082      0.008588378 0.004549082
#> 3:    A=0  2500     1  0.001410709      0.004755201 0.001410709
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001664046     2.8655944           Inf  2.87 relative
#> 2:                0  0.008588378     0.5296789           Inf  0.53 relative
#> 3:                0  0.004755201     0.2966666           Inf  0.30 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.005070701
#> Starting step 5 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.004159953      0.001663394 0.004159953
#> 2:    A=1  2500    -1  0.004573814      0.008588377 0.004573814
#> 3:    A=0  2500    -1 -0.001881159      0.009387628 0.001881159
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001663394     2.5008827           Inf  2.50 relative
#> 2:                0  0.008588377     0.5325586           Inf  0.53 relative
#> 3:                0  0.009387628     0.2003871           Inf  0.20 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004464216
#> Starting step 6 with update epsilon = 0.025 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.003611158      0.001662857 0.003611158
#> 2:    A=1  2500    -1  0.003306361      0.008587958 0.003306361
#> 3:    A=0  2500     2 -0.002259257      0.009071464 0.002259257
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662857     2.1716591           Inf  2.17 relative
#> 2:                0  0.008587958     0.3849997           Inf  0.38 relative
#> 3:                0  0.009071464     0.2490509           Inf  0.25 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.004303761
#> Starting step 7 with update epsilon = 0.025 
#> Update increased ||PnEIC||, halving OneStepEps
#> Starting step 7 with update epsilon = 0.0125 
#>       Trt  Time Event        PnEIC RelativeCriteria    AbsPnEIC
#>    <char> <num> <num>        <num>            <num>       <num>
#> 1:    A=1  2500     1 -0.003356527      0.001662621 0.003356527
#> 2:    A=1  2500    -1  0.003485137      0.008588013 0.003485137
#> 3:    A=0  2500     2  0.000789178      0.009070979 0.000789178
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662621    2.01881634           Inf  2.02 relative
#> 2:                0  0.008588013    0.40581415           Inf  0.41 relative
#> 3:                0  0.009070979    0.08700032           Inf  0.09 relative
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
#> 2:    A=1  2500    -1  0.0031197857      0.008587900 0.0031197857
#> 3:    A=0  2500     2 -0.0005031201      0.009071152 0.0005031201
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662377    1.85237506           Inf  1.85 relative
#> 2:                0  0.008587900    0.36327689           Inf  0.36 relative
#> 3:                0  0.009071152    0.05546375           Inf  0.06 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.003124308
#> Starting step 9 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0028026244      0.001662145 0.0028026244
#> 2:    A=1  2500    -1  0.0028738540      0.008587827 0.0028738540
#> 3:    A=0  2500     2  0.0004397499      0.009071021 0.0004397499
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001662145    1.68614904           Inf  1.69 relative
#> 2:                0  0.008587827    0.33464273           Inf  0.33 relative
#> 3:                0  0.009071021    0.04847855           Inf  0.05 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002837864
#> Starting step 10 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0025329366      0.001661930 0.0025329366
#> 2:    A=1  2500    -1  0.0025729624      0.008587741 0.0025729624
#> 3:    A=0  2500     2 -0.0004428983      0.009071143 0.0004428983
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661930    1.52409338           Inf  1.52 relative
#> 2:                0  0.008587741    0.29960876           Inf  0.30 relative
#> 3:                0  0.009071143    0.04882497           Inf  0.05 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002573131
#> Starting step 11 with update epsilon = 0.0125 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0022695259      0.001661731 0.0022695259
#> 2:    A=1  2500    -1  0.0023229218      0.008587672 0.0023229218
#> 3:    A=0  2500     2  0.0005528891      0.009071007 0.0005528891
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661731    1.36576048           Inf  1.37 relative
#> 2:                0  0.008587672    0.27049494           Inf  0.27 relative
#> 3:                0  0.009071007    0.06095124           Inf  0.06 relative
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
#> 2:    A=1  2500    -1  0.0020438893      0.008587597 0.0020438893
#> 3:    A=0  2500     2 -0.0008071269      0.009071199 0.0008071269
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661549    1.21355837           Inf  1.21 relative
#> 2:                0  0.008587597    0.23800479           Inf  0.24 relative
#> 3:                0  0.009071199    0.08897687           Inf  0.09 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.002175417
#> Starting step 13 with update epsilon = 0.0125 
#> Update increased ||PnEIC||, halving OneStepEps
#> Starting step 13 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0018972094      0.001661467 0.0018972094
#> 2:    A=1  2500    -1  0.0019340127      0.008587569 0.0019340127
#> 3:    A=0  2500     2  0.0002610421      0.009071044 0.0002610421
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661467    1.14188845           Inf  1.14 relative
#> 2:                0  0.008587569    0.22521073           Inf  0.23 relative
#> 3:                0  0.009071044    0.02877751           Inf  0.03 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001915757
#> Starting step 14 with update epsilon = 0.00625 
#>       Trt  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>    <char> <num> <num>         <num>            <num>        <num>
#> 1:    A=1  2500     1 -0.0017715752      0.001661382 0.0017715752
#> 2:    A=1  2500    -1  0.0018037687      0.008587536 0.0018037687
#> 3:    A=0  2500     2 -0.0001308206      0.009071097 0.0001308206
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661382     1.0663261           Inf  1.07 relative
#> 2:                0  0.008587536     0.2100450           Inf  0.21 relative
#> 3:                0  0.009071097     0.0144217           Inf  0.01 relative
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
#> 2:    A=1  2500    -1  1.676581e-03      0.008587505 1.676581e-03
#> 3:    A=0  2500     2  8.137837e-05      0.009071068 8.137837e-05
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio ratio StopRule
#>               <num>        <num>         <num>         <num> <num>   <char>
#> 1:                0  0.001661301     0.9911082           Inf  0.99 relative
#> 2:                0  0.008587505     0.1952349           Inf  0.20 relative
#> 3:                0  0.009071068     0.0089712           Inf  0.01 relative
#>    StopAbsTol
#>         <num>
#> 1:          0
#> 2:          0
#> 3:          0
#> Norm PnEIC = 0.001648847
```
