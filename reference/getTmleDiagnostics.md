# Extract TMLE convergence diagnostics

`getTmleDiagnostics()` returns the empirical efficient influence curve
(EIC) convergence diagnostics stored on a fitted `"ConcreteEst"` object.
Use `type = "components"` to inspect the final component-wise stopping
criteria, `type = "trace"` to inspect the update history, or
`type = "norm"` to inspect the norm of the empirical mean EIC across
update steps.

## Usage

``` r
getTmleDiagnostics(ConcreteEst, type = c("components", "trace", "norm"))
```

## Arguments

- ConcreteEst:

  `"ConcreteEst"` object returned by
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- type:

  character; one of `"components"`, `"trace"`, or `"norm"`.

## Value

A `data.table` containing the requested diagnostics.

## Examples

``` r
library(data.table)
library(concrete)

data <- as.data.table(survival::pbc)
data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
data[, trt := sample(0:1, nrow(data), TRUE)]
#>         id  time status   trt      age    sex
#>      <int> <int>  <int> <int>    <num> <fctr>
#>   1:     1   400      2     1 58.76523      f
#>   2:     2  4500      0     1 56.44627      f
#>   3:     3  1012      2     1 70.07255      m
#>   4:     4  1925      2     1 54.74059      f
#>   5:     5  1504      1     1 38.10541      f
#>  ---                                         
#> 196:   196  2363      0     1 57.04038      f
#> 197:   197  2365      0     1 44.62697      f
#> 198:   198  2357      0     1 35.79740      f
#> 199:   199  1592      0     0 40.71732      f
#> 200:   200  2318      0     0 32.23272      f

concrete.args <- formatArguments(DataTable = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 TargetTime = 2500,
                                 TargetEvent = c(1, 2),
                                 Intervention = makeITT(),
                                 CVArg = list(V = 2),
                                 MaxUpdateIter = 2,
                                 Verbose = FALSE)

# \donttest{
concrete.est <- doConcrete(concrete.args)
#> 
#> Estimating Treatment Propensity:
#> 
#> Estimating Hazards:
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1,3 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> 
#> Starting TMLE Update:
#> Warning: TMLE has not converged by step 2 - Estimates may not have the desired asymptotic properties
getTmleDiagnostics(concrete.est, type = "components")
#>    Intervention  Time Event         PnEIC RelativeCriteria     AbsPnEIC
#>          <char> <num> <num>         <num>            <num>        <num>
#> 1:          A=1  2500     1 -0.0009685427      0.003404555 0.0009685427
#> 2:          A=1  2500     2  0.0138214071      0.008307993 0.0138214071
#> 3:          A=0  2500     1  0.0042990963      0.003185392 0.0042990963
#> 4:          A=0  2500     2  0.0126429279      0.009243972 0.0126429279
#>    AbsoluteCriteria StopCriteria RelativeRatio AbsoluteRatio     ratio  check
#>               <num>        <num>         <num>         <num>     <num> <lgcl>
#> 1:                0  0.003404555     0.2844844           Inf 0.2844844   TRUE
#> 2:                0  0.008307993     1.6636277           Inf 1.6636277  FALSE
#> 3:                0  0.003185392     1.3496284           Inf 1.3496284  FALSE
#> 4:                0  0.009243972     1.3676942           Inf 1.3676942  FALSE
#>    StopRule StopAbsTol Converged ConvergenceStep
#>      <char>      <num>    <lgcl>           <num>
#> 1: relative          0     FALSE               3
#> 2: relative          0     FALSE               3
#> 3: relative          0     FALSE               3
#> 4: relative          0     FALSE               3
# }
```
