class07
================
Tianhao Qiu
4/23/2019

Function Revisited
==================

1. Source file online with function last time
---------------------------------------------

``` r
source("http://tinyurl.com/rescale-R")
```

### Try out the rescale() function last time

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

### Try rescale2() that catches invalid form

``` r
#rescale2(c(1:10, "string"))
```

2. Find missing NA values in two vectors
========================================

Start with simple example
-------------------------

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

Combine them
------------

``` r
naVal <- is.na(x) & is.na(y)
sum(naVal)
```

    ## [1] 1

Now write the function for the task
-----------------------------------

``` r
both_na <- function(x, y) {
  if(length(x) != length(y)) {
    stop("Input x and y should be the same length")
  }
  naVal <- is.na(x) & is.na(y)
  naNum <- (sum(naVal))
  naWhich <- which(naVal)
  message("Found ", naNum, " NA's at position(s):",
          paste(naWhich, collapse = ", "))
  return( list(number = naNum, which = naWhich) )
}
```

``` r
both_na(x, c(NA, NA, NA, NA, NA))
```

    ## Found 2 NA's at position(s):3, 5

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 3 5

!When lengths are inconsistent, will do cycling by default
----------------------------------------------------------

Intersect function
==================

``` r
x <- df1$IDs
y <- df2$IDs
intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
which(x %in% y)
```

    ## [1] 2 3

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
which(y %in% x)
```

    ## [1] 1 3

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],
        y[y %in% x])
}
```

Update intersect function
-------------------------

``` r
gene_intersect2 <- function(df1, df2) { 
  cbind( df1[ df1$IDs %in% df2$IDs, ],
         df2[ df2$IDs %in% df1$IDs, "exp"] )
}
gene_intersect2 (df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

``` r
gene_intersect3 <- function(df1, df2, gene.colname="IDs") {
  cbind( df1[ df1[,gene.colname] %in% df2[,gene.colname], ],
        exp2 = df2[ df2[,gene.colname] %in% df1[,gene.colname], "exp"] )
}
```

``` r
gene_intersect3(df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

``` r
merge(df1, df2, by = "IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

``` r
x <- c(100, 100, 100, 100, 100, 100, 100, 90)
y <- c(100, 90, 90, 90, 90 ,90, 97,80)

averageGrade <- function(x, na.rm = TRUE) {
  if (!is.numeric(x)) {
    stop("Invalid Input")
  }
  initialSum = sum(x)
  minGrade = min(x)
  sumAftDrop = sum(x) - min(x)
  average = sumAftDrop / (length(x) - 1)
  message("Studednt's average score is ", average)
}
```

``` r
averageGrade(x)
```

    ## Studednt's average score is 100

``` r
averageGrade(y)
```

    ## Studednt's average score is 92.4285714285714
