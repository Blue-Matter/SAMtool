# Methods for retro object

plot and summary functions for retro object.

## Usage

``` r
# S4 method for class 'retro,missing'
plot(x, color = NULL)

# S4 method for class 'retro'
summary(object)
```

## Arguments

- x:

  An object of class
  [retro](https://samtool.openmse.com/reference/retro-class.md).

- color:

  An optional character vector of colors for plotting.

- object:

  An object of class
  [retro](https://samtool.openmse.com/reference/retro-class.md).

## Value

A series of plots showing retrospective patterns in fishing mortality,
spawning biomass, recruitment, etc.

## Author

Q. Huynh

## Examples

``` r
# \donttest{
res <- SP(Data = swordfish)
ret <- retrospective(res, figure = FALSE)

summary(ret)
#>                   Mohn's rho
#> Fishing mortality      0.130
#> F/F[MSY]              -0.145
#> Biomass               -0.091
#> B/B[MSY]               0.152
#> B/B[0]                 0.152
plot(ret)





# }
```
