# Find the production parameter based on depletion that produces MSY

For surplus production models, this function returns the production
exponent n corresponding to BMSY/K (Fletcher 1978).

## Usage

``` r
SP_production(depletion, figure = TRUE)
```

## Arguments

- depletion:

  The hypothesized depletion that produces MSY.

- figure:

  Local, plots figure of production function as a function of depletion
  (B/K)

## Value

The production function exponent n (numeric).

## Note

May be useful for parameterizing `n` in
[SP](https://samtool.openmse.com/reference/SP.md) and
[SP_SS](https://samtool.openmse.com/reference/SP.md).

## References

Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson
system. Fishery Bulletin 76:515:521.

## See also

[SP](https://samtool.openmse.com/reference/SP.md)
[SP_SS](https://samtool.openmse.com/reference/SP.md)

## Author

Q. Huynh

## Examples

``` r
SP_production(0.5)

#> [1] 2
SP_production(0.5)
#> [1] 2
```
