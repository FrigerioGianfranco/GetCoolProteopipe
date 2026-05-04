# Plotting an overall Q-Q plot

It crates a Q-Q plot considering all the intensities of protein from all
samples.

## Usage

``` r
GCP_QQplotIntensities(GCPlist, Title = "QQ plot intensities")
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- Title:

  NULL or character of length 1. The title you want to ad on the top of
  the plot.

## Value

a ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{

Fig02_QQ_plot_before_processing <- GCP_QQplotIntensities(GCPlist = GCPlist05,
                                                         Title = "QQ plot - before processing")
export_figures(Fig02_QQ_plot_before_processing)

} # }


```
