# Plotting an overall density plot

It crates an histogram and density curve considering all the intensities
of protein from all samples.

## Usage

``` r
GCP_DensityplotIntensities(GCPlist, Title = "Distribution of intensities")
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

Fig01_IntensityDistribution_before_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist05,
                                                                            Title = "Distribution of intensities, before processing")
export_figures(Fig01_IntensityDistribution_before_processing)

} # }


```
