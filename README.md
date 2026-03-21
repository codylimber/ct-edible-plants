# Connecticut Wild Edibles Phenology

Best times to forage wild edibles in Connecticut

## About

Interactive phenology calendar showing the best times to observe wild edible
species in Connecticut, based on iNaturalist research-grade observations
(2015--2025).

**Live site**: https://codylimber.github.io/ct-edible-plants/

## Setup

1. Install R (>= 4.4) and [Quarto](https://quarto.org/docs/get-started/)
2. Restore R packages:
   ```r
   install.packages("renv")
   renv::restore()
   ```
3. Render the site:
   ```bash
   quarto render
   ```

The rendered site will be in `docs/`.

## Configuration

Edit `params.yml` to customize:
- Taxon and place (iNaturalist IDs)
- Date range, minimum observation threshold
- Rarity thresholds
- Terminology (flights, broods, blooms, etc.)

## Updating data

Re-render to fetch the latest iNaturalist observations. Cached data is stored
in `outs/cache/` — delete this directory to force a full re-fetch.

## Built with

- [phenology-site](https://github.com/codylimber/iNat) template engine
- [iNaturalist API](https://api.inaturalist.org/)
- [Quarto](https://quarto.org/)
- R packages: rinat, tidyverse, plotly, sf, tigris
