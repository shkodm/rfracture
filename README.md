# Seam geometry generator

Set of R scripts for generating rock/coal seam geometries.

## Install
```bash
wget https://github.com/llaniewski/rfracture/archive/master.tar.gz -O rfracture.tar.gz
R CMD INSTALL rfracture.tar.gz
rm rfracture.tar.gz
```

## Use

```R
library(rfracture)

alpha = 4.5
power.iso = function(f) 0.000001*ifelse(f<5,0,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)
refine = 20
closed = 0.1
frac = fracture_geom(refine=refine, closed = closed, corr.profile = corr.profile, power.iso = power.iso, seed=125)
frac = slice(frac,  value="above")
write.stl(frac, "fracture.stl")
```

## Visualize
```R
library(rgl)
fracture3d(frac)
```
