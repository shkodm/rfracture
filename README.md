# Seam geometry generator

Set of R scripts for generating rock/coal seam geometries.

## Install
You can install the package in R, with:
```R
remotes::install_github("llaniewski/rfracture")
```

Or through bash:
```bash
wget https://github.com/llaniewski/rfracture/archive/master.tar.gz -O rfracture.tar.gz
R CMD INSTALL rfracture.tar.gz
rm rfracture.tar.gz
```

## Use

```R
library(rfracture)

# Construct the power spectrum and correlation between upper and lower surface:
alpha = 4.5
power.iso = function(f) 0.000001*ifelse(f<5,0,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)

# Set the refinement level and aperture (gap)
refine = 20
gap = 0.015

# Construct the fracture geometry
frac = fracture_geom(refine=refine, gap=gap, corr.profile = corr.profile, power.iso = power.iso, seed=125)

# Cut it to a square
frac = cut(frac)

# Cut the parts where upper and lower surface overlap
frac = slice(frac, by = "h", level = 0, type = "above")

# Plot fracture in 3D
library(rgl)
fracture3d(frac)

# Save to STL
write.stl(frac, "fracture.stl")
```

You can also generate a fracture as matrices of upper and lower heigth fields.
```R
frac_m = fracture_matrix(dims = c(128,128), gap=gap, corr.profile = corr.profile, power.iso = power.iso, seed=125)

image(frac_m$f1 - frac_m$f2, asp=1)
```
