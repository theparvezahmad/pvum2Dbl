Development of a 2D cartesian code based on PVU-M scheme for compressible boundary layer.
Conservative form of the NS equations used.
Implemented existing PVUM+ scheme for the case of compressible boundary layer.
A flux based, two-step predictor-corrector scheme.
The total flux is separated into the ‘convective’ and ‘non-convective’ parts.
The gradient of the convective fluxes is constructed from the cell face values.
The non-convective fluxes are computed on the grid nodes.
