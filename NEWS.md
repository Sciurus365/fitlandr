# fitlandr 0.1.1.9000

## Important changes

- The core algorithm has been substantially redesigned:
  - A cross-validation procedure is now used to select the kernel width in the drift–diffusion (MVKE) estimation.
  - A finite volume method is introduced to compute the steady-state distribution from the drift–diffusion function.
  - The previous simulation-based approach is now discouraged in favor of this new method.
  - A bootstrapping pipeline has been added to support inference based on the potential landscape.
- As a result of these changes, several previously central functions are now deprecated.

- Dimension notation has been standardized:
  - The dimension of the system is now defined solely by the input data (state variables).
  - This dimension is used consistently for both the vector field and the landscape.
  - The potential function \(U\) is no longer treated as an additional dimension.

## New functionality

- Added one-dimensional vector-field, cross-validation, landscape, and
  bootstrap-inference workflows corresponding to the new two-dimensional
  APIs.
- Added `make_2d_stream()` to estimate a stream function from a
  two-dimensional probability flow using sparse finite-difference least
  squares, together with an `autoplot()` method and residual diagnostics.
- Added `fit_individual_dynamics()` for the complete single-dataset workflow
  from a cross-validated vector field through its landscape, probability flow,
  and stream function. `fit_group_dynamics()` reuses this workflow over
  multiple datasets and provides faceted plots of individual landscapes and
  stream functions. Group results can be passed separately to the landscape-
  or stream-clustering APIs.
- Added `evaluate_stream_clusters()` and `cluster_streams()` for K-means
  clustering of mean-centered stream functions, with elbow and cluster-center
  `autoplot()` methods. `evaluate_landscape_stream_clusters()` and
  `cluster_landscape_streams()` jointly cluster density-space landscapes and
  mean-centered streams after scaling each block to unit mean pairwise squared
  distance, with optional modality weights.
- Added a prototype landscape-clustering API. Use
  `evaluate_landscape_clusters()` to compare candidate K-means solutions with
  an elbow plot and `cluster_landscapes()` to obtain assignments and
  density-space cluster centers transformed back to potential landscapes.
- Standardized plotting APIs: use `autoplot()` for ggplot output and
  `plotly_ld()` for interactive 3D landscapes. Existing `plot()` methods remain
  available with soft-deprecation warnings.

## Other changes

- Fixed `fit_2d_ld(..., vector_position = "middle")` so it works correctly.
- Improved consistency and clarity of user-facing messages and error reporting across core workflows.
- Updated examples to use the current recommended landscape workflow (`make_2d_ld()`).

# fitlandr 0.1.1

- Added the `fit_2d_ld()` function to fit the landscape for 1D data using the MVKE algorithm and simple integration.
- Applied the log trick for the kernel calculation for better handling of small values.
- Debug: 
	- Changed the function form of `find_eqs()` according to the new setting of the `MVKE()` function; added `linear_interp` to `sim_vf_options`.
	- Fixed a typo in the `fit_2d_vf()` function.
	- The parameter `na_action` in `fit_2d_vf()` was not effective for `method = "MVKE"` in the previous version. Now it is fixed. The `vector_position` parameter is now also effective for `method = "MVKE"`.
	- For `MVKE()`, the Gaussian kernel was used despite the user specifying an exponential kernel. This has been fixed and the default was changed to Gaussian.

# fitlandr 0.1.0

- Initial release.
