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

## Other changes

- Standardized plotting APIs: use `autoplot()` for ggplot output and
  `plotly_ld()` for interactive 3D landscapes. Existing `plot()` methods remain
  available with soft-deprecation warnings.
- Fixed a plotting bug in `plot.cv_vectorfield()`.
- Fixed `fit_2d_ld(..., vector_position = "middle")` so it works correctly.
- Improved consistency and clarity of user-facing messages and error reporting across core workflows.
- Updated examples to use the current recommended landscape workflow (`make_2d_ld()`).
- Fixed a runtime error in `summary.bootstrap_2d_ld()` when filtering minor minima.
- Fixed additional `summary.bootstrap_2d_ld()` runtime errors introduced by NSE handling changes.

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
