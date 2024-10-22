# Development version

- Added the `fit_2d_ld()` function to fit the landscape for 1D data using the MVKE algorithm and simple integration.
- Debug: 
	- Changed the function form of `find_eqs()` according to the new setting of the `MVKE()` function; added `linear_interp` to `sim_vf_options`.
	- Fixed a typo in the `fit_2d_vf()` function.
	- The parameter `na_action` in `fit_2d_vf()` was not effective for `method = "MVKE"` in the previous version. Now it is fixed. The `vector_position` parameter is now also effective for `method = "MVKE"`.

# fitlandr 0.1.0

- Initial release.
