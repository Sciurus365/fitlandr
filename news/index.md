# Changelog

## fitlandr 0.1.1

- Added the
  [`fit_2d_ld()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_ld.md)
  function to fit the landscape for 1D data using the MVKE algorithm and
  simple integration.
- Applied the log trick for the kernel calculation for better handling
  of small values.
- Debug:
  - Changed the function form of
    [`find_eqs()`](https://sciurus365.github.io/fitlandr/reference/find_eqs.md)
    according to the new setting of the
    [`MVKE()`](https://sciurus365.github.io/fitlandr/reference/mvke.md)
    function; added `linear_interp` to `sim_vf_options`.
  - Fixed a typo in the
    [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md)
    function.
  - The parameter `na_action` in
    [`fit_2d_vf()`](https://sciurus365.github.io/fitlandr/reference/fit_2d_vf.md)
    was not effective for `method = "MVKE"` in the previous version. Now
    it is fixed. The `vector_position` parameter is now also effective
    for `method = "MVKE"`.
  - For
    [`MVKE()`](https://sciurus365.github.io/fitlandr/reference/mvke.md),
    the Gaussian kernel was used despite the user specifying an
    exponential kernel. This has been fixed and the default was changed
    to Gaussian.

## fitlandr 0.1.0

CRAN release: 2023-02-10

- Initial release.
