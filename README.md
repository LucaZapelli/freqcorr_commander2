# commander2.5
Generelization of Commander2 code for frequency-correlated noise:

 1. Each processor handles a single main band, although importing every noise covariance and computing residual/signal maps for all bands when a N^-1 d operation has to be performed. Other quantities are expanded to every band, such as beams, foreground spectral responses and templates.
2. Ft N^-1 F computation is simplified, limiting the Ft N^-1 F solving to a vector-scalar product rather then a matrix-vector one.
3. Chi^2 computation formula has been changed, replacing the sqrt{N^-1} multiplication with the standard N^-1 one. The reason is that the new routine treats differently the rhs e and the lhs map terms, so that the multiplication of the two no longer equals lhs^2.
4. The code is generalized so that each frequency channel can be associated to a specific mask. However, the format of the imported covariance matrices is left untouched, having their size related to the number of observed pixels. Therefore, each frequency mask is currently forced to be the same mask in order to avoid rectangular off-diagonal block matrices. 

For more information check out the QUALITATIVE_LOG.txt and QUANTITATIVE_LOG.txt files (given in order of how detailed the list of modifications is).


##########################################################################################################################################

# List of new global variables:

- integer(i4b) :: band_iter, map_id_fcn
- logical(lgt) :: freq_corr_noise
- real(dp),     allocatable, dimension(:,:,:)           :: beams_fcn
- real(dp),     allocatable, dimension(:),      private :: reg_noises_fcn
- real(dp),     allocatable, dimension(:,:,:)           :: cmbmaps_fcn, residuals_fcn
- real(dp),     allocatable, dimension(:,:,:)           :: mask_calibs_fcn
- real(dp),     allocatable, dimension(:,:,:,:)         :: fg_pix_spec_responses_fcn
- real(dp),     allocatable, dimension(:,:,:,:)         :: fg_temps_fcn
- real(dp),     allocatable, dimension(:,:,:,:)         :: inv_Ns_lowres_fcn, inv_Ns_scaled_fcn
- real(dp),     allocatable, dimension(:,:),    private :: my_inv_Ns_fcn
- real(dp),     allocatable, dimension(:,:,:),  private :: all_inv_Ns_fcn
- real(dp),     allocatable, dimension(:,:,:),  private :: invN_regs_fcn


##########################################################################################################################################

# New parameters/files needed in the 'parameter_file.txt'

- FREQ_CORR_NOISE          = .true. # or .false. (needed in any case)
- NOISE_RMS0101            = '../..' # instead of standard NOISE_RMS01 (needed only if 'freq_corr_noise' is true)      > it's a .fits map
- INV_N_MAT0101            = '../..' # instead of standard INV_N_MAT01 (needed only if 'freq_corr_noise' is true)      > it's a .unf file
- SQRT_INV_N_MAT0101       = '../..' # instead of standard SQRT_INV_N_MAT01 (needed only if 'freq_corr_noise' is true) > it's a .unf file
