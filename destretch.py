# ---------------------------- DESCRIPTION --------------------------------------- #
"""

Module with destretching related functions for the alignment of the observation modes. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

import numpy as np
import time

import torchmfbd
import torch
import gc

# Own functions
from demodulation import demodulate
from demodulation import mod_matrices_david_ct as mod_matrices
from process_data_utils import balance
from alignment import compute_pseudo_images_leastsq
import logging

def destretch(data, ngrid=8, ngrid_mod = 2, lr=0.50, reference_frame=0, border=6,
            n_iterations=200, lambda_tt=0.01,aling_cam='partial', filter = filter, align_modulations = False):
    """
    Aligns modulations and camera data using the `torchmfbd` package.

    This function utilizes the `torchmfbd` distribution from:
    https://github.com/aasensio/torchmfbd

    It inherits parameters from the original `destretch.py` program and performs
    preprocessing such as Fourier filtering and tip-tilt alignment.

    Parameters
    ----------
    data : np.ndarray
        Input data array with shape (Ncams, Nlambda, Nmods, Nx, Ny). If there is one dimension less, it is assumed that the input is two cameras and single wave
    ngrid : int, optional
        Grid size for tip-tilt estimation. Default is 2.
    lr : float, optional
        Learning rate for the optimizer. Default is 0.50.
    reference_frame : int, optional
        Index of the reference frame to align all other frames to. Default is 0.
    border : int, optional
        Border size excluded from loss computation. Default is 6. IMPORTANT. If we add the field stop, border > 300
    n_iterations : int, optional
        Number of optimization iterations. Default is 200.
    lambda_tt : float, optional
        Regularization weight for tip-tilt smoothness. Default is 0.01.
    aling_cam : string, optional
        Align the two cameras using the first modulation 'partial' or the four modulations 'full' or do nothing 'none'. Default is 'partial'

    Returns
    -------
    np.ndarray
        Filtered and aligned data array with the same shape as the input.

    Notes
    -----
    - Currently, the function does not save the alignment matrices. This should be implemented.
    - For installation and detailed setup of the `destretching_update` branch, refer to:
      https://github.com/PabloSGN/TuMags_Reduction_Pipeline/blob/destretching_update/Documents/Installation.md
    """

    tic = time.time() # Get the time to measure execution time.
    destretch_pars = None

    if len(data.shape) == 4:
        data = data[:, np.newaxis] # To allow for only one lamdba.
    elif len(data.shape) == 5:
        logging.info('data shape is correct')
    elif len(data.shape) == 3:
        logging.info('data shape is for serie option')
    else:
        raise ValueError("Data must be of shape (Ncams, Nlambda, Nmods, Nx, Ny) or (Ncams, Nmods, Nx, Ny)")
    
    if aling_cam == 'full':
        for lambd in range(data.shape[1]):
            logging.info('aligning lambda', lambd)

            #align modulations
            # first thing to do is move data to gpu memory

            cam0_frames =  torch.tensor(data[0,lambd,:,:,:].astype('float32')) 
            cam1_frames =  torch.tensor(data[1,lambd,:,:,:].astype('float32')) 
            # add tensor dimensions
            cam0_frames = cam0_frames.unsqueeze(0).unsqueeze(0)  # (1, 1, nmod, x, y)
            cam1_frames = cam1_frames.unsqueeze(0).unsqueeze(0)  # (1, 1, nmod, x, y)

            # run the destretching
            warped_cam0_frames, shift_0 = torchmfbd.destretch(
                cam0_frames,
                ngrid=ngrid_mod,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )
            warped_cam1_frames, shift_1 = torchmfbd.destretch(
                cam1_frames,
                ngrid=ngrid_mod,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )

            for i in range(4):
                for j in range(2):
                    print('shifts pol= ',i,'(x,y) = (0,1) ',j,' cam0 ',shift_0[i,j,data.shape[-1]//2,data.shape[-1]//2],' cam1 ',shift_1[i,j,data.shape[-1]//2,data.shape[-1]//2])
                    #detach and assoc

            data[0,lambd] = warped_cam0_frames[0, 0].detach().cpu().numpy()
            data[1,lambd] = warped_cam1_frames[0, 0].detach().cpu().numpy()


            for m in range(4):
                dm0 = torch.tensor(data[0, lambd, m].astype('float32'))
                dm1 = torch.tensor(data[1, lambd, m].astype('float32'))
                dm_stack = torch.stack([dm0, dm1], dim=0).unsqueeze(0).unsqueeze(0)

                warped_mod, _ = torchmfbd.destretch(
                    dm_stack, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
                    border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
                )

                data[:, lambd, m] = warped_mod[0, 0].detach().cpu().numpy()

    elif aling_cam=='same':
        for lambd in range(data.shape[1]):
            logging.info(f'aligning lambda {lambd}')

            # =============================
            # (1) DEMODULACIÓN DE LAS CAMARAS
            # =============================
            _, demod = demodulate(data[:, lambd], 
                                filt=filter, 
                                # dmod_matrices='demod_matrices_david',
                                onelambda=True, 
                                BothCams=True)

            # demod tiene shape (2, 4, x, y)
            I_cam1 = demod[0, 0]
            I_cam2 = demod[1, 0]

            scale, gamma = balance(I_cam1, I_cam2)
            I_cam2 = I_cam2 * scale
            logging.info(f"Balance (2D): scale={scale:.6f}, gamma={gamma:.6f}")

            # =============================
            # (2) ALINEAR CAM1 <--> CAM2 CON STOKES I
            # =============================
            logging.info("Global alignment between cameras using Stokes I")

            dm0 = torch.tensor(I_cam1.astype('float32'))
            dm1 = torch.tensor(I_cam2.astype('float32'))
            dm_stack = torch.stack([dm0, dm1], dim=0).unsqueeze(0).unsqueeze(0)

            _, distortion_map = torchmfbd.destretch(
                dm_stack, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
                border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
            )

            for i in range(2):
                for j in range(2):
                    print('shifts cam=', i, '(x,y)=', j,
                            distortion_map[i,j,data.shape[-1]//2,data.shape[-1]//2])

            for mod in range(4):
                mod_data = torch.tensor(data[:, lambd, mod].astype('float32')).unsqueeze(0).unsqueeze(0)
                warped = torchmfbd.apply_destretch(mod_data, distortion_map, mode='bilinear')
                data[:, lambd, mod] = warped[0, 0].detach().cpu().numpy()

            if align_modulations:
                # =============================
                # (3) PSEUDO-IMÁGENES POR MODULACIÓN j
                # =============================
                imgs_cam1 = data[0, lambd]    # (4,x,y)
                imgs_cam2 = data[1, lambd]    # (4,x,y)

                M1 = mod_matrices[filter][0]
                M2 = mod_matrices[filter][1]

                pseudo_imgs, alphas, betas, gammas, residuals = compute_pseudo_images_leastsq(
                    imgs_cam1, imgs_cam2, M1, M2
                )

                for j in range(4):
                    print(f"[λ={lambd}] mod {j}: alpha={alphas[j]:.4f}, beta={betas[j]:.4f}, "
                        f"gamma={gammas[j]:.4f}, residual_norm={residuals[j]:.4e}")

                # =============================
                # (4) ALINEACIÓN ENTRE MODULACIONES
                # =============================
                logging.info("Aligning modulations using pseudo-images...")

                pseudo_imgs_tensor = torch.tensor(pseudo_imgs.astype('float32')).unsqueeze(0).unsqueeze(0)

                _, distortion_map_pseudo_imgs = torchmfbd.destretch(
                    pseudo_imgs_tensor, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
                    border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
                )

                for i in range(4):
                    for j in range(2):
                        print('shifts mod=', i, '(x,y)=', j,
                                distortion_map_pseudo_imgs[i,j,data.shape[-1]//2,data.shape[-1]//2])


                # =============================
                # (5) APLICAR ESTE JITTER A AMBAS CAMARAS
                # =============================
                for cam in range(2):
                    mod_data = torch.tensor(data[cam, lambd, :].astype('float32')).unsqueeze(0).unsqueeze(0)
                    warped = torchmfbd.apply_destretch(mod_data, distortion_map_pseudo_imgs, mode='bilinear')
                    data[cam, lambd, :] = warped[0, 0].detach().cpu().numpy()

    elif aling_cam == 'all':
        for lambd in range(data.shape[1]):
            logging.info('aligning lambda', lambd)

            frames =  torch.tensor(np.reshape(data[:,lambd,:,:,:],(2*data.shape[-3],data.shape[-2],data.shape[-1])).astype('float32')) 
            # add tensor dimensions
            frames = frames.unsqueeze(0).unsqueeze(0)  # (1, 1, nmod, x, y)

            warped_frames, destretch_pars = torchmfbd.destretch(
                frames, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
                border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
            )
            data[:,lambd,:,:,:] = np.reshape(
                warped_frames[0, 0].detach().cpu().numpy(),(2,data.shape[-3],data.shape[-2],data.shape[-1])
                )
            
    elif aling_cam == '0s':
        # # first thing to do is move data to gpu memory
        # frames =  torch.tensor(np.reshape(data,(2*data.shape[-3]*data.shape[-4],data.shape[-2],data.shape[-1])).astype('float32')) 
        # # add tensor dimensions
        # frames = frames.unsqueeze(0).unsqueeze(0)  # (1, 1, nmod, x, y)

        # # run the destretching
        # warped_frames, shifts = torchmfbd.destretch(
        #     frames,
        #     ngrid=ngrid,
        #     lr=lr,
        #     reference_frame=reference_frame,
        #     border=border,
        #     n_iterations=n_iterations,
        #     lambda_tt=lambda_tt,
        # )

        data = np.reshape(warped_frames[0, 0].detach().cpu().numpy(),(2,data.shape[-4],data.shape[-3],data.shape[-2],data.shape[-1]))

        frames = torch.tensor(
            data.reshape(-1, data.shape[-2], data.shape[-1]).astype('float32')
            ).unsqueeze(0).unsqueeze(0)

        warped_frames, _ = torchmfbd.destretch(
            frames, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
            border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
        )

        data = warped_frames[0, 0].detach().cpu().numpy().reshape(data.shape)

    elif aling_cam == 'serie':
        # first thing to do is move data to gpu memory
        # print(data.shape)
        # frames =  torch.tensor(data.astype('float32')) 
        # # add tensor dimensions
        # print(frames.shape)
        # frames = frames.unsqueeze(0).unsqueeze(0)  # (1, 1, nmod, x, y)
        # print(frames.shape)
        # # run the destretching
        # warped_frames, shifts = torchmfbd.destretch(
        #     frames,
        #     ngrid=ngrid,
        #     lr=lr,
        #     reference_frame=reference_frame,
        #     border=border,
        #     n_iterations=n_iterations,
        #     lambda_tt=lambda_tt,
        # )

        data = warped_frames[0, 0].detach().cpu().numpy()

        print(data.shape)
        frames = torch.tensor(data.astype('float32')).unsqueeze(0).unsqueeze(0)
        print(frames.shape)

        warped_frames, _ = torchmfbd.destretch(
            frames, ngrid=ngrid, lr=lr, reference_frame=reference_frame,
            border=border, n_iterations=n_iterations, lambda_tt=lambda_tt,
        )

        data = warped_frames[0, 0].detach().cpu().numpy()
    
    try:
        destretch_pars = destretch_pars.detach().cpu().numpy()
    except:
        print('no destrectch pars')
    torch.cuda.empty_cache()
    gc.collect()

    tac = time.time()

    logging.info(f"Alignment finished in {round(tac - tic, 3)} s.")

    return data, destretch_pars