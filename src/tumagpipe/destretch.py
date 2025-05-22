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


def destretch(data, ngrid=2, lr=0.50, reference_frame=0, border=6, n_iterations=200, lambda_tt=0.01,aling_cam='partial'):
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

    if len(data.shape) == 4:
        data = data[:, np.newaxis] # To allow for only one lamdba.
    elif len(data.shape) == 5:
        print('data shape is correct')
    else:
        raise ValueError("Data must be of shape (Ncams, Nlambda, Nmods, Nx, Ny) or (Ncams, Nmods, Nx, Ny)")
    
    # if filterflag:
    #     data = filter_frecuencies(data, verbose=verbose)

    # for the aligment between cameras, it is best to balance the intensities before
# camera_balance = np.median(d1[0,0,:,size:-size,size:-size]) / np.median(d1[1,0,:,size:-size,size:-size])
# print('camera balance', camera_balance)

    # for lambd in range(1):
    for lambd in range(data.shape[1]):
        print('aligning lambda', lambd)
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
            ngrid=2,
            lr=lr,
            reference_frame=reference_frame,
            border=border,
            n_iterations=n_iterations,
            lambda_tt=lambda_tt,
        )
        warped_cam1_frames, shift_1 = torchmfbd.destretch(
            cam1_frames,
            ngrid=2,
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
        data[0,lambd,:,:,:] = warped_cam0_frames[0, 0].detach().cpu().numpy()
        data[1,lambd,:,:,:] = warped_cam1_frames[0, 0].detach().cpu().numpy()

        # now between the cameras

        # two optians, all four (2 to 2) or just m1 and apply to the rest 
        if aling_cam=='full':
            dm0m_0 = torch.tensor(data[0,lambd,0,:,:].astype('float32'))  # (x, y)
            dm1m_0 = torch.tensor(data[1,lambd,0,:,:].astype('float32'))  # (x, y)
            dm0m_1 = torch.tensor(data[0,lambd,1,:,:].astype('float32'))  # (x, y)
            dm1m_1 = torch.tensor(data[1,lambd,1,:,:].astype('float32'))  # (x, y)
            dm0m_2 = torch.tensor(data[0,lambd,2,:,:].astype('float32'))  # (x, y)
            dm1m_2 = torch.tensor(data[1,lambd,2,:,:].astype('float32'))  # (x, y)
            dm0m_3 = torch.tensor(data[0,lambd,3,:,:].astype('float32'))  # (x, y)
            dm1m_3 = torch.tensor(data[1,lambd,3,:,:].astype('float32'))  # (x, y)

            dm_stack_0 = torch.stack([dm0m_0, dm1m_0], dim=0).unsqueeze(0).unsqueeze(0)
            dm_stack_1 = torch.stack([dm0m_1, dm1m_1], dim=0).unsqueeze(0).unsqueeze(0)
            dm_stack_2 = torch.stack([dm0m_2, dm1m_2], dim=0).unsqueeze(0).unsqueeze(0)
            dm_stack_3 = torch.stack([dm0m_3, dm1m_3], dim=0).unsqueeze(0).unsqueeze(0)

            warped_mod0, _ = torchmfbd.destretch(
                dm_stack_0,
                ngrid=ngrid,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )
            warped_mod1, _ = torchmfbd.destretch(
                dm_stack_1,
                ngrid=ngrid,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )
            warped_mod2, _ = torchmfbd.destretch(
                dm_stack_2,
                ngrid=ngrid,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )
            warped_mod3, _ = torchmfbd.destretch(
                dm_stack_3,
                ngrid=ngrid,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )

            data[:,lambd,0,:,:] = warped_mod0.detach().cpu().numpy()
            data[:,lambd,1,:,:] = warped_mod1.detach().cpu().numpy()
            data[:,lambd,2,:,:] = warped_mod2.detach().cpu().numpy()
            data[:,lambd,3,:,:] = warped_mod3.detach().cpu().numpy()

        if aling_cam=='partial':

            dm0m_0 = torch.tensor(data[0,lambd,0,:,:].astype('float32'))  # (x, y)
            dm1m_0 = torch.tensor(data[1,lambd,0,:,:].astype('float32'))  # (x, y)

            dm_stack_0 = torch.stack([dm0m_0, dm1m_0], dim=0).unsqueeze(0).unsqueeze(0)

            warped_mod0, distortion_map = torchmfbd.destretch(
                dm_stack_0,
                ngrid=ngrid,
                lr=lr,
                reference_frame=reference_frame,
                border=border,
                n_iterations=n_iterations,
                lambda_tt=lambda_tt,
            )

            for i in range(2):
                for j in range(2):
                    print('shifts cam= ',i,'(x,y) = (0,1) ',j,' ',distortion_map[i,j,data.shape[-1]//2,data.shape[-1]//2],distortion_map[i,j,data.shape[-1]//2,data.shape[-1]//2])
                    
            data[:,lambd,0,:,:] = warped_mod0.detach().cpu().numpy()

            #apply to QUV
            data[:,lambd,1,:,:] = torchmfbd.apply_destretch(torch.tensor(data[:,lambd,1,:,:].astype('float32')).unsqueeze(0).unsqueeze(0), distortion_map, mode='bilinear').detach().cpu().numpy()
            data[:,lambd,2,:,:] = torchmfbd.apply_destretch(torch.tensor(data[:,lambd,2,:,:].astype('float32')).unsqueeze(0).unsqueeze(0), distortion_map, mode='bilinear').detach().cpu().numpy()
            data[:,lambd,3,:,:] = torchmfbd.apply_destretch(torch.tensor(data[:,lambd,3,:,:].astype('float32')).unsqueeze(0).unsqueeze(0), distortion_map, mode='bilinear').detach().cpu().numpy()



    tac = time.time()

    print(f"Alignment finished in {round(tac - tic, 3)} s.")

    return data