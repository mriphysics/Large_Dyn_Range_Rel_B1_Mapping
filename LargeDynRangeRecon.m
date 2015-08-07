function [relTx,sig_recon,N_samps] = LargeDynRangeRecon(input_ims,mask,v,cmplx_noi)
%
%   FUNCTION: LargeDynRangeRecon - reconstructs relative B1+ Maps from
%   multi-voltage multi-transmit channel SPGR data.
%
%   INPUTS:     input_ims   Nx x Ny x Nz x Nv x Nt complex image matrix
%                           Nx, Ny & Nz = number of voxels in X, Y & Z directions
%                           Nv = number of voltages
%                           Nt = number of transmit channels
%
%               mask        Nx x Ny x Nz binary matrix of signal regions
%
%               v           Nv x 1 vector of drive levels of each image
%
%               cmplx_noi   Matrix of complex IMAGE-DOMAIN noise samples
%
%   OUTPUTS:    relTx       Nx x Ny x Nz x Nt complex matrix of relative B1+ maps
%
%               sig_recon   Nx x Ny x Nz x Nt complex matrix of
%                           resonstructed signal intensities
%
%               N_samps     Nx x Ny x Nz x Nt matrix of number of samples
%                           used
%
%   Created by: Dr Francesco Padormo 
%   Email: francesco.padormo@kcl.ac.uk
%   Date: 7th August 2015
%
%   This code is free under the terms of the MIT license. Please cite
%   Padormo et al. MRM (2015) doi:10.1002/mrm.25884 in any published work.
%
%%

sz_i = size(input_ims);
v = v(:);

% Normalise Data to max intensity of 1
sc_fact = max(abs(input_ims(:)));
input_ims = input_ims/sc_fact;
cmplx_noi = cmplx_noi/sc_fact;

% Calulcate standard deviation of noise (assume its the same on each channel)
std_dev = std([real(cmplx_noi(:));imag(cmplx_noi(:))]);

% Pre-allocate variables
sig_recon = zeros([sz_i(1:3) sz_i(5)]);
N_samps = zeros([sz_i(1:3) sz_i(5)]);

for tt = 1:sz_i(5)
    for zz = 1:sz_i(3)
        for yy = 1:sz_i(2)
            for xx = 1:sz_i(1)
                if mask(xx,yy,zz)==1
                    
                    % Normalised signal vector
                    dv = squeeze(input_ims(xx,yy,zz,:,tt))./v;
                    % Normalised noise level
                    std_norm = (std_dev./v)';                    
                    
                    % Calculate best estimate for each possible 'number of
                    % estimates', and their corresponding likelihood
                    for kk = 1:sz_i(4)
                        
                        best_est_r(kk) = sum(real(dv(1:kk))./(std_norm(1:kk)'.^2))./sum(1./(std_norm(1:kk)'.^2));
                        best_est_i(kk) = sum(imag(dv(1:kk))./(std_norm(1:kk)'.^2))./sum(1./(std_norm(1:kk)'.^2));
                        
                        t1(kk) = -kk*log(2*pi);
                        t2(kk) = -2*sum(log(std_norm(1:kk)));
                        t3(kk) = -0.5*sum(((real(dv(1:kk))-best_est_r(kk))./std_norm(1:kk)').^2);
                        t4(kk) = -0.5*sum(((imag(dv(1:kk))-best_est_i(kk))./std_norm(1:kk)').^2);
                        
                        likelihood_best_est(kk) = t1(kk) + t2(kk) + t3(kk) + t4(kk);
                        
                    end
                    
                    % Find which kk has highest likelihood
                    [~,best_kk] = max(likelihood_best_est);
                    
                    % Calculate corresponding complex signal
                    sig_recon(xx,yy,zz,tt) = best_est_r(best_kk) + 1i*best_est_i(best_kk);
                    
                    % Create map of 'number of samples used' 
                    N_samps(xx,yy,zz,tt) = best_kk;
                    
                end
                
            end
        end
    end
end

relTx = sig_recon./repmat(sum(abs(sig_recon),4),[1 1 1 sz_i(5)]);

end