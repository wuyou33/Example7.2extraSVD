function [X,DP,QP] = svd2017(z,R,H,X,DP,QP)
% ------------------------------------------------------------------- 
% Riccati_KF_SVD: SVD factorization-based implementation: Copyright(C) 2017 
% Authors: Maria Kulikova:       maria dot kulikova at ist dot utl dot pt 
% License: GNU General Public License Version 2 
%
% References:
%   Kulikova M.V., Tsyganova J.V. (2017) "Improved discrete-time Kalman 
%   filtering within singular value decomposition". 
%   IET Control Theory & Applications, 11(15): 2412-2418
%   DOI:   10.1049/iet-cta.2016.1282 
% ------------------------------------------------------------------- 
      [m,n]     = size(H);
      [~,DR,RQ] = svd(R);                        % SVD for measurement noise covariance
    PreArray    = [DR.^(1/2)*RQ'; DP.^(1/2)*QP'*H';];
     [~,S,QRek] = svd(PreArray);                 % SVD factors of R_{e,k}  
           DRek = S(1:m,1:end)^2;                % SVD factors of R_{e,k}

    norm_residual =  QRek'*(z-H*X);              % normalized residual
       norm_Gain  = (QP*DP*QP')*H'*QRek;         % normalized gain
        inv_DRek  =  DRek\eye(size(DRek));       % inverse once
                X = X + norm_Gain*inv_DRek*norm_residual; 

     Gain  = norm_Gain*inv_DRek*QRek';            % Gain
       PreArray  = [DP.^(1/2)*QP'*(eye(n) - Gain*H)'; DR.^(1/2)*RQ'*Gain';];
      [~,S,QP] = svd(PreArray);                   % Filtered SVD factors of P        
            DP = S(1:n,1:end).^2;                 % Filtered SVD factors of P

return;
