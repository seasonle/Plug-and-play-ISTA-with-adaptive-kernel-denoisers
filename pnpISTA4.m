function [x_curr,err,iters,e,P] = ...
                    pnpISTA4(x0,gradf,W,delta,I,tol,maxiters,H2)
%PNPISTA Plug-and-play ISTA (Proximal gradient method)
% x0 = Initial point
% gradf = Handle to gradient function of the data fidelity term
% H = Handle to denoiser
% delta = Step-size for the descent step
% I = Ground-truth image to compare PSNR (optional), can be empty ([])
% tol = Tolerance (optional)
% maxiters = Max. no. of iterations (optional)
% x_curr = Output signal/image
% converged = Boolean flag indicating if algorithm converged
% iters = No. of iterations completed
% e = Iteration-wise errors
% P = Iteration-wise PSNRs
%

if(~exist('I','var') || isempty(I))
    calcPSNR = false;
else
    calcPSNR = true;
end
if(~exist('tol','var') || isempty(tol))
    tol = 1/255;
end
if(~exist('maxiters','var') || isempty(maxiters))
    maxiters = 30;
end

iters = 1;
x_prev = x0;
e = nan(1,maxiters);
if(calcPSNR)
    P = nan(1,maxiters);
%     P(1) = psnr(x0,I,1);
    P(1) = psnr(x0,I);
end
e(1)=euclNorm(x0 - I);
while(true)
    
    
    v_curr = x_prev - delta*gradf(x_prev);
    x_curr = W(v_curr);
    fprintf('Iteration %d completed',iters);
    if(calcPSNR)
%         P(iters+1) = psnr(x_curr,I,1);
        P(iters+1) = psnr(x_curr,I);
        fprintf(', PSNR = %f\n',P(iters+1));
    end
    err = euclNorm(x_curr - x_prev);
    e(iters+1) = err;
    
    if(iters==maxiters)
   
        break;
    end
    iters = iters+1;
    x_prev = x_curr;
end

e(iters+1:end) = [];
if(calcPSNR)
    P(iters+1:end) = [];
else
    P = [];
end

end


function n = euclNorm(y)
% Euclidean norm of vector or matrix

n = sqrt(sum(y(:).*y(:)));

end
