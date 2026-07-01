function k = dpim_gaussian_kernel(x, h)
%DPIM_GAUSSIAN_KERNEL Gaussian kernel K_h(x).
k = exp(-0.5*(x./h).^2) ./ (sqrt(2*pi)*h);
end
