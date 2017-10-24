function [x,y] = snakedeform_Endocardium(x,y,alpha,beta,gamma,kappa,mu,fx,fy,ITER)
%     alpha:   elasticity parameter
%     beta:    rigidity parameter
%     gamma:   viscosity parameter
%     kappa:   external force weight
%     mu:      shape weight
%     fx,fy:   external force field


N = length(x);

alpha = alpha* ones(1,N); 
beta = beta*ones(1,N);


alpham1 = [alpha(2:N) alpha(1)];    
alphap1 = [alpha(N) alpha(1:N-1)];
betam1 = [beta(2:N) beta(1)];
betap1 = [beta(N) beta(1:N-1)];

a = betam1;                                 % a=beta(i-1)
b = -alpha - 2*beta - 2*betam1;             % b=-alpha(i)-2beta(i)-2beta(i-1)
c = alpha + alphap1 +betam1 + 4*beta + betap1;% c=alpha(i)+alpha(i+1)+beta(i-1)+4beta(i)+beta(i+1)
d = -alphap1 - 2*beta - 2*betap1;            % d=-alpha(i+1)-2beta(i)-2beta(i+1)
e = betap1;                                  % e=beta(i+1)

% generate the parameters matrix  A  
A = diag(a(1:N-2),-2) + diag(a(N-1:N),N-2);
A = A + diag(b(1:N-1),-1) + diag(b(N), N-1);
A = A + diag(c);
A = A + diag(d(1:N-1),1) + diag(d(N),-(N-1));
A = A + diag(e(1:N-2),2) + diag(e(N-1:N),-(N-2));

invAI = inv(A + (mu+gamma) * diag(ones(1,N)));



for count = 1:ITER
   vfx = interp2(fx,x,y,'*linear');
   vfy = interp2(fy,x,y,'*linear');
   idx = isnan(vfx) | isnan(vfy);
    vfx(idx) = 0.0;
    vfy(idx) = 0.0;
    norm_f = sqrt(vfx.*vfx + vfy.*vfy);
    min_f = min(min(norm_f));
    max_f = max(max(norm_f));
    if min_f == max_f
        norm_fn = 0.5;
    else
        k = 1 / (max_f-min_f);
        b = 1 - k*max_f;
        norm_fn = k*norm_f + b; 
    end
    vfx = (vfx ./ (norm_f+1e-5)) .* norm_fn;
    vfy = (vfy ./ (norm_f+1e-5)) .* norm_fn;


% add the circle-shape constraint

xb = mean(x); yb = mean(y);
xbar=repmat(xb,N,1);
ybar=repmat(yb,N,1);
R = sqrt( (x - xb).*(x - xb) + (y- yb).*(y- yb) );
Rbar = mean(R);
Rbar = repmat(Rbar,N, 1);
indx = 2*pi/N*([1:1:N]' - 1);
C = cos(indx); S = sin(indx);
sfx = xbar + Rbar.*C ;
sfy = ybar + Rbar.*S ;

     x = invAI * (gamma* x +kappa*vfx+mu*sfx);
     y = invAI * (gamma* y +kappa*vfy+mu*sfy);
end               
  




