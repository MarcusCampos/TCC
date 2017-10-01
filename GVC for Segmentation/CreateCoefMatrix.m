function A = CreateCoefMatrix(N, alpha, beta)
% Generates the parameters for snake
    alpha = repmat(alpha, N, 1);
    beta = repmat(beta, N, 1);
    % Produce the five diagnal vectors
    alpham1 = [alpha(2:N); alpha(1)];
    alphap1 = [alpha(N); alpha(1:N-1)];
    betam1 = [beta(2:N); beta(1)];
    betap1 = [beta(N); beta(1:N-1)];

    a = betam1;
    b = -alpha - 2*beta - 2*betam1;
    c = alpha + alphap1 +betam1 + 4*beta + betap1;
    d = -alphap1 - 2*beta - 2*betap1;
    e = betap1;

    % Generate the parameters matrix
    A = diag(a(1:N-2), -2) + diag(a(N-1:N), N-2);
    A = A + diag(b(1:N-1), -1) + diag(b(N), N-1);
    A = A + diag(c);
    A = A + diag(d(1:N-1), 1) + diag(d(N), -(N-1));
    A = A + diag(e(1:N-2), 2) + diag(e(N-1:N), -(N-2));
