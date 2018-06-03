function [P] = unmix(data, endmembers)

warning('off', 'all');
options = optimset('Display', 'off', 'LargeScale', 'off');

%endmembers should be column vectors
X = data;

%number of endmembers
M = size(endmembers, 2);
%number of pixels
N = size(X, 2);

%Equation constraint Aeq*x = beq
%All values must sum to 1 (X1+X2+...+XM = 1)
Aeq = ones([1, M]);
beq = 1;

%Boundary Constraints lb >= x >= ub
%All values must be greater than 0 (0 ? X1,0 ? X2,...,0 ? XM)
lb = zeros([M, 1]);
ub = ones([M,1]);

H = 2*(endmembers'*endmembers);
P = zeros(N,M);

for i = 1:N
    F = ((-2*X(:,i)'*endmembers))';
    qp_ans = quadprog(H, F, [], [], Aeq, beq, lb, ub, 0, options);
    P(i,:) = qp_ans;
end

P(P<0) = 0;

