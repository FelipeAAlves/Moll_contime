%INITIAL GUESS
v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;

% z_ave = la2/(la1+la2)*z(1) + la1/(la1+la2)*z(2);
% v0(:,1) = (z_ave + r.*a).^(1-s)/(1-s)/rho;
% v0(:,2) = (z_ave + r.*a).^(1-s)/(1-s)/rho;

v = v0;
tic;
for n=1:1000
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (z + r.*amin).^(-s); %state constraint boundary condition

    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)

    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    ssf = zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    ssb = zz + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = zz + r.*aa;
    dV0 = c0.^(-s);

    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically

    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);

    %CONSTRUCT MATRIX
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;

    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B = (rho + 1/Delta)*speye(2*I) - A;

    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];

    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS

    V = [V_stacked(1:I),V_stacked(I+1:2*I)];

    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    % if dist(n)<crit
    %     disp('Value Function Converged, Iteration = ')
    %     disp(n)
    %     break
    % end
end
toc;
