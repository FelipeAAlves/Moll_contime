%Written by SeHyoun Ahn - credit crunch extension by Gustavo Mellior
%NEEDS INPUT FROM huggett_initial_creditcrunch.m AND huggett_terminal_creditcrunch.m
clear all; close all; clc;

load huggett_initial_creditcrunch.mat %Equilibrium with lax debt limit
g0 = sparse(g);
gg0 = sparse(gg);
r00 = r;
A0 = A;
AT0 = A';

load huggett_terminal_creditcrunch.mat %Terminal condition
g_st = sparse(g); 
v_st = v;
%Initial and terminal distributions of wealth
plot(a,g0(:,1),'-.','linewidth',3,'color',[0 158/255 1])
hold on
plot(a,g0(:,2),'-.','linewidth',3,'color',[1 0.4 0])
hold on
plot(a,g(:,1),'linewidth',3,'color',[0 0 158/255])
hold on
plot(a,g(:,2),'linewidth',3,'color',[1 0 0])
hold off
legend('g_1(a,t_0)','g_2(a,t_T)','g_1(a,t_T)','g_2(a,t_T)','interpreter','latex', 'location','northeast')
legend('boxoff')
xlim([amin-num*da 0.1]);
ylim([0 3.5]);
xlabel('$a$','interpreter','latex','fontsize',25)
ylabel('$g_i(a)$','interpreter','latex','fontsize',25)
alpha(0.15)
grid on

r_st = r;
clear r Delta;
T = 50;
N = 400;
dt = T/N;
N1=N;


%initial guess of interest rate sequence
r0 = r_st*ones(N,1);

S = zeros(N,1);
SS = zeros(N,1);
SS_tm = zeros(N,1);
dSS = zeros(N,1);

v = zeros(I,2,N);
gg = cell(N+1,1);
v(:,:,N)= v_st;

rnew = r0;

A_t=cell(N,1);

maxit = 600;
r_it=zeros(N,maxit);
Sdist=zeros(maxit,1);
SS_it=zeros(N,maxit);
convergence_criterion = 10^(-5);

%If you want to color code plots to check algorithm is converging this will
%create a matrix of RGBs. The algorithm will run faster if you shut off the
%plots
figure(1)
itr = maxit;
CM = jet(itr);

%Initialize the speed of updateing the interest rate.
xi=5; 
for it=1:maxit
    disp('ITERATION = ')
    disp(it)
    r_t = rnew;
    r_it(:,it)=r_t;
    
    V = v_st;
    
    for n=N1:-1:1
        v(:,:,n)=V;
        % forward difference
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVf(I,:) = (z + r_t(n).*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
        % backward difference      
        dVb(num+2:I,:) = (V(num+2:I,:)-V(num+1:I-1,:))/da;
        dVb(num+1,:) = (z + r_t(n).*amin).^(-s); %state constraint boundary condition as before
            
        %Consumption, derivative of value function at steady state and new
        %state constraints
        c0 = zz + r_t(n).*aa;
        for j=num:-1:1
            dVb(num+1-j,:) = (z + r_t(n).*(amin-j*da)-1*da).^(-s); %new state constraints
            c0(num+1-j) = c0(num+1-j)-1*da;
        end
        %Households in the newly created inadmissible region
        %save "delta a" every period until they reach the new debt limit.
        dV0 = c0.^(-s);
        
        I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
                
        %consumption and savings with forward difference
        cf = dVf.^(-1/s);
        ssf = zz + r_t(n).*aa - cf;
        %consumption and savings with backward difference
        cb = dVb.^(-1/s);
        ssb = zz + r_t(n).*aa - cb;
        
        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift
        If = ssf > 0; %positive drift --> forward difference
        Ib = ssb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state
        I0(I0<0)=0;
        %make sure backward difference is used at amax
        Ib(I,:) = 1; If(I,:) = 0;
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
        
        %%Note the syntax for the cell array
        A_t{n} = A; %Remember we are still going from N to 1, backwards
        B = (1/dt + rho)*speye(2*I) - A;
        
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
        
        V = [V_stacked(1:I),V_stacked(I+1:2*I)];
        ss = zz + r_t(n).*aa - c;
    end
    
    
    gg{1}=gg0;
    for n=1:N 
        AT=A_t{n}';
        %Implicit method in Updating Distribution.
        gg{n+1}= (speye(2*I) - AT*dt)\gg{n};
        %gg{n+1}=gg{n}+AT*gg{n}*dt; %This is the explicit method.
        %check(n) = gg(:,n)'*ones(2*I,1)*da;
        SS(n) = gg{n}(1:I)'*a*da + gg{n}(I+1:2*I)'*a*da;
    end
    
    SS_it(:,it)=SS;
    
    %Compute the change in the aggregate savings.
    SS(1)=0;    %This was done to increase the accuracy by a little (Otherwise dSS has a bigger noise at t=0).
    
    %Forward and backward approximations are used alternately because once
    %the excess savings level becomes small enough, the "rounding error" can
    %lead the update in the interest rate, and repeatedly adding the same
    %rounding error can lead to a propagating error.
    if mod(it,2)==0
        dSS(1:N-1)=SS(2:N)-SS(1:N-1);
        dSS(N)=SS(N)-SS(N-1);
    else
        dSS(2:N)=SS(2:N)-SS(1:N-1);
        dSS(1)=dSS(2);
    end
    
    %Update the interest rate to reduce aggregate savings amount.
    rnew = r_t - xi'.*dSS;
    
    rnew(N-5:N)=rnew(N-5); %This was done to minimize the "rounding error" at the end. Otherwise, the end point keeps decreasing

    %To improve speed, for the first few updates, the update will be done
    %fast, but to reduce wave pattern that gets created, updated interest
    %rate will be smoothed (since the wave patern is due to how things are
    %being updated). After getting "good" initial starting r_t, standard update
    %with declining update weight will be used for convergence. Note that the
    %smoothing will be stopped before the convergence.
    if it<20
        %choose is just a parameter so that you can select how far to the
        %left do you want to smooth the interest rate series.
        choose = 100;
        rnew(N-100:N)=r_st;                 
        rnew(choose:N)=smooth(rnew(choose:N),50);
    elseif it<60
        rnew(N-10:N)=r_st;
        rnew(choose:N)=smooth(rnew(choose:N),10);
    elseif it<100
        rnew(N-5:N)=r_st;
        rnew(choose:N)=smooth(rnew(choose:N),5);
    elseif it==100
        xi=10*exp(-0.0006*(1:N));
    elseif it==150
        xi=10*exp(-0.00006*(1:N)); %Give bigger update weight to later times.
    end
    
    %This can be removed for speed.
    if mod(it,50)==0
        %figure(1);
        %clf;
        subplot(2,2,1);
        plot(r_st*ones(N,1),'r--');
        hold on;
        title('r_t');
        plot(r_t,'color',CM(it,:));
        subplot(2,2,2);
        plot(SS);
        hold on
        plot(zeros(1,400),'r--')
        title('SS');
        subplot(2,2,3);
        plot(zeros(1,20),'r--');
        xlim([0 20])
        hold on;
        plot(SS(1:20));
        title('SS(1:20)');
        subplot(2,2,4);
        plot(SS(N-20:N));
        xlim([0 20])
        hold on;
        title('SS(N-20:N)');
        hold on
        pause(0.1);
    end
    
    Sdist(it) = max(abs(SS));
    if Sdist(it)<convergence_criterion
        break
    end
    
end

time = (1:N)'*dt;

save huggett_transition_creditcrunch.mat

figure(1);
clf;
plot(Sdist(1:it))

figure(2);
plot(1:N,SS_it(:,1),1:N,SS_it(:,it))
legend('Iteration 1','Last Iteration')
ylabel('Excess Supply')
xlabel('Time Period')

T1 = -0.5;
N1 = -T1/dt;
time1 = T1 + (1:N1)'*dt;
time2 = [time1;time];
r_t2 = [r00*ones(N1,1);r_t];


amax1 = 0.5;
gmax = 3;
set(gcf,'PaperPosition',[0 0 15 10])
n=2;
figure(3)
subplot(2,2,1)
set(gca,'FontSize',16)
h1 = plot(a,gg{n}(1:I),'b',a,gg{n}(I+1:2*I),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
hold on
plot(a,g0(:,1),'b--',a,g0(:,2),'r--','LineWidth',2)
xlim([amin amax1])
ylim([0 gmax])
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
title('t = 0.1')

t = 2;
n = t/dt;
subplot(2,2,2)
set(gca,'FontSize',16)
h1 = plot(a,gg{n}(1:I),'b',a,gg{n}(I+1:2*I),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
hold on
plot(a,g0(:,1),'b--',a,g0(:,2),'r--','LineWidth',2)
xlim([amin amax1])
ylim([0 gmax])
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
title('t = 2')

t = 5;
n = t/dt;
subplot(2,2,3)
set(gca,'FontSize',16)
h1 = plot(a,gg{n}(1:I),'b',a,gg{n}(I+1:2*I),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
hold on
plot(a,g0(:,1),'b--',a,g0(:,2),'r--','LineWidth',2)
xlim([amin amax1])
ylim([0 gmax])
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
title('t = 5')


t = T;
n = t/dt;
subplot(2,2,4)
set(gca,'FontSize',16)
h1 = plot(a,gg{n}(1:I),'b',a,gg{n}(I+1:2*I),'r','LineWidth',2)
legend(h1,'g_1(a)','g_2(a)')
hold on
plot(a,g0(:,1),'b--',a,g0(:,2),'r--','LineWidth',2)
xlim([amin amax1])
ylim([0 gmax])
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
title('t = \infty')
print -depsc transition_distribution.eps


figure(4)
plot(time2,r_t2,time2,r_st*ones(N1+N,1),'k--','LineWidth',3)
xlim([T1 10])
ylim([0 0.035])
ylabel('$r$','interpreter','latex','fontsize', 25)
grid on
xlabel('Year','interpreter','latex','fontsize', 25)
alpha(0.15)

