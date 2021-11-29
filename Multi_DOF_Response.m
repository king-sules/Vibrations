
%Code for multi-degree of freedom
clear all
format short
close all
clc
syms t tau

%Set Plotting Time and interval
tMax=100;
inter=0.1;

%input mass matrix
m = 10; %we can change this accordingly 
M = [1,0,0;
    0,2,0;
    0,0,3] * m;
n=length(M);
M_inv= inv(M);

%input stiffness matrix
k= 20; %we can also change this
K= [7,-3,0;
    -3,5,-2;
    0,-2,2] *k;

%Damping Ratio Coeddicients
a=0.403953871462269;
b=0.0886325897013889;
%% 
%Now for eigenvalue problem
MK = M\K;

[V,D] = eig(MK);

[eig_MK,ind] = sort(diag(D));
V=V(:,ind);

%% 
% Lets also set initial conditions
%for any external forces
Fg=[5;0;0];
wf=[2;0;0];
funct='cos'; %Is force in sin or cos function?

% IC and in principle space
x0= [0;0;0];
v0 = [0;0;0];

qx0= V'* M * x0;
qv0 = V'* M * v0;

%% 
%Now for any degree in parellel with M matrix

if strcmp(funct,'sin')
    for i=1:n
        w(i) = sqrt(eig_MK(i));

        % This is damping ratio
        z(i)= (a+(b*(w(i)^2)))/(2*w(i));
        wd(i)= sqrt(1-z(i)^2)*w(i);
        Ft(i) = Fg(i)*sin(wf(i)*t);

    end

elseif strcmp(funct,'cos')
    for i=1:n
        w(i) = sqrt(eig_MK(i));

        % This is damping ratio
        z(i)= (a+(b*(w(i)^2)))/(2*w(i));
        wd(i)= sqrt(1-z(i)^2)*w(i);
        Ft(i) = Fg(i)*cos(wf(i)*t);

    end

else
    'Please input only "sin" or "cos" functions'
    return
    
end

z=z';
wd=wd';
w=w';
Ft=Ft';
%% 
%Computinng principle force
Qt= V'*Ft;

%This section changes finds q(t) for all modes;
for i=1:n
    Qf(t) = Qt(i);

    expr = Qf(tau)*exp(-z(i)*w(i)*(t-tau))*sin(wd(i)*(t-tau));
    F=matlabFunction(int(expr,tau,0,t));
    fun{i} =@(t)  exp(-z(i)*w(i)*t) * (qx0(i)*cos(wd(i)*t) + ((z(i)*qx0(i))/(sqrt(1-z(i)^2)))*sin(wd(i)*t) + (qv0(i)/wd(i))*sin(wd(i)*t));
    tempFunc=fun{i};
    qt{i}=@(t) tempFunc(t) + F(t)/wd(i);

end


qt=qt';

%% 
%Finally get the response with increasing time
tspan=0:inter:tMax;

for i=1:n
    for j=1:numel(tspan)
        xt(i,j)=0; 
        for k=1:n
              xt(i,j)=V(i,k)*qt{k}(inter*j)+xt(i,j);
        end
    end
end


%% 
%Last we plot our response
i=1;
for i=1:n
    figure(i); 
    plot(tspan,xt(i,:));
    xlabel('Time (s)')
    ylabel('Displacement')
    title(['Mass ',num2str(i), ' Response'])
end