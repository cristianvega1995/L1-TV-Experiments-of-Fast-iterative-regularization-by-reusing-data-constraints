
function [BD,Feas,Dist,ValP,EST,ESI,Time_Tol,Iter_Tol,Time,error]  = PROJ_primalandaptative(xk,pk,p_old,uk,lambda,MaxIt,A,b_delta,b,x_exact,mu,normA,tol,factor)
ESI=0;% Iteration number
EST=0;% Execution time
Time_Tol=0; %Time until reach the tolerance
Iter_Tol=0; %Iteration until reach the tolerance
Time=0;%Total time
error=norm(xk-x_exact);
BD=zeros(MaxIt,1);%Bregman distance % BD_ERG=zeros(MaxIt,1);
Feas=zeros(MaxIt,1);%Feasibility % Feas_ERG=zeros(MaxIt,1);
Dist=zeros(MaxIt,1);%Norm error % Dist_ERG=zeros(MaxIt,1);
ValP=zeros(MaxIt,1);%Function value
alphita1=0;%Initial step sizes for landweber
N=length(A(1,:)); M=length(A(:,1));%Matrix size

p_bar=pk+xk-p_old; % Initialization
tic;
uk=A'*uk;
tau=lambda;
for k=1:MaxIt
        x_aux=pk;
        for i=1:N
            if x_aux(i,1)<-lambda
                xk(i,1)=x_aux(i)+lambda;
            elseif x_aux(i,1)>lambda
                xk(i,1)=x_aux(i)-lambda;
            else
                xk(i,1)=0;
            end
        end
        if norm(xk-x_exact)>tol
            Time_Tol=toc;
            Iter_Tol=k;
        end
        if norm(xk-x_exact)<error
            EST=toc;
            ESI=k;
            error=norm(xk-x_exact);
        end
        ValP(k,1)=abs(sum(abs(xk))-sum(abs(x_exact)));%Function value
        BD(k,1)=abs(sum(abs(xk))-sum(abs(x_exact))+ mu'*(A*xk-b_delta)); %Bregman distance
        Feas(k,1)=norm(A*xk-b); % Feasibility
        Dist(k,1)=norm(xk-x_exact); % Reconstruction error
        pk=xk;
        p1=A*pk-b_delta;
        cof1=norm(p1)^2;
        p2=lambda*p1;
        p2=tau*(A'*p2);
        uk=uk+p2;
        p3=A'*p1;
        cof2=norm(p3)^2;
        alphita=cof1/cof2;%Adaptive step
        pk = pk-factor*alphita*p3;%Landweber step
        p4=A*pk-b_delta;
        p4=lambda*p4;
        p4=tau*A'*p4;
        pk=pk-p4-uk;
end
 Time=toc;%Execution time