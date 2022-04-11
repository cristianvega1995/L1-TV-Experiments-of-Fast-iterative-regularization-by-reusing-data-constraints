
function [BD,Feas,Dist,ValP,EST,ESI,Time_Tol,Iter_Tol,Time,error] = PROJ_primalparallel(xk,pk,p_old,uk,lambda,MaxIt,A,b_delta,b,x_exact,mu,normA,Af,tol)
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
alphita=1/norm(A,'fro')^2;%Step-sizes for parallel projectios
N=length(A(1,:)); M=length(A(:,1));%Matrix size
p_bar=pk+xk-p_old;  % Initialization
tic;
uk=A'*uk;
tau=lambda;%Step sizes
for k=1:MaxIt
         x_aux=pk;      
%     
%         uk=uk+lambda*(A*p_bar-b_delta);
%         
%         x_aux=pk-lambda*A'*uk;
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
        ValP(k,1)=abs(sum(abs(xk))-sum(abs(x_exact)));
        BD(k,1)=abs(sum(abs(xk))-sum(abs(x_exact))+sum(A'*mu.*(xk-x_exact)));% BregDivP_ERG(k,1)=sum(abs(PPk))-sum(abs(x))+sum(A'*mu.*(PPk-x));
        Feas(k,1)=norm(A*xk-b); % FeasP_ERG(k,1)=norm(A*XbarPk-b);
        Dist(k,1)=norm(xk-x_exact); % DistP_ERG(k,1)=norm(PPk-x);
        pk=xk;
        p1=A*pk-b_delta;
        p2=lambda*p1;
        p2=tau*(A'*p2);
        uk=uk+p2;
        p3=A'*p1;
        pk = pk-alphita*p3; %Parallel projections
        p4=A*pk-b_delta;
        p4=lambda*p4;
        p4=tau*A'*p4;
        pk=pk-p4-uk;      
end
Time=toc;%Execution time