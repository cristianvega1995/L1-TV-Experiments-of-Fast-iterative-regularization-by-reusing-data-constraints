function [norma] = NormA(a,b,N,M)
K=sparse(3*M*N,3*M*N);
a=1;
b=80;
M=256;
N=256;
D=a*eye(N,N);
for i=1:(N-1)
    D(i,i+1)=b;
    D(i+1,i)=b;
end
D1=(1/(a+3*b))*eye(M,M);
D2=(1/(a+4*b))*eye(M,M);
D1(1,1)=1/(a+2*b);
D2(1,1)=1/(a+3*b);
D1(M,M)=1/(a+2*b);
D2(M,M)=1/(a+3*b);
K1=D1*D;
K2=b*D1;
K3=D2*D;
K4=b*D2;
DX=-eye(N,N);
DX(N,N)=0;
for i=1:(N-1)
    DX(i,i+1)=1;
end
for i=1:N
    if i==1
     K(1:M,1:M)=K1;
     K(1:M,(M+1):2*M)=K2;
     
     K(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-eye(M,M);
     K(((N+i-1)*M+1):M*(i+N), 1:M)=DX;
     
     K(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-eye(M,M);
     K(((2*N+i-1)*M+1):M*(i+2*N),1:M)=-eye(M,M);
     K(((2*N+i-1)*M+1):M*(i+2*N),(M+1):2*M)=eye(M,M);
    elseif i==N
     K(((N-1)*M+1):M*N,((N-1)*M+1):M*N)=K1;
     K(((N-1)*M+1):M*N,((N-2)*M+1):M*(N-1))=K2;
     
     K(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-eye(M,M);
     K(((N+i-1)*M+1):M*(i+N), ((N-1)*M+1):M*N)=DX;
     
      K(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-eye(M,M);
    else
     K(((i-1)*M+1):M*i,((i-1)*M+1):M*i)=K3;
     K(((i-1)*M+1):M*i,((i-2)*M+1):M*(i-1))=K4;
     K(((i-1)*M+1):M*i,(i*M+1):M*(i+1))=K4;
     
     K(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-eye(M,M);
     K(((N+i-1)*M+1):M*(i+N),((i-1)*M+1):M*i)=DX;
     
     K(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-eye(M,M);
     K(((2*N+i-1)*M+1):M*(i+2*N),((i-1)*M+1):M*i)=-eye(M,M);
     K(((2*N+i-1)*M+1):M*(i+2*N),(i*M+1):M*(i+1))=eye(M,M);
    end 
end
norma=normest(K);
end