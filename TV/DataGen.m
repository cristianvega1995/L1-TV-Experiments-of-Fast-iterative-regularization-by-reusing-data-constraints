M=256;
N=256;
K=sparse(M*N,M*N);
A=sparse(3*M*N,3*M*N);
x0 = load_image('boat',M);
x0 = rescale( sum(x0,3) );
B=9;
a=1;
b=80;
M=256;
N=256;
D=a*speye(N,N);
for i=1:(N-1)
    D(i,i+1)=b;
    D(i+1,i)=b;
end
D1=(1/(a+3*b))*speye(M,M);
D2=(1/(a+4*b))*speye(M,M);
D1(1,1)=1/(a+2*b);
D2(1,1)=1/(a+3*b);
D1(M,M)=1/(a+2*b);
D2(M,M)=1/(a+3*b);
K1=D1*D;
K2=b*D1;
K3=D2*D;
K4=b*D2;
DY=-speye(N,N);
DY(N,N)=0;
for i=1:(N-1)
    DY(i,i+1)=1;
end
for i=1:N
    if i==1
     K(1:M,1:M)=K1;
     K(1:M,(M+1):2*M)=K2;
    elseif i==N
     K(((N-1)*M+1):M*N,((N-1)*M+1):M*N)=K1;
     K(((N-1)*M+1):M*N,((N-2)*M+1):M*(N-1))=K2;
    else
     K(((i-1)*M+1):M*i,((i-1)*M+1):M*i)=K3;
     K(((i-1)*M+1):M*i,((i-2)*M+1):M*(i-1))=K4;
     K(((i-1)*M+1):M*i,(i*M+1):M*(i+1))=K4;
    end 
end
K5=K;
for i=1:(B-1)
    K5=K*K5;
end
K=K5;
A(1:M*N,1:M*N)=K;
for i=1:N
    if i==1
     A(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-speye(M,M);
     A(((N+i-1)*M+1):M*(i+N), 1:M)=DY;
     
     A(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-speye(M,M);
     A(((2*N+i-1)*M+1):M*(i+2*N),1:M)=-speye(M,M);
     A(((2*N+i-1)*M+1):M*(i+2*N),(M+1):2*M)=speye(M,M);
    elseif i==N
     A(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-speye(M,M);
     A(((N+i-1)*M+1):M*(i+N), ((N-1)*M+1):M*N)=DY;
     
     A(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-speye(M,M);
    else
     A(((N+i-1)*M+1):M*(i+N),((N+i-1)*M+1):M*(i+N))=-speye(M,M);
     A(((N+i-1)*M+1):M*(i+N),((i-1)*M+1):M*i)=DY;
     
     A(((2*N+i-1)*M+1):M*(i+2*N),((2*N+i-1)*M+1):M*(i+2*N))=-speye(M,M);
     A(((2*N+i-1)*M+1):M*(i+2*N),((i-1)*M+1):M*i)=-speye(M,M);
     A(((2*N+i-1)*M+1):M*(i+2*N),(i*M+1):M*(i+1))=speye(M,M);
    end 
end
x=zeros(3*M*N,1);
x(1:M*N)=reshape(x0,M*N,1);
% b=A*x;
% b_delta=b;
% imageplot(reshape(b(1:M*N),M,N))
norma=normest(A)
save('data_gauss')
