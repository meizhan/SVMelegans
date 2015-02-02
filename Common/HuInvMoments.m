function H=HuInvMoments(I)
% H=HuInvMoments(I) calculates the seven Hu invariant moments of input image
% I in an output vector H.

[M, N]=size(I);
[x,y]=meshgrid(1:N, 1:M);
x=x(:);
y=y(:);
I=double(I(:));

mu00=sum(I);
xbar=sum(I.*x)/mu00;
ybar=sum(I.*y)/mu00;
powers=[1 1; 1 2; 2 1; 0 2; 2 0; 0 3; 3 0];


%Calculate Normalized Central Moments
for i=1:length(powers)
    p=powers(i,1);
    q=powers(i,2);
    eval(['mu' num2str(p) num2str(q) '=sum(I.*((x-xbar).^p.*(y-ybar).^q));']);
    eval(['eta' num2str(p) num2str(q) '=' 'mu' num2str(p) num2str(q) '/mu00^((p+q)/2+1);']);
end
%Calculate Hu Invariants
H(1)=eta20+eta02;
H(2)=(eta20-eta02)^2+4*eta11^2;
H(3)=(eta30-3*eta12)^2+(3*eta21-eta03)^2;
H(4)=(eta30+eta12)^2+(eta21+eta03)^2;
H(5)=(eta30-3*eta12)*(eta30+eta12)*((eta30+eta12)^2-3*(eta21+eta03)^2)...
    +(3*eta21-eta03)*(eta21+eta03)*(3*(eta30+eta12)^2-(eta21+eta03)^2);
H(6)=(eta20-eta02)*((eta30+eta12)^2-(eta21+eta03)^2)...
    +4*eta11*(eta30+eta12)*(eta21+eta03);
H(7)=(3*eta21-eta03)*(eta30+eta12)*((eta30+eta12)^2-3*(eta21+eta03)^2)...
    -(eta30-3*eta12)*(eta21+eta03)*(3*(eta30+eta12)^2-(eta21+eta03)^2);
end



