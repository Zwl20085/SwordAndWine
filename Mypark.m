%   Mypark
%psiA M x N




%%MultiState
if (Mode == 1||Mode == 3)
AngleStep = 2;
AngleFinal = 120;
end
if(Mode == 2)
AngleStep = 6;
AngleFinal = 120;
end
if(Mode == 4)
AngleStep = 1;
AngleFinal = 120;
end


%%
P = @(x)[cosd(x) cosd(x-120) cosd(x+120);
    -sind(x) -sind(x-120) -sind(x+120);
    0.5 0.5 0.5];
Angle = 0:AngleStep:AngleFinal;

[M N] = size(psiA);
for j = 1:1:N
    for i = 1:1:M
        %[psiD(i,j) psiQ(i,j) psi0(i,j)]= (P((i-1)*AngleStep)*([psiA(i,j) psiB(i,j) psiC(i,j)]'))';
        K= (P((i-1)*AngleStep)*([psiA(i,j) psiB(i,j) psiC(i,j)]'))';
        psiD(i,j) = -K(1);
        psiQ(i,j) = -K(2);
        psiO(i,j) = -K(3);
    end
end


