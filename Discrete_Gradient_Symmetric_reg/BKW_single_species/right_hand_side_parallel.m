function [U_x,U_y,dissipation] = right_hand_side_parallel(W,Vx,Vy,Vx_new,Vy_new,Vrx,Vry,dv,epsilon,C_gamma,gamma,Np)

%midpoint calculation
Vmid_x = 0.5*(Vx_new + Vx);
Vmid_y = 0.5*(Vy_new + Vy);
Vdif_x = Vx_new - Vx;
Vdif_y = Vy_new - Vy;

%Gauss quadrature nodes and weights
[gv,gw] = lgwt(4,0,1);
Ng = length(gv);

inside = zeros(Np,Ng);
gF_x = zeros(Np,1); gF_y = zeros(Np,1);

parfor i = 1:Np
    for j = 1:Ng
        inside(i,j) = ...
            sum(W.*psi_2d(Vrx(i)-(Vx+gv(j)*Vdif_x),Vry(i)-(Vy+gv(j)*(Vdif_y)),epsilon));
    end
end
inside = log(inside);


term1_x = zeros(Np,Ng);
term1_y = zeros(Np,Ng);

parfor i = 1:Np
    for j = 1:Ng
        [A,B] = gpsi_2d((Vx(i)+gv(j)*Vdif_x(i))-Vrx,(Vy(i)+gv(j)*Vdif_y(i))-Vry,epsilon);
        term1_x(i,j) = dv^2*sum(A.*inside(:,j));
        term1_y(i,j) = dv^2*sum(B.*inside(:,j));
    end
end

parfor i = 1:Np
    gF_x(i) = sum(gw'.*term1_x(i,:));
    gF_y(i) = sum(gw'.*term1_y(i,:));
end


U_x = zeros(Np,1);
U_y = zeros(Np,1);
% parfor for i
D_term = zeros(Np,1);
parfor i = 1:Np
    % A is a 2x2 symmertric matrix
    len = sqrt((Vmid_x(i)-Vmid_x).^2 +(Vmid_y(i)-Vmid_y).^2);
    
    A11 = C_gamma*len.^gamma.*(len.^2-(Vmid_x(i)-Vmid_x).^2);
    A12 = -C_gamma*len.^gamma.*(Vmid_x(i)-Vmid_x).*(Vmid_y(i)-Vmid_y);
    
    A21 = A12;
    A22 = C_gamma*len.^gamma.*(len.^2-(Vmid_y(i)-Vmid_y).^2);

    termA = gF_x(i)-gF_x;
    termB = gF_y(i)-gF_y;
    
    U_x(i) = -sum(W.*(A11.*termA+A12.*termB));
    U_y(i) = -sum(W.*(A21.*termA+A22.*termB));

    D_term(i) = sum(W.*(termA.*(A11.*termA + A12.*termB)+...
    termB.*(A21.*termA + A22.*termB))) 
end

dissipation = 0.5*sum(W.*D_term);
