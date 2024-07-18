function gF = right_hand_side(w,x,x_new,xr,dx,epsilon,n)

%difference calculation
x_dif = x_new - x;

%Gauss quadrature nodes and weights
[gx,gw] = lgwt(4,0,1);
ng = length(gx);

inside  = zeros(ng,n);
gF = zeros(1,n);

for i = 1:n
    for j = 1:ng
        inside(j,i) = ...
            sum(w.*psi_1d(xr(i)-(x+gx(j)*x_dif),epsilon));
    end
end

term1 = zeros(ng,n);
inside = log(inside);

for i = 1:n
    for j = 1:ng
        A = gpsi_1d((x(i)+gx(j)*x_dif(i))-xr,epsilon);
        term1(j,i) = dx*sum(A.*inside(j,:));
    end
end

for i = 1:n
    gF(i) = sum(gw.*term1(:,i))+.5*sum(w.*(x(i)-x+x_new(i) - x_new));
end





