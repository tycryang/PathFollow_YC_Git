function [kappa,psi] = PJcurvature(x,y)
% 利用相邻三个点用二次曲线拟合后，利用曲率公式求曲率
    x = reshape(x,3,1);
    y = reshape(y,3,1);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x;
    b = M\y;

    kappa = 2.*(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    psi = atan2(b(2),a(2));% 限定在-pi~pi中
%     psi = mod(atan2(b(2),a(2)),2*pi);% 限定在0~2pi中
end
