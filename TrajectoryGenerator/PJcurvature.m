function [kappa,psi] = PJcurvature(x,y)
% ���������������ö���������Ϻ��������ʹ�ʽ������
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
    psi = atan2(b(2),a(2));% �޶���-pi~pi��
%     psi = mod(atan2(b(2),a(2)),2*pi);% �޶���0~2pi��
end
