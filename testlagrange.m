n = 10;

% x 값 범위 설정
x_range = linspace(0, 2*pi, n);

% 주어진 함수 위의 점 생성
x = x_range;
y = sin(2*x) + cos(x);

% 보간할 x 값 범위 설정
x_interp = linspace(0, 2*pi, 100);

% 1차원 Lagrange 보간 함수
lagrange_interp = zeros(size(x_interp));
for i = 1:length(x_interp)
    for j = 1:length(x)
        % 각 점에서의 Lagrange 다항식 계산
        L = 1;
        for k = 1:length(x)
            if k ~= j
                L = L * (x_interp(i) - x(k)) / (x(j) - x(k));
            end
        end
        % 각 점에서의 보간값 계산
        lagrange_interp(i) = lagrange_interp(i) + y(j) * L;
    end
end

% 결과 그래프 출력
plot(x, y, 'o', x_interp, sin(2*x_interp) + cos(x_interp), '-', x_interp, lagrange_interp, '--');
legend('주어진 함수 위의 점', '원본 함수', '1차원 Lagrange 보간');
xlabel('x');
ylabel('y');
title('1차원 Lagrange 보간');
hold off
plot(x_interp, lagrange_interp, 'x')