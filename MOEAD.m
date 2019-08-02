function MOEAD(Problem,M)
clc;format compact;tic;


%参数设定
Generations = 700;
delta = 0.9;
nr = 2;
if M == 2
    N = 100;
    H = 99;
else M == 3
    N = 105;
    H = 13;
end

    %初始化向量
    Evaluations = Generations*N;
    [N,W] = EqualWeight(H,M);
    W(W==0) = 0.000001;
    T = floor(N/10);
    Generations = floor(Evaluations/N);

    %邻居判断
    B = zeros(N);
    for i = 1 : N-1
        for j = i+1 : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %初始化种群
    [Population,Boundary] = Objective(0,Problem,M,N);
    FunctionValue = Objective(1,Problem,M,Population);
    Z = min(FunctionValue);

    %开始迭代
    for Gene = 1 : Generations
        %对每个个体执行操作
        for i = 1 : N
            %选出父母
            if rand < delta
                P = B(i,:);
            else
                P = 1:N;
            end
            k = randperm(length(P));
            
            %产生子代
            Offspring = Gen(Population(i,:),Population(P(k(1)),:),Population(P(k(2)),:),Boundary);
            OffFunValue = Objective(1,Problem,M,Offspring);

            %更新最优理想点
            Z = min(Z,OffFunValue);
            
            %更新P中的个体
            c = 0;
            for j = randperm(length(P))
                if c >= nr
                    break;
                end
                g_old = max(abs(FunctionValue(P(j),:)-Z).*W(P(j),:));
                g_new = max(abs(OffFunValue-Z).*W(P(j),:));              
                if g_new < g_old
                    %更新当前向量的个体
                    Population(P(j),:) = Offspring;
                    FunctionValue(P(j),:) = OffFunValue;
                    c = c+1;
                end
            end

        end
        cla;
        DrawGraph(FunctionValue);
        hold on;
        switch Problem
            case 'DTLZ1'
                if M == 2
                    pareto_x = linspace(0,0.5);
                    pareto_y = 0.5 - pareto_x;
                    plot(pareto_x, pareto_y, 'r');
                elseif M == 3
                    [pareto_x,pareto_y]  = meshgrid(linspace(0,0.5));
                    pareto_z = 0.5 - pareto_x - pareto_y;
                    axis([0,1,0,1,0,1]);
                    mesh(pareto_x, pareto_y, pareto_z);
                end
            otherwise
                if M == 2
                    pareto_x = linspace(0,1);
                    pareto_y = sqrt(1-pareto_x.^2);
                    plot(pareto_x, pareto_y, 'r');
                elseif M == 3
                    [pareto_x,pareto_y,pareto_z] =sphere(50);
                    axis([0,1,0,1,0,1]);
                    mesh(1*pareto_x,1*pareto_y,1*pareto_z);
                end
        end
        pause(0.01);
        %clc;
    end
end