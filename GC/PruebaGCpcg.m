C = char('lehmer', 'moler', 'poisson');
m = size(C, 1);
display('  GC                                  |    pcg ');
display('------------------------------------------------------------------------------------------------------');
display('  tiempo        iter   ||A * x - b||  |    tiempo         iter   ||A * x - b||   |   Matriz        n');
display('------------------------------------------------------------------------------------------------------');

%display('  GC  & & &                                  pcg & & & & \\');
%display('\hline');
%display('  tiempo   &     iter  & $||A * x - b||_2$  &    tiempo    &     iter  & $||A * x - b||_2$   &   Matriz   &     n\\');
%display('\hline');

n1 = zeros(10, 1);

for i = 1:10
    [~, n1(i)] = min(abs(i * 1000 - (1:100).^2));
end

for j = 1:m
    c = deblank(C(j, :));
    for k = 1:10
        if strcmp(c, 'poisson')
            nj = n1(k);
            n = nj^2;
        else
            nj = k * 1000;
            n = nj;
        end
        x0 = zeros(n, 1);
        A = gallery(c, nj);
        b = A * ones(n, 1);
        %tic;
        [x1, iter1, t1] = gradiente_conjugado(A, b, x0, 1.0e-8, n);
        %t1 = toc;
        tic;
        [x2, ~, ~, iter2] = pcg(A, b, 1.0e-8, n);
        t2 = toc;
        fprintf('%8.4f s    %5i     %1.4e    |  %8.4f s     %5i     %1.4e     |  %7s     %5i\n',...
            t1, iter1, norm(A * x1 - b), t2, iter2, norm(A * x2 - b),...
            c, n);
        %fprintf('$%8.4f s$ &   $%5i$  &   $%1.4e$    &  $%8.4f s$  &  $ %5i $  &  $%1.4e$     &  %7s &    $%5i$\\\\\n',...
        %    t1, iter1, norm(A * x1 - b), t2, iter2, norm(A * x2 - b),...
        %    c, n);
    end
end

clear;
