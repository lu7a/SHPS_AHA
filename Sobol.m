%% Sobol平稳序列
function best_x= Sobol(D,pop,Lb,Ub)
p = sobolset(D);
R=[];
for i=1:pop
    r=p(i,:);
    r=Lb+r.*(Ub-Lb);
    R=[R; r];
end
best_x=R;
end