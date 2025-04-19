
clear all
clc
tic;


MaxIt=15000;
nPop=30;   
Function_name='F3'; 

[lb,ub,dim,fobj] = Get_Functions_details(Function_name);
[Best_pos,Best_score,HisBestFit,VisitTable]=SHPSAHA(MaxIt,nPop,lb,ub,dim,fobj);

          
