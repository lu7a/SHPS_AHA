function [BestX,BestF,HisBestFit,VisitTable]=SHPSAHA(MaxIt,nPop,Low,Up,Dim,fobj)


    PopPos=zeros(nPop,Dim);
    PopFit=zeros(1,nPop);
    PopPos= Sobol(Dim,nPop,Low,Up);
    for i=1:nPop
        PopFit(i)=fobj(PopPos(i,:));
    end

    BestF=inf;
    BestX=[]; 
    HisBestFit=zeros(1,MaxIt);
    for i=1:nPop
        if PopFit(i)<=BestF
            BestF=PopFit(i);
            BestX=PopPos(i,:);
        end
        HisBestFit(i)=BestF;
    end
    %% Initialize visit table
    VisitTable=zeros(nPop) ;
    VisitTable(logical(eye(nPop)))=NaN;    

    It = 0;
    SF=[];
    a = 0.5;
    MaxIt = MaxIt-nPop;
    while It<MaxIt
        DirectVector=zeros(nPop,Dim);% Direction vector/matrix

        for i=1:nPop

            r=rand;
            if r<1/3     % Diagonal flight 
                RandDim=randperm(Dim);
                if Dim>=3

                    RandNum=ceil(rand*(Dim-2)+1);
                else
                    RandNum=ceil(rand*(Dim-1)+1);
                end
 
                DirectVector(i,RandDim(1:RandNum))=1;
            else
                if r>2/3  % Omnidirectional flight 
                    DirectVector(i,:)=1;
                else  % Axial flight 
                    RandNum=ceil(rand*Dim);
                    DirectVector(i,RandNum)=1;
                end
            end
 %%
            if rand<0.5   % Guided foraging 
            
                [MaxUnvisitedTime,TargetFoodIndex]=max(VisitTable(i,:));
     
                MUT_Index=find(VisitTable(i,:)==MaxUnvisitedTime);
                
                if length(MUT_Index)>1
                    [~,Ind]= min(PopFit(MUT_Index));
                    TargetFoodIndex=MUT_Index(Ind);
                end
    
                newPopPos=PopPos(TargetFoodIndex,:)+randn*DirectVector(i,:).* ...
                    (PopPos(i,:)-PopPos(TargetFoodIndex,:));
              
                newPopPos=SpaceBound(newPopPos,Up,Low);
%            
                newPopFit=fobj(newPopPos);
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                end
                if PopFit(i)<=BestF
                    BestF=PopFit(i);
                    BestX=PopPos(i,:);
                end
                It=It+1;
                if It>MaxIt
                    break;
                end
                HisBestFit(It+nPop)=BestF;
            else
             
                r1 = 2 * randn() - 1;
                newPopPos= PopPos(i,:)+r1*DirectVector(i,:).*PopPos(i,:);
                newPopPos=SpaceBound(newPopPos,Up,Low);
                newPopFit=fobj(newPopPos);
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                end
                if PopFit(i)<=BestF
                    BestF=PopFit(i);
                    BestX=PopPos(i,:);
                end
                It=It+1;
                if It>MaxIt
                    break;
                end
                HisBestFit(It+nPop)=BestF;
            end
        end
        Fa(i) = a*randn();
        if mod(It,2*nPop)==0 % Migration foraging 
            


            [~, MigrationIndex]=max(PopFit);
            newPopPos =PopPos(MigrationIndex,:)+Fa(i)*DirectVector(MigrationIndex,:).*(BestX-PopPos(MigrationIndex,:));
            newPopFit=fobj(newPopPos);
            if newPopFit<PopFit(MigrationIndex)
                SF=[SF;Fa(i)];
                PopFit(MigrationIndex)=newPopFit;
                PopPos(MigrationIndex,:)=newPopPos;
                VisitTable(MigrationIndex,:)=VisitTable(MigrationIndex,:)+1;
                VisitTable(:,MigrationIndex)=max(VisitTable,[],2)+1;
                VisitTable(MigrationIndex,MigrationIndex)=NaN; 
            else
                VisitTable(MigrationIndex,:)=VisitTable(MigrationIndex,:)+1;
            end
           
            if PopFit(i)<=BestF
                    BestF=PopFit(i);
                    BestX=PopPos(i,:);
            end
            It=It+1;
            if It>MaxIt
                    break;
            end
            HisBestFit(It+nPop)=BestF;
            c = 0.2;
            LA = mean(SF);
            a = (1-c)*a+c*LA;
        end

    end


