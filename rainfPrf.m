function f=rainfPrf(NCic)
%if no cicle limit is specified, it is set to 1e7
if exist('NCic')==0
    NCic=1e7;
end
%directory where the files are
dirPath="D:\\Google Drive\\Faculdade\\TCC\\Parafusos Transiente\\frmtstr\\";
%file list with all dirs and files at the origin directory
list=dir(dirPath);
%File size : Number of variables extracted x node number (inf) in order to keep adding data until it is full
sizeF = [5,inf];
%initialize variable to calculate damage as zeros
Dano=zeros(size(list,1),9);
%loop to calculate normal stress then shear stress
for k=1:2
    %loop to calculate the damage in the time history data of each screw
    for j=1:size(list,1)
        %if it is not a file, skip
        if list(j).isdir==1
            continue;
        end
        
        vk=strfind(list(j).name,'rstFatigueAnalysis');
        if isempty(vk)~=1
        if list(j).isdir==1||(vk>0)
            continue;
        end
        end
        %open the file, read and fill the A matrix with the data then close the fiile
        arquivo=fopen(dirPath+list(j).name,'r');
        A=fscanf(arquivo,'%f',sizeF);
        A=A';
        fclose(arquivo);
        Name=erase(list(j).name,'.txt');
        Dano(j,1)=str2double(Name);
        %call sig2ext to find inflexion points in the time history data
        extB=sig2ext(A(:,2*k));
        extT=sig2ext(A(:,2*k+1));
        %call rainflow, returns a 3xn matrix, where :
        %1: Stress Amplitude
        %2: Mean Stress
        %3: Number of Cycles Counted
        rfB=rainflow(extB);
        rfT=rainflow(extT);
        ampB=rfB(1,:);
        ampT=rfT(1,:);
        meanB=rfB(2,:);
        meanT=rfT(2,:);
        nB=rfB(3,:);
        nT=rfT(3,:);
        %Start the miner Sum at 0
        minerSumB=0;
        minerSumT=0;
        %sums from 1 to n cycles, divided by the cycle limit
        for i=1:size(nT)
            if meanT(i)+ampT(i)>0
                minerSumT=minerSumT+nT(i)/FLife(ampT(i)*2,50*k);
            else
            end
        end
        for i=1:size(nB)
            if meanB(i)+ampB(i)>0
                minerSumB=minerSumB+nB(i)/FLife(ampB(i)*2,50*k);
            else
            end
        end
        Dano(j,2*k)=minerSumB;
        Dano(j,2*k+1)=minerSumT;
    end
    %loop to calculate the damage under normal operation
    for j=1:size(list,1)
        %if it is not a file, skip
        if list(j).isdir==1
            continue;
        end
        vk=strfind(list(j).name,'rstFatigueAnalysis');
        if isempty(vk)~=1
        if list(j).isdir==1||(vk>0)
            continue;
        end
        end
        %open the file, read and fill the A matrix with the data then close the fiile
        arquivo=fopen(dirPath+list(j).name,'r');
        A=fscanf(arquivo,'%f',sizeF);
        A=A';
        %creates an auxiliary vector to check if time is >= 1
        StartTime=4;
        vet=(A(:,1)>=StartTime);
        A(:,k+1)=vet.*A(:,k+1);
        fclose(arquivo);
        extB=sig2ext(A(:,2*k));
        extT=sig2ext(A(:,2*k+1));
        %call rainflow, returns a 3xn matrix, where :
        %1: Stress Amplitude
        %2: Mean Stress
        %3: Number of Cycles Counted
        rfB=rainflow(extB);
        rfT=rainflow(extT);
        ampB=rfB(1,:);
        ampT=rfT(1,:);
        meanB=rfB(2,:);
        meanT=rfT(2,:);
        nB=rfB(3,:);
        nT=rfT(3,:);
        %Start the miner Sum at 0
        minerSumB=0;
        minerSumT=0;
        for i=1:size(nT)
            %Only take into account if it is under tensile stress
            if meanT(i)+ampT(i)>0
                % Only sums damage if the stress amplitude is over the stress limit
                if ampT(i)>0.1*50*k
                    minerSumT=minerSumT+nT(i)*(NCic)/4.166667/FLife(ampT(i)*2,50*k);
                end
            else
            end
        end
        for i=1:size(nB)
            if meanB(i)+ampB(i)>0
                % Only sums damage if the stress amplitude is over the stress limit
                if ampB(i)>0.1*50*k
                    minerSumB=minerSumB+nB(i)*(NCic)/4.166667/FLife(ampB(i)*2,50*k);
                end
            else
            end
        end
        Dano(j,2*k+4)=minerSumB;
        Dano(j,2*k+5)=minerSumT;
    end
end
f=max(Dano(:,2:9));
fileID=fopen([dirPath "rstFatigueAnalysis.txt",'w']);
Con=formatDummy("D:\\Google Drive\\Faculdade\\TCC\\teste calculo paraf\\");
[a,b]=ismember(Dano(:,1),Con(:,1));
fprintf(fileID,"Elem. DanoTop/Acion DanoBot/Acion CisDanoTop/Acion CisDanoBot/Acion DanoTop DanoBot CisDanoTop CisDanoBot\n");
for i=1:size(Dano,1)
    if Dano(i,1)==0
        continue
    end
fprintf(fileID,"%f %e %e %e %e %e %e %e %e %f\n",string(Dano(i,1)),Dano(i,2),Dano(i,3),Dano(i,4),Dano(i,5),Dano(i,6),Dano(i,7),Dano(i,8),Dano(i,9),string(Con(b(i),2)));
end
fclose(fileID);
end
