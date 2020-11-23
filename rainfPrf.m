function f=rainfPrf(NCic)
% NCic=4*100;
%Diret�rio onde est�o os arquivos
if exist('NCic')==0
    NCic=1e7;
end
dirPath="D:\\Google Drive\\Faculdade\\TCC\\Parafusos Transiente\\frmtstr\\";
%Lista com todos diret�rios e arquivos no diret�rio Pai
list=dir(dirPath);
%Tamanho dos arquivos : numero de vari�veis x n�mero de n�s (inf)
sizeF = [5,inf];
%Inicializa a Vari�vel Dano
Dano=zeros(size(list,1),9);
%La�o para calcular tens�o normal depois cisalhante
for k=1:2
    %La�o para calcular o dano do hist�rico de cada parafuso
    for j=1:size(list,1)
        %Se n�o for arquivo pula para o pr�ximo j
        if list(j).isdir==1
            continue;
        end
        
        vk=strfind(list(j).name,'rstFatigueAnalysis');
        if isempty(vk)~=1
        if list(j).isdir==1||(vk>0)
            continue;
        end
        end
        %abre o arquivo l� e joga as informa��es na matriz A e fecha o arq
        arquivo=fopen(dirPath+list(j).name,'r');
        A=fscanf(arquivo,'%f',sizeF);
        A=A';
        fclose(arquivo);
        Name=erase(list(j).name,'.txt');
        Dano(j,1)=str2double(Name);
        %Chama a fun��o sig2ext que extrai os pontos de inflex�o do
        %hist�rico para usar a fun��o rainflow
        extB=sig2ext(A(:,2*k));
        extT=sig2ext(A(:,2*k+1));
        %Chama a fun��o rainflow que retorna uma matriz 3xn, onde :
        %1: Amplitude
        %2: M�dia
        %3: N�mero de Ciclos no Hist�rico
        rfB=rainflow(extB);
        rfT=rainflow(extT);
        ampB=rfB(1,:);
        ampT=rfT(1,:);
        meanB=rfB(2,:);
        meanT=rfT(2,:);
        nB=rfB(3,:);
        nT=rfT(3,:);
        %Inicializa a Soma de Miner
        minerSumB=0;
        minerSumT=0;
        %Soma de um a n do n�mero de ciclos dividido pelo limite
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
    %La�o para calcular o Dano em funcionamento
    for j=1:size(list,1)
        %Se n�o for arquivo pula para o pr�ximo j
        if list(j).isdir==1
            continue;
        end
        vk=strfind(list(j).name,'rstFatigueAnalysis');
        if isempty(vk)~=1
        if list(j).isdir==1||(vk>0)
            continue;
        end
        end
        %abre o arquivo l� e joga as informa��es na matriz A e fecha o arq
        arquivo=fopen(dirPath+list(j).name,'r');
        A=fscanf(arquivo,'%f',sizeF);
        A=A';
        %Cria um vetor auxiliar que verifica se o tempo � superior a 1
        StartTime=4;
        vet=(A(:,1)>=StartTime);
        A(:,k+1)=vet.*A(:,k+1);
        fclose(arquivo);
        extB=sig2ext(A(:,2*k));
        extT=sig2ext(A(:,2*k+1));
        %Chama a fun��o rainflow que retorna uma matriz 3xn, onde :
        %1: Amplitude
        %2: M�dia
        %3: N�mero de Ciclos no Hist�rico
        rfB=rainflow(extB);
        rfT=rainflow(extT);
        ampB=rfB(1,:);
        ampT=rfT(1,:);
        meanB=rfB(2,:);
        meanT=rfT(2,:);
        nB=rfB(3,:);
        nT=rfT(3,:);
        %Inicializa a Soma de Miner
        minerSumB=0;
        minerSumT=0;
        for i=1:size(nT)
            %Somente leva em conta o dano caso seja trativo
            if meanT(i)+ampT(i)>0
                %Somente leva em conta o dano caso exceda 30% do limite
                if ampT(i)>0.1*50*k
                    minerSumT=minerSumT+nT(i)*(NCic)/4.166667/FLife(ampT(i)*2,50*k);
                end
            else
            end
        end
        for i=1:size(nB)
            if meanB(i)+ampB(i)>0
                %Somente leva em conta o dano caso exceda 30% do limite
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