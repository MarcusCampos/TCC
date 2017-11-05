path = 'C:\Users\marcu\Documents\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Automação';
DIRs = dir(path);
ListagemdeDiretorios = struct2cell(DIRs);

%[filename,pathname] = uigetfile('*.jpg','open file');
%fname = fullfile(pathname,filename);
%file_name = table(f_name);
%file = dir(fullfile(pathname,'*.jpg'));
arquivoSemSegmentacao = dir(strcat(path,'\',ListagemdeDiretorios{1,3}));
arquivoSegmentacao = dir(strcat(path,'\',ListagemdeDiretorios{1,4}));
ListagemdeDiretoriosSEMSEG = struct2cell(arquivoSemSegmentacao);%
ListagemdeDiretoriosComSEG = struct2cell(arquivoSegmentacao);%uso
ResultadosSeg = strcat(path,'\Resultados'); % uso
Seed = strcat(path,'\Seed');


for fileSize = 3:size(ListagemdeDiretoriosSEMSEG)
    DiretorioAtualSemSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}));
    DiretorioAtualComSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,4},'\',ListagemdeDiretoriosComSEG{1,fileSize}));
   
    DiretorioAtualSemSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}));
    ListagemdeSemSeg = struct2cell(DiretorioAtualSemSeg);
    ListagemdeComSeg = struct2cell(DiretorioAtualComSeg);
    %Trata Ordenação
    % Convert to a matrix
    ListagemdeSemSegfields = fieldnames(DiretorioAtualSemSeg);
    ListagemdeComSegfields = fieldnames(DiretorioAtualComSeg);
    sz = size(ListagemdeSemSeg)  
    sz1 = size(ListagemdeComSeg)  
    ListagemdeSemSeg = reshape(ListagemdeSemSeg, sz(1), []);      % Px(MxN)
    ListagemdeComSeg = reshape(ListagemdeComSeg, sz1(1), []);

    ListagemdeSemSeg = ListagemdeSemSeg';
    ListagemdeComSeg = ListagemdeComSeg'; % (MxN)xP
    ListagemdeSemSeg = sortrows(ListagemdeSemSeg, 1);
    ListagemdeComSeg = sortrows(ListagemdeComSeg, 1);
    
    ListagemdeSemSeg = reshape(ListagemdeSemSeg', sz);
    ListagemdeComSeg = reshape(ListagemdeComSeg', sz1);
    ListagemdeSemSeg = cell2struct(ListagemdeSemSeg, ListagemdeSemSegfields, 1);
    ListagemdeComSeg = cell2struct(ListagemdeComSeg, ListagemdeComSegfields, 1);

    
    ResultadosNSegString = strcat(ResultadosSeg,'\',ListagemdeDiretoriosSEMSEG{1,fileSize},'\',ListagemdeDiretorios{1,3});
    ResultadosSegString = strcat(ResultadosSeg,'\', ListagemdeDiretoriosComSEG{1,fileSize},'\',ListagemdeDiretorios{1,4});
       
    NF = length(ListagemdeSemSeg)-2;
    images = cell(NF,1);
    [~,index] = sortrows({ListagemdeComSeg.name}.'); ListagemdeComSeg = ListagemdeComSeg(index(end:-1:1)); clear index
    [~,index] = sortrows({ListagemdeComSeg.name}.'); ListagemdeComSeg = ListagemdeComSeg(index); clear index

    for k = 2 : NF
        [P1,P2,P3] = getPontosIniciaisBySelecao((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}, '\',ListagemdeSemSeg(k+1,1).name)));
        Imagem = strsplit(ListagemdeSemSeg(k+1,1).name, '.');
        NomeArquivoResultados = strcat(ResultadosSeg ,'\', 'Pontos Iniciais' , '\', ListagemdeDiretoriosSEMSEG{1,fileSize});
        OkFolderResultadosNSeg = mkdir(NomeArquivoResultados);
        
        ArquivoSaida = strcat(NomeArquivoResultados,'\', 'Pontos' , '.txt'); 
        
        fid = fopen(ArquivoSaida, 'at' );
         fprintf( fid, ' Pontos de %s', ListagemdeSemSeg(k+1,1).name );
        fprintf( fid, ' Pontos P1: [%f,%f]\n', P1(1), P1(2));
        fprintf( fid, ' Pontos P2: [%f,%f]\n', P2(1), P2(2));
        fprintf( fid, ' Pontos P3: [%f,%f]\n', P3(1), P3(2));
        fprintf( fid, '-----------------------------');
        fclose(fid);
    end

end 