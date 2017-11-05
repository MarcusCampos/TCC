
%--read the input image
clear;
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
OkFolder = mkdir(ResultadosSeg)
OkFolder = mkdir(Seed)
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P21', 'P21');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P22', 'P22');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P23', 'P23');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P24', 'P24');
catch
    lastwarn('');
end

try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P25', 'P25');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P26', 'P26');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P27', 'P27');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P28', 'P28');
catch
    lastwarn('');
end
try
  Run_Segmentation_byPasta(Seed, ResultadosSeg, ListagemdeDiretorios,path, 'P29', 'P29');
catch
    lastwarn('');
end
