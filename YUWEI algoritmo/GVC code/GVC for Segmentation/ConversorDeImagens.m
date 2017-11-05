[filename,pathname] = uigetfile('*.jpg','open file');
fname = fullfile(pathname,filename)
    G = imread(fname);
    G = imcrop(G, [714.5 0.5 452 456]);
   imshow(G)
    imwrite(G, fname);
