function [ P1,P2,P3 ] = getPontosIniciaisBySelecao( stringFilepath )
    I = imread(stringFilepath);
    G = I;

    [row,col,t] = size(G);
    if t>1
    G = G(:,:,1);
    end

    %---preprocessing (filtering,dilating,erosion,etc)--------
    G = double(G);             
    I = 1 - G/255;          
    Ima = Gaus_filter(I,.5);

    %-----------edge map------------
    [Ima_x, Ima_y] = gradient(Ima);
    f = sqrt(Ima_x.^2 + Ima_y.^2);
    %figure(p); 
    imdisp(G);
    title(stringFilepath);
     
    %------- Compute the external force filed  
    [u,v] = MyGVC(f, 30,30,.8,4.6); 
    %Nomalizing the GVF external force
    mag = sqrt(u.*u+v.*v);
    px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error

    pontos = cell(3,2);

    [x1, y1] = ginput(1);
    [x2,y2]  = ginput(1);
    [x3, y3]  = ginput(1);
%     P1 = [pontos(2,1), pontos(2,2)];
%     P2 = [pontos(3,1),pontos(3,2)]  ;
%    P3 = [pontos(1,1), pontos(1,2)];
   
     P1 = [x1, y1];
    P2 = [x2,y2]  ;
   P3 = [x3, y3];
end

