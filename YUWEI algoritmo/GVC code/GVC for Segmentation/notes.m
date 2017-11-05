    %  if k == 1
    %      endoxFirst = tempx;
    %      endoyFirst = tempy; 
    %      epixFirst = x; 
    %      epiyFirst =y;
    %      centerXFirst = centerX;
    %      centerYFirst = centerY;
    %  end 
    %  if k == 2
    %     figureNew = figure; 
    %     hold on, 
    %     %imshow(figureNew);
    %     %plot([x;x(1,1)],[y;y(1,1)],'g','LineWidth',1);
    %   %  plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
    %     plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'k','LineWidth',1);
    %     %plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'b','LineWidth',1);
    %     %title(['Final Segmentation Result - Selection Compare,  iter = ' num2str(i*5)]);
    %     axis off
    %     hold off;
    %     strgeral = 'Salvame';
    %     stringnumer = num2str(k);
    %     fileNameFormat = '.jpg';
    %     stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
    %     fullname = fullfile(pathname, stringNameFile);
    %     saveas(figureNew, fullfile(pathname, stringNameFile));
    %     figureNew2 = figure; 
    % 
    %     hold on, 
    %     plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'k','LineWidth',1);
    %      axis off 
    %     hold off;
    %     strgeral = 'Salvame2';
    %     stringnumer = num2str(k);
    %     fileNameFormat = '.jpg';
    %     stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
    %     fullname2 = fullfile(pathname, stringNameFile);
    %     saveas(figureNew2, fullfile(pathname, stringNameFile));
    %      N = verifyIntersection(fullname,fullname2);
    %     hold off;
    %  end