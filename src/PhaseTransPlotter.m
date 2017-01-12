% PhaseTransPlotter

function [] = PhaseTransPlotter(bc,sigVec,rhoVec,x,bc_iso, bc_nem,f_iso,f_nem,...
    f_rec,nFrames,FileNameMov,FileNameDist, FileNameRhoSig,FPS,GotPhaseTrans)
% Plot rho,sig vs concentration
h1 = figure;

subplot(2,1,1)
plot(bc,sigVec,bc,rhoVec)
title('\rho and \sigma vs concentration')

legend('\sigma','\rho')
xlabel('scaled concentration bc')
ylabel('\rho and \sigma')
axis([min(bc) max(bc) 0  (max(max(rhoVec,sigVec)) + max(max(rhoVec,sigVec)) / 10)] )

subplot(2,1,2)
plot(bc, log(bc), bc,sigVec,bc,bc.*rhoVec,bc, (log(bc) + sigVec + bc.*rhoVec) )
title('Terms in FE vs concentration')

legend('ln(bc)','\sigma','bc*\rho','sum of terms')
xlabel('scaled concentration bc')
ylabel('\rho and \sigma')

saveas(h1,FileNameRhoSig,'jpg')
% Plot the distribution at the phase transition

if GotPhaseTrans
    h2 = figure;
    
    subplot(2,1,1)
    plot(x,f_iso)
    xlabel('x')
    ylabel('f(x)')
    TitStr = sprintf('Equilbrium distribution at dimenionless concentration c'' = %f',bc_iso);
    title(TitStr)
    
    subplot(2,1,2)
    plot(x,f_nem)
    xlabel('x')
    ylabel('f(x)')
    TitStr = sprintf('Equilbrium distribution at dimenionless concentration c'' = %f',bc_nem);
    title(TitStr)
    saveas(h2,FileNameDist,'jpg')
end
% Make a movie

% Use VideoWriter
vidObj = VideoWriter(FileNameMov);    %Create a structure for the AVI file with this filename
vidObj.FrameRate = FPS;
open(vidObj);                       %Open it so I can write to it

% Set up the figure
figure

xlabel('x');ylabel('f');title('distribution changing with concentration')
set(gca,'XLim',[min(x) max(x)],'YLim',[min(min(f_rec)) max(max(f_rec))]);  %Set axis
set(gca,'NextPlot','replaceChildren'); %This makes it so only data, not axis, changes
% Sometimes (for 2D and 3D plots) you need this
set(gcf,'renderer','zbuffer')    %This is what fixed the issue for me. Fixes problem with getframe

% keyboard
% Textbox stuff
mTextBox = uicontrol('style','text');
set(mTextBox,'Position', [250 250 80 20])

for ii = 1:nFrames
    TBStr = sprintf( 'c'' = %f ',bc(ii) );
    set(mTextBox,'String',TBStr)
    plot(x, f_rec(ii,:));
    frametemp   = getframe(gcf);    %Current frame
    writeVideo(vidObj,frametemp) %Write to open structure
end

close(vidObj) %We are done writing the structure. File already saved

% Move it to your results directory
movefile('*.avi', 'C:\Users\MwS\Documents\Research\BG\Figure\INtrans\')
movefile('*.jpg', 'C:\Users\MwS\Documents\Research\BG\Figure\INtrans\')