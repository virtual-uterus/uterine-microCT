%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot entropy measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heart 1 Original
InputPath{1} = '../DataHeart1/STBinary/';
Prefix{1} = 'Heart1Original';

% Heart 2 Original 
InputPath{2} = '../Data/STBinary/';
Prefix{2} = 'Heart2Original';

% Heart 2 Air Dry
InputPath{3} = '../DataHeartDry/STBinary/';
Prefix{3} = 'Heart2AirDry';

% Other data
Levels = [2:5];

% Load in entropy data
EntropyData = zeros(5,1,3,4);
for h=1:3
    for l=1:4
       S=load(sprintf('%s%s_ShannonLogEnergyEntropyData_L%1d_v2.mat',InputPath{h},Prefix{h},Levels(l)));
       EntropyData(:,:,h,l) = S.InfoContent;
    end
end

figure(1);
bar3(squeeze(EntropyData(3,:,:,:))-min(min(squeeze(EntropyData(3,:,:,:)))));
xlabel('Smoothing level','FontName','Arial','FontSize',10);
set(gca,'XTickLabel',{'L2','L3','L4','L5'},'FontName','Arial','FontSize',10);
ylabel('Heart','FontName','Arial','FontSize',10);
set(gca,'YTickLabel',{'H1 Orig','H2 Orig','H2 Dry'},'FontName','Arial','FontSize',10);
zlabel('\Delta Shannon Entropy (bits)','FontName','Arial','FontSize',10);
box on;
camproj perspective;
view(59.1193,33.6);
print(gcf,'../Reports/ChangeShannonEntropyHeartsLevelsHelixAngle.png','-dpng','-r1200');

figure(2);
bar3(squeeze(EntropyData(4,:,:,:))-min(min(squeeze(EntropyData(4,:,:,:)))));
xlabel('Smoothing level','FontName','Arial','FontSize',10);
set(gca,'XTickLabel',{'L2','L3','L4','L5'},'FontName','Arial','FontSize',10);
ylabel('Heart','FontName','Arial','FontSize',10);
set(gca,'YTickLabel',{'H1 Orig','H2 Orig','H2 Dry'},'FontName','Arial','FontSize',10);
zlabel('\Delta Shannon Entropy (bits)','FontName','Arial','FontSize',10);
box on;
camproj perspective;
view(59.1193,33.6);
print(gcf,'../Reports/ChangeShannonEntropyHeartsLevelsFA.png','-dpng','-r1200');

return;

% For each data type (row) plot column 2 (log energy vs Shannon entropy)
figure(1); clf;
CS = {'rs','ro','rd';'gs','go','gd';'bs','bo','bd';'cs','co','cd'};
for t=1:5
    subplot(1,5,t); hold on;
    for h=1:3
        for l=1:4
          scatter(EntropyData(t,1,h,l),EntropyData(t,2,h,l),50,CS{l,h},'filled'); 
          box on; grid on;
          set(gca,'xscale','log','yscale','log');
          xlabel('Shannon entropy','FontName','Arial','FontSize',10); 
          ylabel('Log energy entropy','FontName','Arial','FontSize',10); 
          title(S.InfoLabels{t},'FontName','Arial','FontSize',9,'FontWeight','normal');
          set(gca,'FontName','Arial','FontSize',10);
        end
    end
    hold off;
end


