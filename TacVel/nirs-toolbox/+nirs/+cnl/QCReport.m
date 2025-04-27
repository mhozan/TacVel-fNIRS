function QCReport(data)

% figur(mfilename);clf;
figure;

n=length(data);
axpos = subplotpos(5,1,'xgap',0.01,'ygap',0.01);

for i=1:n
    figur([data(i).demographics.values{48},'_',data(i).demographics.values{49}])
    ax1=axes();
    ax1.Position = axpos(1,:);
    data(i).draw;
    legend off;
    ttle = [data(i).demographics.values{end-1},' ',data(i).demographics.values{end}];
    title(ax1,ttle,'Interpreter','none')

    sni=nirs.math.structnoiseindex(data(i).data);
    types=unique(data(i).probe.link.type);
    lst=ismember(data(i).probe.link.type,types(1));
    [colors,lineStyles]=nirs.util.values2lines(sni(lst),[-10 10]);
    ax2=axes();
    ax2.Position = axpos(2,:);
    data(i).probe.draw(colors,lineStyles,ax2);
    text(0,.38,num2str(types(1)),'FontSize',20)



    lst=ismember(data(i).probe.link.type,types(2));
    [colors,lineStyles]=nirs.util.values2lines(sni(lst),[-10 10]);

    ax3=axes();
    ax3.Position = axpos(3,:);
    data(i).probe.draw(colors,lineStyles,ax3);
    text(0,.38,num2str(types(2)),'FontSize',20)

    ax4=axes();
    ax4.Position = axpos(4,:);
    ax4.Box = 'off';
    hist(sni);
    ax5=axes();
    ax5.Position = axpos(5,:);
    ax5.NextPlot ='add';
    for aa=1:data(i).auxillary.count
        aux = data(i).auxillary.values{aa}.data;
        thresh = prctile(aux,99);
        aux(aux>thresh)=thresh;
        aux = 1*aux/thresh;
        thresh = prctile(aux,1);
        aux(aux<thresh)=thresh;
        % aux = 1*aux/thresh;

        lh = plot(ax5,data(i).auxillary.values{aa}.time,...
            2*(aa-1)+aux);
        lh.DisplayName = data(i).auxillary.keys{aa};

    end
    legend('Interpreter','none')
    linkaxes([ax1,ax5],'x')

end
AdjustFontSize('maxfontsize',7,'minfontsize',2)
