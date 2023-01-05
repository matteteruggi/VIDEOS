load rates.mat
load date.mat
data=data(2:end);
% store concrentration rates and standard deviations
CRATES=[];
STDS=[];
% store level of spread
Spread10minus2=round((rates(2:end,8)-rates(2:end,4))*100,2);
rates=diff(rates);
v = VideoWriter('volatilities.avi');
h = VideoWriter('loadfactors.avi');
r = VideoWriter('correlations.avi');
shock = VideoWriter('1sigmashock.avi');
% attenzione sono secondi alla men uno
duratainHz=2;
v.FrameRate=duratainHz;
h.FrameRate=duratainHz;
r.FrameRate=duratainHz;
shock.FrameRate=duratainHz;
open(v);
open(h);
open(r);
open(shock);
X=[0.25 0.5 1 2 3 5 7 10 20 30];
Xtilda=linspace(1,10,10);
Volatlityelasticiy=[];
pd = makedist('Gamma','a',20/2,'b',2);
begindate=[];
for i=1:10:(length(rates)-250)
    pippo=rates(i:i+250,:);
    Spread=Spread10minus2(i);
    a=mahal(pippo,pippo);
    qqplot(a,pd);
    title(["QQ plot of Sqaured Malhanobis distance vs Chi-square with 10 DoF"," 1 year daily observations"])
    share=(sum(a>18.307)/length(a));
    subtitle(["bgn data: "+string(data(i))+" 10 minus 2 spread the first day= "+string(Spread)," Multivariate normal assumption violation: "+string(share)])
    saveas(gcf,"QQ"+string(i)+".png")
    R=corr(pippo);
    heatmap(X,X,R);
    colormap(turbo);
    clim([0 1])
    title("Correlation matrix of 1 year daily changes, bgn data: "+string(data(i)))
    frame = getframe(gcf);
    writeVideo(r,frame);
    V=cov(pippo);
    [coeffs,latents,explaineds] = pcacov(V);
    % check by V-coeff*diag(latents)*coeff'
    volatilities=diag(V).^0.5;
    %caution: is already multiplied by desired shock
    pcsigma=3*latents.^0.5;
    STDS=[STDS;volatilities'];
    CRATES=[CRATES; cumsum(explaineds(1:3))']
    Volelast=polyfit(X,volatilities,1);
    Volatlityelasticiy=[Volatlityelasticiy;Volelast];
    scatter(Xtilda,volatilities)
    xlim([1 10])
    xticklabels(X)
    title("volatility term structure, 1 yr daily observations")
    subtitle(["bgn data: "+string(data(i))+" 10 minus 2 spread the first day= "+string(Spread),"elasticity= "+string(Volelast(1))])
    xlabel('time to maturity')
    ylabel('volatility')
    saveas(gcf,"VOLATILITY"+string(i)+".png")
    frame = getframe(gcf);
    writeVideo(v,frame);
    hold off
    plot(Xtilda,coeffs(:,1:3));
    legend('first component','second component','third component')
    title('PCA Load factors - first three components, 1 yr daily observations')
    subtitle("bgn data:"+string(data(i))+" 10 minus 2 spread the first day= "+string(Spread))
    xlabel('time to maturity')
    xticklabels(X)
    xlim([1 10])
    ylim([-0.6 1])
    ylabel('load component')
    saveas(gcf,"LOADF"+string(i)+".png")
    frame = getframe(gcf);
    writeVideo(h,frame);
    sigmashock=coeffs(:,1:3).*repmat(pcsigma(1:3)',10,1);
    sigmashock=100*sigmashock;
    plot(Xtilda,sigmashock);
    legend('first component','second component','third component')
    lgd.FontSize=6;
    title('Shift size - first three components - 3sigma shock, 1 yr daily observations')
    subtitle("bgn data:"+string(data(i))+" 10 minus 2 spread the first day= "+string(Spread))
    xlabel('time to maturity')
    xticklabels(X)
    xlim([1 10])
    ylabel('shift size - bps')
    if year(data(i))<1990
      ylim([-20 80])
    else
      ylim([-10 20])
    end
    saveas(gcf,"sigmashock"+string(i)+".png")
    frame = getframe(gcf);
    writeVideo(shock,frame);
    begindate=[begindate data(i)];
end
close (v)
close(h)
close(r)
close(shock)
plot(CRATES)
xlim([1,1010])
assex=round(linspace(1,length(CRATES),10),0);
legend('First component','CR2','CR3')
title('Concentration Rates')
subtitle('how much is explained by the first N components?')
xlabel('bgn date-1 year daily observations')
xticks(assex)
xticklabels(begindate(assex))
ylabel('Cumulative sums up to first, second and third component')
saveas(gcf,'CRATES.fig')
surf(STDS)
title('volatility surface')
ylabel('begin date- 1 year daily observations')
yticks(assex)
yticklabels(begindate(assex))
xlabel('term structure - time to maturity')
xticks(Xtilda)
xticklabels(X)
zlabel('volatility')
saveas(gcf,'vol.fig')
scatter(Spread10minus2(1:47),CRATES(:,1))
saveas(gcf,'CRATES1REG.fig')
scatter(Spread10minus2(1:47),CRATES(:,2))
saveas(gcf,'CRATES2REG.fig')
scatter(Spread10minus2(1:47),CRATES(:,3))
saveas(gcf,'CRATES3REG.fig')
scatter(Spread10minus2(1:47),Volatlityelasticiy(:,1))
saveas(gcf,'volest.fig')
