% corrleation vs information for two gaussian variables
C = linspace(-1,1,200);
NC = length(C);

I = zeros(size(C));
for i=1:NC
    
    thsI = log(2*pi*exp(1)) - 1 - log(2*pi) - 0.5*log(det([ 1 C(i); C(i) 1]));
    I(i) = thsI;
end


%%
figure
plot(C,I,'b','LineWidth',3)
xlabel('Correlation','FontSize',10)
ylabel('Mutual Information (bits)','FontSize',10)
%title('Two Gaussian Variables','FontSize',16)
xlim([-1.1 1.1])

xline(-1,'m:','LineWidth',5)
xline(1,'m:','LineWidth',5)
axis square
set(gca,'FontSize',10)
set(gca,'TickLength',[0.02 0.025])
set(gca,'TickDir','out')
box off
