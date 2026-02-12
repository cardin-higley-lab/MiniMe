function scatterPlot(vector, minVal, maxVal)


x = 1:size(vector,1);
coefficients = polyfit(x, vector(:,1), 1);
xFit = linspace(min(x), max(x), 1000);
P = polyfit(x, vector(:,1), 1);
yFit = polyval(coefficients , xFit);
mdl = fitlm(x,vector(:,1));
r2 = [' R2 = ', num2str(mdl.Rsquared.Adjusted)];
pval = [' p = ', num2str(mdl.ModelFitVsNullModel.Pvalue)];
eqn = string(" Linear: y = " + P(1)) + "x + " + string(P(2));

hold on
scatter(x, vector(:,1), 'ko')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
ylim([minVal maxVal])
text(min(x),maxVal,eqn,"HorizontalAlignment","left","VerticalAlignment","top")
text(min(x),maxVal*0.9,r2,"HorizontalAlignment","left","VerticalAlignment","top")
text(min(x),maxVal*0.8,pval,"HorizontalAlignment","left","VerticalAlignment","top")



end 