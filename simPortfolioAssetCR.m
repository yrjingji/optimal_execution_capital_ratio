function [capRatioSeq] = simPortfolioAssetCR(initialPrice, permImpact,temImpact,tau,sigma,initPosition,liability,riskWeight,tradingStrats)
    xiSeq = randn(2,5);
    marketPrice =[initialPrice];
    execPrice = [[0;0]];
    remainPosition = initPosition - cumsum(tradingStrats);
    riskWeightMatrix = diag(riskWeight);
    revenueSeq = [0];
    for i =2:6
        marketPrice(:,i) = marketPrice(:,i-1) + sqrt(tau)*sigma*xiSeq(:,i-1) - permImpact*tradingStrats(:,i-1);
        execPrice(:,i) = marketPrice(:,i-1) - temImpact*tradingStrats(:,i-1);
        denominator(i-1) = marketPrice(:,i)'*riskWeightMatrix*remainPosition(:,i-1);
        revenueSeq(i) = revenueSeq(i-1) + execPrice(:,i)'*tradingStrats(:,i-1);
        numerator(i-1) = max(revenueSeq(i) + marketPrice(:,i)'*remainPosition(:,i-1) - liability,0);
        capRatioSeq(i-1) = numerator(i-1)/denominator(i-1); 
    end
