function [capRatioSeq] = simSingleAssetCR(initialPrice, permImpact,temImpact,tau,sigma,initPosition,liability,riskWeight,tradingStrats)
    xiSeq = randn(1,5);
    marketPrice =[initialPrice];
    execPrice = [0];
    remainPosition = initPosition - cumsum(tradingStrats);
    revenueSeq =[0];
    for i =2:6
        marketPrice(i) = marketPrice(i-1) + sqrt(tau)*sigma*xiSeq(i-1) - permImpact*tradingStrats(i-1);
        execPrice(i) = marketPrice(i-1) - temImpact*tradingStrats(i-1);
        numerator(i-1) = max(execPrice(2:i)*tradingStrats(1:i-1)' + marketPrice(i)*remainPosition(i-1) - liability,0);
        denominator(i-1) = riskWeight*marketPrice(i)*remainPosition(i-1);
        capRatioSeq(i-1) = numerator(i-1)/denominator(i-1);
    end

    

