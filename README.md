# Replication of the paper: "Stock returns, real activity, inflation, and money" by Eugene Fama, 1981.
Fama's analysis makes the case that the relationship between inflation and the stock market work through real activity. He does this with three sets of regressions: (1) inflation and indicators of real activity; (2) relations among the real activity indicators; (3) real stock returns and future real activity. For (2) and (3), Fama incrementally added the expected and unexpected components of inflation, so that we can assess their marginal contribution in explaining the development of the real variables and real stock returns.

I extended Fama's sample to include the data from 1977-2019. In a sharp contrast to the simple inflation-stock market regressions in the first set of regressions, the estimates of the coefficients in the sets of regressions (2)-(3), as well as the overall explanatory power of the regressors, remain relatively stable over the sample. The sign of the coefficients is stable throughout the sample: the positive relationship between different indicators of real activity, and between those indicators and real stock returns; the insignificant marginal contribution of the different components of inflation. 

## Code
This project is a collection of Matlab scripts I used to replicate the paper. Some of the files pertain to U.S. data, and some to Israeli data (notably those with a suffix _il).

# Files
## replication scripts
tab1a.m  # Table 1a, Fama (1981), US
tab1b.m  # Table 1b, Fama (1981), US

tab2.m  # Table 2, Fama (1981), US
tab3.m  # Table 3, Fama (1981), US
tab4.m  # Table 4, Fama (1981), US
tab5.m  # Table 5, Fama (1981), US
tab5xr.m  # Table 5 with excess returns, Fama (1981), US
tab6.m  # Table 6, Fama (1981), US

tab1a_il.m  # Table 1a, Fama (1981), Israel

tab2_il.m  # Table 2, Fama (1981), Israel
tab3_il.m  # Table 3, Fama (1981), Israel
tab4_il.m  # Table 4, Fama (1981), Israel
tab5_il.m  # Table 5, Fama (1981), Israel
tab6_il.m  # Table 6, Fama (1981), Israel



# other papers: Nelson (1976), and Fama and Gibbons (1982)

## Nelson, 1976. "Inflation and rates of return on common stocks." - JF
nelson76_tab1_2.m  # Tables 1-2, Nelson 1976, US
nelson76_EI_UI.m  # Tables 1-2 with a breakdown to expected and unexpected, Nelson 1976, US
nelson76_EI_UI_OOS.m  # Tables 1-2 with a breakdown to expected and unexpected (out of sample), Nelson 1976, US

nelson76_EI_UI_il.m  # Tables 1-2 with a breakdown to expected and unexpected, Nelson 1976, Israel
nelson76_tab1_2_il.m  # Tables 1-2, Nelson 1976, Israel
nelson76_il_struct_break.m  # working sciprt to check for structural breaks in Nelson's model


## Fama, 1975. "Short-Term Interest Rates as Predictors of Inflation." - AER
regs_01.m # Estimating expected inflation from interest rate models.


## Fama & Gibbons, 1982. "Inflation, real returns and capital investment." - JME
wandering_intercept_script.m


## scripts to plot and study the time series of variable used in Fama (1981) regressions
cx.m # capital expenditure
ns.m # net capital stock
profits.m
tax.m # component in deriving profits after tax
money_supply.m
interest_paid.m # component in deriving profits after tax
fixed_investment.m
capital_consumption.m

cx_il.m # capital expenditure, Israel
ns_il.m # net capital stock, Israel


## scripts to generate expected inflation
gen_EITB.m # interest rate, US
gen_EITB_OOS.m # interest, US, out of sample
gen_EIMD.m # money demand, US

gen_EITB_il.m # interest rate, Israel
gen_EIMD_il.m # money demand, Israel


gen_EITB_il_ZRD.m # interest rate, Israel, using also the real rate (ZRD, yield curve estimate) -- essentially BEI
gen_EITB_il_GALIL.m # interest rate, Israel, using also the real rate (GALIL, observed real rate) -- essentially BEI

gen_EITB_il_ZRD_nss.m # interest rate, Israel, using also the real rate (ZRD_nss, yield curve estimate from Nelson-Siegel-Svensson) -- essentially BEI
gen_EITB_il_ZRD_oos.m # interest rate, Israel, using also the real rate (ZRD, yield curve estimate) -- essentially BEI. out of sample model

HAZ_PI.m # Israeli professional forecasters expected inflation

## script to study how various variables can predict the residual from the interest rate model.
predict_UITB_il.m



## basic regressions of the stock market index (LHS) on expected inflation (RHS)
expected_inflation_stock_markets.m
expected_inflation_stock_markets_monthly.m  # monthly version of the file before
expected_inflation_stock_markets_monthly_Shiller.m  # variation using R. Shiller data

TB_model.m # Conditional model to predict next month inflation using one-month interest rates. 



## scripts for plotting various data
plot_UITB.m
plot_fear_index.m
plot_EITB_il.m
plot_UITB_il.m
plot_1980s_il.m
plot_ROC_CX_NS.m
plotting_EI_UI.m
plot_UIAR_UIPF.m
plot_UITB_UIAR.m
plot_UITB_UIAR_UIHA.m







