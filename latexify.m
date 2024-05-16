function [outputLatex] = latexify(startStr,mainStr)
% Converts to Latex format suitable for report using find and replace
temp = latex(mainStr);
temp = replace(temp, "\mathrm{mdot}", "\dot{m}");
temp = replace(temp, "\mathrm{cH2O}", "c_{\mathrm{H}_2\mathrm{O}}");
temp = replace(temp, "\mathrm{Qdotger}", "\dot{Q}_\mathrm{ger}");
temp = replace(temp, "\Delta _{P,\mathrm{H2O}}", "\Delta P_{\mathrm{H}_2\mathrm{O}}");
temp = replace(temp, "\Delta _{P,\mathrm{R}}", "\Delta P_{R}");
temp = replace(temp, "\Delta _{P,R}", "\Delta P_{R}");

% Only use if you know what to expect
%temp = replace(temp, "\left(\begin{array}", "\left\{\begin{matrix}");
%temp = replace(temp, "\end{array}\right)", "\end{matrix}\right.");

outputLatex = append('\begin{equation}',startStr,temp,'\end{equation}')
disp("Latex prepared above ^^")
end