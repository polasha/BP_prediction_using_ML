SevendataNan = readtable('seven2polash.csv');
%cleandData = readtable('cleaned_further.csv');

SevendataNanM = table2array(SevendataNan);
SevendataNanM(isinf(SevendataNanM)|isnan(SevendataNanM)) = 0;


csvwrite('SevendataNanM.csv', SevendataNanM);