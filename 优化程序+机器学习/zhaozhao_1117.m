%PV
[num_pv,txt_pv,raw_pv]=xlsread('pv_2018_filled.xlsx');
pv_2018 = num_pv(:,6);
%prcie
[num_price,txt_price,raw_price]=xlsread('elspot-prices_2018_hourly_dkk.xls');
dk2 = txt_price(4:8764,10);
price = zeros(8761,1);
for i = 1:8761
    avail = dk2(i);
    number = str2num(avail{1});
    if isempty(number)
        continue
    end
    price(i) = number(1)+number(2)*0.01;
end
price_2018 = zeros(8760,1);
price_2018(1:1994,:) = price(1:1994,:);
price_2018(1995:8760,:) = price(1996:8761,:);

%ten with &without
[num_ten,txt_ten,raw_ten]=xlsread('20180301_ten with and without heatpumps.ods');
load_withhp = num_ten(:,11);
load_withouthp = num_ten(:,22);

