function [  ] = save_kymos( kymos, fold, i  )
%

st =find(~isnan(mean(kymos)),1,'first');
en = find(~isnan(mean(kymos)),1,'last');
if st < en
  name = strcat([fold 'mol' num2str(i) '.tif']);
  imwrite(uint16(kymos(:,st:en)), name);
end

end

