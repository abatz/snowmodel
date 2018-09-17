
function ew=satvap(temp)
a1=(16.78 * temp(:) - 116.9)./(temp(:) + 273.3);
ew=exp(a1) ./((273.15 + temp(:)) * 0.4615);