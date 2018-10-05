function ew=satvap(temp)
a1=0.611*exp((17.3* temp(:))./(temp(:) + 237.2));
ew=(a1)./((273.15 + temp(:)) * 0.4615);