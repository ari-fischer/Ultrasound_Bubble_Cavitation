% data processing

out_compare = [Ts,OHs,Ars]

SS_compare = SS_data'
SS_compare=SS_compare([1:7,9,10,12:end],:)
rel_abs_error = abs(out_compare-SS_compare)./SS_compare

mean(rel_abs_error,1)