function [ Dataout ] = reorderFiles(Datain, order_TCD, order_Monash)
erp_TCD = single(zeros(size(Datain.data)));
erp_TCD(order_TCD,:,:) = Datain.data(order_Monash,:,:);
Dataout.data = erp_TCD(1:64,:,:);

erp_TCD = single(zeros(size(Datain.dataCSD)));
erp_TCD(order_TCD,:,:) = Datain.dataCSD(order_Monash,:,:);
Dataout.dataCSD = erp_TCD(1:64,:,:);

end

