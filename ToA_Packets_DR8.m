function [Payload_CRC_ToA_DR8,Payload_CRC_ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M)    
   
    for PL=1:length(Payload)
    Payload_CRC_ToA_DR8(PL) = Header_ToA_DR8  + ceil((Payload(PL) + 2)/M)*(102/1000); 
    Payload_CRC_ToA_DR8_WH(PL) = ceil((Payload(PL) + 2)/M)*(102/1000); 
    end
end