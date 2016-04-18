function [fxh1,fPh1,flike] = kalman(fy,fH,fxp1,fxp3,fV,fV2,fV1,fPp1)
innov = fy - (fH * fxp1) - fxp3;  
S = fH  * fPp1 * ctranspose(fH) + fV1 + (2*fV)*eye(4) + fV2;
fK = fPp1*ctranspose(fH)*inv(S); 
fxh1= fxp1 + fK*innov;
fPh1 = [eye(4)-fK*fH]*fPp1;
flike = exp(-.5*ctranspose(innov)*inv(S)*innov);