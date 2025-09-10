function a = idctMfg(a)
    % This function computes the Inverse Discrete Cosine Transform 
    % of an image a
    % Input: a = DCT coefficients
    % Output: a = image
    % code Reference https://www.nlb.gov.sg/Libraries/teaching/computer-science-resources/transforming-images/idct/idct-en.html
    


persistent siz ww 
a = a.';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if the variable size has changed and we need to
%%% precompute weights and indicies
precompute=0;
if  ~exist('siz','var')
    precompute=1;
elseif abs(numel(siz)-ndims(a))>0
    precompute=1;
elseif sum(abs(siz-size(a)),2)>0
    precompute=1;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precompute weights and indicies
if precompute
    siz=size(a);
    
    ww = cos( pi/(4*siz(1)+2) * (1:2:2*siz(1)-1)' * (2*siz(1)-1:-2:1) );
    ww = ww * (2/sqrt(2*siz(1)+1));
    
    
end

a = (ww*a).';


end