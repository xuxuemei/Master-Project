% To make N*M classification
% Idx1 is labels of N
% Idx2 is labels of M
% re_Idx is labels of final
function re_Idx = classify_xu(Idx1,Idx2,N,M)
    re_Idx = zeros(size(Idx1));
    for i = 1:N
        for j = 1:M
            for t = 1:365
                if Idx1(t) == i && Idx2(t) == j
                    re_Idx(t) = (i-1)*M+j;
                end
            end
        end
    end           
end
