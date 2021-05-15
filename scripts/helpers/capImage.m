function[O] = capImage(I,val,type)

    if type == "abs"
        cap = val;
    elseif type == "prc"
        cap = prctile(reshape(I,[],1),val);
    end
    
    O = I;
    O(O>cap) = cap;

end