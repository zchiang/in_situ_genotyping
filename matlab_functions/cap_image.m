function[O] = cap_image(I,cap)

    O = I;
    O(O>cap) = cap;

end