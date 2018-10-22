function in_circ = ptb_is_in_circ( p1, p2, dist )
%ptb_is_in_circ Check if point in circle
if ((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2) < dist^2
    in_circ = 1;
else
    in_circ = 0;
end

