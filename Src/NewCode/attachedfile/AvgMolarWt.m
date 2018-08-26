function Mavg=AvgMolarWt(yw1,yw2,yw3,yw4,yw5,yw6,yw7)
 ns=7;
 M(1:ns)=[69 32 28 18 2 16 44];
Mavg=M(1)*yw1 +M(2)*yw2 +M(3)*yw3 +M(4)*yw4 +M(5)*yw5 +M(6)*yw6 +M(7)*yw7;
end