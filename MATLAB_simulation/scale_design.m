function g=scale_design(kind,g,t)
if kind==2,g(1:4)=g(1:4).*[t^2 t^3 t t];else,g(1:2)=g(1:2)*t;end
end
