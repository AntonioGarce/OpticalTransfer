function [out] = circle_window(x,y,r)
    Nx = length(x);
    Ny = length(y);
    out = zeros(Nx,Ny);
    for i =1:Nx
        for j=1:Ny
            if((x(i)^2+y(j)^2)<r^2)
                out(i,j)=1;
            end
        end
    end
end
