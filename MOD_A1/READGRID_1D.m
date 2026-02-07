% Grid reader 1D

function grid = READGRID_1D(gridfile)

fid=fopen(char(gridfile),"r");

Nf=fscanf(fid,"%d",1);
xf=zeros(Nf,1);


    for i= 1:Nf
    
        xf(i) = fscanf(fid,"%f",1);
   
    
   end



grid.Nf = Nf;
grid.xf = xf;


fclose(fid);