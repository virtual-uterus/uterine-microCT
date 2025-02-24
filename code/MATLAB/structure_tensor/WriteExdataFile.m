function WriteExdataFile(N,freq,X,Y,Z,E1,E2,E3,CO,COFS,COFN,COSN,L1,L2,L3,FA1,FA2,Angles,GNorm,IFG,HelixF,HelixS,HelixN,fname,gpname,DataN)

  fid = fopen(fname,'w');
  fprintf(fid,' Group name : %s\n #Fields=19\n',gpname);
  fprintf(fid,' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  Value index= 1, #Derivatives= 0\n   y.  Value index= 2, #Derivatives= 0\n   z.  Value index= 3, #Derivatives= 0\n');
  fprintf(fid,' 2) Fiber, field, rectangular cartesian, #Components=3\n   x.  Value index= 4, #Derivatives= 0\n   y.  Value index= 5, #Derivatives= 0\n   z.  Value index= 6, #Derivatives= 0\n');
  fprintf(fid,' 3) Sheet, field, rectangular cartesian, #Components=3\n   x.  Value index= 7, #Derivatives= 0\n   y.  Value index= 8, #Derivatives= 0\n   z.  Value index= 9, #Derivatives= 0\n');
  fprintf(fid,' 4) Normal, field, rectangular cartesian, #Components=3\n   x.  Value index= 10, #Derivatives= 0\n   y.  Value index= 11, #Derivatives= 0\n   z.  Value index= 12, #Derivatives= 0\n');
  fprintf(fid,' 5) CO, field, rectangular cartesian, #Components=1\n   1.  Value index= 13, #Derivatives= 0\n');
  fprintf(fid,' 6) COfs, field, rectangular cartesian, #Components=1\n   1.  Value index= 14, #Derivatives= 0\n');
  fprintf(fid,' 7) COfn, field, rectangular cartesian, #Components=1\n   1.  Value index= 15, #Derivatives= 0\n');
  fprintf(fid,' 8) COsn, field, rectangular cartesian, #Components=1\n   1.  Value index= 16, #Derivatives= 0\n');
  fprintf(fid,' 9) Lf, field, rectangular cartesian, #Components=1\n   1.  Value index= 17, #Derivatives= 0\n');
  fprintf(fid,' 10) Ls, field, rectangular cartesian, #Components=1\n   1.  Value index= 18, #Derivatives= 0\n');
  fprintf(fid,' 11) Ln, field, rectangular cartesian, #Components=1\n   1.  Value index= 19, #Derivatives= 0\n');
  fprintf(fid,' 12) FA1, field, rectangular cartesian, #Components=1\n   1.  Value index= 20, #Derivatives= 0\n');
  fprintf(fid,' 13) FA2, field, rectangular cartesian, #Components=1\n   1.  Value index= 21, #Derivatives= 0\n');
  fprintf(fid,' 14) Angle, field, rectangular cartesian, #Components=1\n   1.  Value index= 22, #Derivatives= 0\n');
  fprintf(fid,' 15) GNorm, field, rectangular cartesian, #Components=1\n   1.  Value index= 23, #Derivatives= 0\n');
  fprintf(fid,' 16) IFGrad, field, rectangular cartesian, #Components=1\n   1.  Value index= 24, #Derivatives= 0\n');
  fprintf(fid,' 17) HelixF, field, rectangular cartesian, #Components=1\n   1.  Value index= 25, #Derivatives= 0\n');
  fprintf(fid,' 18) HelixS, field, rectangular cartesian, #Components=1\n   1.  Value index= 26, #Derivatives= 0\n');
  fprintf(fid,' 19) HelixN, field, rectangular cartesian, #Components=1\n   1.  Value index= 27, #Derivatives= 0\n');
%DataN = 200000; % 200001

for g = 1:length(X),

%for k=1:freq:N(3),
%  for j=1:freq:N(2),
%    for i=1:freq:N(1),
%      g = ((k-1)*N(2)+j-1)*N(1)+i;
%      xx = Origin(1)+(i-1)*Dx;
%      yy = Origin(2)+(j-1)*Dy;
%      zz = Origin(3)+(k-1)*Dz;
        DataN = DataN+1;
        fprintf(fid,' Node: %d\n',DataN);
        fprintf(fid,'   %f %f %f\n',X(g),Y(g),Z(g));
        fprintf(fid,'   %f %f %f\n',E3(g,1),E3(g,2),E3(g,3));
        fprintf(fid,'   %f %f %f\n',E2(g,1),E2(g,2),E2(g,3));
        fprintf(fid,'   %f %f %f\n',E1(g,1),E1(g,2),E1(g,3));
        fprintf(fid,'   %f %f %f %f\n',CO(g),COFS(g),COFN(g),COSN(g));   
        fprintf(fid,'   %f %f %f\n',L3(g),L2(g),L1(g));   
        fprintf(fid,'   %f %f %f %f %f\n',FA1(g),FA2(g),Angles(g),GNorm(g),IFG(g));   
        fprintf(fid,'   %f %f %f\n',HelixF(g),HelixS(g),HelixN(g));   
%    end;
%  end;
end;
fclose(fid);

return;
