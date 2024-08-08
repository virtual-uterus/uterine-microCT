function WriteGeneralExdataFile(X,Y,Z,DNum,DataS,DataV,Filename,GroupName,DataSNames,DataVNames)

% DataS - scalar data in N points by M data (i.e. labels) array
% DataV  - vector data in cell array of M data (i.e. labels)

  fid = fopen(sprintf('%s.exdata',Filename),'w');
  fprintf(fid,' Group name: %s\n #Fields=%d\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  Value index= 1, #Derivatives= 0\n   y.  Value index= 2, #Derivatives= 0\n   z.  Value index= 3, #Derivatives= 0\n',GroupName,1+size(DataS,2)+size(DataV,2));
  Index = 4;
  for i=1:size(DataS,2)
    fprintf(fid,' %d) %s, field, rectangular cartesian, #Components=1\n   1.  Value index= %d, #Derivatives= 0\n',i+1,DataSNames{i},Index);
    Index = Index+1;
  end
  for i=1:length(DataV)
    fprintf(fid,' %d) %s, field, rectangular cartesian, #Components=3\n   x.  Value index= %d, #Derivatives= 0\n   y.  Value index= %d, #Derivatives= 0\n   z.  Value index= %d, #Derivatives= 0\n',size(DataS,2)+i+1,DataVNames{i},Index,Index+1,Index+2);
    Index = Index+3;
  end

  for n=1:size(X,1)
    fprintf(fid,' Node: %d\n',DNum(n));
    fprintf(fid,' %f %f %f',X(n),Y(n),Z(n));
    for i=1:size(DataS,2)
      fprintf(fid,' %f',DataS(n,i));
    end
    fprintf(fid,'\n');

    for i=1:length(DataV)
        fprintf(fid,' %f %f %f\n',DataV{i}(n,:));
    end
  end

  fclose(fid);
return;
