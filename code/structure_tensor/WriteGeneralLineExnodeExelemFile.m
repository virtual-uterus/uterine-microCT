function WriteGeneralLineExnodeExelemFile(X,Y,Z,EVert,NNum,ENum,Data,Filename,GroupName,DataNames)
 
  % Write exnode file
  fid = fopen(sprintf('%s.exnode',Filename),'w');
  fprintf(fid,' Group name: %s\n #Fields=%d\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  Value index= 1, #Derivatives= 0\n   y.  Value index= 2, #Derivatives= 0\n   z.  Value index= 3, #Derivatives= 0\n',GroupName,size(Data,2)+1);
  for i=1:size(Data,2)
    fprintf(fid,' %d) %s, field, rectangular cartesian, #Components=1\n   1.  Value index= %d, #Derivatives= 0\n',i+1,DataNames{i},3+i);
  end
  for n=1:size(X,1)
    fprintf(fid,' Node: %d\n',NNum(n));
    fprintf(fid,' %f %f %f',X(n),Y(n),Z(n));
    for i=1:size(Data,2)
      fprintf(fid,' %f',Data(n,i));
    end
    fprintf(fid,'\n');
  end
  fclose(fid);

  % Write exelem file
  fid = fopen(sprintf('%s.exelem',Filename),'w');
  fprintf(fid,' Group name: %s\n',GroupName);
  fprintf(fid,' Shape.  Dimension=1\n #Scale factor sets= 1\n   l.Lagrange, #Scale factors= 2\n #Nodes=           2\n #Fields=%d\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n  z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n',size(Data,2)+1);
  for i=1:size(Data,2)
      fprintf(fid,' %d) %s, field, rectangular cartesian, #Components=1\n',i, DataNames{i});
      fprintf(fid,'    value. l.Lagrange, no modify, standard node based.\n');
      fprintf(fid,'      #Nodes=2\n');
      fprintf(fid,'      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n');
  end
  for e=1:size(EVert,1)
      fprintf(fid,' Element:            %d 0 0\n',ENum(e));
      fprintf(fid,'   Nodes:\n     %d %d\n',NNum(EVert(e,:)));
      fprintf(fid,'   Scale factors:\n     1 1\n');
  end
  fclose(fid);
return;
