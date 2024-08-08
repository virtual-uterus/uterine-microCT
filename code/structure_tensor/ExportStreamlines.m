function ExportStreamlines(Paths,EField,Fname,GroupName, Region,NodeOffset,ElementOffset)

% Write exnode and exelem file headers
fidn = fopen(sprintf('%s.exnode',Fname),'w');
fprintf(fidn,' Group name: %s\n Region: %s\n #Fields=3\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  Value index= 1, #Derivatives= 0\n   y.  Value index= 2, #Derivatives= 0\n   z.  Value index= 3, #Derivatives= 0\n2) FA1, field, rectangular cartesian, #Components=1\n  1.  Value index= 4, #Derivatives= 0\n',GroupName, Region);
fprintf(fidn,' 3) angle, field, rectangular cartesian, #Components=1\n  1.  Value index= 5, #Derivatives= 0\n');

fide = fopen(sprintf('%s.exelem',Fname),'w');
fprintf(fide,' Group name: %s\n Region: %s\n Shape.  Dimension=1\n #Scale factor sets= 1\n   l.Lagrange, #Scale factors= 2\n #Nodes=           2\n #Fields=4\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n      y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n      z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n       Value indices:     1\n       Scale factor indices:   1\n      2.  #Values=1\n       Value indices:     1\n       Scale factor indices:   2\n  2) EField, field, rectangular cartesian, #Components=1\n   potential.  l.Lagrange, no modify, grid based.\n   #xi1=  1\n 3) FA1, field, rectangular cartesian, #Components=1\n   value.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n        Value indices: 1\n        Scale factor indices: 1\n      2.  #Values=1\n        Value indices: 1\n        Scale factor indices: 2\n',GroupName, Region);
fprintf(fide,' 4) angle, field, rectangular cartesian, #Components=1\n   value.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n      1.  #Values=1\n        Value indices: 1\n        Scale factor indices: 1\n      2.  #Values=1\n        Value indices: 1\n        Scale factor indices: 2\n');
NodeNumber = NodeOffset;
ElementNumber = ElementOffset;

% Loop over paths
for p=1:length(Paths)

    if ~isempty(Paths{p})
        % Loop over foward and back in path
        for d = 1:1 % just go one way with streamlines for now.

            % Loop over path nodes - treat forward and backward path
            % independently
            L3 = Paths{p}{d}(1,1); L2 = Paths{p}{d}(1,2); L1 = Paths{p}{d}(1,3);
            Trace = (L1+L2+L3)/3;
            Denom = sqrt(L1.^2+L2.^2+L3.^3+1e-8);
            FA1 = sqrt(3/2)*(sqrt((L1-Trace).^2+(L2-Trace).^2+(L3-Trace).^2))./Denom;
            angle = Paths{p}{3}(1, 4);

            fprintf(fidn,' Node: %d\n',NodeNumber);
            fprintf(fidn,'  %e %e %e %e %e\n',-Paths{p}{d}(1,1),-Paths{p}{d}(1,2),Paths{p}{d}(1,3),FA1, angle);
            NodeNumber = NodeNumber+1;
            for n=2:size(Paths{p}{d},1)
                L3 = Paths{p}{d}(n,1); L2 = Paths{p}{d}(n,2); L1 = Paths{p}{d}(n,3);
                Trace = (L1+L2+L3)/3;
                Denom = sqrt(L1.^2+L2.^2+L3.^3+1e-8);
                FA1 = sqrt(3/2)*(sqrt((L1-Trace).^2+(L2-Trace).^2+(L3-Trace).^2))./Denom;
                angle = Paths{p}{3}(n, 4);

                fprintf(fidn,' Node: %d\n',NodeNumber);
                fprintf(fidn,'  %e %e %e %e %e\n',-Paths{p}{d}(n,1),-Paths{p}{d}(n,2),Paths{p}{d}(n,3),FA1, angle);
                NodeNumber = NodeNumber+1;

                fprintf(fide,' Element: %d 0 0\n',ElementNumber);
                fprintf(fide,'  Values:\n   %f %f\n',1,1);
                fprintf(fide,'  Nodes:\n   %d %d\n',NodeNumber-2,NodeNumber-1);
                fprintf(fide,'  Scale factors:\n   1 1 \n');
                ElementNumber = ElementNumber+1;
            end

        end
    end
end

fclose(fidn);
fclose(fide);

return;
