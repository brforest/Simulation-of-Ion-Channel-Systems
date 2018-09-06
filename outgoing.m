function [] = outgoing(rows, columns, time, data, filename)
arrToSave = zeros(rows,columns+1);
for i=1:rows
   arrToSave(i,1) = time(i,1); 
end
for i=1:rows
   for j=2:columns+1
       arrToSave(i,j) = data(i,j-1);
   end
end
dlmwrite(filename,arrToSave,'\t');
end
