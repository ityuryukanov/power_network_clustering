function out=mat2latex(M, n, Tlim)
% X
M = M(:,1:end-1); % remove Theur for being unnecessary
ncol = size(M,2);
nrow = size(M,1);
out  = cell(nrow,1);
for i = 1:1:nrow
  row = sprintf('\\small{%.0f} & \\small{%.0f} & ', n, i+1);
  for j = 2:1:ncol
    if j==2
      sigdig = 3;
      a = round(M(i,2),sigdig,"significant");
      d = dig_after_dot(a);
      if d>=sigdig, a=round(a,sigdig-1); d=sigdig-1; end
      colj = sprintf(['\\small{%.',num2str(d),'f} & '], a); 
    end
    if j==3
      sigdig = 2;
      if M(i,j)>=Tlim
        a = round(M(i,4),sigdig,"significant");
        d = dig_after_dot(a);
        if d>=sigdig, a=round(a,sigdig-1); d=sigdig-1; end
        colj = sprintf(['\\small{%.',num2str(d),'f~\\%%} & '], a);        
      else
        a = round(M(i,3),sigdig,"significant");
        d = dig_after_dot(a);
        if d>=sigdig, a=round(a,sigdig-1); d=sigdig-1; end
        colj = sprintf(['\\small{%.',num2str(d),'f~s} & '], a);        
      end
    end
    if j==4
      continue
    end
    if j>=5 && j<=ncol-1
      sigdig = 3;
      a = round(M(i,j),sigdig,"significant");
      d = dig_after_dot(a);
      if d>=sigdig, a=round(a,sigdig-1); d=sigdig-1; end
      colj = sprintf(['\\small{%.',num2str(d),'f} & '], a);      
    end    
    if j==ncol
      sigdig = 3;
      a = round(M(i,j),sigdig,"significant");
      d = dig_after_dot(a);      
      if d>=sigdig, a=round(a,sigdig-1); d=sigdig-1; end
      colj = sprintf(['\\small{%.',num2str(d),'f}\\\\'], a);
    end
    row = strcat(row,colj);
  end
  if i==nrow
    row = strcat(row,'\hline');
  end
  out{i} = row;
end

end


function n = dig_after_dot(x)
x = abs(x);  %in case of negative numbers
n=0;
while abs(floor(x*10^n)-x*10^n)>1e-14
  n=n+1;
end
end