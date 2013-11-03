function [x_avg,x_sd,x_count] = averageXbyYbin(x,y,y_edges)
 
[~,idx] = histc(y,y_edges); % idx returns which bin each point in y goes into

x_avg = zeros(length(y_edges)-1,1);
x_sd = zeros(length(y_edges)-1,1);
x_count = zeros(length(y_edges)-1,1);

x_avg = zeros(length(y_edges)-1,1);
x_sd = zeros(length(y_edges)-1,1);

for iBin = length(y_edges):-1:1 % for each bin...
    
   if sum(idx == iBin) ~= 0 % at least one sample in this bin
      x_avg(iBin) = nanmean(x(idx == iBin)); % compute average of those x's that go in it
      x_sd(iBin) = nanstd(x(idx == iBin)); % compute sd of those x's that go in it
      x_count(iBin) = numel(x(idx == iBin));
   end
 
end