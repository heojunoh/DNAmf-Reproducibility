
x = dlmread('generate_text/temp_to_matlab.txt', ',', 1, 0);
y = [];
CPU_time = [];
for i=1:size(x,1)
    a1 = x(i,1);
    a2 = x(i,2);
    a3 = x(i,3);
    t = x(i,4);
    tic; 
    [freq] = Plate(a1,a2,a3,t);
    CPU_time(i) = toc;    
    y(i) = freq;
end
final_data = [x,y',CPU_time'];
dlmwrite('generate_text/temp_to_r.txt',final_data,'delimiter',',','precision', 32);