function [flush_z, flush_t] = flushing(i, z, time, flush_z, flush_t, flushing_time)

i=i-1;
time_i = time(i);
time_end = time(i)-flushing_time;
while time_i>=time_end
    flush_z(i)=z(i);
    flush_t(i)=time(i);
    i=i-1;
    time_i=time(i);
end

end

