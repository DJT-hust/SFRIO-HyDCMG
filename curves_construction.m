arrive_curve=[];
min_departure_curve=[];
arrive_temp=0;
departure_temp=0;
N_b=length(delay_torlant_workload.at);
for t=1:T
    for k=1:N_b
        if delay_torlant_workload.at(k)==t
            arrive_temp=arrive_temp+delay_torlant_workload.rt(k);
        end
        if delay_torlant_workload.en(k)==t
            departure_temp=departure_temp+delay_torlant_workload.rt(k);
        end
    end
    arrive_curve=[arrive_curve arrive_temp];
    min_departure_curve=[min_departure_curve departure_temp];
end
