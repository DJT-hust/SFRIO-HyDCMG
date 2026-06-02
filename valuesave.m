Theta_real = aT*temperature_forecast_ontime_24h(1:T) ...
           + aH*RH_forecast ...
           + aP*(pref-Press_forecast);

n_el_value = value(n_el);
P_elz_value = value(P_elz);
n_fc_value = value(n_fc);
P_fc_value = value(P_fc);
n_H2_store_value = value(n_H2_store);
P_grid_value = value(P_grid);
P_cool_value = value(P_cool);
P_re_value = value(P_re);
f_RT_value = value(f_RT);
f_DT_value = value(f_DT);
P_DC_value = value(P_DC);
H_DC_value = value(H_DC);
W_RT_value = value(W_RT_sdp);
W_DT_value = value(W_DT);
W_DT_s_value = value(W_DT_s);
real_time_workload_arrive_p_value = value(real_time_workload_arrive_p);
windpower_uncertainty_p_value = value(windpower_uncertainty_p);
real_time_workload_arrive_o_value = value(real_time_workload_arrive_o);
windpower_uncertainty_o_value = value(windpower_uncertainty_o);
if exist('Theta_err_p', 'var')
    Theta_err_p_value = value(Theta_err_p);
else
    Theta_err_p_value = NaN(T, 1);
end
if exist('Theta_err_o', 'var')
    Theta_err_o_value = value(Theta_err_o);
else
    Theta_err_o_value = NaN(T, 1);
end

obj1_value = price_grid * sum(P_grid_value) + price_el * sum(P_elz_value) + price_fc * sum(P_fc_value) + price_re * sum(P_re_value) + price_cool * sum(P_cool_value);
obj2_p_value = sum(windpower_forecast'+windpower_uncertainty_p_value - P_re_value);
obj2_o_value = sum(windpower_forecast'+windpower_uncertainty_o_value - P_re_value);
objmean_value = 0.5 * (obj2_p_value + obj2_o_value);
objdev_value = 0.5 * (obj2_p_value - obj2_o_value);

obj_1_list=[obj_1_list; obj1_value];
    obj_2_o_list=[obj_2_o_list; obj2_o_value];
    obj_2_p_list=[obj_2_p_list; obj2_p_value];
    obj_2_mean_list=[obj_2_mean_list; objmean_value];
    obj_2_dev_list=[obj_2_dev_list; objdev_value];


temperature_res = zeros(T,1);
temperature_res(1) = temperature_1;

for t = 2:T
    temperature_res(t) = ...
        varkappa^(t-1)*temperature_1 ...
      + (1-varkappa)*varkappa^(t-1) * ...
        sum(kappa_vec(1:t-1) .* ...
        (Theta_real(1:t-1) ...
        + heat2teperature_para * H_DC_value(1:t-1) ...
        - eta_ele * P_cool_value(1:t-1)));
end

temperature_over_time_ratio = sum(temperature_res(2:T)>temperature_set)/(T-1);

if exist('Theta_err_samples', 'var')
    Theta_err_validation = Theta_err_samples(1:T,:);
else
    Theta_err_validation = aT*temperature_forecast_err(1:T,:) ...
        + aH*humidity_forecast_err(1:T,:) ...
        + aP*(-pressure_forecast_err(1:T,:));
end
temperature_scenario_res = zeros(T, size(Theta_err_validation, 2));
temperature_scenario_res(1, :) = temperature_1;
for t = 2:T
    thermal_input = Theta_forecast(1:t-1) ...
        + heat2teperature_para * H_DC_value(1:t-1) ...
        - eta_ele * P_cool_value(1:t-1);
    thermal_input = thermal_input(:);
    temperature_scenario_res(t, :) = ...
        varkappa^(t-1)*temperature_1 ...
      + (1-varkappa)*varkappa^(t-1) * ...
        sum(kappa_vec(1:t-1) .* (thermal_input + Theta_err_validation(1:t-1,:)), 1);
end
temperature_scenario_violate = any(temperature_scenario_res(2:T,:) > temperature_set + 1e-8, 1);
temperature_trajectory_violation_prob = mean(temperature_scenario_violate);
temperature_over_avg_time_ratio = mean(mean(temperature_scenario_res(2:T,:) > temperature_set + 1e-8, 1));
temperature_over_prob = temperature_over_avg_time_ratio;
temperature_over_max_violation = max(temperature_scenario_res(:) - temperature_set);
if exist('omega', 'var') && exist('mu', 'var') && exist('rho', 'var')
    temperature_sget_out_of_set_prob = ...
        mean(any(omega * (Theta_err_validation - mu') > rho + 1e-8, 1));
else
    temperature_sget_out_of_set_prob = NaN;
end

zone_count = size(f_RT_value, 2);
P_zone_value = P_idle + (P_peak - P_idle) / f_max * (f_RT_value + f_DT_value);
H_zone_value = kappa_heat * P_zone_value;
P_cool_zone_value = zeros(T, zone_count);
zone_power_sum = sum(P_zone_value, 2);
positive_power_idx = zone_power_sum > 1e-8;
if any(positive_power_idx)
    P_cool_zone_value(positive_power_idx, :) = ...
        bsxfun(@times, P_cool_value(positive_power_idx), ...
        bsxfun(@rdivide, P_zone_value(positive_power_idx, :), ...
        zone_power_sum(positive_power_idx)));
end
if any(~positive_power_idx)
    P_cool_zone_value(~positive_power_idx, :) = ...
        repmat(P_cool_value(~positive_power_idx) ./ zone_count, 1, zone_count);
end

temperature_zone_scenario_res = zeros(T, size(Theta_err_validation, 2), zone_count);
temperature_zone_scenario_res(1, :, :) = temperature_1;
temperature_zone_over_prob = zeros(zone_count, 1);
temperature_zone_trajectory_violation_prob = zeros(zone_count, 1);
temperature_zone_over_max_violation = zeros(zone_count, 1);
for z = 1:zone_count
    for t = 2:T
        thermal_input_zone = Theta_forecast(1:t-1) ...
            + heat2teperature_para * H_zone_value(1:t-1, z) ...
            - eta_ele * P_cool_zone_value(1:t-1, z);
        thermal_input_zone = thermal_input_zone(:);
        temperature_zone_scenario_res(t, :, z) = ...
            varkappa^(t-1)*temperature_1 ...
          + (1-varkappa)*varkappa^(t-1) * ...
            sum(kappa_vec(1:t-1) .* ...
            (thermal_input_zone + Theta_err_validation(1:t-1,:)), 1);
    end
    zone_scenario = temperature_zone_scenario_res(2:T, :, z);
    zone_violate = zone_scenario > temperature_set + 1e-8;
    temperature_zone_over_prob(z) = mean(zone_violate(:));
    temperature_zone_trajectory_violation_prob(z) = mean(any(zone_violate, 1));
    temperature_zone_over_max_violation(z) = max(zone_scenario(:) - temperature_set);
end

if exist('wind_err_samples', 'var')
    re_err_validation = wind_err_samples(1:T,:);
else
    re_err_validation = 0.5 * windpowerforuncertaintyset(1:T,:);
end
re_available_scenario = windpower_forecast(:) + re_err_validation;
re_scenario_violate = P_re_value(:) > re_available_scenario + 1e-8;
re_over_prob = mean(re_scenario_violate(:));
re_trajectory_violation_prob = mean(any(re_scenario_violate, 1));
re_violation_amount = P_re_value(:) - re_available_scenario;
re_over_max_violation = max(re_violation_amount(:));
if exist('omega_w', 'var') && exist('mu_w', 'var') && exist('rho_w', 'var')
    re_sget_out_of_set_prob = ...
        mean(any(omega_w * (re_err_validation - mu_w') > rho_w + 1e-8, 1));
else
    re_sget_out_of_set_prob = NaN;
end
constraint_over_prob = max(temperature_over_prob, re_over_prob);

if ~exist('temperature_over_time_ratio_list', 'var')
    temperature_over_time_ratio_list = [];
end
if ~exist('temperature_trajectory_violation_prob_list', 'var')
    temperature_trajectory_violation_prob_list = [];
end
if ~exist('temperature_over_avg_time_ratio_list', 'var')
    temperature_over_avg_time_ratio_list = [];
end
if ~exist('temperature_over_max_violation_list', 'var')
    temperature_over_max_violation_list = [];
end
if ~exist('temperature_validation_sample_count_list', 'var')
    temperature_validation_sample_count_list = [];
end
if ~exist('temperature_sget_out_of_set_prob_list', 'var')
    temperature_sget_out_of_set_prob_list = [];
end
if ~exist('temperature_zone_over_prob_mat', 'var')
    temperature_zone_over_prob_mat = [];
end
if ~exist('temperature_zone_trajectory_violation_prob_mat', 'var')
    temperature_zone_trajectory_violation_prob_mat = [];
end
if ~exist('temperature_zone_over_max_violation_mat', 'var')
    temperature_zone_over_max_violation_mat = [];
end
if ~exist('re_over_prob_list', 'var')
    re_over_prob_list = [];
end
if ~exist('re_trajectory_violation_prob_list', 'var')
    re_trajectory_violation_prob_list = [];
end
if ~exist('re_over_max_violation_list', 'var')
    re_over_max_violation_list = [];
end
if ~exist('re_sget_out_of_set_prob_list', 'var')
    re_sget_out_of_set_prob_list = [];
end
if ~exist('constraint_over_prob_list', 'var')
    constraint_over_prob_list = [];
end
if ~exist('H_DC_solution_mat', 'var')
    H_DC_solution_mat = [];
end
if ~exist('P_cool_solution_mat', 'var')
    P_cool_solution_mat = [];
end
if ~exist('P_DC_solution_mat', 'var')
    P_DC_solution_mat = [];
end
if ~exist('P_re_solution_mat', 'var')
    P_re_solution_mat = [];
end
if ~exist('P_zone_solution_tensor', 'var')
    P_zone_solution_tensor = [];
end
if ~exist('P_cool_zone_solution_tensor', 'var')
    P_cool_zone_solution_tensor = [];
end
temperature_res_list=[temperature_res_list temperature_res];
temperature_over_prob_list=[temperature_over_prob_list temperature_over_prob];
temperature_over_time_ratio_list=[temperature_over_time_ratio_list temperature_over_time_ratio];
temperature_trajectory_violation_prob_list=[temperature_trajectory_violation_prob_list temperature_trajectory_violation_prob];
temperature_over_avg_time_ratio_list=[temperature_over_avg_time_ratio_list temperature_over_avg_time_ratio];
temperature_over_max_violation_list=[temperature_over_max_violation_list temperature_over_max_violation];
temperature_validation_sample_count_list=[temperature_validation_sample_count_list size(Theta_err_validation, 2)];
temperature_sget_out_of_set_prob_list=[temperature_sget_out_of_set_prob_list temperature_sget_out_of_set_prob];
temperature_zone_over_prob_mat=[temperature_zone_over_prob_mat temperature_zone_over_prob(:)];
temperature_zone_trajectory_violation_prob_mat=[temperature_zone_trajectory_violation_prob_mat temperature_zone_trajectory_violation_prob(:)];
temperature_zone_over_max_violation_mat=[temperature_zone_over_max_violation_mat temperature_zone_over_max_violation(:)];
re_over_prob_list=[re_over_prob_list re_over_prob];
re_trajectory_violation_prob_list=[re_trajectory_violation_prob_list re_trajectory_violation_prob];
re_over_max_violation_list=[re_over_max_violation_list re_over_max_violation];
re_sget_out_of_set_prob_list=[re_sget_out_of_set_prob_list re_sget_out_of_set_prob];
constraint_over_prob_list=[constraint_over_prob_list constraint_over_prob];
H_DC_solution_mat=[H_DC_solution_mat H_DC_value(:)];
P_cool_solution_mat=[P_cool_solution_mat P_cool_value(:)];
P_DC_solution_mat=[P_DC_solution_mat P_DC_value(:)];
P_re_solution_mat=[P_re_solution_mat P_re_value(:)];
if isempty(P_zone_solution_tensor)
    P_zone_solution_tensor = P_zone_value;
else
    P_zone_solution_tensor = cat(3, P_zone_solution_tensor, P_zone_value);
end
if isempty(P_cool_zone_solution_tensor)
    P_cool_zone_solution_tensor = P_cool_zone_value;
else
    P_cool_zone_solution_tensor = cat(3, P_cool_zone_solution_tensor, P_cool_zone_value);
end
