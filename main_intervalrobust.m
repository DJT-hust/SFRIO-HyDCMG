clear
close all
warning off

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
addpath(genpath(scriptDir));
dataDir = fullfile(scriptDir, 'data');

T = 12;
delta_list = [0.4];
epsilon_list = [0.4];
paradef

load(fullfile(dataDir, 'temperatureforuncertaintyset1.mat'))
load(fullfile(dataDir, 'temperature_ontime.mat'))
load(fullfile(dataDir, 'windpower_24h.mat'))
load(fullfile(dataDir, 'windpowerforuncertaintyset.mat'))
load(fullfile(dataDir, 'realtimeworkload_forecast_tencent.mat'))
load(fullfile(dataDir, 'humidity_24h.mat'))
load(fullfile(dataDir, 'pressure_24h.mat'))
load(fullfile(dataDir, 'humidity_ontime.mat'))
load(fullfile(dataDir, 'pressure_ontime.mat'))

RH_forecast   = RH_forecast_ontime_24h(1:T);
Press_forecast= pressure_ontime_24h(1:T);

aT = 1;
aH = 0.3;
aP = 0.02;

pref = 1013;

temperature_forecast_ontime_24h = temperature_forecast_ontime_24h';

ops = sdpsettings('solver', 'gurobi+', 'allowmilp', 1, 'gurobi.QCPDual', 1);

n_0 = n_H2_store_min + (n_H2_store_max - n_H2_store_min) / 2;
SoC_0 = 0.22;
temperature_1 = 32.5;
accumulate_dt_workload = 0;
curves_construction

Theta_forecast = aT*temperature_forecast_ontime_24h(1:T) ...
    + aH*RH_forecast ...
    + aP*(pref-Press_forecast);




for epsit = 1:length(epsilon_list)
    epsilon = epsilon_list(epsit);
    disp(['epsilon = ', num2str(epsilon)]);

    for delit = 1:length(delta_list)
        delta = delta_list(delit);
        disp(['delta = ', num2str(delta)]);
        cost_mat=[];
        wind_mat=[];
        interval_para=0.2;
        chance_target_epsilon = epsilon;
        sf_confidence_delta = delta;
        for bigint=1:1
            options = sdpsettings('solver', 'gurobi', 'verbose', 1);
            Theta_err_samples = ...
                aT*temperature_forecast_err(1:T,:) ...
                + aH*humidity_forecast_err(1:T,:) ...
                + aP*(-pressure_forecast_err(1:T,:));
            [omega, mu, rho] = Construcrt_Box_uncertainty_set(Theta_err_samples,epsilon,delta);

            wind_err_samples = 0.5 * windpowerforuncertaintyset(1:T,:);
            [omega_w, mu_w, rho_w]=Construcrt_Box_uncertainty_set(wind_err_samples,epsilon,delta);
            windpower_forecast=windpower_forecast(1:T);
            [windpower_max, windpower_min] = uncertainty_set_construct_wind_power(windpower_forecast,interval_para);
            obj_1_list=[];
            obj_2_o_list=[];
            obj_2_p_list=[];
            obj_2_mean_list=[];
            obj_2_dev_list=[];
            temperature_res_list=[];


            constraint = [];
            constraint_p = [];
            constraint_o = [];
            vardef
            constraint = [constraint n_el == eta_elz / H_HV_H2 * P_elz];
            constraint = [constraint P_el_min * ones(T, 1).*I_elz <= P_elz <= P_el_max * ones(T, 1).*I_elz];
            constraint = [constraint P_fc' == eta_fc * L_HV_H2 * n_fc'];
            constraint=[constraint I_fc+I_elz<=ones(T,1)];
            constraint = [constraint P_fc_min * ones(T, 1).*I_fc <= P_fc <= P_fc_max * ones(T, 1).*I_fc];
            n_pre = [n_0; n_H2_store(1:T - 1)];
            constraint = [constraint n_H2_store == n_pre + n_el - n_fc];
            constraint = [constraint n_H2_store_min * ones(T, 1) <= n_H2_store <= n_H2_store_max * ones(T, 1)];
            constraint_p = [constraint_p zeros(T, 1) <= P_re, P_re <= windpower_forecast'+windpower_uncertainty_p, windpower_forecast'+windpower_uncertainty_p>=zeros(T, 1)];
            constraint_o = [constraint_o zeros(T, 1) <= P_re, P_re <= windpower_forecast'+windpower_uncertainty_o, windpower_forecast'+windpower_uncertainty_o>=zeros(T, 1)];
            constraint_p = [constraint_p windpower_min<=windpower_uncertainty_p'<=windpower_max];
            constraint_o = [constraint_o windpower_min<=windpower_uncertainty_o'<=windpower_max];
            re_sf_margin = zeros(T, 1);
            re_sf_order = max(1, ceil(chance_target_epsilon * size(wind_err_samples, 2)));
            re_sf_order = min(size(wind_err_samples, 2), re_sf_order);
            for t = 1:T
                re_err_effect = sort(wind_err_samples(t, :).');
                re_sf_margin(t) = re_err_effect(re_sf_order);
                constraint = [constraint ...
                    P_re(t) <= windpower_forecast(t) + re_sf_margin(t)];
            end
            constraint = [constraint 1 / ddl_real_time_workload * ones(T, N_s) <= f_RT * serve_rate - W_RT_sdp];
            constraint=[constraint f_DT>=zeros(size(f_DT))];
            constraint=[constraint f_RT>=zeros(size(f_RT))];
            constraint=[constraint W_DT_s>=zeros(size(W_DT_s))];
            for t = 1:T
                constraint=[constraint sum(W_RT_sdp(t, :)) == realtimeworkload_forecast(t)];
                constraint = [constraint sum(W_DT_s(t, :)) == W_DT(t)];
                constraint = [constraint min_departure_curve(t) <= sum(W_DT(1:t)) <= arrive_curve(t)];
                constraint = [constraint f_DT(t, :) == para_dt_workload * W_DT_s(t, :)];
                constraint = [constraint f_min * ones(1, N_s) <= f_DT(t, :) + f_RT(t, :) <= f_max * ones(1, N_s)];
                constraint = [constraint P_DC(t) == sum(P_idle * ones(1, N_s) + (P_peak - P_idle) / f_max * (f_RT(t, :) + f_DT(t, :)))];
            end











            constraint = [constraint H_DC == kappa_heat * P_DC];

            kappa_vec = [];
            for t = 1:T
                kappa_vec = [kappa_vec; varkappa^(-t)];
            end

            Theta_err_p = sdpvar(T,1);
            Theta_err_o = sdpvar(T,1);

            for t = 2:T
                constraint_p = [constraint_p ...
                    varkappa^(t-1)*temperature_1 + ...
                    (1-varkappa)*varkappa^(t-1) * ...
                    sum(kappa_vec(1:t-1) .* ...
                    (Theta_forecast(1:t-1) + Theta_err_p(1:t-1) ...
                    + heat2teperature_para * H_DC(1:t-1) ...
                    - eta_ele * P_cool(1:t-1))) <= temperature_set];

                constraint_o = [constraint_o ...
                    varkappa^(t-1)*temperature_1 + ...
                    (1-varkappa)*varkappa^(t-1) * ...
                    sum(kappa_vec(1:t-1) .* ...
                    (Theta_forecast(1:t-1) + Theta_err_o(1:t-1) ...
                    + heat2teperature_para * H_DC(1:t-1) ...
                    - eta_ele * P_cool(1:t-1))) <= temperature_set];
            end

            constraint_p = [constraint_p omega*(Theta_err_p-mu')<=rho];
            constraint_o = [constraint_o omega*(Theta_err_o-mu')<=rho];

            temperature_sf_margin = zeros(T, 1);
            temperature_sf_order = max(1, ceil((1 - chance_target_epsilon) * size(Theta_err_samples, 2)));
            temperature_sf_order = min(size(Theta_err_samples, 2), temperature_sf_order);
            for t = 2:T
                theta_err_effect = ...
                    (1-varkappa)*varkappa^(t-1) * ...
                    sum(kappa_vec(1:t-1) .* Theta_err_samples(1:t-1,:), 1);
                theta_err_effect = sort(theta_err_effect(:));
                temperature_sf_margin(t) = theta_err_effect(temperature_sf_order);
                temperature_nominal_sf = ...
                    varkappa^(t-1)*temperature_1 + ...
                    (1-varkappa)*varkappa^(t-1) * ...
                    sum(kappa_vec(1:t-1) .* ...
                    (Theta_forecast(1:t-1) ...
                    + heat2teperature_para * H_DC(1:t-1) ...
                    - eta_ele * P_cool(1:t-1)));
                constraint = [constraint ...
                    temperature_nominal_sf + temperature_sf_margin(t) <= temperature_set];
            end

            constraint = [constraint P_cool_min * ones(T, 1) <= P_cool <= P_cool_max * ones(T, 1)];
            constraint = [constraint P_grid + P_fc + P_re == P_DC + P_cool + P_elz];
            constraint = [constraint P_grid_min * ones(T, 1) <= P_grid <= P_grid_max * ones(T, 1)];
            obj1 = price_grid * sum(P_grid) + price_el * sum(P_elz) + price_fc * sum(P_fc) + ...
                price_re * sum(P_re) + price_cool * sum(P_cool);
            obj2_p = sum(windpower_forecast'+windpower_uncertainty_p - P_re);
            obj2_o = sum(windpower_forecast'+windpower_uncertainty_o - P_re);
            objmean = 0.5 * (obj2_p + obj2_o);
            objdev = 0.5 * (obj2_p - obj2_o);
            ops = sdpsettings('solver', 'gurobi+', 'allowmilp', 1, 'gurobi.QCPDual', 1,'verbose', 1,'gurobi.MIPGap',0.003,'gurobi.timelimit',1800);
            [KKT_Constraints_p, details_p] = kkt(constraint_p, -obj2_p, [P_re H_DC P_cool], ops);
            [KKT_Constraints_o, details_o] = kkt(constraint_o, obj2_o, [P_re H_DC P_cool], ops);
            dual_p=sdpvar(size(details_p.A,1),1);
            dual_o=sdpvar(size(details_p.A,1),1);
            KKT_cons_p = [];
            KKT_cons_p = [KKT_cons_p details_p.c + (details_p.A)' * dual_p == 0];
            KKT_cons_p = [KKT_cons_p dual_p .* ((details_p.A) * details_p.primal - details_p.b) == 0];
            KKT_cons_p = [KKT_cons_p dual_p >= 0, dual_p <= 2000, dual_p <= details_p.dualbounds];

            KKT_cons_o = [];
            KKT_cons_o = [KKT_cons_o details_o.c + (details_o.A)' * dual_o == 0];
            KKT_cons_o = [KKT_cons_o dual_o .* ((details_o.A) * details_o.primal - details_o.b) == 0];
            KKT_cons_o = [KKT_cons_o dual_o >= 0, dual_o <= 2000, dual_o <= details_o.dualbounds];
            disp('KKT done')
            Constraints = [constraint constraint_p constraint_o KKT_cons_p KKT_cons_o];
            weight_price = 1:80;
            weight_mean = 1:80;
            weight_dev = 1:80;
            forleagth=size(weight_price,2);

            obj1v=3.3424;
            devv=5.2853;
            meanv=2.0131;


            ops = sdpsettings('solver', 'gurobi+', 'allowmilp', 1, 'gurobi.QCPDual', 1,'verbose', 0,'gurobi.MIPGap',0.003);
            Payoff_mat=zeros(2);
            temperature_over_prob_list=[];
            temperature_over_time_ratio_list=[];
            temperature_trajectory_violation_prob_list=[];
            temperature_over_avg_time_ratio_list=[];
            temperature_over_max_violation_list=[];
            temperature_validation_sample_count_list=[];
            temperature_sget_out_of_set_prob_list=[];
            temperature_zone_over_prob_mat=[];
            temperature_zone_trajectory_violation_prob_mat=[];
            temperature_zone_over_max_violation_mat=[];
            re_over_prob_list=[];
            re_trajectory_violation_prob_list=[];
            re_over_max_violation_list=[];
            re_sget_out_of_set_prob_list=[];
            constraint_over_prob_list=[];
            H_DC_solution_mat=[];
            P_cool_solution_mat=[];
            P_DC_solution_mat=[];
            P_zone_solution_tensor=[];
            P_cool_zone_solution_tensor=[];
            solvingtime_list=[];
            t1=clock;
            sol1 = optimize(Constraints, obj1, ops);
            t2=clock;
            deltatime=etime(t2,t1);
            solvingtime_list=[solvingtime_list deltatime];
            valuesave
            Payoff_mat(1,1)=obj1_value;
            Payoff_mat(1,2)=objmean_value;
            t1=clock;
            sol2 = optimize(Constraints, objmean, ops);
            t2=clock;
            deltatime=etime(t2,t1);
            solvingtime_list=[solvingtime_list deltatime];
            valuesave
            Payoff_mat(2,1)=obj1_value;
            Payoff_mat(2,2)=objmean_value;

            opt_f=min(Payoff_mat);
            pes_f=max(Payoff_mat);

            K=100;
            L=20;
            m=0;

            obj_1_list=[];
            obj_2_o_list=[];
            obj_2_p_list=[];
            obj_2_mean_list=[];
            obj_2_dev_list=[];
            temperature_over_prob_list=[];
            temperature_over_time_ratio_list=[];
            temperature_trajectory_violation_prob_list=[];
            temperature_over_avg_time_ratio_list=[];
            temperature_over_max_violation_list=[];
            temperature_validation_sample_count_list=[];
            temperature_sget_out_of_set_prob_list=[];
            temperature_zone_over_prob_mat=[];
            temperature_zone_trajectory_violation_prob_mat=[];
            temperature_zone_over_max_violation_mat=[];
            re_over_prob_list=[];
            re_trajectory_violation_prob_list=[];
            re_over_max_violation_list=[];
            re_sget_out_of_set_prob_list=[];
            constraint_over_prob_list=[];
            H_DC_solution_mat=[];
            P_cool_solution_mat=[];
            P_DC_solution_mat=[];
            P_zone_solution_tensor=[];
            P_cool_zone_solution_tensor=[];
            temperature_res_list=[];

            ops = sdpsettings('solver', 'gurobi+', 'allowmilp', 1, 'gurobi.QCPDual', 1,'verbose', 0,'gurobi.MIPGap',0.003);
            for k=1:K-1
                eps_k=pes_f(2)-(pes_f(2)-opt_f(2))*K^0.75*k^0.25/K;
                t1=clock;
                sol = optimize([Constraints objmean<=eps_k], obj1, ops);
                t2=clock;
                deltatime=etime(t2,t1);
                solvingtime_list=[solvingtime_list deltatime];
                valuesave
                m=m+1;
                if mod(m,10)==0
                    disp(m)
                end
            end

            ress=[obj_1_list obj_2_o_list obj_2_p_list obj_2_mean_list];
            C=[ress(:,1),ress(:,4)];
            objsssss=[Payoff_mat' C'];
            cost_mat=[cost_mat objsssss(1,:)];
            wind_mat=[wind_mat objsssss(2,:)];
        end
        temperature_zone_over_prob_list = temperature_zone_over_prob_mat;
        temperature_zone_trajectory_violation_prob_list = temperature_zone_trajectory_violation_prob_mat;
        temperature_zone_over_max_violation_list = temperature_zone_over_max_violation_mat;
        resDir = fullfile(scriptDir, 'results');
        if ~exist(resDir, 'dir'); mkdir(resDir); end
        filename = fullfile(resDir, ['MultiDCMG3_DataInter_NoGuess_e', int2str(10 * epsilon), '_d', int2str(10 * delta), '.mat']);
        save(filename, 'cost_mat', 'wind_mat', 'Payoff_mat', ...
            'temperature_over_prob_list', ...
            'temperature_over_time_ratio_list', ...
            'temperature_trajectory_violation_prob_list', ...
            'temperature_over_avg_time_ratio_list', ...
            'temperature_over_max_violation_list', ...
            'temperature_validation_sample_count_list', ...
            'temperature_sget_out_of_set_prob_list', ...
            'temperature_zone_over_prob_list', ...
            'temperature_zone_trajectory_violation_prob_list', ...
            'temperature_zone_over_max_violation_list', ...
            'temperature_zone_over_prob_mat', ...
            'temperature_zone_trajectory_violation_prob_mat', ...
            'temperature_zone_over_max_violation_mat', ...
            're_over_prob_list', 're_trajectory_violation_prob_list', ...
            're_over_max_violation_list', 're_sget_out_of_set_prob_list', ...
            'constraint_over_prob_list', ...
            'H_DC_solution_mat', 'P_cool_solution_mat', 'P_DC_solution_mat', ...
            'P_re_solution_mat', ...
            'P_zone_solution_tensor', 'P_cool_zone_solution_tensor', ...
            'Theta_err_samples', 'temperature_sf_margin', ...
            'temperature_sf_order', 'wind_err_samples', ...
            're_sf_margin', 're_sf_order', ...
            'chance_target_epsilon', 'sf_confidence_delta', ...
            'solvingtime_list');
        disp(['Saved: ', filename]);

    end
end
