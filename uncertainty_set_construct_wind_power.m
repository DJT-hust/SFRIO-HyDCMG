function [wind_power_max,wind_power_min]=uncertainty_set_construct_wind_power(wind_power_forecast,a)

wind_power_max=a*wind_power_forecast;
wind_power_min=-wind_power_max;
