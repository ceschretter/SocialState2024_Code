function [bin_vect,t_vect] = BinTimeseries(timeseries,samplerate,bin_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



t = [1:length(timeseries)]./samplerate;
edges = (0:bin_size:max(t));
bins = discretize(t, edges);
t_vect = 0.5*(edges(1:end-1)+edges(2:end));
bin_vect = grpstats(timeseries,bins);



end