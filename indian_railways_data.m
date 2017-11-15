clear all
clc
power = readtable("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/irn_power.csv");

grid.index =table2array(power(:,1));
grid.id = table2cell(power(:,2));
grid.line = table2cell(power(:,3)); 
clear power;

natural = readtable("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/irn_tsunami.csv");
tsunami.index =table2array(natural(:,1));
tsunami.id = table2cell(natural(:,2));
tsunami.line = table2cell(natural(:,3)); 
clear natural;

target  = readtable("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/irn_cyber.csv");
cyber.index =table2array(target(:,1));
cyber.id = table2cell(target(:,2));
cyber.line = table2cell(target(:,3));
clear target;

adj_IRN = csvread("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/adjacencyindianrailways.csv");
adj_IRN_weighted = csvread("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/weightedadjacencyirn.csv");

stations = readtable("datawithlatlong.csv");
nodes_irn.index =table2array(stations(:,1));
nodes_irn.id = table2array(stations(:,2));
nodes_irn.lat = table2array(stations(:,3));
nodes_irn.lon = table2array(stations(:,4));
clear stations
save("/Users/udit/Desktop/mbta/mbta matlab codes June 2017/Indian_datafiles");
clear filename
%clear all
close all



