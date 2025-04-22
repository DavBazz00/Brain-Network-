function [A, CoordTable, Coord] = AdjacencyMatrix(edgeFile, nodeFile)
%BUILDADJACENCYMATRIX Loads and processes the adjacency matrix and node coordinate table.
%
%   [A, CoordTable, Coord] = buildAdjacencyMatrix(edgeFile, nodeFile)
%
%   INPUTS:
%       edgeFile - String, full path to the CSV file containing the adjacency matrix.
%       nodeFile - String, full path to the CSV file containing node coordinates and region info.
%
%   OUTPUTS:
%       A         - The adjacency matrix with the first row removed.
%       CoordTable- The table of node coordinates updated with a new 'Label' column.
%       Coord     - A numeric matrix extracted from the first three columns of CoordTable.
%
%   The function performs the following steps:
%       1. Loads the adjacency matrix from edgeFile and removes the first row.
%       2. Reads the node coordinate table from nodeFile.
%       3. Extracts the numeric coordinates (first three columns) from the table.


    % Load the adjacency matrix and remove the first row
    A = load(edgeFile);
    A(1, :) = []; 
    
    % Read the node coordinate table
    CoordTable = readtable(nodeFile);
    
    % Extract numeric coordinates (assumed to be in the first three columns)
    Coord = table2array(CoordTable(:, 1:3));
end
