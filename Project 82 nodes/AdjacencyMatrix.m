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
%       4. Defines a set of sequences and creates labels for node pairs in the format 'nL' and 'nR'.
%       5. Updates the coordinate table with the new 'Label' column.
%

    % Load the adjacency matrix and remove the first row
    A = load(edgeFile);
    A(1, :) = []; 
    
    % Read the node coordinate table
    CoordTable = readtable(nodeFile);
    
    % Extract numeric coordinates (assumed to be in the first three columns)
    Coord = table2array(CoordTable(:, 1:3));
    
    % Define the desired sequences
    seqA = 1:11;      % 1, 2, ..., 11
    seqB = 17:30;     % 17, 18, ..., 30
    seqC = 32;        % 32
    secD = 34:48;
    
    % Concatenate all sequences into one vector
    allNumbers = [seqA, seqB, seqC, secD];  % This gives the desired order
    
    % Calculate the number of pairs
    numPairs = length(allNumbers);
    
    % Preallocate cell array for the new labels (each pair has two labels)
    Labels = cell(numPairs * 2, 1);
    
    % Generate labels in the format "nL" and "nR" for each number in allNumbers
    for i = 1:numPairs
        Labels{2*i - 1} = sprintf('%dL', allNumbers(i));
        Labels{2*i}     = sprintf('%dR', allNumbers(i));
    end
    
    % Update the Label column in the table
    CoordTable.("Label") = Labels;
end
