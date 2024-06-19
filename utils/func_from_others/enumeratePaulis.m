% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
%
%  Copyright © 2019 James R. Seddon
%  This file is part of project: Channel Magic
%
function [pauli_array,pauli_labels,Z_paulis,Z_indices] = enumeratePaulis(n);
%%% Create array containing all Paulis (+ phase) for n-qubit system.
%%% pauli_labels outputs the names of the Paulis eg 'IIXI' 'ZXZX' and so
%%% on. Z_paulis outputs a vector that flags which matrices in the array
%%% are composed from only Zs and Is, and contain at least 1 Z. eg. 'IIIZ'
%%% and 'ZZZZ' would be flagged, but 'IIII', 'YIYI' and 'XZZZ' would not.
%%% Z_indices gives the indices of the rows flagged in Z_paulis.

pauliSingle(:,:,1) = eye(2); % I

pauliSingle(:,:,2) = [ 0 1; % X
                       1 0];
pauliSingle(:,:,3) = 1i*[0 -1; % Y
                         1 0];
pauliSingle(:,:,4)  = [1 0;    % Z
                       0 -1];

pauliName{1} = 'I';

pauliName{2} = 'X';
                    
pauliName{3} = 'Y';
pauliName{4}  = 'Z';
                   
                 
number_paulis = 4^n;

Z_paulis = zeros(number_paulis,1);
pauli_labels = {};
Z_indices = [];  

for mm = 1:number_paulis
    this_pauli_encoding = fliplr(de2bi(mm-1,n,4));
    current_matrix = 1;
    current_string = '';
    for pp = 1:size(this_pauli_encoding,2);
        pauli_index = this_pauli_encoding(pp)+1;
        current_matrix = kron(current_matrix,pauliSingle(:,:,pauli_index));
        current_string = [current_string pauliName{pauli_index}];
    end
    Z_paulis(mm) = ismember(3,this_pauli_encoding) && ...
                ~ismember(1,this_pauli_encoding) && ...
                    ~ismember(2,this_pauli_encoding);
    if Z_paulis(mm)
        Z_indices = [Z_indices; mm];
    end
    pauli_array(:,:,mm) = current_matrix;
    pauli_labels = [pauli_labels; current_string];
end

end
