clear
clc

patient_ids = [1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];

selected_numbers = zeros(1, 57);  % Initialize an array to store the selected numbers

for i = 1:57
    index = randi(71 - i + 1);  % Generate a random index in the remaining numbers
    selected_numbers(i) = patient_ids(index);  % Store the selected number
    patient_ids(index) = [];  % Remove the selected number from the list
end

selected_numbers = sort(selected_numbers);  % Sort the selected numbers in ascending order

disp(selected_numbers)
