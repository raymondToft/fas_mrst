diary 'TestCaseFCycle_4a_to_4e_FAS.txt'

%%
fprintf('\n FASTestCase_4a_12_1to5_F_ho \n');
for i=1:5
   run FASTestCase_4a_12_1to5_F_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_4b_16_1to5_F_ho \n');
for i=1:5
   run FASTestCase_4b_16_1to5_F_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\n FASTestCase_4c_20_1to5_F_ho \n');
for i=1:5
   run FASTestCase_4c_20_1to5_F_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_4d_24_1to5_F_ho \n');
for i=1:5
   run FASTestCase_4d_24_1to5_F_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_4e_28_1to5_F_ho \n');
for i=1:5
   run FASTestCase_4e_28_1to5_F_ho
end
%%
figure; % Lazy way of getting notice of the end of the simulations
diary off
